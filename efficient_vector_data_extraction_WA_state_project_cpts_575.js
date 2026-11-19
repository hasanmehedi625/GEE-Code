// ===================== STUDY AREA =====================
var counties = ee.FeatureCollection('TIGER/2018/Counties');
var waCounties = counties.filter(ee.Filter.eq('STATEFP', '53'));
var waGeometry = waCounties.geometry();

Map.centerObject(waCounties, 6);
Map.addLayer(waCounties, {}, 'WA Counties');

var startDate = '2023-05-01';
var endDate = '2023-06-30';

print('Study area (km²):', waGeometry.area(1).divide(1e6));

// ===================== ENVIRONMENTAL LAYER FUNCTIONS =====================
// NDVI (green space)
function computeNDVI(feature, startDate, endDate) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
    .select(['B4', 'B8'])
    .median();

  return s2.normalizedDifference(['B8', 'B4'])
    .rename('NDVI')
    .clip(feature.geometry());
}

// NO2 (air quality)
function computeNO2(feature, startDate, endDate) {
  return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('tropospheric_NO2_column_number_density')
    .mean()
    .clip(feature.geometry());
}

// Land surface temperature (°C)
function computeLST(feature, startDate, endDate) {
  return ee.ImageCollection('MODIS/061/MOD11A2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('LST_Day_1km')
    .mean()
    .multiply(0.02)      // scale factor
    .subtract(273.15)    // Kelvin to Celsius
    .clip(feature.geometry());
}

// PM2.5 (NASA GEOS-CF)
function computePM25(feature, startDate, endDate) {
  return ee.ImageCollection('NASA/GEOS-CF/v1/rpl/tavg1hr')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('PM25_RH35_GCC')
    .mean()
    .clip(feature.geometry());
}

// ===================== COMPUTE LAYERS =====================
print('Computing environmental layers...');

var ndviLayer = computeNDVI(ee.Feature(waGeometry), startDate, endDate);
var no2Layer  = computeNO2(ee.Feature(waGeometry), startDate, endDate);
var lstLayer  = computeLST(ee.Feature(waGeometry), startDate, endDate);
var pm25Layer = computePM25(ee.Feature(waGeometry), startDate, endDate);

// ===================== VISUALIZATION =====================
var ndviVis = {min: 0.2, max: 0.8, palette: ['white', 'yellow', 'green']};
var no2Vis  = {min: 0, max: 0.0002, palette: ['blue', 'yellow', 'red']};
var lstVis  = {min: 10, max: 30, palette: ['blue', 'white', 'red']};
var pm25Vis = {min: 0, max: 30, palette: ['green', 'yellow', 'orange', 'red', 'purple']};

Map.addLayer(ndviLayer, ndviVis, 'NDVI (Green Space)', false);
Map.addLayer(no2Layer, no2Vis, 'NO2 (Air Pollution)');
Map.addLayer(lstLayer, lstVis, 'Temperature (°C)', false);
Map.addLayer(pm25Layer, pm25Vis, 'PM2.5 (μg/m³)', false);

print('Environmental layers computed successfully!');
print('Toggle layers on/off in the Layers panel on the map.');

// ===================== EXPORTS =====================
Export.image.toDrive({
  image: ndviLayer,
  description: 'WA_NDVI_2023_MayJune',
  scale: 30,
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: no2Layer,
  description: 'WA_NO2_2023_MayJune',
  scale: 1000,
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: lstLayer,
  description: 'WA_LST_2023_MayJune',
  scale: 1000,   // MODIS 1 km
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: pm25Layer,
  description: 'WA_PM25_2023_MayJune',
  scale: 250,
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// ===================== LEGENDS =====================
function createColorLegend(title, palette, min, max, position) {
  var legendTitle = ui.Label(title, {fontWeight: 'bold', fontSize: '14px'});
  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0)
      .multiply((max - min) / 100.0).add(min)
      .visualize({min: min, max: max, palette: palette}),
    params: {bbox: [0, 0, 100, 10], dimensions: '100x10'},
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '20px'}
  });
  var legendLabels = ui.Panel([
    ui.Label(min.toFixed(2)),
    ui.Label(((min + max) / 2).toFixed(2)),
    ui.Label(max.toFixed(2))
  ], ui.Panel.Layout.flow('horizontal'));

  return ui.Panel({
    widgets: [legendTitle, colorBar, legendLabels],
    style: {position: position || 'bottom-left', padding: '8px 15px'}
  });
}

var ndviLegend = createColorLegend('NDVI (Green Space)', ndviVis.palette, ndviVis.min, ndviVis.max);
var no2Legend  = createColorLegend('NO2 (mol/m²)', no2Vis.palette, no2Vis.min, no2Vis.max);
var lstLegend  = createColorLegend('Temperature (°C)', lstVis.palette, lstVis.min, lstVis.max);
var pm25Legend = createColorLegend('PM2.5 (μg/m³)', pm25Vis.palette, pm25Vis.min, pm25Vis.max);

var allLegends = [ndviLegend, no2Legend, lstLegend, pm25Legend];

// Remove all legends, then show the one requested (if any)
function showLegend(legend) {
  allLegends.forEach(function(l) { Map.remove(l); });
  if (legend) Map.add(legend);
}

var controlPanel = ui.Panel({
  widgets: [
    ui.Label('LEGEND CONTROLS', {fontWeight: 'bold', fontSize: '16px', margin: '0 0 10px 0'}),
    ui.Button({label: 'Show NDVI Legend', onClick: function() { showLegend(ndviLegend); }}),
    ui.Button({label: 'Show NO2 Legend', onClick: function() { showLegend(no2Legend); }}),
    ui.Button({label: 'Show Temperature Legend', onClick: function() { showLegend(lstLegend); }}),
    ui.Button({label: 'Show PM2.5 Legend', onClick: function() { showLegend(pm25Legend); }}),
    ui.Button({label: 'Hide All Legends', onClick: function() { showLegend(null); }, style: {color: 'red'}})
  ],
  style: {position: 'top-right', padding: '8px', backgroundColor: 'white'}
});

Map.add(controlPanel);

// ===================== LAYER STATISTICS (WA) =====================
function layerStats(image) {
  return image.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.minMax(),
      sharedInputs: true
    }),
    geometry: waGeometry,
    scale: 1000,
    maxPixels: 1e13
  });
}

print('ENVIRONMENTAL LAYER STATISTICS');
print('NDVI (Green Space) Stats:', layerStats(ndviLayer));
print('NO2 (Air Pollution) Stats:', layerStats(no2Layer));
print('Temperature (°C) Stats:', layerStats(lstLayer));
print('PM2.5 Stats:', layerStats(pm25Layer));

// ===================== TRIP =====================
/*
 == TRIP LOGIC SECTION ==
 To run a different trip, comment out the active one and uncomment another.
*/

/*
var trip = ee.Feature(
  ee.Geometry.LineString([
    [-117.1540, 46.7324],   // Pullman
    [-117.4260, 47.6588]    // Spokane
  ]),
  {'trip_id': 'Pullman_to_Spokane_2_vert'}
);
*/

/*
var trip = ee.Feature(
  ee.Geometry.LineString([
    [-117.1863, 46.7483],   // PHS
    [-117.1628, 46.7313]    // Terrell Library
  ]),
  {'trip_id': 'Pullman_High_to_WSU'}
);
*/

// US 195 approx
var trip = ee.Feature(
  ee.Geometry.LineString([
    [-117.1818, 46.7313],   // Pullman
    [-117.2952, 46.7815],   // Halfway Pullman-Colfax
    [-117.3644, 46.8802],   // Colfax
    [-117.3810, 47.0630],   // Steptoe
    [-117.3660, 47.1300],   // Thornton (approx)
    [-117.3802, 47.4289],   // Spangle
    [-117.4260, 47.6588]    // Spokane
  ]),
  {'trip_id': 'Pullman_to_Spokane_7_vert'}
);

print('Distance (km):', trip.geometry().length({'maxError': 1}).divide(1000));
Map.addLayer(trip, {color: 'black', width: 3}, 'Sample Trip');

// ===================== TRIP EXPOSURE FUNCTIONS =====================
function createPointsAlongLine(lineGeometry, intervalMeters) {
  var length = lineGeometry.length();
  var distances = ee.List.sequence(0, length, intervalMeters);

  var points = distances.map(function(dist) {
    var pt = lineGeometry.cutLines([dist]).geometries().get(0);
    return ee.Feature(ee.Geometry(pt));
  });

  return ee.FeatureCollection(points);
}

function calculateDetailedExposure(tripFeature, envLayer, samplingInterval) {
  var points = createPointsAlongLine(tripFeature.geometry(), samplingInterval);

  var sampledPoints = points.map(function(point) {
    var value = envLayer.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry().buffer(25),
      scale: 250,
      maxPixels: 1e9
    });
    return point.set('exposure', value.values().get(0));
  });

  return tripFeature.set({
    'mean_exposure': sampledPoints.aggregate_mean('exposure'),
    'max_exposure': sampledPoints.aggregate_max('exposure'),
    'min_exposure': sampledPoints.aggregate_min('exposure'),
    'sampling_points_count': sampledPoints.size(),
    'sampling_points_fc': sampledPoints
  });
}

function calculateMultiExposure(tripFeature, layersDict, samplingInterval) {
  var result = tripFeature;

  Object.keys(layersDict).forEach(function(name) {
    var exposure = calculateDetailedExposure(tripFeature, layersDict[name], samplingInterval);

    result = result
      .set('mean_' + name, exposure.get('mean_exposure'))
      .set('max_' + name, exposure.get('max_exposure'))
      .set('min_' + name, exposure.get('min_exposure'))
      .set('sampling_points_count_' + name, exposure.get('sampling_points_count'))
      .set('sampling_points_fc_' + name, exposure.get('sampling_points_fc'));
  });

  return result;
}

// ===================== MULTI-EXPOSURE CALCULATION =====================
var layers = {
  'NDVI': ndviLayer,
  'NO2': no2Layer,
  'Temperature': lstLayer,
  'PM25': pm25Layer
};

var samplingInterval = 1000; // meters

var multiExposure = calculateMultiExposure(trip, layers, samplingInterval);
print('Trip Exposure Summary (All Layers):', multiExposure);

// ===================== TRIP STATS TABLE + CHARTS =====================
var layerNames = ['NDVI', 'NO2', 'Temperature', 'PM25'];

// One feature per layer/stat: Min, Mean, Max
var statFeatures = [];
layerNames.forEach(function(name) {
  statFeatures.push(ee.Feature(null, {layer: name, stat: 'Min',  value: multiExposure.get('min_' + name)}));
  statFeatures.push(ee.Feature(null, {layer: name, stat: 'Mean', value: multiExposure.get('mean_' + name)}));
  statFeatures.push(ee.Feature(null, {layer: name, stat: 'Max',  value: multiExposure.get('max_' + name)}));
});
var statsFc = ee.FeatureCollection(statFeatures);

print('Trip stats table (layer, stat, value):', statsFc);

// Small chart per layer
layerNames.forEach(function(name) {
  var layerFilter = statsFc.filter(ee.Filter.eq('layer', name));
  var smallChart = ui.Chart.feature.byFeature(layerFilter, 'stat', 'value')
    .setChartType('ColumnChart')
    .setOptions({
      title: name + ' — Min / Mean / Max',
      legend: {position: 'none'},
      hAxis: {title: ''},
      vAxis: {title: 'Value'}
    });
  print(smallChart);
});
