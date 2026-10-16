// ===================== CONFIG =====================
var STATE_FIPS = '53'; // Washington
var DATE_RANGE = ['2023-05-01', '2023-06-30'];
var SAMPLING_INTERVAL_M = 1000;

var TRIPS = {
  pullman_spokane_2vert: ee.Feature(
    ee.Geometry.LineString([
      [-117.1540, 46.7324], // Pullman
      [-117.4260, 47.6588]  // Spokane
    ]),
    {trip_id: 'Pullman_to_Spokane_2_vert'}
  ),
  pullman_high_to_wsu: ee.Feature(
    ee.Geometry.LineString([
      [-117.1863, 46.7483], // PHS
      [-117.1628, 46.7313]  // Terrell Library
    ]),
    {trip_id: 'Pullman_High_to_WSU'}
  ),
  pullman_spokane_us195: ee.Feature(
    ee.Geometry.LineString([
      [-117.1818, 46.7313], // Pullman
      [-117.2952, 46.7815], // Halfway Pullman-Colfax
      [-117.3644, 46.8802], // Colfax
      [-117.3810, 47.0630], // Steptoe
      [-117.3660, 47.1300], // Thornton (approx)
      [-117.3802, 47.4289], // Spangle
      [-117.4260, 47.6588]  // Spokane
    ]),
    {trip_id: 'Pullman_to_Spokane_7_vert'}
  )
};
var ACTIVE_TRIP_KEY = 'pullman_spokane_us195'; // change this to switch trips

// ===================== SPATIAL CHUNKING UTILITIES (currently unused below) =====================
// Kept here in case large-region processing is needed later; not called
// in the current workflow, which processes all of WA in one pass.

function createSpatialGrid(geometry, cellSize) {
  // cellSize in degrees (e.g., 0.1 = ~11km at this latitude)
  var bounds = geometry.bounds();
  var coords = ee.List(bounds.coordinates().get(0));
  var xMin = ee.List(coords.get(0)).get(0);
  var yMin = ee.List(coords.get(0)).get(1);
  var xMax = ee.List(coords.get(2)).get(0);
  var yMax = ee.List(coords.get(2)).get(1);

  var xSteps = ee.Number(xMax).subtract(xMin).divide(cellSize).ceil();
  var ySteps = ee.Number(yMax).subtract(yMin).divide(cellSize).ceil();

  var cells = ee.FeatureCollection(
    ee.List.sequence(0, xSteps.subtract(1)).map(function(x) {
      return ee.List.sequence(0, ySteps.subtract(1)).map(function(y) {
        var xStart = ee.Number(xMin).add(ee.Number(x).multiply(cellSize));
        var yStart = ee.Number(yMin).add(ee.Number(y).multiply(cellSize));
        var cell = ee.Geometry.Rectangle([
          xStart, yStart,
          xStart.add(cellSize), yStart.add(cellSize)
        ]);
        return ee.Feature(cell.intersection(geometry, 1), {
          grid_x: x,
          grid_y: y
        });
      });
    }).flatten()
  );

  return cells.filter(ee.Filter.neq('area', 0));
}

function processCountyInChunks(county, envFunction) {
  var grid = createSpatialGrid(county.geometry(), 0.05); // ~5.5km cells
  return grid.map(function(cell) { return envFunction(cell); });
}

function getAdaptiveScale(geometry) {
  var area = geometry.area(1); // sq meters
  return ee.Number(area).sqrt().divide(1000).max(30).min(500);
}

// ===================== ENVIRONMENTAL LAYER FUNCTIONS =====================
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

function computeNO2(feature, startDate, endDate) {
  return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('tropospheric_NO2_column_number_density')
    .mean()
    .clip(feature.geometry());
}

function computeLST(feature, startDate, endDate) {
  return ee.ImageCollection('MODIS/061/MOD11A2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('LST_Day_1km')
    .mean()
    .multiply(0.02)        // scale factor
    .subtract(273.15)      // Kelvin -> Celsius
    .clip(feature.geometry());
}

// PM2.5: multiple source options, pick one via computePM25()
function computePM25_GEOSCF(feature, startDate, endDate) {
  return ee.ImageCollection('NASA/GEOS-CF/v1/rpl/tavg1hr')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('PM25_RH35_GCC')
    .mean()
    .clip(feature.geometry());
}

// NOTE: different signature (takes a single `year`, not a date range) —
// not drop-in compatible with the others; keep that in mind if you swap it in.
function computePM25_WashU(feature, year) {
  return ee.ImageCollection('projects/sat-io/open-datasets/GLOBAL_SATELLITE_PM25')
    .filter(ee.Filter.eq('year', year))
    .first()
    .clip(feature.geometry());
}

function computePM25_GHAP(feature, startDate, endDate) {
  return ee.ImageCollection('projects/sat-io/open-datasets/GHAP/GHAP_D1K_PM25')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .mean()
    .clip(feature.geometry());
}

function computePM25(feature, startDate, endDate) {
  return computePM25_GEOSCF(feature, startDate, endDate); // default source
}

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

function calculateDetailedExposure(tripFeature, envLayer, samplingInterval, sampleScale) {
  var points = createPointsAlongLine(tripFeature.geometry(), samplingInterval);

  var sampledPoints = points.map(function(point) {
    var value = envLayer.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry().buffer(25),
      scale: sampleScale,
      maxPixels: 1e9
    });
    return point.set('exposure', value.values().get(0));
  });

  return tripFeature.set({
    mean_exposure: sampledPoints.aggregate_mean('exposure'),
    max_exposure: sampledPoints.aggregate_max('exposure'),
    min_exposure: sampledPoints.aggregate_min('exposure'),
    sampling_points_count: sampledPoints.size(),
    sampling_points_fc: sampledPoints
  });
}

function calculateMultiExposure(tripFeature, layersDict, samplingInterval, sampleScale) {
  var tripWithExposures = tripFeature;

  Object.keys(layersDict).forEach(function(name) {
    var safeName = name.replace(/\./g, '');
    var layer = layersDict[name];
    var exposure = calculateDetailedExposure(tripWithExposures, layer, samplingInterval, sampleScale);

    tripWithExposures = tripWithExposures
      .set('mean_' + safeName, exposure.get('mean_exposure'))
      .set('max_' + safeName, exposure.get('max_exposure'))
      .set('min_' + safeName, exposure.get('min_exposure'))
      .set('sampling_points_count_' + safeName, exposure.get('sampling_points_count'))
      .set('sampling_points_fc_' + safeName, exposure.get('sampling_points_fc'));
  });

  return tripWithExposures;
}

// ===================== LEGEND / UI HELPERS =====================
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

// ===================== MAIN: LOAD REGION =====================
var counties = ee.FeatureCollection('TIGER/2018/Counties');
var waCounties = counties.filter(ee.Filter.eq('STATEFP', STATE_FIPS));
var waGeometry = waCounties.geometry();

Map.centerObject(waCounties, 6);
Map.addLayer(waCounties, {}, 'WA Counties');

print('Study area (km^2):', waGeometry.area(1).divide(1e6));

// ===================== MAIN: ENVIRONMENTAL LAYERS =====================
var ENV_LAYERS = [
  {
    name: 'NDVI', label: 'NDVI (Green Space)',
    image: computeNDVI(ee.Feature(waGeometry), DATE_RANGE[0], DATE_RANGE[1]),
    vis: {min: 0.2, max: 0.8, palette: ['white', 'yellow', 'green']},
    shown: false,
    exportScale: 30,
    exportDescription: 'WA_NDVI_2023_MayJune'
  },
  {
    name: 'NO2', label: 'NO2 (Air Pollution)',
    image: computeNO2(ee.Feature(waGeometry), DATE_RANGE[0], DATE_RANGE[1]),
    vis: {min: 0, max: 0.0002, palette: ['blue', 'yellow', 'red']},
    shown: true,
    exportScale: 1000,
    exportDescription: 'WA_NO2_2023_MayJune'
  },
  {
    name: 'Temperature', label: 'Temperature (°C)',
    image: computeLST(ee.Feature(waGeometry), DATE_RANGE[0], DATE_RANGE[1]),
    vis: {min: 10, max: 30, palette: ['blue', 'white', 'red']},
    shown: false,
    exportScale: 1000,
    exportDescription: 'WA_LST_2023_MayJune'
  },
  {
    name: 'PM25', label: 'PM2.5 (ug/m3)',
    image: computePM25(ee.Feature(waGeometry), DATE_RANGE[0], DATE_RANGE[1]),
    vis: {min: 0, max: 30, palette: ['green', 'yellow', 'orange', 'red', 'purple']},
    shown: false,
    exportScale: 250, // native GEOS-CF ~25km; saved finer per original script
    exportDescription: 'WA_PM25_2023_MayJune'
  }
];

print('Computing environmental layers...');
ENV_LAYERS.forEach(function(layer) {
  Map.addLayer(layer.image, layer.vis, layer.label, layer.shown);
});
print('Environmental layers computed successfully! Toggle layers in the Layers panel.');

// ===================== MAIN: EXPORT ENVIRONMENTAL LAYERS =====================
ENV_LAYERS.forEach(function(layer) {
  Export.image.toDrive({
    image: layer.image,
    description: layer.exportDescription,
    scale: layer.exportScale,
    region: waGeometry,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
});

// ===================== MAIN: LEGENDS + CONTROL PANEL =====================
var legends = {};
ENV_LAYERS.forEach(function(layer) {
  legends[layer.name] = createColorLegend(
    layer.label, layer.vis.palette, layer.vis.min, layer.vis.max, 'bottom-left'
  );
});

function showOnlyLegend(name) {
  Object.keys(legends).forEach(function(key) { Map.remove(legends[key]); });
  if (name) Map.add(legends[name]);
}

var controlPanel = ui.Panel({
  widgets: [
    ui.Label('LEGEND CONTROLS', {fontWeight: 'bold', fontSize: '16px', margin: '0 0 10px 0'})
  ].concat(ENV_LAYERS.map(function(layer) {
    return ui.Button({
      label: 'Show ' + layer.label + ' Legend',
      onClick: function() { showOnlyLegend(layer.name); }
    });
  })).concat([
    ui.Button({
      label: 'Hide All Legends',
      onClick: function() { showOnlyLegend(null); },
      style: {color: 'red'}
    })
  ]),
  style: {position: 'top-right', padding: '8px', backgroundColor: 'white'}
});

Map.add(controlPanel);

// ===================== MAIN: REGION-WIDE STATISTICS =====================
var layerStats = {};
ENV_LAYERS.forEach(function(layer) {
  layerStats[layer.name] = layer.image.reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.minMax(), sharedInputs: true}),
    geometry: waGeometry,
    scale: 1000,
    maxPixels: 1e13
  });
});

print('ENVIRONMENTAL LAYER STATISTICS');
ENV_LAYERS.forEach(function(layer) {
  print(layer.label + ' Stats:', layerStats[layer.name]);
});

// ===================== MAIN: TRIP EXPOSURE =====================
var trip = TRIPS[ACTIVE_TRIP_KEY];
print('Distance (km):', trip.geometry().length({maxError: 1}).divide(1000));
Map.addLayer(trip, {color: 'black', width: 3}, 'Sample Trip');

var layersForExposure = {
  NDVI: ENV_LAYERS[0].image,
  NO2: ENV_LAYERS[1].image,
  Temperature: ENV_LAYERS[2].image,
  PM25: ENV_LAYERS[3].image
};

var EXPOSURE_SAMPLE_SCALE = 250; // matches PM25 native resolution

var multiExposure = calculateMultiExposure(trip, layersForExposure, SAMPLING_INTERVAL_M, EXPOSURE_SAMPLE_SCALE);
print('Trip Exposure Summary (All Layers):', multiExposure);

// ===================== MAIN: TRIP STATS TABLE + CHARTS =====================
var layerNames = Object.keys(layersForExposure);

var statsFc = ee.FeatureCollection(
  layerNames.map(function(name) {
    var meanVal = multiExposure.get('mean_' + name);
    var minVal = multiExposure.get('min_' + name);
    var maxVal = multiExposure.get('max_' + name);
    return ee.List([
      ee.Feature(null, {layer: name, stat: 'Min', value: minVal}),
      ee.Feature(null, {layer: name, stat: 'Mean', value: meanVal}),
      ee.Feature(null, {layer: name, stat: 'Max', value: maxVal})
    ]);
  }).reduce(function(a, b) { return ee.List(a).cat(b); }, ee.List([]))
);

print('Trip stats table (layer, stat, value):', statsFc);

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
