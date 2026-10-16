// ===================== CONFIG =====================
var REGION = Sylhet;
var EXPORT_FOLDER = 'GEE Sylhet Flood';
var SCALE = 30;
var MAX_PIXELS = 1E8;

var DATE_RANGES = {
  before: ['2022-04-10', '2022-04-20'],
  after:  ['2022-05-20', '2022-05-30']
};

var FLOOD_THRESHOLD_DB = -18; // VH backscatter threshold for water/flood

// ===================== SAR PREPROCESSING FUNCTIONS =====================
function toNatural(image) {
  return ee.Image(10).pow(image.divide(10));
}

function toDB(image) {
  return image.log10().multiply(10);
}

// Refined Lee speckle filter
function RefinedLee(img) {
  var weights = ee.Kernel.square({ radius: 7, units: 'pixels', normalize: true });

  var mean = img.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: weights
  });
  var variance = img.reduceNeighborhood({
    reducer: ee.Reducer.variance(),
    kernel: weights
  });

  var ci = variance.sqrt().divide(mean); // coefficient of variation
  var cu = 0.523; // noise threshold, uniform regions
  var cmax = 1.73; // max threshold, highly textured areas (unused currently)

  var alpha = ee.Image(1.0)
    .subtract(ee.Image.constant(cu).multiply(cu).divide(ci.multiply(ci)))
    .max(0).min(1);

  return img.multiply(alpha).add(mean.multiply(ee.Image(1).subtract(alpha)));
}

function speckleFilterDB(image) {
  return ee.Image(toDB(RefinedLee(toNatural(image))));
}

// ===================== LOAD SENTINEL-1 =====================
var collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.or(
    ee.Filter.eq('orbitProperties_pass', 'DESCENDING'),
    ee.Filter.eq('orbitProperties_pass', 'ASCENDING')
  ))
  .filterBounds(REGION);

var before = collection.filter(ee.Filter.date(DATE_RANGES.before[0], DATE_RANGES.before[1]));
var after  = collection.filter(ee.Filter.date(DATE_RANGES.after[0], DATE_RANGES.after[1]));

print('Before collection:', before);
print('After collection:', after);

var before_flood = before.select('VH').mosaic().clip(REGION);
var after_flood  = after.select('VH').mosaic().clip(REGION);

// ===================== FILTER + CLASSIFY =====================
var before_filter = speckleFilterDB(before_flood);
var after_filter  = speckleFilterDB(after_flood);

// Flood: was dry before, is water-like after
var flood = before_filter.gt(FLOOD_THRESHOLD_DB).and(after_filter.lt(FLOOD_THRESHOLD_DB));
var flood_mask = flood.updateMask(flood.eq(1));

// Permanent water: water-like in both periods
var water = before_filter.lt(FLOOD_THRESHOLD_DB).and(after_filter.lt(FLOOD_THRESHOLD_DB));
var water_mask = water.updateMask(water.eq(1));

// ===================== VISUALIZATION =====================
Map.centerObject(REGION, 10);
Map.addLayer(before_filter, {min: -25, max: 0}, 'Before Flood Filtered');
Map.addLayer(after_filter, {min: -25, max: 0}, 'After Flood Filtered');
Map.addLayer(flood_mask, {palette: ['Red']}, 'Flood Inundation');
Map.addLayer(water_mask, {palette: ['Blue']}, 'Permanent Waterbodies');

// ===================== EXPORT =====================
var exportsList = [
  { image: water_mask, description: 'Permanent_Water' },
  { image: flood,       description: 'Sylhet_Flood' },     // unfiltered boolean (0/1 everywhere)
  { image: water,       description: 'Water' },             // unfiltered boolean (0/1 everywhere)
  { image: flood_mask,  description: 'Flood_Mask' }         // masked (only flood pixels)
];

exportsList.forEach(function(item) {
  Export.image.toDrive({
    image: item.image,
    description: item.description,
    folder: EXPORT_FOLDER,
    region: REGION,
    scale: SCALE,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });
});