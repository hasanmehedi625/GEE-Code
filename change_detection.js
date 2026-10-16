// ===================== CONFIG =====================
var REGION = Dhaka;
var SCALE = 30;
var BANDS = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];

var YEARS = {
  2020: {
    dateRange: ['2020-03-01', '2020-03-31'],
    samples: {
      Waterbodies: Waterbodies_2020,
      Vegetation: Vegetation_2020,
      Bareland: Bareland_2020,
      Agricultural_land: Agricultural_land_2020,
      Built_up: Built_up_2020
    },
    exportFolder: 'MS Thesis GEE Map 2020'
  },
  2025: {
    dateRange: ['2025-03-01', '2025-03-31'],
    samples: {
      Waterbodies: Waterbodies_2025,
      Vegetation: Vegetation_2025,
      Bareland: Bareland_2025,
      Agricultural_land: Agricultural_land_2025,
      Built_up: Built_up_2025
    },
    exportFolder: 'MS Thesis GEE Map 2025'
  }
};

// class codes: 0 water, 1 vegetation, 2 bareland, 3 agricultural land, 4 built-up
var CLASS_CODES = {
  Waterbodies: 0,
  Vegetation: 1,
  Bareland: 2,
  Agricultural_land: 3,
  Built_up: 4
};

var LANDCOVER_PALETTE = [
  '419BDF',  // 0: Water
  '397D49',  // 1: Vegetation
  'CBAF07',  // 2: Bareland
  '067317',  // 3: Agricultural land
  'CC0013'   // 4: Built-up
];

// ===================== CLASSIFY ONE YEAR =====================
function classifyYear(year) {
  var cfg = YEARS[year];

  var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(REGION)
    .filterMetadata('CLOUD_COVER', 'less_than', 1)
    .filterDate(cfg.dateRange[0], cfg.dateRange[1])
    .median()
    .clip(REGION);

  var labeled = Object.keys(cfg.samples).map(function(key) {
    return cfg.samples[key].map(function(f) {
      return f.set('class', CLASS_CODES[key]);
    });
  });
  var samples = labeled.reduce(function(acc, fc) { return acc.merge(fc); });
  print('Samples ' + year + ':', samples);

  var training = image.select(BANDS).reduceRegions({
    collection: samples,
    reducer: ee.Reducer.mean(),
    scale: SCALE
  }).randomColumn();

  var trainSet = training.filter(ee.Filter.lt('random', 0.8));
  var testSet = training.filter(ee.Filter.gte('random', 0.8));

  var svm = ee.Classifier.libsvm().train({
    features: trainSet,
    classProperty: 'class',
    inputProperties: BANDS
  });

  var classified = image.select(BANDS).classify(svm);

  Map.addLayer(classified, {min: 0, max: 4, palette: LANDCOVER_PALETTE},
    'SVM Land Cover Classification ' + year);

  var validated = testSet.classify(svm);
  var matrix = validated.errorMatrix('class', 'classification');
  print('SVM Error Matrix ' + year + ':', matrix);
  print('SVM Accuracy ' + year + ':', matrix.accuracy());

  Export.image.toDrive({
    image: classified,
    description: 'LULC_SVM_L8_Dhaka_' + year,
    folder: cfg.exportFolder,
    region: REGION,
    scale: SCALE,
    maxPixels: 1e10
  });

  return classified;
}

Map.addLayer(REGION);
Map.centerObject(REGION, 10);

var class_2020 = classifyYear(2020);
var class_2025 = classifyYear(2025);

// ===================== CHANGE / TRANSITION MAP =====================
// Code = class_2020 * 10 + class_2025 (each class is 0-4)
var changeMap = class_2020.multiply(10).add(class_2025);

// Lookup table of transitions worth labeling; anything not listed -> 0 (unlabeled/no-change)
// NOTE: this only labels a subset of possible transitions (the ones you had).
// Many real transitions (e.g. water<->anything, vegetation->agri, etc.) are
// NOT in this table and will fall into bucket 0 along with true "no change" —
// see caveat below.
var TRANSITIONS = {
  11: 0, // Vegetation -> Vegetation (no change)
  12: 4, // Vegetation -> Bareland
  14: 1, // Vegetation -> Built-up
  22: 0, // Bareland -> Bareland (no change)
  24: 2, // Bareland -> Built-up
  33: 0, // Agri -> Agri (no change)
  34: 3, // Agri -> Built-up
  32: 5  // Agri -> Bareland
};

var transitionExpr = Object.keys(TRANSITIONS).map(function(code) {
  return 'transition == ' + code + ' ? ' + TRANSITIONS[code];
}).join(' : ') + ' : 0';

var transitions = changeMap.expression(transitionExpr, { transition: changeMap });

var CHANGE_PALETTE = [
  'ffffff', // 0 - No change / unlabeled
  'ff0000', // 1 - Vegetation -> Built-up
  'ff9900', // 2 - Bareland -> Built-up
  'cc00ff', // 3 - Agricultural land -> Built-up
  '00ffff', // 4 - Vegetation -> Bareland
  '00ff00'  // 5 - Agricultural land -> Bareland
];

var clippedTransitions = transitions.clip(REGION);
Map.addLayer(clippedTransitions, {min: 0, max: 5, palette: CHANGE_PALETTE}, 'Change Detection Detailed');

Export.image.toDrive({
  image: clippedTransitions,
  description: 'LandCoverChange_Dhaka_L8',
  folder: 'Change Detection 2020 to 2025 L8',
  fileNamePrefix: 'Change_2020_2025_Dhaka',
  region: REGION,
  scale: SCALE,
  maxPixels: 1e13,
  crs: 'EPSG:4326'
});

// ===================== AREA STATISTICS =====================
function calculateArea(image, region, scale) {
  var pixelArea = ee.Image.pixelArea().divide(1e6); // sq km
  var areaImage = pixelArea.addBands(image.rename('class'));

  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum().group({
      groupField: 1,
      groupName: 'class_type'
    }),
    geometry: region,
    scale: scale,
    maxPixels: 1e13
  });

  var statsList = ee.List(stats.get('groups')).map(function(item) {
    var dict = ee.Dictionary(item);
    var classValue = dict.get('class_type');
    var label = ee.String('Class_').cat(ee.Number(classValue).format());
    return ee.Feature(null, {
      'Class': label,
      'Area_sq_km': dict.get('sum')
    });
  });

  return ee.FeatureCollection(statsList);
}

print('Area per class (2020) in sq.km:', calculateArea(class_2020, REGION, SCALE));
print('Area per class (2025) in sq.km:', calculateArea(class_2025, REGION, SCALE));

var changeStats = calculateArea(transitions.rename('change_type'), REGION, SCALE);
print('Land Cover Change Statistics (2020–2025) in sq.km:', changeStats);