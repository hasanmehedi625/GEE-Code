// ===================== CONFIG =====================
var REGION = Dhaka;
var DATE_RANGE = ['2025-03-01', '2025-03-31'];
var EXPORT_FOLDER = 'MS Thesis GEE Map 2025';

var PALETTE = [
  '419BDF',  // 0: Water
  '397D49',  // 1: Vegetation
  'CBAF07',  // 2: Bareland
  '067317',  // 3: Agricultural land
  'CC0013'   // 4: Built-up
];

var SENSORS = {
  L8: {
    bands: ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],
    scale: 30,
    getImage: function() {
      return ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(REGION)
        .filterMetadata('CLOUD_COVER', 'less_than', 1)
        .filterDate(DATE_RANGE[0], DATE_RANGE[1])
        .median()
        .clip(REGION);
    }
  },
  S2: {
    bands: ['B2','B3','B4','B8','B11','B12'],
    scale: 10,
    getImage: function() {
      return ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(REGION)
        .filterDate(DATE_RANGE[0], DATE_RANGE[1])
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1))
        .median()
        .clip(REGION);
    }
  },
  S1: {
    bands: ['VV', 'VH', 'VV_div_VH'],
    scale: 10,
    getImage: function() {
      var s1 = ee.ImageCollection("COPERNICUS/S1_GRD")
        .filterBounds(REGION)
        .filterDate(DATE_RANGE[0], DATE_RANGE[1])
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.eq('resolution_meters', 10))
        .median()
        .clip(REGION);

      var ratio = s1.select('VV').divide(s1.select('VH')).rename('VV_div_VH');
      return s1.addBands(ratio);
    }
  }
};

// ===================== SAMPLES =====================
// Build once, reuse for every sensor/classifier
var classMap = [
  [Waterbodies, 0],
  [Vegetation, 1],
  [Bareland, 2],
  [Agricultural_land, 3],
  [Built_up, 4]
];

var labeled = classMap.map(function(pair) {
  return pair[0].map(function(f) { return f.set('class', pair[1]); });
});

var samples = labeled.reduce(function(acc, fc) { return acc.merge(fc); });
print('Samples:', samples);

// ===================== CLASSIFIERS =====================
function getClassifier(name) {
  if (name === 'CART') return ee.Classifier.smileCart();
  if (name === 'RF') return ee.Classifier.smileRandomForest(50);
  if (name === 'SVM') return ee.Classifier.libsvm({
    kernelType: 'RBF', gamma: 0.05, cost: 100
  });
  throw 'Unknown classifier: ' + name;
}

// ===================== CORE PIPELINE =====================
function classifyAndEvaluate(sensorName, classifierName) {
  var sensor = SENSORS[sensorName];
  var bands = sensor.bands;
  var scale = sensor.scale;
  var image = sensor.getImage().select(bands);

  // SVM needs normalized reflectance bands; skip normalization for radar (S1)
  var workingImage = (classifierName === 'SVM' && sensorName !== 'S1')
    ? image.divide(10000)
    : image;

  var training = workingImage.reduceRegions({
    collection: samples,
    reducer: ee.Reducer.mean(),
    scale: scale
  });

  if (classifierName === 'SVM') {
    training = training.filter(ee.Filter.notNull([bands[0]]));
  }

  var withRandom = training.randomColumn();
  var trainSet = withRandom.filter(ee.Filter.lt('random', 0.8));
  var testSet = withRandom.filter(ee.Filter.gte('random', 0.8));

  var classifier = getClassifier(classifierName).train({
    features: trainSet,
    classProperty: 'class',
    inputProperties: bands
  });

  var classified = workingImage.classify(classifier);

  var label = sensorName + ' ' + classifierName;
  Map.addLayer(classified, {min: 0, max: 4, palette: PALETTE}, label + ' Classification');

  var validated = testSet.classify(classifier);
  var errorMatrix = validated.errorMatrix('class', 'classification');
  print(label + ' Error Matrix:', errorMatrix);
  print(label + ' Accuracy:', errorMatrix.accuracy());
  print(label + ' Kappa:', errorMatrix.kappa());

  Export.image.toDrive({
    image: classified,
    description: 'LULC_' + classifierName + '_' + sensorName + '_Dhaka_2025',
    folder: EXPORT_FOLDER,
    region: REGION,
    scale: scale,
    maxPixels: 1e10
  });

  return classified;
}

// ===================== RUN =====================
Map.centerObject(REGION, 10);
Map.addLayer(REGION);

var RUNS = {
  L8: ['CART', 'RF', 'SVM'],
  S2: ['CART', 'RF', 'SVM'],
  S1: ['RF']   // SAR typically only run with RF in your original script
};

Object.keys(RUNS).forEach(function(sensorName) {
  RUNS[sensorName].forEach(function(classifierName) {
    classifyAndEvaluate(sensorName, classifierName);
  });
});