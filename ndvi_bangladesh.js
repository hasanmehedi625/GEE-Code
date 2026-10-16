var Dataset = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');


var image1 = ee.ImageCollection(Dataset
        .filterDate('2024-01-01', '2024-12-31')
        .filterBounds(bd))
        .filterMetadata('CLOUD_COVER_LAND','less_than',3) // cloud cover <3%
        .sort('CLOUD_COVER');

        // .sort('CLOUD_COVER').first();  

var image = image1.mosaic().clip(bd);
// var image = image2.clip(bd);
print(image)

/* Select and rename bands */
image = image.select(['SR_B5','SR_B4', 'SR_B3','SR_B2'],['nir','red', 'green', 'blue']);
print(image)


// NDVI
// Compute the Normalized Difference Vegetation Index (NDVI).
// NDVI = (NIR - Red) / (NIR + Red)
// var ndvi = 'B5'.subtract('B4').divide('B5'.add('B4')).rename('NDVI');

var ndvi = image.expression('(nir - red)/(nir + red)',{
  'nir': image.select('nir'),
      'red': image.select('red'),
} )

// var ndvi = image.normalizedDifference(['nir','red']);
print(ndvi)

var evi = image.expression('2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))',{
  'nir': image.select('nir'),
      'red': image.select('red'),
      'blue': image.select('blue')
} )



// Display results.
Map.addLayer(
  ndvi,
  {
    min: -.4,
    max: .5,
    palette: ['blue', 'white', 'green']
  },
  'NDVI image'
  );
Map.addLayer(evi, {min: 0, max: 2, palette: ['blue', 'white', 'green']}, 'EVI image');



  
  
Map.addLayer(
  image, 
  {
    bands : ['nir','red','blue'],
    min: 8000,
    max: 17000
  },
  'False Color ()');

// Map.addLayer(bd)


// Export
Export.image.toDrive({
  image: ndvi,
  description: 'NDVI_Clodu_Cover_3',
  folder: '507lab',
  crs: "EPSG:4326",
  scale: 1000,
  maxPixels: 1000,
  region: bd,
  fileFormat: 'GeoTIFF'
})
