function Shoreline_Detection(start_date, end_date, geometry){

//INPUT date filter
var start = ee.Date(start_date);
var end = ee.Date(end_date);

// Specify an equals filter for image timestamps for cfmask~image join. 
var filterDateEq = ee.Filter.equals({
    leftField:'system:time_start',
    rightField: 'system:time_start'
});


//Merging Landsat 4, 5, 7, 8
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//Dictionaries and lists for use in image collection filtering.
var possibleSensors = ee.List(['L4','L5','L7','L8']);
var bandNames = ee.List(['blue','green','red','nir','swir1','temp','swir2']);
var bandNumbers = [0,1,2,3,4,5,6];
var sensor_band_dict =ee.Dictionary({L8 : ee.List([1,2,3,4,5,9,6]),
                      L7 : ee.List([0,1,2,3,4,5,7]),
                      L5 : ee.List([0,1,2,3,4,5,6]),
                      L4 : ee.List([0,1,2,3,4,5,6])
});
var collection_dict = {L8: 'LC8',
                         L7: 'LE7',
                         L5: 'LT5',
                         L4: 'LT4'
  };
var spacecraft_dict = {'Landsat4': 'L4','Landsat5': 'L5','Landsat7': 'L7', 'LANDSAT_8': 'L8'};
///////////////////////////////////////////////////////////////////////////
//Function for creating a collection.
function getCollection(sensor,startDate,endDate){
  var collectionName = collection_dict[sensor];
  
  //Start with an un-date-confined collection of iamges
  var WOD = ee.ImageCollection(collectionName)
          .filterBounds(geometry)
          .map(ee.Algorithms.Landsat.TOA);
          
  //Pop off an image to serve as a template if there are no images in the date range
  var dummy = ee.Image(WOD.first());
  
  //Filter by the dates
  var ls = WOD.filterDate(startDate,endDate);
  ls = ls.select(sensor_band_dict.get(sensor),bandNames)
  return ls;
 }

//Create the collections
var l4s = ee.ImageCollection(ee.Algorithms.If(possibleSensors.contains('L4'),getCollection('L4',start, end),getCollection('L4',ee.Date('1000-01-01'),ee.Date('1001-01-01'))));
var l5s = ee.ImageCollection(ee.Algorithms.If(possibleSensors.contains('L5'),getCollection('L5',start, end),getCollection('L5',ee.Date('1000-01-01'),ee.Date('1001-01-01'))));
var l7s = ee.ImageCollection(ee.Algorithms.If(possibleSensors.contains('L7'),getCollection('L7',start, end),getCollection('L7',ee.Date('1000-01-01'),ee.Date('1001-01-01'))));
var l8s = ee.ImageCollection(ee.Algorithms.If(possibleSensors.contains('L8'),getCollection('L8',start, end),getCollection('L8',ee.Date('1000-01-01'),ee.Date('1001-01-01'))));
//Merge the collections
var imagery = ee.ImageCollection(l4s.merge(l5s).merge(l7s).merge(l8s));
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//Determine Projection
var proj =  ee.Image(imagery.first()).select([1]);
proj = proj.projection().crs();
 
//cfmask band from Surface Reflectance collections:
var L8FMASK = ee.ImageCollection('LANDSAT/LC8_SR').select('cfmask').filterBounds(geometry);
var L7FMASK = ee.ImageCollection('LANDSAT/LE7_SR').select('cfmask').filterBounds(geometry);
var L5FMASK = ee.ImageCollection('LANDSAT/LT5_SR').select('cfmask').filterBounds(geometry);
var L4FMASK = ee.ImageCollection('LANDSAT/LT4_SR').select('cfmask').filterBounds(geometry);
var FMASK = ee.ImageCollection(L4FMASK.merge(L5FMASK).merge(L7FMASK).merge(L8FMASK));

// INNER JOIN Imagery & Masks
//Inner join created 
var innerJoin = ee.Join.inner();
// Apply the join.
var innerJoined = innerJoin.apply(imagery, FMASK, filterDateEq);
// Map a function to merge the image and mask in the output FeatureCollection if dates match.   
var joined_ImageandMask = innerJoined.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
});
 
 //Masking Function
function applyMask(image) {
  image = ee.Image(image);
  // Mask out SHADOW, SNOW, and CLOUD classes, leaving WATER and CLEAR classes.
  return image.updateMask(image.select('cfmask').lt(2));
}
var filteredCollection = ee.ImageCollection(joined_ImageandMask.map(applyMask)).sort('CLOUD_COVER');

// (FUNCTION) Clipper Function, clips by Area Of Interest
var clipper = function(image) {
  return image.clip(geometry);
};
filteredCollection=filteredCollection.map(clipper);
Map.addLayer(filteredCollection);
 
//Create a mean of the capture date timeframe
var mean = filteredCollection.mean();
//Reproject back into original coordinate system (bug fix for alteration of projection when mean filtering)
mean = mean.reproject(proj, null, 30)
Map.addLayer(mean);
 
//Rename bands for compatability between Landsat series. Only applied if <L7 collections are used.
var na = [""]
na = mean.bandNames();
var mixed_LS = na.contains('blue');
var divbands = ['B7', 'B3'];
if(mixed_LS){
  divbands = ['swir2', 'green'];
}
//Apply Normalised Difference to isolate land and water.
var ratio = mean.normalizedDifference(divbands).rename('RATIO');

 
//Image Smoothing Attempt (boxcar ~ convolve)
// Define a boxcar or low-pass kernel.
var boxcar = ee.Kernel.circle({
  radius: 4, units: 'pixels', normalize: true
});
 
// Smooth the image by convolving with the boxcar kernel.
//var ratio = ratio.convolve(boxcar);
Map.addLayer(ratio)
 
// Compute the histogram of the band. The mean and variance are only FYI.
var histogram = ratio.reduceRegion({
  reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true),
  geometry: geometry,
  scale: 30,
  bestEffort: true
});
      //print(histogram);
// Chart the histogram
      //print(Chart.image.histogram(ratio, geometry, 30));
 
 
//////////// OTSU
// Return the DN that maximizes interclass variance (in the region).
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
 
  var indices = ee.List.sequence(1, size);
 
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });
 
      //print(ui.Chart.array.values(ee.Array(bss), 0, means));
 
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};
////////////

//Filter imagery by threshold
var threshold = otsu(histogram.get('RATIO_histogram')); //Pulls threshold from Otsu.
      //print('threshold', threshold); //Displays threshold float.
var classA = ratio.gt(threshold); //Selects all pixels less than threshold.

//Applying Mask and Adding to Map
var Areamask = classA.mask(classA);
//Remove Border Issues
var buffer = geometry.buffer(-100);
Areamask = Areamask.clip(buffer)
Map.addLayer(Areamask, {palette: 'blue'}, start_date + '_' + end_date);//Displays selection on the map as a layer.

//My Vectorize
var vectors = Areamask.mask(Areamask).reduceToVectors({
  geometry: geometry, //extent works for the whole map but gives offset results
  crs: ratio.projection(),
  scale: 30,
  maxPixels: 1e13
});
      Map.addLayer(vectors, {color: '50E52B'})

//Export to drive
/*
Export.table.toDrive({
  collection: vectors,
  description: start_date + '_' + end_date,
  fileFormat: 'KML'
   });

   //Export to asset (Optional)
Export.image.toAsset({
  image: Areamask,
  description: start_date + '_' + end_date,
  scale: 30,
  maxPixels: 1e13
   });
*/ 
}



// User Input: (start) YYYY-MM-DD, (end) YYYY-MM-DD, (AOI) geometry;
Shoreline_Detection('2000-01-01', '2001-01-01', geometry)
Map.centerObject(geometry, 13);
