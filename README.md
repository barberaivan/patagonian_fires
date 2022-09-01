# patagonian_fires

Landsat-based Fire database for northwestern Patagonia, from 1999 to 2022.  

patagonian_fires.shp has the filtered fires, excluding those smaller than 10 ha, not intersecting the study area or occurred before 1999.  
patagonian_fires_extended.shp includes a few more fires meeting those conditions, with a few occurred before 1999.
patagonian_fires_raw.shp shows the raw result from the burned area segmentation performed in Google Earth Engine. As these polygons were not checked, many false positives can be found. However, it may be usefuel to find fires that were discarded but the user is certain about.
