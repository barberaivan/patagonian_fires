# patagonian_fires

Landsat-based fire database for northwestern Patagonia, from 1999 to 2022.  

patagonian_fires.shp has the filtered fires, excluding those smaller than 10 ha, not intersecting the study area or occurred before 1999.  
patagonian_fires_extended.shp includes a few more fires meeting those conditions, with a few occurred before 1999.  
patagonian_fires_raw.shp shows the raw result from the burned area segmentation performed in Google Earth Engine. As these polygons were not checked, many false positives can be found. However, it may be useful to find fires that were discarded but the user is certain about.  
  
patagonian_fires_eastern.shp includes the patagonian_fires files and around 70 more fires mapped towards the East. It was produced mainly by Juan Paritsis to be used in a pine plantations study. These extra fires may show a lower-quality mapping, because the mapping procedure was tuned using only a few steppe fires which occurred in more productive steppes. For example, a few fires show many isolated tiny polygons that are probably noise, not burned area.  
  
All files are in WGS84 coordinate reference system.  
  
The fire date was obtained mainly from Landsat imagery, but also from MODIS for large fires. This is indicated in 'datesrc'. If MODIS was used, 'date_l' (for lower) is the first hotspot found, and 'date_u' (for upper) is the latest. If Landsat was used, 'date_l' is the latest date before the NBR dropped, and 'date_u' is the earliest date when the NBR was lower than expected in the absence of fire. 'collec' means nothing useful for users.  
  
In a few cases, date_l is in late summer or autumn (fire season ending) and date_u is in the spring of the next fire season (beginning), which makes hard to identify in which season the fire occurred. These cases are labelled as 'uncertain_year' in the 'obs' field.  
  
In the kml file, the field 'fire_id' was renamed as 'name', and the original 'name' was renamed 'name2' for them to display the fire_id on the left menu on Google Earth.
