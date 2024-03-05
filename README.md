# (North-western) Patagonian fires

Data and code for the article "Spatio-temporal fire patterns in north-western Patagonia".

The ```data_analyses``` folder contains all data and code to reproduce the results shown in the paper.  

The final fire database is ```patagonian_fires.shp```. It has the filtered fires, excluding those < 10 ha, not intersecting the study area or occurred before July 1998. Its fields are:  
- year: fire-year, running from July to June, to be centred in summer.  
- fire_id: unique fire event identifier.  
- area_ha: burned area (ha).  
- date: point estimate of the date of burning.  
- date_l: lower bound for date, set from images.  
- date_u: upper bound for date, set from images.  
- datesrc: source from which the date was determined.  
- name: human-readable name, available for known fires.  
- obs: observations.  
- collec: collection, internal identifier.  
- date_lr: here, "r" stands for "records". Date lower bound based on fire records.  
- date_ur: date upper bound based on records.
- date_info: data used to obtain date.  
  
```patagonian_fires_extended.shp``` includes a few more fires meeting those conditions, with a few occurred before 1999.  
```patagonian_fires_raw.shp``` shows the raw result from the burned area segmentation performed in Google Earth Engine. As these polygons were not checked, many false positives can be found. However, it may be useful to find fires that were discarded but the user is certain about.  
```patagonian_fires_eastern.shp``` includes the patagonian_fires files and around 70 more fires mapped towards the east. It was produced mainly by Juan Paritsis and Guadalupe Franco to be used in a study about fires in pine plantations. These extra fires may show a lower-quality mapping, because the mapping procedure was tuned using only a few steppe fires which occurred in more productive steppes. For example, a few fires show many isolated tiny polygons that are probably noise, not burned area.  

In the kml file, the field 'fire_id' was renamed as 'name', and the original 'name' was renamed 'name2' for them to display the fire_id on the left menu on Google Earth.
  
The automatic phase of fire mapping might be reproduced by adapting the following Google Earth Engine codes. As the process required too much memory to be performed in one step, we divided the mapping in two codes: first, we created the annual images for the whole study area (NDVI- and NBR- summaries), and then we used those images to map fires year by year.  
Code to produce annual images:  
https://code.earthengine.google.com/ae2e04af6cdbd8ac097015ac8217598a?noload=true  
Code to map fires:  
https://code.earthengine.google.com/85c2ee42ffbc1ae5037c2986533762c9?noload=true  
Code to check fires visually (by inspecting the NBR time series):  
https://code.earthengine.google.com/7c013e4f16690c9621acfe1192e0ebf0?noload=true  
Code to extract most of the data used for analyses, and more:  
https://code.earthengine.google.com/b012ccfe54574a8610c60924ea77f3e3?noload=true
  
For details about fire mapping methods, see the supplementary information of the article (Appendix 2).
