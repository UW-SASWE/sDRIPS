# Importing required libraries 
import ee 
import logging
import colorlog
from tqdm import tqdm
# import geemap
import os, datetime, requests, zipfile, time, math, glob
import urllib.request
import datetime,math
import numpy as np
import pandas as pd
import subprocess as sub
import rasterio as rio
import geopandas as gpd
from pyproj import Proj, Transformer, CRS
from shapely.geometry import Polygon
from rasterio.mask import mask as riomask
import matplotlib.pyplot as plt
import configparser
import sys
import traceback
import cartopy.crs as ccrs
from contextily import add_basemap  
from matplotlib.colors import Normalize
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrow
from rasterio.coords import BoundingBox
import warnings
warnings.filterwarnings("ignore")

# Initialise Earth Engine
# ee.Authenticate() # If you are running the earth engine first time in your machine, you need to run the ee.Authenticate() command first.
ee.Initialize()


#########*******Module (1) USER INPUT (ROI and Intereset Dates and Planting Date)********########
script_config = configparser.ConfigParser()

script_config.read('..\Config_files\Script_Config.ini')

# Accessing the information from 'Irrigation_Canals_shapefile' section
irrigation_canals_path = script_config.get('Irrigation_Canals_shapefile', 'path')
feature_name = script_config.get('Irrigation_Canals_shapefile', 'feature_name')
bounds_leftlon = float(script_config.get('Irrigation_Canals_shapefile_Bounds', 'leftlon'))
bounds_rightlon = float(script_config.get('Irrigation_Canals_shapefile_Bounds', 'rightlon'))
bounds_toplat = float(script_config.get('Irrigation_Canals_shapefile_Bounds', 'toplat'))
bounds_bottomlat = float(script_config.get('Irrigation_Canals_shapefile_Bounds', 'bottomlat'))
save_data_loc = script_config.get('Save_Data_Location', 'save_data_loc')
gee_asset_id = script_config.get('GEE_Asset_ID', 'id')
# Accessing the information from 'Date_Running' section
start_date = script_config.get('Date_Running', 'start_date')
default_run_week = script_config.getboolean('Date_Running', 'default_run_week')
if default_run_week:
    run_week = ['lastweek', 'currentweek']
else:
    run_week = eval(script_config.get('Date_Running', 'run_week'))

clear_directory_condition = script_config.getboolean('Clean_Directory', 'clear_directory_condition')

ET_estimation = script_config.getboolean('Run_ET_Estimation', 'ET_estimation')

GLCC_mask_condition = script_config.getboolean('GLCC_Mask', 'glcc_mask')
if GLCC_mask_condition:
    glcc_mask = True
else:
    glcc_mask = False

precipitation_condition = script_config.getboolean('Precipitation_Config', 'consider_preciptation')

weather_condition = script_config.getboolean('Weather_Config', 'consider_forecasted_weather')

percolation_condition = script_config.getboolean('Percolation_Config', 'consider_percolation')

region_stats_condition = script_config.getboolean('Region_stats', 'estimate_region_stats')

Tiff2PNGs_condition = script_config.getboolean('Tiff2PNGs', 'tiff2_pngs')
if Tiff2PNGs_condition:
    Tiff2PNGs = True
else:
    Tiff2PNGs = False
secrets = configparser.ConfigParser()
secrets.read('..\Config_files\Secrets.ini')
imerg_username = secrets.get('IMERG_Account', 'username')
imerg_password = secrets.get('IMERG_Account', 'password')

## SETUP LOGGING
os.makedirs(f'{save_data_loc}',exist_ok = True)
os.makedirs(f'{save_data_loc}/logs/',exist_ok = True)
dt_fmt = '%Y%m%d_%H%M%S'
logging.basicConfig(filename=f'{save_data_loc}/logs/{datetime.datetime.today().strftime(dt_fmt)}.log', format='%(levelname)s - %(asctime)s :%(message)s', level=logging.DEBUG)
logging.info("")
logging.info("Location of the Data Folder: "+save_data_loc)
logging.info("")
logging.critical('Start Date:'+str(start_date) )
logging.critical('Run Week:'+ str(run_week) )

# Accessing the information from 'OSGEO_Path' section
osgeo_path = script_config.get('OSGEO_Path', 'path')
canal_config_path = script_config.get('Canal_Config', 'path')
irrigation_canals= ee.FeatureCollection(gee_asset_id) 
canal_list = irrigation_canals.reduceColumns(ee.Reducer.toList(1), [feature_name]).get('list').getInfo()
glcc = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")

# """ Reading the canal config file to parse the planting date, 
# crop type and soil coefficient for each command region """

config = configparser.ConfigParser()

# Read the INI file
# config.read('../Config_files/New_Canal_Config.ini')
config.read(canal_config_path)

# Function to read settings from a section
def read_canal_settings(section):
    return {
        'planting_date': config.get(section, 'planting_date', fallback=config.get('DEFAULT', 'default_planting_date')),
        'crop_type': config.get(section, 'crop_type', fallback=config.get('DEFAULT', 'default_crop_type')),
        'soil_coef': config.getfloat(section, 'soil_coef', fallback=config.getfloat('DEFAULT', 'defualt_soil_coef'))
    }

# Check if use_default is True
if config.getboolean('DEFAULT', 'use_default'):
    canal_settings = {'DEFAULT': read_canal_settings('DEFAULT')}
else:
    canal_settings = {}
    for section in config.sections():
        if section != 'DEFAULT':  # Exclude the DEFAULT section
            canal_settings[section] = read_canal_settings(section)

def read_crop_coefficients(config_path, crop):
    crop_config = configparser.ConfigParser()
    crop_config.read(config_path)
    coefficients = {}
    for key in crop_config[crop]:
        values = crop_config[crop][key].split(', ')
        if values[0] == 'linear':
            # Format: 'linear, start_value, end_value, days'
            coefficients[(int(key.split('-')[0]), int(key.split('-')[1]))] = ('linear', float(values[1]), float(values[2]), int(values[3]))
        else:
            coefficients[(int(key.split('-')[0]), int(key.split('-')[1]))] = float(values[0])
    return coefficients

def get_growth_kc(coefficients, num_days):
    for day_range, kc in coefficients.items():
        if day_range[0] <= num_days <= day_range[1]:
            if isinstance(kc, tuple) and kc[0] == 'linear':
                # Linear interpolation
                start_value, end_value, days = kc[1:]
                delta_per_day = (end_value - start_value) / days
                days_into_range = num_days - day_range[0]
                return start_value + delta_per_day * days_into_range
            return kc
    return None  # Default case if num_days doesn't fit any range


penman_mean_values = []
sebal_mean_values = []
irr_mean_values = []
def sebal_eto_Landsat():
  total_iterations_sebal = len(run_week) * len(canal_list)
  with tqdm(total=total_iterations_sebal, desc="Estimating Penman ET and SEBAL ET", unit=" Command Area") as pbar:
    for wktime in run_week:
      with open(rf"{save_data_loc}/landsat/stats_" + wktime + ".txt", 'w') as txt:
        txt.write("Region,Penman_ET,Sebal_ET,Irrigation\n")
      if wktime == "lastweek":
        startdate = start_date
        dayVal = 14
      if wktime == "currentweek":
        startdate = start_date
        dayVal = 7
      enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 8+2)
      startDate=ee.Date(startdate)
      endDate=ee.Date(enddate)
      logging.critical('Running Week:'+str(wktime))
      logging.critical("Running Week's Start Date:"+str(startdate))
      logging.critical("Running Week's End Date:"+str(enddate))
      for canal in canal_list:
        regionid = canal[0]
        #logging.info(regionid)
        regionn = canal[0]
        regionid = regionid.replace(" ", "-")
        if os.path.exists(rf"{save_data_loc}/landsat/sebal/" + wktime + r"/sebal_eto_" + regionid + ".zip") ==False:
          try:
            table = irrigation_canals.filter(ee.Filter.equals(feature_name, regionn))
            glccmask = glcc.select("discrete_classification").clip(table)
            glcc_1 = glccmask.gt(20)
            glcc_2 = glccmask.lte(40)
            glcc_crop = glccmask.updateMask(glcc_1).updateMask(glcc_2)
            ROI=table
            NDVIhot_low = 0.03               # Lower NDVI treshold for hot pixels
            NDVIhot_high = 0.25              # Higher NDVI treshold for hot pixels
            Cold_Pixel_Constant=0.50
            Hot_Pixel_Constant=2.00
            selscale = 30
            #########*******Module (2) Satellite Image********########/
            # Reference a Landsat scene.
            l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA").filterBounds(ROI) \
              .filterMetadata('CLOUD_COVER', 'less_than', 90)
            l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA").filterBounds(ROI) \
                .filterMetadata('CLOUD_COVER', 'less_than', 90)
            img =  ee.ImageCollection(l8.merge(l9)).filterDate(startDate, endDate) \
                  .filterBounds(ROI) \
                  .filterMetadata('CLOUD_COVER', 'less_than', 90)#.sort('system:time_start').filterDate(date_range).sort('CLOUD_COVER').first()
            
            collectionSize = img.size()

            # Function to get the latest two images if the collection size is greater than 2
            def get_latest_two_or_all(imageCollection):
                return ee.ImageCollection(ee.Algorithms.If(
                    collectionSize.gte(2),
                    imageCollection.sort('system:time_start', False).limit(2),
                    imageCollection
                ))

            # Apply the function to get the desired image collection
            filteredCollection = get_latest_two_or_all(img)

            # Create a composite image using the median of the filtered collection
            composite = filteredCollection.median()

            # composite = img.median() ### So that Landsat Image completely covers the ROI
            # composite_date = datetime.datetime.fromtimestamp(composite.get('system:time_start').getInfo() / 1000).date()
            sample = ee.Image(img.first())
            img_date = ee.Date(sample.get('system:time_start'))
            my_datetime = datetime.datetime.fromtimestamp(sample.get('system:time_start').getInfo() / 1000).date()  # Apply fromtimestamp function
            if list(canal_settings.keys())[0] != 'DEFAULT':
              planting_YY_MM_DD = canal_settings[canal[0]]['planting_date']
              crop_type = canal_settings[canal[0]]['crop_type']
              soil_coef = canal_settings[canal[0]]['soil_coef']
              logging.info("Planting Date For:"+str(canal_settings[canal[0]])+' is :'+ planting_YY_MM_DD)
              logging.info("Crop Type For:"+str(canal_settings[canal[0]])+' is :'+ crop_type)
            else:
              planting_YY_MM_DD = canal_settings['DEFAULT']['planting_date']
              crop_type = canal_settings['DEFAULT']['crop_type']
              soil_coef = canal_settings['DEFAULT']['soil_coef']
            
            planting_date = datetime.datetime.strptime(planting_YY_MM_DD, '%Y-%m-%d').date()
            num_days = abs((my_datetime - planting_date).days)
            # print(f'Planting Date for command area {canal[0]}:',planting_date)
            # print('Number of days = ',abs((my_datetime - planting_date).days))
            # print((f'Crop grown in command area {canal[0]}:',crop_type))
            logging.info("Planting Date for the command area {area}:{date}".format(area = canal[0],date =planting_date ) )
            logging.info("Number of days = %s",str(abs((my_datetime - planting_date).days)))
            logging.info('Crop grown in command area {canalname}: {crop_type}'.format(canalname = canal[0],crop_type=crop_type))

            crop_config_path = '../Config_files/Crop_Config.ini'
            coefficients = read_crop_coefficients(crop_config_path, crop_type)
            growth_kc = get_growth_kc(coefficients, num_days)
            # print(f"Growth Kc for {crop_type} at {num_days} days: {growth_kc}")  
            logging.info("Growth Kc for {crop_type} at {num_days} days: {growth_kc}".format(crop_type=crop_type,num_days=num_days,growth_kc=growth_kc))
      
            #########*******Module (1.b crop coeffcient)********########/  
            
            #########*******END of Module (2)********########/
            
            #########*******Module (3) Forcing Data From GFS********########/
            gfs_forcings = ee.ImageCollection("NOAA/GFS0P25") \
              .filterDate(ee.Date(my_datetime.strftime("%Y-%m-%d")))\
              .filterBounds(ROI)
            gfs_forcings_composite=gfs_forcings.first()
            temp = gfs_forcings_composite.select('temperature_2m_above_ground').add(273.15)
            max_temp_dict = temp.reduceRegion(
                reducer=ee.Reducer.max(),
                geometry=ROI,
                scale=selscale  # You may need to adjust the scale according to your data
            )
            min_temp_dict = temp.reduceRegion(
                reducer=ee.Reducer.min(),
                geometry=ROI,
                scale=selscale  # You may need to adjust the scale according to your data
            )
            maxtemp= ee.Image.constant(max_temp_dict.get('temperature_2m_above_ground').getInfo())
            mintemp= ee.Image.constant(min_temp_dict.get('temperature_2m_above_ground').getInfo())
            wind_u = gfs_forcings_composite.select('u_component_of_wind_10m_above_ground') 
            wind_v = gfs_forcings_composite.select('v_component_of_wind_10m_above_ground')     
            # winwndd=gfs_forcings_composite.select('Wind_f_inst')   # at 10 m height
            wind = ee.Image().expression(
                'sqrt(a**2+b**2)', {'a': wind_u, 'b': wind_v}
            ).rename('wind')
            pressure_Pa = ee.Image.constant(101325).rename('pressure')
            gfs_forcings_composite.addBands(pressure_Pa)
            QH=gfs_forcings_composite.select('specific_humidity_2m_above_ground')
            # PS=gfs_forcings_composite.select('pressure_Pa')

            #### To calculate the surface pressure, we use hypsometric equation
            seaLevelPressure = 101325; ##// Standard sea level pressure in Pascals
            Rd = 287.05; ##// Gas constant for dry air in J/(kg*K)
            g = 9.80665;## // Gravitational acceleration in m/s^2
            srtm= ee.Image("USGS/SRTMGL1_003")
            exponent = srtm.multiply(g).divide(temp.multiply(Rd)).multiply(-1)

            # Calculate surface pressure using the hypsometric equation
            PS = ee.Image(seaLevelPressure).multiply(exponent.exp())
            
            ## equations for saturated vapor pressure has temperature in celsius
            tmp1=ee.Image.constant(17.67).multiply(temp.subtract(ee.Image.constant(273.15)))
            tmp2=temp.add(ee.Image.constant(243.5-273.15))
            tmp3=(tmp1.divide(tmp2)).exp()
            es=ee.Image.constant(6.112).multiply(tmp3)    # Unit is millibar (mb)
            ## finished estimating saturated vapour pressure

            ## Started estimating Relative Humidity
            tmp4=QH.multiply(PS)
            # tmp4=QH
            tmp5=(ee.Image.constant(0.378).multiply(QH)).add(ee.Image.constant(0.622))
            e=tmp4.divide(tmp5)
            # Unit is Pascal; divide by 100 to convert to millibar (mb)
            e=e.divide(ee.Image.constant(100))
            RH=e.divide(es).multiply(ee.Image.constant(100))
            RH=RH.rename('RH')
            ## Finished estimating Relative Humidity

            # Reproject Data to 30m pixel spacing (a bilinear resampled image)
            temp_30m = temp.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            maxtemp_30m = maxtemp.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            mintemp_30m = mintemp.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            
            wind_30m = wind.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
          
            e_30m = e.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            es_30m = es.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            srtm= ee.Image("USGS/SRTMGL1_003")
                      #.filterBounds(Delta)
            #Map.addLayer(srtm.clip(ROI),{},'srtm')
            # slope=ee.Terrain.slope(srtm)
            # aspect=ee.Terrain.aspect(srtm)

            DEM_30m = srtm.clip(ROI).reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            
            #Map.addLayer(es_30m.clip(ROI),{},'dem30')
            
            # proj = DEM_30m.projection()
            latlon = ee.Image.pixelLonLat().reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            LAT=latlon.select('latitude').reproject(
              crs= sample.select('B1').projection().crs(),
              scale= selscale
            )
            # LON=latlon.select('longitude').reproject(
            #   crs= sample.select('B1').projection().crs(),
            #   scale= selscale
            # )
            
            #########*******END of Module (3) Forcing Data********########/
            
            #########*******Module (4) Albedo********########/
            Surf_albedo=composite.select('B1').multiply(0.130).add(composite.select('B2').multiply(0.115)).add(composite.select('B3').multiply(0.143)).add(composite.select('B4').multiply(0.180)).add(composite.select('B5').multiply(0.281)).add(composite.select('B6').multiply(0.108)).add(composite.select('B7').multiply(0.042)).subtract(0.03).divide(ee.Number(0.89).pow(ee.Number(0.89)))
            Surf_albedo=Surf_albedo.rename('Albedo')
            maskLow=Surf_albedo.gt(0.00)
            maskHigh=Surf_albedo.lt(0.60)
            Surf_albedo_mask = Surf_albedo.updateMask(maskLow)
            Surf_albedo_mask=Surf_albedo_mask.unmask(0.00)
            Surf_albedo_mask = Surf_albedo_mask.updateMask(maskHigh)
            Surf_albedo_mask=Surf_albedo_mask.unmask(0.60)
            #########*******END of Module (4) Albedo********########/
            
            #########*******Module (5) Surface Temperature********########/
            NDVI=composite.normalizedDifference(['B5', 'B4'])
            tmp1=ee.Image.constant(0.8).subtract(NDVI.select('nd'))
            tmp2=ee.Image.constant(0.8).subtract(ee.Image.constant(0.125))
            vegt_cover = ee.Image.constant(1.0).subtract((tmp1.divide(tmp2)).pow(ee.Image.constant(0.70)))
            vegt_cover=vegt_cover.rename('vegt_cover')
            LAI_1 = (((vegt_cover.subtract(ee.Image.constant(1)))).multiply(ee.Image.constant(-1))).log().divide(ee.Image.constant(-0.45))
            LAI_1=LAI_1.rename('LAI')
            LAIHigh=LAI_1.lt(8.0)
            LAI_1=LAI_1.updateMask(LAIHigh)
            LAI_1=LAI_1.unmask(8.0)
            tmp1=(NDVI.pow(ee.Image.constant(3))).multiply(ee.Image.constant(9.519))
            tmp2=(NDVI.pow(ee.Image.constant(2))).multiply(ee.Image.constant(0.104))
            tmp3=NDVI.multiply(ee.Image.constant(1.236))
            tmp4=ee.Image.constant(0.257)
            LAI_2 =tmp1.add(tmp2).add(tmp3).add(tmp4)
            tmp=LAI_1.add(LAI_2)
            LAI=tmp.divide(ee.Image.constant(2))
            # tir_emis = ee.Image.constant(1.009).add(ee.Image.constant(0.047).multiply(NDVI.log()))
            b10_emissivity = LAI
            emissLow=LAI.lte(3.0)
            emissValue1=ee.Image.constant(0.95).add(ee.Image.constant(0.01).multiply(LAI))
            b10_emissivity=b10_emissivity.updateMask(emissLow).unmask(emissValue1)
            emissHigh=LAI.gt(3.0)
            emissValue2=ee.Image.constant(0.98)
            b10_emissivity=b10_emissivity.updateMask(emissHigh).unmask(emissValue2)
            Temp_TOA_10=composite.select('B10')
            Temp_TOA_11=composite.select('B11')
            tt1=ee.Image.constant(1.378).multiply((Temp_TOA_10).subtract(Temp_TOA_11))
            tt2=ee.Image.constant(0.183).multiply((Temp_TOA_10.subtract(Temp_TOA_11)).pow(ee.Image.constant(2)))
            tt3=(ee.Image.constant(54.30).subtract(ee.Image.constant(2.238).multiply(e_30m))).multiply(ee.Image.constant(1).subtract(b10_emissivity))
            Surface_temp = Temp_TOA_10.add(tt1).add(tt2).add(tt3).subtract(ee.Image.constant(0.268))
            
            #########*******END of Module (5) Surface Temperature********########/
            
            #########*******Module (6) Daily Radiation (mm/day)********########/
            # Daily 24 hr radiation - For flat terrain only !
            # 1) Shortwave Radiation
            featureROI=ee.Algorithms.Feature(ROI)
            centroid = featureROI.centroid(maxError= 3)
            ##logging.info(centroid)
            cenLAT=centroid.geometry().coordinates().get(1).getInfo()
            # cenLON=centroid.geometry().coordinates().get(0).getInfo()
            deg2rad=ee.Image.constant(3.14).divide(ee.Image.constant(180))
            phi=LAT.multiply(deg2rad)
            sortImg = img.sort('system:time_start', False)
            listOfImages = sortImg.toList(sortImg.size())
            recentImg=ee.Image(listOfImages.get(0))
            DOY=recentImg.date().getRelative('day', 'year')
            
            ## finding solar declination
            tmp = ee.Number(2 * math.pi * DOY.getInfo() / 365 - 1.39)
            delta=ee.Number(0.409).multiply(ee.Number(tmp).sin())
            ## finished finding solar declination
            
            ## finding sunset hour angle 
            tmp = ee.Number(DOY.getInfo() * 2 * math.pi / 365).cos()
            tmp1=delta.tan()
            ##logging.info(tmp1.getInfo())
            tmp2=ee.Number(-1).multiply(ee.Number(cenLAT*math.pi/180).tan())
            ##logging.info(cenLAT)
            ws=(tmp1.multiply(tmp2)).acos()
            ## finished finding sunset hour angle
            
            ## finding earth to sun distance 
            dr=ee.Number(1).add(ee.Number(0.33).multiply(tmp))
            ## finished finding earth to sun distance 

            ## Finding Ra,Rns
            Gsc=1367
            tmp1=phi.cos().multiply(ee.Image.constant(delta).cos()).multiply(ee.Image.constant(ws).sin())
            tmp2=phi.sin().multiply(ee.Image.constant(delta).sin()).multiply(ee.Image.constant(ws))
            Ra=ee.Image.constant(Gsc).multiply(ee.Image.constant(dr)).divide(3.14).multiply(tmp1.add(tmp2))
            tmp=ee.Number(2).multiply(ee.Number(10).pow(-5))
            Rs=(ee.Image.constant(0.75).add(ee.Image.constant(tmp).multiply(DEM_30m))).multiply(Ra)
            Rns=(ee.Image.constant(1).subtract(Surf_albedo_mask)).multiply(Rs)
            ## Finished Finding Ra,Rns
            
            # 2) Longwave Radiation
            ea=e_30m.multiply(0.10) # convert vapour pressure to kpa from milli bars
            sigma=ee.Number(10).pow(-9).multiply(ee.Number(4.89802))
            tmp1=(mintemp_30m.pow(4).add(maxtemp_30m.pow(4))).divide(ee.Image.constant(2))
            tmp2=ee.Image.constant(0.34).subtract((ea.sqrt()).multiply(0.14))
            tmp3 = ee.Image.constant(1.35).multiply(Ra).divide(Rs)
            tmp4 = tmp3.subtract(ee.Image.constant(1.35))
            Rnl=ee.Image.constant(sigma).multiply(tmp1).multiply(tmp2).multiply(tmp3).multiply(tmp4)
            # 3) Net Radiation
            Rn=Rns.subtract(Rnl)
            #########*******END of Module (6) Daily Radiation (mm/day)********########/
            ##logging.info(Rn)
            #Map.addLayer(Rns,{},'Rns')
            #########*******Module (7) Soil Heat Flux (mm/day)********########/
            # Soil Heat Flux Radiation (G)
            tmp1=ee.Image.constant(0.0038).add(ee.Image.constant(0.0074).multiply(Surf_albedo_mask))
            tmp2=ee.Image.constant(1).subtract(ee.Image.constant(0.978).multiply((NDVI).pow(4)))
            tmp3=Surface_temp.subtract(ee.Image.constant(273.15))
            G= tmp1.multiply(tmp2).multiply(tmp3).multiply(Rn)
            #########*******END of Module (7) Soil Heat Flux (mm/day)********########/
            #########*******Module (8) Selection of Cold/Hot Pixels********########/
            maxNDVI= NDVI.reduceRegion(
              reducer= ee.Reducer.max(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
            )
            stdNDVI= NDVI.reduceRegion(
              reducer= ee.Reducer.stdDev(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11)
            maxNDVINumber = ee.Number(maxNDVI.get('nd'))
            stdNDVINumber = ee.Number(stdNDVI.get('nd'))
            zeroCold = maxNDVINumber.subtract(stdNDVINumber.multiply(Cold_Pixel_Constant))
            maskNDVI=NDVI.select('nd').gt(ee.Number(zeroCold))
            cold_pixels_vegetation =Surface_temp.select('B10').rename('cold')
            cold_pixels_vegetation=cold_pixels_vegetation.updateMask(maskNDVI)
            # cold_pixels_vegetation=cold_pixels_vegetation.unmask(0)
            tempCold= cold_pixels_vegetation.reduceRegion(
              reducer= ee.Reducer.min(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11)
            hot_pixels =Surface_temp.select('B10').rename('hot')
            hot_pixels_check= hot_pixels.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
            )
            hotLow=hot_pixels.gt(ee.Number(NDVIhot_low))
            hotHigh=hot_pixels.lt(ee.Number(NDVIhot_high))
            hot_pixels_mask = hot_pixels.updateMask(hotLow)
            hot_pixels_mask=hot_pixels_mask.unmask(ee.Number(NDVIhot_low))
            hot_pixels_mask =hot_pixels_mask.updateMask(hotHigh)
            hot_pixels_mask=hot_pixels_mask.unmask(ee.Number(NDVIhot_high))
            avgtempHot= hot_pixels.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
            )
            stdtempHot= hot_pixels.reduceRegion(
              reducer= ee.Reducer.stdDev(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            tempHot=ee.Number(avgtempHot.get('hot').getInfo()).add(ee.Number(Hot_Pixel_Constant).multiply(stdtempHot.get('hot').getInfo()))

            #########*******END of Module (8) Selection of Cold/Hot Pixels********########/
            ##logging.info('tempHot',tempHot)
            #########*******Module (9) Sensible Heat Flux********########/
            
            # calculate the windspeed and friction by using the Raupach or NDVI model
            # constants
            k_vk = ee.Number(0.41)      # Von Karman constant
            h_grass = ee.Number(0.12)   # Grass height (m)
            # cd = ee.Number(53)          # Free parameter for displacement height, default = 20.6
            zx=ee.Number(10)    # wind speed height in meters
            # Surface roughness using NDVI Model (other option is Raupach Model)
            # a, b need to be determined by fitting relation between ln(zom) vs NDVI/ α. Using zom=0.12h (h=vegetation height)
            # LAI method: zom=0.018 × LAI
            tmp=ee.Image.constant(1.096).multiply(NDVI.select('nd')).divide(Surf_albedo_mask)
            zom_NDVI =(tmp.subtract(ee.Image.constant(5.307))).exp()
            #zom_NDVI[water_mask == 1.0] = 0.001
            #Map.addLayer(zom_NDVI)
            maxzomNDVI= zom_NDVI.lt(10.0)
            zom_NDVI=zom_NDVI.updateMask(maxzomNDVI)
            zom_NDVI=zom_NDVI.unmask(10.0)
            Surf_roughness = zom_NDVI
            zom_grass = ee.Number(0.123).multiply(h_grass)
            # Friction velocity for grass (m/s):
            tmp1=ee.Image.constant(k_vk).multiply(wind_30m)
            tmp2=(zx.divide(zom_grass)).log()
            ustar_grass = tmp1.divide(ee.Image.constant(tmp2))
            # Wind speed (m/s) at the "blending height" (200m):
            tmp=(ee.Number(200).divide(zom_grass)).log()
            u_200 = ustar_grass.multiply(tmp).divide(ee.Image.constant(k_vk))
            tmp=(ee.Image.constant(200).divide(Surf_roughness)).log()
            ustar = ee.Image.constant(k_vk).multiply(u_200).divide(tmp)
            # areodynamic rah (at stable conditions)
            tmp1=(ee.Image.constant(2).divide(ee.Image.constant(0.01))).log()
            tmp2=ee.Image.constant(k_vk).multiply(ustar)
            rah= tmp1.divide(tmp2)
            # Generally, air temperature decreases by about 0.65 celsius when elevation increases by 100 m under neutral stability conditions.
            Temp_lapse_rate = 0.0065   # or 0.01199   # Temperature lapse rate (°K/m)
            tmp=(ee.Image.constant(293).subtract(ee.Image.constant(Temp_lapse_rate).multiply(DEM_30m))).divide(ee.Image.constant(293))
            Pair = ee.Image.constant(101.3).multiply(tmp.pow(5.26))   #units':KPa
            #Map.addLayer(Pair)
            # Air denisty using Ideal gas law
            air_dens = Pair.multiply(ee.Image.constant(1000)).divide(ee.Image.constant(1.01).multiply(Surface_temp).multiply(ee.Image.constant(287)))
            
            srtmclip=srtm.clip(ROI)
            #Map.addLayer(srtmclip,{},'srtmclip')
            DEM_300m = srtmclip.resample('bilinear')
            
            # Giving the first guess for the stability (LE = 0 Therefore Rn-G = H)
            dT_init = (Rn.subtract(G)).multiply(rah).divide(air_dens.multiply(ee.Image.constant(1004)))
            dT_init=dT_init.select('constant').rename('dt')
            medDT= dT_init.reduceRegion(
              reducer= ee.Reducer.median(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            tmp=tempHot.getInfo()-(tempCold.get('cold').getInfo())
            ##logging.info(tmp)
            slope_dt = dT_init.divide(ee.Image.constant(tmp))
            offset_dt = dT_init.subtract(slope_dt.multiply(ee.Image.constant(tempHot.getInfo())))
            dT=offset_dt.add(slope_dt.multiply(Surface_temp))
            H_initial = air_dens.multiply(ee.Image.constant(1004)).multiply(dT).divide(rah)
            initial_H_max = H_initial.reduceRegion(
                reducer=ee.Reducer.max(),
                geometry=ROI,
                scale=selscale,
                maxPixels=1e11
            ).get('constant')
            H = H_initial
            H_prev = H_initial
            
            # Sensible Heat

            #### New approach for iteration
            H = air_dens.multiply(ee.Image.constant(1004)).multiply(dT).divide(rah)
            
            
            # Iterative process is required here for correcting ustar & rah
          
            for iter in range(3):
              tmp1=ee.Image.constant(-1004).multiply(air_dens).multiply(ustar.pow(3)).multiply(Surface_temp)
              tmp2=ee.Image.constant(k_vk).multiply(ee.Image.constant(9.81)).multiply(H)
              L_MO = tmp1.divide(tmp2)
              # L_MO < 0 Unstable, L_MO > 0 Stable, L_MO = 0 Neutral.
              # Stability Condition
              psi_m200_stable = ee.Image.constant(-10).divide(L_MO)
              psi_h2_stable = ee.Image.constant(-10).divide(L_MO)
              psi_h001_stable = ee.Image.constant(-0.05).divide(L_MO)
              # Neutral Condition
              # psi_m200_neutral = ee.Image.constant(0)
              # psi_h2_neutral = ee.Image.constant(0)
              # psi_h001_neutral = ee.Image.constant(0)
              # UnStability Condition
              tmp=ee.Image.constant(1).subtract(ee.Image.constant(16).multiply(ee.Image.constant(2)).divide(L_MO))
              x2 = tmp.pow(0.25)  # x at 2m
              tmp=ee.Image.constant(1).subtract(ee.Image.constant(16).multiply(ee.Image.constant(200)).divide(L_MO))
              x200 = tmp.pow(0.25)  # x at 200m
              tmp=ee.Image.constant(1).subtract(ee.Image.constant(16).multiply(ee.Image.constant(0.01)).divide(L_MO))
              x001 = tmp.pow(0.25)  # x at 0.01m
              tmp=(ee.Image.constant(1).add(x2.pow(2))).divide(ee.Image.constant(2)).log()
              psi_h2_unstable= ee.Image.constant(2).multiply(tmp)
              tmp=(ee.Image.constant(1).add(x001.pow(2))).divide(ee.Image.constant(2)).log()
              psi_h001_unstable= ee.Image.constant(2).multiply(tmp)
              
              tmp1=(ee.Image.constant(1).add(x200)).divide(ee.Image.constant(2)).log()
              tmp2=(ee.Image.constant(1).add(x200.pow(2))).divide(ee.Image.constant(2)).log()
              tmp3=ee.Image.constant(2).multiply(x200.atan())
              psi_m200_unstable= (ee.Image.constant(2).multiply(tmp1)).add(tmp2).subtract(tmp3).add(ee.Image.constant(0.5).multiply(math.pi))
              tmp1=ee.Image.constant(k_vk).multiply(u_200)
              tmp2=((ee.Image.constant(200).divide(Surf_roughness)).log()).subtract(psi_m200_unstable)
              ustar_corr_unstable = tmp1.divide(tmp2)
              tmp2=((ee.Image.constant(200).divide(Surf_roughness)).log()).subtract(psi_m200_stable)
              ustar_corr_stable = tmp1.divide(tmp2)
              ustar_corr=ustar_corr_unstable
              L_unstable=L_MO.lt(1.0)
              ustar_corr_mask=ustar_corr.updateMask(L_unstable)  #masking stable pixels
              ustar_corr=ustar_corr_mask.unmask(ustar_corr_stable)
              tmp1=((ee.Image.constant(2).divide(ee.Image.constant(0.01))).log()).subtract(psi_h2_stable).add(psi_h001_stable)
              tmp2=ee.Image.constant(k_vk).multiply(ustar_corr)
              rah_corr_stable=tmp1.divide(tmp2)
              tmp1=((ee.Image.constant(2).divide(ee.Image.constant(0.01))).log()).subtract(psi_h2_unstable).add(psi_h001_unstable)
              rah_corr_unstable=tmp1.divide(tmp2)
              rah_corr=rah_corr_unstable
              L_unstable=L_MO.lt(1.0)
              rah_corr_mask=rah_corr.updateMask(L_unstable)  #masking stable pixels
              rah_corr=rah_corr_mask.unmask(rah_corr_stable)
              dT_corr = (Rn.subtract(G)).multiply(rah_corr).divide(air_dens.multiply(ee.Image.constant(1004)))
              dT_corr=dT_corr.select('constant').rename('dt')
              medDT= dT_corr.reduceRegion(
                reducer= ee.Reducer.median(),
                geometry= ROI,
                scale= selscale,
                maxPixels= 1e11
                )
              slope_dt= ee.Number(medDT.get('dt').getInfo()).divide(tempHot.subtract(tempCold.get('cold')))
              offset_dt = ee.Number(medDT.get('dt').getInfo()).subtract(slope_dt.multiply(tempHot))
              dT=ee.Image.constant(offset_dt).add(ee.Image.constant(slope_dt).multiply(Surface_temp))
              # Sensible Heat
              H = air_dens.multiply(ee.Image.constant(1004)).multiply(dT).divide(rah_corr)
             
            
            H_mask=H.updateMask(H.gte(0))
            H_mask=H_mask.unmask(0)
            
            
            #########*******END of Module (9) Sensible Heat Flux********########/
            #########*******Module (10) Reference ET********########/
            # Reference Evapotranspiration (Penman-Monteith)
            # Effective leaf area index involved, see Allen et al. (2006):
            
            tmp=(ee.Image.constant(0.30).multiply(LAI)).add(ee.Image.constant(1.2))
            LAI_eff = LAI.divide(tmp)
            rl = 130 # Bulk stomatal resistance of the well-illuminated leaf (s/m) [See box 5 in FAO56]
            rs_min = ee.Image.constant(rl).divide(LAI_eff)  # Min (Bulk) surface resistance (s/m)
            # Latent heat of vaporization (J/kg) (calculated above)
            # Reference evapotranspiration- grass
            # Penman-Monteith of the combination equation (eq 3 FAO 56) (J/s/m2)
            # For reference ETo, the albedo is 0.23
            Rns_ref=(ee.Image.constant(1).subtract(ee.Image.constant(0.23))).multiply(Rs)
            Rnl_ref=Rnl
            Rn_ref=Rns_ref.subtract(Rnl_ref)
            
            # convert units of vapour pressure from millibar to KPa
            es_30m=es_30m.multiply(ee.Image.constant(0.10))
            e_30m=e_30m.multiply(ee.Image.constant(0.10))
            tmp=(temp_30m.subtract(ee.Image.constant(273.15)).add(ee.Image.constant(237.3))).pow(2)
            es_slope =ee.Image.constant(4098).multiply(es_30m).divide(tmp)  #unit is KPa/°C
            rah_grass = ee.Image.constant(208.0).divide(wind_30m)
            # Psychrometric constant (kPa / °C), FAO 56, eq 8.:
            #Temp_lapse_rate = 0.0065   # or 0.01199   # Temperature lapse rate (°K/m)
            #tmp=(ee.Image.constant(293).subtract(ee.Image.constant(Temp_lapse_rate).multiply(DEM_30m))).divide(ee.Image.constant(293))
            #Pair = ee.Image.constant(101.3).multiply(tmp.pow(5.26))   #units:KPa
            Psychro_c = ee.Image.constant(0.665).multiply(ee.Image.constant(10).pow(-3)).multiply(Pair)
            tmp1=ee.Image.constant(1).add(ee.Image.constant(70).divide(rah_grass))
            tmp2=(tmp1.multiply(Psychro_c)).add(es_slope)
            tmp3=air_dens.multiply(ee.Image.constant(1004)).multiply(es_30m.subtract(e_30m)).divide(rah_grass)
            LET_ref_24 =(es_slope.multiply(Rn_ref).add(tmp3)).divide(tmp2)
            tmp=ee.Image.constant(0.002361).multiply(Surface_temp.subtract(ee.Image.constant(273.15)))
            # calculate lamda or latent heat vaporization for latent heat
            Lhv = (ee.Image.constant(2.501).subtract(tmp)).multiply(ee.Image.constant(10).pow(6))
            # Reference evaportranspiration (mm/day):
            ETref_24 = LET_ref_24.divide(Lhv.multiply(ee.Image.constant(1000))).multiply(ee.Image.constant(86400000))
            # Potential Evapotranspiration mm/day)
            # Penman-Monteith of the combination equation (eq 3 FAO 56) (J/s/m2)
            #rah_pm_act=((np.log((2.0-0.0)/(Surf_roughness*0.1))*np.log((2.0-0.0)/(Surf_roughness)))/(k_vk*1.5**2))*((1-5*(-9.82*dT*(2.0-0.0))/((273.15+Temp_inst)*1.5**2))**(-0.75))
            #rah_pm_act[rah_pm_act<25]=25
            #LETpot_24 = ((sl_es_24 * (Rn_24 - Refl_rad_water) + air_dens * 1004 *(esat_24 - eact_24)/rah_pm_pot) / (sl_es_24 + Psychro_c * (1 + rs_min/rah_pm_pot)))
            #ETpot_24 = LETpot_24 / (Lhv * 1000) * 86400000
            #ETpot_24[ETpot_24 > 15.0] = 15.0
            ETref_24=ETref_24.select('constant').rename('etr')
            if glcc_mask == GLCC_mask_condition:
                ETref_24 = ETref_24.updateMask(glcc_crop)  
            else:
                ETref_24 = ETref_24.clip(ROI) 
            # #logging.info(ETref_24.getInfo())
            penmanET = ee.Image.constant(soil_coef).multiply(ETref_24.multiply(ee.Image.constant(growth_kc))).multiply(ee.Image.constant(dayVal))
            # #logging.info(penmanET_24.getInfo())
            medET_ref= penmanET.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            # #logging.info(medET_ref.getInfo())
            median_etref= medET_ref.get('constant').getInfo()
            
            #########*******END of Module (10) Reference ET********########/
            proj = ee.Projection('EPSG:4326')
            penmanET = penmanET.clip(ROI).reproject(
              crs= proj,
              scale= selscale
            )
            params={}
            params = {'name': "penman_eto_" + regionid, 'filePerBand': "false", "scale": penmanET.projection().nominalScale(), 'region': ROI.geometry()}
            url = penmanET.getDownloadURL(params)
            # logging.info("Download URL:" + url)
            # params_et = {'name': "penman_et_" + regionid, 'filePerBand': "false", "scale": ETref_24.projection().nominalScale(), 'region': ROI.geometry()}
            # url_penman_et = ETref_24.getDownloadURL(params_et)
            urllib.request.urlretrieve(url, rf"{save_data_loc}/landsat/penman/" + wktime + r"/penman_eto_" + regionid + ".zip")
            # urllib.request.urlretrieve(url_penman_et, r"Data/landsat/penman/" + wktime + r"/penman_et_" + regionid + ".zip")

            #########*******Module (11) Actual ET********########/
            # Evaporative Fraction (crop coefficient)
            
            
            LE= Rn.subtract(G).subtract(H)   # Latent Heat
            
            EF= LE.divide(Rn.subtract(G))    # Evaporative fraction
            tmp=ee.Image.constant(0.002361).multiply(Surface_temp.subtract(ee.Image.constant(273.15)))
            # calculate lamda or latent heat vaporization for latent heat
            Lhv = (ee.Image.constant(2.501).subtract(tmp)).multiply(ee.Image.constant(10).pow(6))
            ETA_24 = EF.multiply(Rn).divide(Lhv.multiply(ee.Image.constant(1000))).multiply(ee.Image.constant(86400000))
            ETA_24=ETA_24.select('constant').rename('eta')
            ETA_24_mask=ETA_24.updateMask(ETA_24.lt(30))
            ETA_24_mask=ETA_24_mask.updateMask(ETA_24.gte(0))
            # ETA_24_mask=ETA_24_mask.unmask(0)
            ETA_24_mask=ETA_24_mask.updateMask(NDVI.gte(0.2)).unmask(0)
            if glcc_mask == GLCC_mask_condition:
                ETA_24_mask = ETA_24_mask.updateMask(glcc_crop)
            else:
                ETA_24_mask = ETA_24_mask.clip(ROI)
            ETA_mask = ETA_24_mask.multiply(ee.Image.constant(dayVal))
            #Map.addLayer(NDVI,{},'NDVI')
            ETA_mask=ETA_mask.clip(ROI).reproject(
              crs= proj,
              scale= selscale
            )
            ETA_mask = ETA_mask.min(penmanET)

            totET= ETA_mask.reduceRegion(
              reducer= ee.Reducer.sum(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            ##logging.info(ETA_24)
            total_etc = totET.get('eta').getInfo()
            params={}
            
            params = {'name': "sebal_eto_" + regionid, 'filePerBand': "false", "scale": ETA_mask.projection().nominalScale(), 'region': ROI.geometry()}
            url = ETA_mask.getDownloadURL(params)
            # #logging.info("Download URL:" + url)
            
            urllib.request.urlretrieve(url, rf"{save_data_loc}/landsat/sebal/" + wktime + r"/sebal_eto_" + regionid + ".zip")
            
            avgET= ETA_mask.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            ##logging.info(ETA_24)
            avg_etc = avgET.get('eta').getInfo()
            overirri = ETA_mask.subtract(penmanET)
            if glcc_mask == GLCC_mask_condition:
                overirri = overirri.updateMask(glcc_crop)
            else:
                overirri = overirri.clip(ROI)
            irristatus= overirri.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            avg_irri = irristatus.get('eta').getInfo()
            params={}
            
            params = {'name': "irrigation_" + regionid, 'filePerBand': "false", "scale": overirri.projection().nominalScale(), 'region': ROI.geometry()}
            url = overirri.getDownloadURL(params)
            urllib.request.urlretrieve(url, rf"{save_data_loc}/landsat/irrigation/" + wktime + r"/irrigation_" + regionid + ".zip")
            logging.info('ET Values For The Processed Command Area - Format: Command Area, Average Penman ET, Average SEBAL ET, Average Overirrigation (SEBAL - Penman)')
            logging.info((canal[0] + ',' + str(median_etref) + "," + str(avg_etc)+ ',' + str(avg_irri)))
            penman_mean_values.append((regionid, median_etref))
            sebal_mean_values.append((regionid, avg_etc))
            irr_mean_values.append((regionid, avg_irri))
            
          except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(f"Error found: {e}")
            logging.error("Week: " + str(wktime) + ", Region: " + str(regionid) + " Showed Error During Processing. Error Line:" + str(exc_tb.tb_lineno)+" See error below.\n")
            logging.error(traceback.format_exc())  # Logging full traceback
            continue
          pbar.update(1)
  return penman_mean_values, sebal_mean_values, irr_mean_values


def ET_Values(penman_mean_values,sebal_mean_values,irr_mean_values):
    formatted_penman_mean_values = []
    formatted_sebal_mean_values = []
    formatted_irr_mean_values = []

    for file, mean_value in penman_mean_values:
        formatted_file = file#.split('penman_eto_')[1].split('.constant.tif')[0]
        formatted_penman_mean_values.append((formatted_file, mean_value))

    for file, mean_value in sebal_mean_values:
        formatted_file = file#.split('penman_eto_')[1].split('.constant.tif')[0]
        formatted_sebal_mean_values.append((formatted_file, mean_value))

    for file, mean_value in irr_mean_values:
        formatted_file = file#.split('penman_eto_')[1].split('.constant.tif')[0]
        formatted_irr_mean_values.append((formatted_file, mean_value))

    PET_df = pd.DataFrame(formatted_penman_mean_values, columns=[feature_name, 'PET'])
    SET_df = pd.DataFrame(formatted_sebal_mean_values, columns=[feature_name, 'SET'])
    Irr_df = pd.DataFrame(formatted_irr_mean_values, columns=[feature_name, 'Irr'])
    ET_df = PET_df.merge(SET_df, on = feature_name,how = 'outer')
    ET_df = ET_df.merge(Irr_df,  on = feature_name,how = 'outer')
    ET_df.to_csv(rf'{save_data_loc}/landsat/ET_Stats.csv',index = False)


def clearData():
    logging.info("Clearing Old Data For Both The Weeks")
    sensors = ['landsat']
    # weeks = ['lastweek', 'currentweek']
    weeks = run_week
    params = ['penman','sebal','irrigation']
    for sensor in sensors:
        for param in params:
            for week in weeks:
                datadir = rf"{save_data_loc}/" + sensor + r'/' + param + r'/' + week + r'/'
                for fn in os.listdir(datadir):
                    if '.zip' in fn or '.tif' in fn:
                        os.remove(datadir + fn)
    datadir = rf"{save_data_loc}/uploads/"
    try:
      for fn in os.listdir(datadir):
          os.remove(datadir + fn)
      logging.info("Successfully Cleared The Directories")
    except Exception as logerror:
        logging.error("Found Error While Executing clearData(): " )
        logging.info("") 
        logging.error(logerror)
        pass
    params = ['precip','avgt','tmax','tmin','ugrd','vgrd','wind','temp']
    for param in params:
        datadir = rf"{save_data_loc}/{param}/"
        try:
            for fn in os.listdir(datadir):
                os.remove(datadir + fn)
            logging.info("Successfully Cleared The Directories")
        except Exception as logerror:
            logging.error("Found Error While Executing clearData(): " )
            logging.info("") 
            logging.error(logerror)
        
    
def makeDirectory():
   logging.info("Created Directory If They Did Not Exist")
   folders_list = ['landsat','uploads','precip','avgt','tmax','tmin','ugrd','vgrd','wind','uploads','logs']
   param = ['penman','sebal','irrigation']
   weeks = run_week
   for folder in folders_list:
      try:
        if folder == 'landsat':
          for params in param:
              for week in weeks:
                os.makedirs(rf"{save_data_loc}/" + folder+r'/'+ params+r'/' + week+r'/', exist_ok=True)
              
        else:
          os.makedirs(rf"{save_data_loc}/"+ folder+r'/', exist_ok=True)
        
      except Exception as logerror:
        logging.error("Found Error While Executing makeDirectory(): " )
        logging.info("") 
        logging.error(logerror)
        pass
      

def unziptiffs():
    logging.critical("Unzipping The TIFF")
    sensors = ['landsat']
    weeks = run_week
    variables = ['sebal', 'penman','irrigation'] 
    for sensor in sensors:
        for var in variables:
            for weeki in weeks:
                indir = rf"{save_data_loc}/" + sensor + r'/' + var + r'/' + weeki + r'/'  
                files_to_process = [fn for fn in os.listdir(indir) if ".zip" in fn]
                for fn in tqdm(files_to_process, desc=f"Processing {var} for {weeki}", unit=" file"):
                    if ".zip" in fn:
                        #logging.info(fn)
                        infile = indir + fn
                        try:
                            with zipfile.ZipFile(infile, "r") as zip_input:
                                zip_input.extractall(indir)
                        except Exception as logerror:
                          logging.error("Found Error While Executing unziptiffs(): " )
                          logging.info("") 
                          logging.error(logerror)
                          continue
    logging.critical("Finished Unzipping The TIFF")  

def convertTiffs():
  logging.critical("Renaming The TIFF")
  for sensor in ['landsat']:
    for weeki in run_week:
    # for weeki in ['lastweek']:
        for param in ['sebal','penman','irrigation']:
            try:
              files_to_process = [fn for fn in os.listdir(rf"{save_data_loc}" + sensor + r'/' + param + r'/' + weeki + r'/') if fn.endswith('.tif')]
              for fn in tqdm(files_to_process, desc=f"Processing {param} for {weeki}", unit=" file"):
                    in_fn= rf"{save_data_loc}"+ r'/'+sensor + r'/' + param + r'/' + weeki + r'/' + fn
                    outfn = rf"{save_data_loc}uploads/"+sensor + '_' + weeki + "_" + fn.replace('.constant','').replace('.eta','').replace('_eto','').replace('.etr','')
                    kwargs = {
                        'format': 'GTiff',
                        'noData': 0
                    }
                    with rio.open(in_fn) as src:
                      profile = src.profile
                      profile.update(driver='GTiff', nodata=0)
                      with rio.open(outfn, 'w', **profile) as dst:
                          dst.write(src.read())
            except Exception as logerror:
              logging.error("Found Error While Executing convertTiffs(): " )
              logging.info("") 
              logging.error(logerror)
              pass
  logging.critical("Finished convertTiffs")
                      
def convertingToETo():
  logging.critical("Saving The TIFF For Uploading Purpose")
  for canal in tqdm(canal_list, desc="Processing Canals TIFF to Readable Format", unit = ' canal'):
      try:
        regionid = canal[0]
        regionn = canal[0]
        regionid = regionid.replace(" ", "-")
        for sensor in ['landsat']:
          for weeki in run_week:
            penman = rf"{save_data_loc}uploads/" +sensor + '_' + weeki + "_" + "penman_" +  regionid+ ".tif"
          with rio.open(penman) as dataset:
            penman = dataset.read(1)  # Assuming you want to display the first band
            metadata = dataset.meta.copy()  # Make a copy of the metadata
          output_tif =  rf"{save_data_loc}uploads/" +sensor + '_' + weeki + "_" + "penman_eto_" +  regionid+ ".tif"
          with rio.open(output_tif, 'w', **metadata) as dst:
            dst.write(penman, 1)
          dst.close()
          dataset.close()
      except Exception as logerror:
        logging.error("Found Error While Executing convertingToETo(): " )
        logging.info("") 
        logging.error(logerror)
        pass
  logging.info("Successfully Saved TIFF In Upload Folder.")

def etinfo_update():
    logging.critical('Working on ET Info')
    try:
      for sensor in ['landsat']:
          for week_time in run_week:
              if week_time == 'currentweek':
                  firstfilepath = open(rf"{save_data_loc}" + sensor + rf'/stats_{week_time}.txt', 'r')
                  lines1 = firstfilepath.readlines()
                  
                  regions1 = []
                  penman1 = []
                  sebal1 = []
                  irri1 = []
                  for line1 in lines1:
                      if 'Region' not in line1:
                          elems = line1.split(',')
                          regions1.append(elems[0])
                          penman1.append(float(elems[1]))
                          sebal1.append(float(elems[2]))
                          irri1.append(float(elems[3]))
                  with open(rf"{save_data_loc}" + sensor + r'_stats_currentweek.txt', 'w') as txt:
                      txt.write("Region,Penman_ET,Sebal_ET,Irrigation\n")
                      for i in range(len(regions1)):
                          txt.write(regions1[i] + ',' + '{0:.3f}'.format(penman1[i]) + "," + '{0:.3f}'.format(sebal1[i]) + ',' + '{0:.3f}'.format(irri1[i]) + '\n')
                      
              elif week_time == 'lastweek':            
                  secondfilepath = open(rf"{save_data_loc}" + sensor + rf'/stats_{week_time}.txt', 'r')
                  lines2 = secondfilepath.readlines()
                  
                  regions2 = []
                  penman2 = []
                  sebal2 = []
                  irri2 = []
                  for line2 in lines2:
                      if 'Region' not in line2:
                          elems = line2.split(',')
                          regions2.append(elems[0])
                          penman2.append(float(elems[1]))
                          sebal2.append(float(elems[2]))
                          irri2.append(float(elems[3]))
                  penmanf = []
                  sebalf = []
                  logging.info((len(regions2), len(regions1)))
                  
                  for i in range(len(regions2)):
                      iiii = 0
                      for j in range(len(regions1)):
                          if regions2[i] == regions1[j]:
                              iiii = 1
                              if sebal2[i]<sebal1[j]:
                                  sebalf.append(sebal1[j]*2.0)
                              else:
                                  sebalf.append(sebal2[i])
                              # logging.info(sebal2[i], sebal1[j], sebalf[])
                              if penman2[i]<penman1[j]:
                                  penmanf.append(penman1[j]*2.0)
                              else:
                                  penmanf.append(penman2[i])
                      if iiii == 0:
                          sebalf.append(sebal2[i])
                          penmanf.append(penman2[i])
              
                  with open(rf"{save_data_loc}" + sensor + r'_stats_lastweek.txt', 'w') as txt: 
                      txt.write("Region,Penman_ET,Sebal_ET,Irrigation\n")
                      for i in range(len(regions2)):
                          txt.write(regions2[i] + ',' + '{0:.3f}'.format(penmanf[i]) + "," + '{0:.3f}'.format(sebalf[i]) + ',' + '{0:.3f}'.format(irri2[i]) + '\n')
                    
    except Exception as logerror:
      logging.error("Found Error While Executing etinfo_update(): " )
      logging.info("") 
      logging.error(logerror)
      pass
    logging.critical("Finished ET Info Successfully.")

def imergprecip():
    logging.critical('IMERG Precipitation Download Started')
    start_date_Script = datetime.datetime.strptime(start_date, '%Y-%m-%d').date()
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
    rasterstr = "{osgeo_path} gdal_calc".format(osgeo_path=osgeo_path)
    letterstr = "ABCDEFG"

    today_date = datetime.datetime.now().date()
    days_difference = (start_date_Script - today_date).days
    if abs(days_difference) <= 7:
      for dayi in range(1,8):
          deltime =  startDate + datetime.timedelta(days=-1*dayi)       
          datestr = deltime.strftime("%Y%m%d")
          yrstr = deltime.strftime("%Y")
          monthstr = deltime.strftime("%m")
          path = rf'https://{imerg_username}:{imerg_password}@jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V07B.1day.tif'    # Goto link and check if it exist, if error comes
          # print('Downloaded from: ',path)
          logging.info((yrstr, monthstr, datestr, path))
          fpath = rf"{save_data_loc}/precip/global.precip.imerg." + datestr + '.tif'
          # print('Trying to log into jsimpsonhttps')
          r = requests.get(r'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V07B.1day.tif', auth=(imerg_username, imerg_password))  # Goto link and check if it exist, if error comes
          with open(fpath, 'wb') as f:
              f.write(r.content)
          os.system(f'{osgeo_path} gdalwarp -of GTiff -tr 0.1 0.1 -te {bounds_leftlon} {bounds_bottomlat} {bounds_rightlon} {bounds_toplat} -overwrite ' + fpath + rf" {save_data_loc}/precip/precip.imerg." + datestr + '.tif')
          #os.system('C:\OSGeo4W64\OSGeo4W.bat gdalwarp -of GTiff -tr 0.1 0.1 -te 89.0 23.0 93.0 25.5 -overwrite ' + fpath + rf" {save_data_loc}precip/precip.imerg."" + datestr + '.tif')
          rasterstr = rasterstr + " -" + letterstr[dayi-1] + " " + rf"{save_data_loc}/precip/precip.imerg." + datestr + '.tif' + " "
    else:
      for dayi in range(1,8):
          deltime =  start_date_Script + datetime.timedelta(days=-1*dayi)       
          datestr = deltime.strftime("%Y%m%d")
          yrstr = deltime.strftime("%Y")
          monthstr = deltime.strftime("%m")
          if deltime.year == 2024:
              path = rf'https://{imerg_username}:{imerg_password}@jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/{deltime.year}/{deltime.month}/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V07B.1day.tif'  
              r = requests.get('https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/'+str(deltime.year)+'/'+str(deltime.month).zfill(2)+'/3B-HHR-L.MS.MRG.3IMERG.'+ datestr + '-S233000-E235959.1410.V07B.1day.tif', auth=(imerg_username, imerg_password))
          elif deltime.year == 2023:
              path = f'https://{imerg_username}:{imerg_password}@jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'+str(deltime.year)+'/'+str(deltime.month)+'/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V07B.1day.tif'  
              r = requests.get('https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'+str(deltime.year)+'/'+str(deltime.month).zfill(2)+'/3B-HHR-E.MS.MRG.3IMERG.'+ datestr + '-S233000-E235959.1410.V07B.1day.tif', auth=(imerg_username, imerg_password))
            #   print('IMERG Request Link: '+'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'+str(start_date_Script.year)+'/'+str(start_date_Script.month).zfill(2)+'/3B-HHR-E.MS.MRG.3IMERG.'+ datestr + '-S233000-E235959.1410.V07B.1day.tif')
          else:
            path = rf'https://{imerg_username}:{imerg_password}@jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/{deltime.year}/{deltime.month}/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V07B.1day.tif'    # Goto link and check if it exist, if error comes
            r = requests.get('https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'+str(deltime.year)+'/'+str(deltime.month).zfill(2)+'/3B-HHR-E.MS.MRG.3IMERG.'+ datestr + '-S233000-E235959.1410.V07B.1day.tif', auth=(imerg_username, imerg_password))
          # print('Downloaded from: ',path)
          logging.info((yrstr, monthstr, datestr, path))
          fpath = rf"{save_data_loc}/precip/global.precip.imerg." + datestr + '.tif'
          # print('Trying to log into jsimpsonhttps')
        #   r = requests.get(rf'https://{imerg_username}:{imerg_password}@jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/{start_date_Script.year}/{start_date_Script.month}/3B-HHR-E.MS.MRG.3IMERG.' + datestr + '-S233000-E235959.1410.V06E.1day.tif', auth=(imerg_username, imerg_password))  # Goto link and check if it exist, if error comes
        #   r = requests.get(path, auth=(imerg_username, imerg_password))  # Goto link and check if it exist, if error comes
          with open(fpath, 'wb') as f:
              f.write(r.content)
          # print('Just wrote',fpath)    
          os.system(f'{osgeo_path} gdalwarp -of GTiff -tr 0.1 0.1 -te {bounds_leftlon} {bounds_bottomlat} {bounds_rightlon} {bounds_toplat} -overwrite ' + fpath + rf" {save_data_loc}/precip/precip.imerg." + datestr + '.tif')
          # print('Just did mycall2') 
          #os.system('C:\OSGeo4W64\OSGeo4W.bat gdalwarp -of GTiff -tr 0.1 0.1 -te 89.0 23.0 93.0 25.5 -overwrite ' + fpath + rf" {save_data_loc}precip/precip.imerg."" + datestr + '.tif')
          rasterstr = rasterstr + " -" + letterstr[dayi-1] + " " + rf"{save_data_loc}/precip/precip.imerg." + datestr + '.tif' + " "
    
    # logging.info(rasterstr + ' --outfile=precip.currentweek1.tif --calc=(A+B+C+D+E+F+G)/10')
    logging.info(rasterstr + '--overwrite  --outfile=precip.currentweek1.tif --calc=(A+B+C+D+E+F+G)/10')
    x = os.popen(rasterstr + f"--overwrite  --outfile={save_data_loc}/precip/precip.imerg.currentweek1.tif --calc=(A+B+C+D+E+F+G)") 
    logging.info(x.read())
    try:
        
        x = os.popen(rf"{osgeo_path} gdal_calc -A {save_data_loc}/precip/precip.imerg.currentweek1.tif --overwrite --outfile={save_data_loc}/precip/precip.currentweek.tif --calc=A/10") 
        logging.info(x.read())

    except Exception as error:
                            logging.error('Error found while creating precip.currentweek.tif')
                            logging.error(error)
                            pass
  
    logging.critical("IMERG Precipitation Data Processed Successfully.")

def download_tif_from_ee(image, path, name):
    """Generate download URL from EE and download the GFS TIFF using urllib."""
    
    params = {
        'name':f'{name}',
        'filePerBand': "false",
        'scale': 25000,
        'crs': 'EPSG:4326',
        'fileFormat': 'GeoTIFF',
        'region': ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat])
    }
    url = image.getDownloadURL(params)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    urllib.request.urlretrieve(url, path)
    logging.info(f'Downloaded Successfully From GEE: {path}')

def gfsdata_ee():
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
    datestr = startDate.strftime("%Y%m%d")
    endDate = startDate + datetime.timedelta(days=1)
    dataset = ee.ImageCollection('NOAA/GFS0P25')
    dataset = dataset.filterDate(startDate.strftime('%Y-%m-%d'), endDate.strftime('%Y-%m-%d')).filterBounds(ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat]))

    # APCP download for next week, 168-hour forecast
    apcp_image = dataset.select('total_precipitation_surface').filter(ee.Filter.eq('forecast_hours', 168)).first()
    apcp_path = os.path.join(save_data_loc, "precip", f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.zip")
    download_tif_from_ee(apcp_image, apcp_path,name =f"precip_gfs_{startDate.strftime('%Y%m%d')}_nextweek" )
    extract_dir = os.path.join(save_data_loc, "precip")
    with zipfile.ZipFile(apcp_path, "r") as zip_input:
        for member in zip_input.namelist():
            extracted_path = zip_input.extract(member, extract_dir)
            # new_path = os.path.join(extract_dir, f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.tif")
            new_path = os.path.join(extract_dir, f"precip.nextweek.tif")
            os.rename(extracted_path, new_path)
    os.remove(apcp_path)
    weeks = ['currentweek', 'nextweek']
    params_ee = ['u_component_of_wind_10m_above_ground', 'v_component_of_wind_10m_above_ground', 'temperature_2m_above_ground']
    total_iterations = len(weeks) * (len(params_ee)+1) * 14  # 14 comes from the range(12, 169, 12)
    with tqdm(total=total_iterations, desc="Downloading GFS Data From GEE", unit=" file") as pbar:
        for week in weeks:
            week_date = startDate + datetime.timedelta(days=-7) if week == 'currentweek' else startDate
            week_date_end = week_date + datetime.timedelta(days=1)
            dataset = ee.ImageCollection('NOAA/GFS0P25')
            dataset = dataset.filterDate(week_date.strftime('%Y-%m-%d'), week_date_end.strftime('%Y-%m-%d')).filterBounds(ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat]))
            for param in params_ee:
                for hours in range(12, 169, 12):
                    image = dataset.select(param).filter(ee.Filter.eq('forecast_hours', hours)).first()
                    folder = 'ugrd' if param == 'u_component_of_wind_10m_above_ground' else ('vgrd' if param == 'v_component_of_wind_10m_above_ground' else 'temp')
                    path = os.path.join(save_data_loc, folder, f"{folder}.{week_date.strftime('%Y%m%d')}.{week}.{str(hours).zfill(3)}.zip")
                    download_tif_from_ee(image, path,name = f"{folder}_{week_date.strftime('%Y%m%d')}_{week}_{str(hours).zfill(3)}")
                    if folder != 'temp':
                        extract_dir = os.path.join(save_data_loc, folder)
                        with zipfile.ZipFile(path, "r") as zip_input:
                            for member in zip_input.namelist():
                                extracted_path = zip_input.extract(member, extract_dir)
                                gtifpath = os.path.join(extract_dir, f"{folder}.gfs.{startDate.strftime('%Y%m%d')}.{week}{str(hours).zfill(3)}.tif")
                                os.rename(extracted_path, gtifpath)
                        ascpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.asc'
                        finalpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                        x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
                        logging.info(x.read())
                        x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
                        logging.info(x.read())
                        logging.info((finalpath))
                        # os.remove(path)
                        pbar.update(1)
                    else:
                        folders=['tmax','tmin']
                        for folder in folders:
                            extract_dir = os.path.join(save_data_loc, folder)
                            with zipfile.ZipFile(path, "r") as zip_input:
                                for member in zip_input.namelist():
                                    extracted_path = zip_input.extract(member, extract_dir)
                                    gtifpath = os.path.join(extract_dir, f"{folder}.gfs.{startDate.strftime('%Y%m%d')}.{week}{str(hours).zfill(3)}.tif")
                                    os.rename(extracted_path, gtifpath)
                            ascpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.asc'
                            finalpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                            x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
                            logging.info(x.read())
                            x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
                            logging.info(x.read())
                            logging.info((finalpath))
                            pbar.update(1)
                    os.remove(path)


def gfsdata_noaa():
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
    logging.info("GFS forecast data is being downloaded from NCEP NOAA server for the date "+str(startDate))
    logging.info("GFS forecast data is being downloaded from NCEP NOAA server ...")
    datestr = (startDate).strftime("%Y%m%d")
    path = fr'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_APCP=on&subregion=&leftlon={bounds_leftlon}&rightlon={bounds_rightlon}&toplat={bounds_toplat}&bottomlat={bounds_bottomlat}&dir=%2Fgfs.' + datestr + '%2F00%2Fatmos'
    # print('Downloaded the data from: ', path)
    fpath = rf"{save_data_loc}/precip/precip.gfs." + datestr + '.nextweek.tif'
    logging.info(('Download Link For GFS Forecast:',path))
    logging.info(('Downloaded Path For GFS Forecast:',fpath))

    args = f'curl.exe --output {fpath} "{path}"'
    sub.run(args, shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    try:
        os.remove(f"{save_data_loc}/precip/precip.gfs.nextweek.asc")
    except OSError:
        pass
    logging.info('Downloaded the nextweek.tiff and deleted existing (if any) asc file')
    x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 2 ' + fpath + f" {save_data_loc}/precip/precip.gfs.nextweek.asc")
    logging.info(x.read())
    logging.info('Created ASC file (precip.gfs.nextweek.asc)')
    try:
        os.remove(f"{save_data_loc}/precip/precip.nextweek.tif")
    except OSError:
        pass
    try:
        x = os.popen(f'{osgeo_path} gdal_translate -of GTiff {save_data_loc}/precip/precip.gfs.nextweek.asc {save_data_loc}/precip/precip.nextweek.tif')
        logging.info(x.read())
    except Exception as error:
                            logging.info('Error found while creating precip.nextweek.tif')
                            logging.error(error)
                            print('Error Here')
                            # pass
    weeks = ['currentweek', 'nextweek']
    params = ['tmax', 'tmin', 'ugrd', 'vgrd']
    paramids = ['TMAX', 'TMIN', 'UGRD', 'VGRD']
    letterstr = "ABCDEFGHIJKLMN"
    # Total iterations for progress calculation
    total_iterations = len(weeks) * len(params) * 14  # 14 comes from the range(12, 169, 12)
    with tqdm(total=total_iterations, desc="Downloading and GFS Data From NOAA Server", unit=" file") as pbar:
        for week in weeks:
            for param in params:
                datestr = startDate.strftime("%Y%m%d")
                if week=='currentweek':
                    datestr = (startDate+datetime.timedelta(days=-7)).strftime("%Y%m%d")
                for hours in range(12, 169, 12):
                    paramid = paramids[params.index(param)]
                    serverpath = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_' + paramid + '=on&subregion=&leftlon={bounds_leftlon}&rightlon={bounds_rightlon}&toplat={bounds_toplat}&bottomlat={bounds_bottomlat}&dir=%2Fgfs.' + datestr + '%2F00%2Fatmos'
                    gtifpath = rf"{save_data_loc}" + param + r'/' + param +'.' + datestr + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                    ascpath = rf"{save_data_loc}" + param + r'/' + param + '.' + week + '.' + str(hours).zfill(3) + '.asc'
                    finalpath = rf"{save_data_loc}" + param + r'/' + param + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                    logging.info(f'GFS Link for {param}: {serverpath}')

                    try:
                        # urllib.request.urlretrieve(serverpath, gtifpath)
                        args = f'curl.exe --output {gtifpath} "{serverpath}"'
                        sub.run(args, shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
                    except:
                        serverpath = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_' + paramid + '=on&subregion=&leftlon={bounds_leftlon}&rightlon={bounds_rightlon}&toplat={bounds_toplat}&bottomlat={bounds_bottomlat}&dir=%2Fgfs.' + datestr + '%2F00'
                        # urllib.request.urlretrieve(serverpath, gtifpath)
                        args = f'curl.exe --output {gtifpath} "{serverpath}"'
                        sub.run(args, shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
                    x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
                    logging.info(x.read())
                    x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
                    logging.info(x.read())
                    logging.info((finalpath))
                    pbar.update(1)
    
def gfsdata():
    logging.critical('GFS Data Download Started')
    # startDate = datetime.datetime.today()+datetime.timedelta(days = -1)   #should be 0
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
    # startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d")##+datetime.timedelta(days=15)
    date_difference = (datetime.datetime.today() - startDate).days
    datestr = startDate.strftime("%Y%m%d")
    if date_difference > 10:
        logging.info("Switching to Google Earth Engine for GFS Forecast Data Download Due To Ten-Day Difference.")
        gfsdata_ee()
    else:
        gfsdata_noaa()      
    weeks = ['currentweek', 'nextweek']
    params = ['tmax', 'tmin', 'ugrd', 'vgrd']
    paramids = ['TMAX', 'TMIN', 'UGRD', 'VGRD']
    letterstr = "ABCDEFGHIJKLMN"
    # Total iterations for progress calculation
    total_iterations = len(weeks) * len(params)
    with tqdm(total=total_iterations, desc="Processing GFS Data", unit=" file") as pbar:
      for week in weeks:
          for param in params:
            datestr = startDate.strftime("%Y%m%d")
            if week=='currentweek':
                datestr = (startDate+datetime.timedelta(days=-7)).strftime("%Y%m%d")
            batchscript = r"{osgeo_path} gdal_calc".format(osgeo_path=osgeo_path)
            for hours in range(12, 169, 12):
                finalpath = rf"{save_data_loc}/"+ param + r'/' + param + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                batchscript = batchscript + " -" + letterstr[int(hours/12)-1] + " " + finalpath
            outpath = rf"{save_data_loc}/"+ param + r'/' + param + '.' + week + '.tif'
            logging.info((batchscript +  outpath + r' --calc=(A+B+C+D+E+F+G+H+I+J+K+L+M+N)/14'))     
            try:
                os.remove(outpath)
            except OSError:
                pass
            command = batchscript + ' --outfile=' + outpath + r' --calc=(A+B+C+D+E+F+G+H+I+J+K+L+M+N)/14'
            try:
                process = sub.Popen(command, stdout=sub.PIPE, stderr=sub.PIPE, shell=False)
                output, error = process.communicate()
            except Exception as e:
                logging.error(("If Error Exist Then See The Error:"+ str(error)))
            pbar.update(1)
        # my_call9 = [r"C:\OSGeo4W\OSGeo4W.bat",
        #     'gdal_calc',outpath,' --calc=(A+B+C+D+E+F+G+H+I+J+K+L+M+N)/14']
        # x= sub.Popen(my_call9,stdout=sub.PIPE, stderr=sub.PIPE)
        # stdout, stderr = x.communicate()
        # logging.info(stdout)
        # logging.info(x.read())
        
        ## Average temperature calculation
          tmaxpath = rf"{save_data_loc}/tmax/tmax." + week + '.tif'
          tminpath = rf"{save_data_loc}/tmin/tmin." + week + '.tif'
          avgtpath = rf"{save_data_loc}/avgt/avgt." + week + '.tif'
          try:
            os.remove(avgtpath)
          except OSError:
            pass

          logging.info((r"{osgeo_path}.format(osgeo_path=osgeo_path) gdal_calc -A " + tmaxpath + ' -B' + tminpath + ' --outfile=' + avgtpath + r' --calc=(A+B)/2'))
          x = os.popen(rf"{osgeo_path} gdal_calc -A " + tmaxpath + ' -B ' + tminpath + ' --outfile=' + avgtpath + r' --calc=(A+B)/2')
            # my_call10 = [r"C:\OSGeo4W\OSGeo4W.bat",
            #         'gdal_calc','-A',tmaxpath,'-B',tminpath,avgtpath, r'--calc=(A+B)/2']
            # x= sub.Popen(my_call10,stdout=sub.PIPE, stderr=sub.PIPE)
            # stdout, stderr = x.communicate()
            # logging.info(stdout)
            
          logging.info(x.read())
                
            ## Average wind speed calculation
          ugrdpath = rf"{save_data_loc}/ugrd/ugrd." + week + '.tif'
          vgrdpath = rf"{save_data_loc}/vgrd/vgrd." + week + '.tif'
          windpath = rf"{save_data_loc}/wind/wind." + week + '.tif'
          try:
              os.remove(windpath)
              os.remove(vgrdpath)
              os.remove(ugrdpath)
          except OSError:
              pass
          logging.info(rf"{osgeo_path}.format(osgeo_path=osgeo_path) gdal_calc -A " + ugrdpath + ' -B ' + vgrdpath + '  --outfile=' + windpath + r' --calc=sqrt(A*A+B*B)')
            # my_call10 = [r"C:\OSGeo4W\OSGeo4W.bat",
            #         'gdal_calc','-A',ugrdpath,'-B',vgrdpath,windpath, r'--calc=sqrt(A*A+B*B)']
            # x= sub.Popen(my_call10,stdout=sub.PIPE, stderr=sub.PIPE)
            # stdout, stderr = x.communicate()
            # logging.info(stdout)
          x = os.popen(rf"{osgeo_path} gdal_calc -A " + ugrdpath + ' -B ' + vgrdpath + '  --outfile=' + windpath + r' --calc=sqrt(A*A+B*B)')
          logging.info(x.read())
    logging.critical('Finished GFS Data Download And Processing')     


def get_mean_value(raster, geometry):
    """Function to get mean value of raster data for the given geometry"""
    out_image, out_transform = riomask(raster, [geometry], crop=True)
    # Assuming that the data of interest is in the first band
    data = out_image[0]
    valid_data = data[data != raster.nodatavals[0]]
    return valid_data.mean()

def get_mean_value_precip(raster, geometry):
    """Function to get mean value of precipitation raster data for the given geometry"""
    masked_data, _ = riomask(raster, geometry, crop=True, all_touched = True)
    nodata_value = raster.nodata
    data_no_nodata = masked_data[masked_data != nodata_value]
    return np.ma.mean(data_no_nodata)

def union_info():
    canal_shp =gpd.read_file(irrigation_canals_path,crs='EPSG:4326')
    if canal_shp.crs is None:
      canal_shp.crs = 'EPSG:4326'
    sensors = ['landsat']
    params = ['precip']
    weeki = ['nextweek']
    logging.critical('Started Command Area Info')
    # rasterpath = rf"{save_data_loc}/landsat"  + r'/irrigation/' + run_week[0] + r'/irrigation_' + canal_list[1][0].replace(' ','-') + '.eta.tif'  
    i = 0
    max_attempts = 100  # or another appropriate limit
    while i < max_attempts:
        rasterpath = rf"{save_data_loc}/landsat" + r'/irrigation/' + run_week[0] + r'/irrigation_' + canal_list[i][0].replace(' ', '-') + '.eta.tif'
        if os.path.exists(rasterpath):
            break
        i += 1
    raster = rio.open(rasterpath)
    raster_crs = raster.crs
    raster.close()
    canal_shp_reprojected = canal_shp.to_crs(raster_crs)
    canal_shp_reprojected['Currentweek PPT'] = np.nan
    canal_shp_reprojected['Nextweek PPT'] = np.nan
    simple_canal_list = [item[0] for item in canal_list]
    canal_shp_reprojected.set_index(feature_name, inplace=True)
    canal_shp_reprojected = canal_shp_reprojected.reindex(simple_canal_list)
    canal_shp_reprojected.reset_index(inplace=True)
    if set(run_week) == {'currentweek', 'lastweek'}:
        total_canals = len(canal_list) + 1
    else:
        total_canals = len(canal_list) * len(run_week)
    with tqdm(total=total_canals, desc="Processing Canal And Creating Command Area Stats", unit=" canal") as pbar:
        if set(run_week) == {'currentweek', 'lastweek'}:
            canal_shp_reprojected['7Day Irrigation'] = np.nan
            canal_shp_reprojected['14Day Irrigation'] = np.nan
            ### Adding Penman and SEBAL ET Below
            canal_shp_reprojected['7Day Penman ET'] = np.nan
            canal_shp_reprojected['14Day Penman ET'] = np.nan
            canal_shp_reprojected['7Day SEBAL ET'] = np.nan
            canal_shp_reprojected['14Day SEBAL ET'] = np.nan
            for week in run_week:
                if week == 'currentweek':
                    for i in range(len(canal_list)):
                        rasterpath = rf"{save_data_loc}/landsat"  + r'/irrigation/' + week + r'/irrigation_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                        ### Adding Penman and SEBAL ET Below
                        rasterpath_penman = rf"{save_data_loc}/landsat"  + r'/penman/' + week + r'/penman_eto_' + canal_list[i][0].replace(' ','-') + '.constant.tif' 
                        rasterpath_sebal = rf"{save_data_loc}/landsat"  + r'/sebal/' + week + r'/sebal_eto_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                        try:

                            ### Adding Penman and SEBAL ET Below
                            ## Penman
                            with rio.open(rasterpath_penman) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['7Day Penman ET'][i] = masked_mean

                            ## Sebal
                            with rio.open(rasterpath_sebal) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['7Day SEBAL ET'][i] = masked_mean
                            ### Finished Adding Penman and SEBAL ET Below

                            pbar.update(1)
                        except Exception as error:
                            logging.error('Error Found for currentweek for',canal_list[i][0], 'Error Message:',error)
                            continue
                    canal_shp_reprojected['7Day Irrigation'] = canal_shp_reprojected['7Day SEBAL ET'] - canal_shp_reprojected['7Day Penman ET']
                    percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
                    canal_shp_reprojected = canal_shp_reprojected.merge(percolation_df, on = [f'{feature_name}'])
                elif week == 'lastweek':
                    for i in range(len(canal_list)):
                        rasterpath = rf"{save_data_loc}/landsat"  + r'/irrigation/' + week + r'/irrigation_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                        ### Adding Penman and SEBAL ET Below
                        rasterpath_penman = rf"{save_data_loc}/landsat"  + r'/penman/' + week + r'/penman_eto_' + canal_list[i][0].replace(' ','-') + '.constant.tif' 
                        rasterpath_sebal = rf"{save_data_loc}/landsat"  + r'/sebal/' + week + r'/sebal_eto_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                        try:

                            ### Adding Penman and SEBAL ET Below
                            ## Penman
                            with rio.open(rasterpath_penman) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['14Day Penman ET'][i] = masked_mean

                            ## Sebal
                            with rio.open(rasterpath_sebal) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['14Day SEBAL ET'][i] = masked_mean
                            ### Finished Adding Penman and SEBAL ET Below

                        except Exception as error:
                            logging.error('Error Found for lastweek for',canal_list[i][0], 'Error Message:',error)
                            continue
                    
                        rasterpath = f"{save_data_loc}/precip/precip.currentweek.tif"
                        raster = rio.open(rasterpath)
                        canal_shp_reprojected['Currentweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])
                        rasterpath = f"{save_data_loc}/precip/precip.nextweek.tif"
                        raster = rio.open(rasterpath)
                        canal_shp_reprojected['Nextweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])
                    canal_shp_reprojected['14Day Irrigation'] = canal_shp_reprojected['14Day SEBAL ET'] - canal_shp_reprojected['14Day Penman ET']
                    percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
                    canal_shp_reprojected = canal_shp_reprojected.merge(percolation_df, on = [f'{feature_name}'])
                    canal_shp_reprojected['net_water_req'] = canal_shp_reprojected['Currentweek PPT'] + canal_shp_reprojected['Nextweek PPT'] + canal_shp_reprojected['14Day Irrigation'] - canal_shp_reprojected['MedianPercolation'] * 7
                    pbar.update(1)
        elif set(run_week) == {'currentweek'}:
            canal_shp_reprojected['7Day Irrigation'] = np.nan
            ### Adding Penman and SEBAL ET Below
            canal_shp_reprojected['7Day Penman ET'] = np.nan
            canal_shp_reprojected['7Day SEBAL ET'] = np.nan
            for week in run_week:
                for i in range(len(canal_list)):
                    rasterpath = rf"{save_data_loc}/landsat"  + r'/irrigation/' + week + r'/irrigation_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                    ### Adding Penman and SEBAL ET Below
                    rasterpath_penman = rf"{save_data_loc}/landsat"  + r'/penman/' + week + r'/penman_eto_' + canal_list[i][0].replace(' ','-') + '.constant.tif' 
                    rasterpath_sebal = rf"{save_data_loc}/landsat"  + r'/sebal/' + week + r'/sebal_eto_' + canal_list[i][0].replace(' ','-') + '.eta.tif'  
                    try:
                        ### Adding Penman and SEBAL ET Below
                        ## Penman
                        with rio.open(rasterpath_penman) as src:
                            raster_crs = src.crs
                            raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                            raster_data = src.read(1, masked=True)
                            masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                        # Calculate the mean value of the masked raster
                        # masked_mean = np.ma.mean(masked_data)
                        masked_mean = np.mean(raster_data)
                        canal_shp_reprojected['7Day Penman ET'][i] = masked_mean

                        ## Sebal
                        with rio.open(rasterpath_sebal) as src:
                            raster_crs = src.crs
                            raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                            raster_data = src.read(1, masked=True)
                            masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                        # Calculate the mean value of the masked raster
                        # masked_mean = np.ma.mean(masked_data)
                        masked_mean = np.mean(raster_data)
                        canal_shp_reprojected['7Day SEBAL ET'][i] = masked_mean
                        ### Finished Adding Penman and SEBAL ET Below


                        rasterpath = f"{save_data_loc}/precip/precip.currentweek.tif"
                        raster = rio.open(rasterpath)
                        canal_shp_reprojected['Currentweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])

                        rasterpath = f"{save_data_loc}/precip/precip.nextweek.tif"
                        raster = rio.open(rasterpath)
                        canal_shp_reprojected['Nextweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])
                        pbar.update(1)
                    except Exception as error:
                            logging.error('Error Found for currentweek for',canal_list[i][0], 'Error Message:',error)
                            continue
                canal_shp_reprojected['7Day Irrigation'] = canal_shp_reprojected['7Day SEBAL ET'] - canal_shp_reprojected['7Day Penman ET']
                percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
                canal_shp_reprojected = canal_shp_reprojected.merge(percolation_df, on = [f'{feature_name}'])
                canal_shp_reprojected['net_water_req'] = canal_shp_reprojected['Currentweek PPT'] + canal_shp_reprojected['Nextweek PPT'] + canal_shp_reprojected['7Day Irrigation'] - canal_shp_reprojected['MedianPercolation'] * 7
        elif set(run_week) == {'lastweek'}:
            canal_shp_reprojected['14Day Irrigation'] = np.nan
            ### Adding Penman and SEBAL ET Below
            canal_shp_reprojected['14Day Penman ET'] = np.nan
            canal_shp_reprojected['14Day SEBAL ET'] = np.nan
            for week in run_week:
                for i in range(len(canal_list)):
                    rasterpath = rf"{save_data_loc}/landsat"  + r'/irrigation/' + week + r'/irrigation_' + canal_list[i][0].replace(' ','-') + '.eta.tif' 
                    ### Adding Penman and SEBAL ET Below
                    rasterpath_penman = rf"{save_data_loc}/landsat"  + r'/penman/' + week + r'/penman_eto_' + canal_list[i][0].replace(' ','-') + '.constant.tif' 
                    rasterpath_sebal = rf"{save_data_loc}/landsat"  + r'/sebal/' + week + r'/sebal_eto_' + canal_list[i][0].replace(' ','-') + '.eta.tif'   
                    try:
                            ### Adding Penman and SEBAL ET Below
                            ## Penman
                            with rio.open(rasterpath_penman) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['14Day Penman ET'][i] = masked_mean

                            ## Sebal
                            with rio.open(rasterpath_sebal) as src:
                                raster_crs = src.crs
                                raster_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                                raster_data = src.read(1, masked=True)
                                masked_data, _ = riomask(src, canal_shp_reprojected.geometry, crop=True)
                            # Calculate the mean value of the masked raster
                            # masked_mean = np.ma.mean(masked_data)
                            masked_mean = np.mean(raster_data)
                            canal_shp_reprojected['14Day SEBAL ET'][i] = masked_mean
                            ### Finished Adding Penman and SEBAL ET Below

                            rasterpath = f"{save_data_loc}/precip/precip.currentweek.tif"
                            raster = rio.open(rasterpath)
                            canal_shp_reprojected['Currentweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])

                            rasterpath = f"{save_data_loc}/precip/precip.nextweek.tif"
                            raster = rio.open(rasterpath)
                            canal_shp_reprojected['Nextweek PPT'][i] = get_mean_value_precip(raster, canal_shp_reprojected['geometry'][i:i+1])
                            pbar.update(1)
                    except Exception as error:
                            logging.error('Error Found for lastweek for',canal_list[i][0], 'Error Message:',error)
                            continue
                canal_shp_reprojected['14Day Irrigation'] = canal_shp_reprojected['14Day SEBAL ET'] - canal_shp_reprojected['14Day Penman ET']
                percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
                canal_shp_reprojected = canal_shp_reprojected.merge(percolation_df, on = [f'{feature_name}'])
                canal_shp_reprojected['net_water_req'] = canal_shp_reprojected['Currentweek PPT'] + canal_shp_reprojected['Nextweek PPT'] + canal_shp_reprojected['14Day Irrigation'] - canal_shp_reprojected['MedianPercolation'] * 7
    # except Exception as error:
    #     logging.error('Error Found:',error)  
    #     pass
    canal_shp_reprojected = canal_shp_reprojected.drop(columns=['geometry'])
    
    # canal_shp_reprojected['net_water_req'] = canal_shp_reprojected['Currentweek PPT'] + canal_shp_reprojected['Nextweek PPT'] + canal_shp_reprojected['7Day Irrigation']
    canal_shp_reprojected['net_water_req'] = canal_shp_reprojected['net_water_req'].apply(lambda x: x if x <= 0 else 0) 
    sorted_canal_shp_reprojected = canal_shp_reprojected.sort_values(by='ID')
    sorted_canal_shp_reprojected.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats.csv",index = False)
    logging.critical('Finished Command Area Info')

def converrt_ETtiff2images(variable, folder_path=f'{save_data_loc}/uploads/'):
    pattern = os.path.join(folder_path, f"landsat_*_{variable}_*.tif*")  
    tiff_files = glob.glob(pattern)
    if not tiff_files:
        return

    # Iterate through all matching TIFF files with tqdm progress bar
    for tiff_file in tqdm(tiff_files, desc="Converting ET TIF Files to PNGs", unit=" file"):
        filename = os.path.basename(tiff_file)
        region_name = filename.split('_')[-1].split('.tif')[0]

        with rio.open(tiff_file) as src:
            band = src.read(1)
            nodata = src.nodata
            bounds = src.bounds
            crs = src.crs

        if nodata is not None:
            band = np.where(band == nodata, np.nan, band)

        # Extended bounds calculation
        margin = 0.1
        width = bounds.right - bounds.left
        height = bounds.top - bounds.bottom

        extended_bounds = [
            bounds.left - margin * width,
            bounds.right + margin * width,
            bounds.bottom - margin * height,
            bounds.top + margin * height,
        ]

        new_bounds = BoundingBox(
            left=extended_bounds[0],
            bottom=extended_bounds[2],
            right=extended_bounds[1],
            top=extended_bounds[3]
        )

        # Plot settings
        vmin = np.nanmin(band)
        vmax = np.nanmax(band)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("Blues")
        rgba_image = cmap(norm(band))
        rgba_image[np.isnan(band)] = [0, 0, 0, 0]

        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = plt.axes(projection=ccrs.PlateCarree())

        img = ax.imshow(rgba_image, extent=[bounds.left, bounds.right, bounds.bottom, bounds.top],
                        origin='upper', transform=ccrs.PlateCarree(), zorder=2, alpha=0.7)

        ax.set_xlim(extended_bounds[0], extended_bounds[1])
        ax.set_ylim(extended_bounds[2], extended_bounds[3])

        add_basemap(ax, crs=ccrs.PlateCarree(), source='https://a.tile.openstreetmap.org/{z}/{x}/{y}.png')

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False

        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.05)
        cbar.set_label('Evapotranspiration (mm)')

        # North arrow
        arrow_x = new_bounds.right - 0.1 * (new_bounds.right - new_bounds.left)
        arrow_y = new_bounds.top - 0.05 * (new_bounds.top - new_bounds.bottom)
        ax.annotate('N', xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y - 0.0008),
                    arrowprops=dict(facecolor='black', arrowstyle='simple', lw=1.5),
                    ha='center', va='center', fontsize=14, fontweight='bold', transform=ccrs.PlateCarree())

        # Scalebar position and length
        latitude_midpoint = (new_bounds.bottom + new_bounds.top) / 2
        meters_per_degree = 111320 * np.cos(np.deg2rad(latitude_midpoint))
        scalebar_length_deg = 100 / meters_per_degree

        x_start = new_bounds.right - 0.15 * (new_bounds.right - new_bounds.left)
        y_start = new_bounds.top - 0.25 * (new_bounds.top - new_bounds.bottom)
        ax.plot([x_start, x_start + scalebar_length_deg], [y_start, y_start],
                color='black', linewidth=3, transform=ccrs.PlateCarree())
        ax.text(x_start + scalebar_length_deg / 2, y_start + 0.0005 * (new_bounds.top - new_bounds.bottom),
                '100 m', ha='center', va='bottom', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
        date_png = datetime.datetime.strptime(start_date, "%Y-%m-%d").strftime("%Y-%m-%d")
        # Set title and save
        if variable == 'penman':
            ax.set_title(f'Penman-Montieth (Demand) ET for {region_name} on {date_png}', size=15)
            png_filename = filename.replace(".tiff", ".png").replace(".tif", ".png")
            plt.tight_layout()
            os.makedirs(f'{save_data_loc}/PNGs/', exist_ok=True)
            plt.savefig(f'{save_data_loc}/PNGs/{region_name}_DemandET_{date_png}.png', bbox_inches='tight')

        elif variable == 'sebal':
            ax.set_title(f'SEBAL (Observed) ET for {region_name} on {date_png}', size=15)
            png_filename = filename.replace(".tiff", ".png").replace(".tif", ".png")
            plt.tight_layout()
            os.makedirs(f'{save_data_loc}/PNGs/', exist_ok=True)
            plt.savefig(f'{save_data_loc}/PNGs/{region_name}_ObservedET_{date_png}.png', bbox_inches='tight')


        
def converrt_Irrigationtiff2images(folder_path=f'{save_data_loc}/uploads/'):
    pattern = os.path.join(folder_path, "landsat_*_irrigation_*.tif*")  
    tiff_files = glob.glob(pattern)
    if not tiff_files:
        return
    
    # Iterate through all matching TIFF files with tqdm progress bar
    for tiff_file in tqdm(tiff_files, desc="Converting Irrigation TIF Files to PNGs", unit=" file"):
        filename = os.path.basename(tiff_file)
        region_name = filename.split('_')[-1].split('.tif')[0]

        # Load the TIFF file
        with rio.open(tiff_file) as src:
            band = src.read(1)
            nodata = src.nodata
            bounds = src.bounds
            crs = src.crs

        # Replace nodata values with NaN
        if nodata is not None:
            band = np.where(band == nodata, np.nan, band)

        # Calculate extended bounds
        margin = 0.1
        width = bounds.right - bounds.left
        height = bounds.top - bounds.bottom
        extended_bounds = [
            bounds.left - margin * width,
            bounds.right + margin * width,
            bounds.bottom - margin * height,
            bounds.top + margin * height,
        ]

        new_bounds = BoundingBox(
            left=extended_bounds[0],
            bottom=extended_bounds[2],
            right=extended_bounds[1],
            top=extended_bounds[3]
        )

        # Processing the data and plotting steps (unchanged)
        abs_max = np.nanmax(np.abs(band))
        norm = TwoSlopeNorm(vcenter=0, vmin=-abs_max, vmax=abs_max)
        cmap = plt.get_cmap("RdBu")
        rgba_image = cmap(norm(band))
        rgba_image[np.isnan(band)] = [0, 0, 0, 0]

        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = plt.axes(projection=ccrs.PlateCarree())
        img = ax.imshow(rgba_image, extent=[bounds.left, bounds.right, bounds.bottom, bounds.top],
                        origin='upper', transform=ccrs.PlateCarree(), zorder=2, alpha=0.7)

        ax.set_xlim(extended_bounds[0], extended_bounds[1])
        ax.set_ylim(extended_bounds[2], extended_bounds[3])

        add_basemap(ax, crs=ccrs.PlateCarree(), source='https://a.tile.openstreetmap.org/{z}/{x}/{y}.png')

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False

        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.05)
        cbar.set_label('Deficit (Red) / Surplus (Blue) Irrigation')

        tick_values = np.linspace(-abs_max, abs_max, num=11)
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels([f"{val:.1f}" for val in tick_values])

        # Add custom North arrow
        arrow_x = new_bounds.right - 0.1 * (new_bounds.right - new_bounds.left)
        arrow_y = new_bounds.top - 0.05 * (new_bounds.top - new_bounds.bottom)
        ax.annotate('N', xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y - 0.0008),
                    arrowprops=dict(facecolor='black', arrowstyle='simple', lw=1.5),
                    ha='center', va='center', fontsize=14, fontweight='bold', transform=ccrs.PlateCarree())

        # Calculate scalebar position and plot
        latitude_midpoint = (new_bounds.bottom + new_bounds.top) / 2
        meters_per_degree = 111320 * np.cos(np.deg2rad(latitude_midpoint))
        scalebar_length_deg = 100 / meters_per_degree

        x_start = new_bounds.right - 0.15 * (new_bounds.right - new_bounds.left)
        y_start = new_bounds.top - 0.25 * (new_bounds.top - new_bounds.bottom)
        ax.plot([x_start, x_start + scalebar_length_deg], [y_start, y_start],
                color='black', linewidth=3, transform=ccrs.PlateCarree())
        ax.text(x_start + scalebar_length_deg / 2, y_start + 0.0005 * (new_bounds.top - new_bounds.bottom),
                '100 m', ha='center', va='bottom', fontsize=10, fontweight='bold', transform=ccrs.PlateCarree())
        date_png = datetime.datetime.strptime(start_date, "%Y-%m-%d").strftime("%Y-%m-%d")
        ax.set_title(f'Deficit or Surplus Irrigation for {region_name} on {date_png}', size=15)

        png_filename = filename.replace(".tiff", ".png").replace(".tif", ".png")
        plt.tight_layout()
        os.makedirs(f'{save_data_loc}/PNGs/', exist_ok=True)
        plt.savefig(f'{save_data_loc}/PNGs/{region_name}_DeficitET_{date_png}.png', bbox_inches='tight')
        
def percolation_estimation():
    logging.critical('Started Percolation Estimation')
    for wktime in run_week:
        try:
            if wktime == "currentweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 8+2)
            elif wktime == "lastweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 8+2+7)
            # start_date = (datetime.today() - datetime.timedelta(days=7)).strftime('%Y-%m-%d')
            startDate=ee.Date(start_date)
            endDate=ee.Date(enddate)
            # Load the shapefile as a feature collection
            percolation_feature_name = f'{feature_name}'
            logging.critical('Running Week:'+str(wktime))
            logging.critical("Running Week's Start Date:"+str(start_date))
            logging.critical("Running Week's End Date:"+str(enddate))

            field_capacity = ee.Image("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01").select('b10').divide(ee.Number(100)).clip(irrigation_canals)

            # Initialize a list to store dataframes for each region
            region_dataframes = []

            # Use tqdm to add a progress bar for the number of regions
            region_list = irrigation_canals.toList(irrigation_canals.size()).getInfo()
            for region in tqdm(region_list, desc="Estimating Soil Mositure", unit=" Command Area"):
                region_feature = ee.Feature(region)
                region_id = region_feature.get(f'{feature_name}').getInfo()
                region_geometry = region_feature.geometry()

                # Load Sentinel-1 Image Collection for the last 7 days and apply speckle filtering
                s1_collection = (ee.ImageCollection('COPERNICUS/S1_GRD')
                                .filterBounds(region_geometry)
                                .filterDate(startDate, endDate)
                                .filter(ee.Filter.eq('instrumentMode', 'IW'))
                                .select(['VV']))

                # Apply smoothing, calculate soil moisture, and create a feature collection for the output
                def process_image(image):
                    # Apply focal median filter for smoothing
                    smoothed = image.addBands(image.focal_max(30, 'circle', 'meters').rename("Smooth"))

                    # Calculate wet and dry indices using smoothed collection
                    wet_index = s1_collection.max().select('VV')
                    dry_index = s1_collection.min().select('VV')
                    sensitivity = wet_index.subtract(dry_index)

                    # Define urban and water masks
                    urban_mask = smoothed.select('Smooth').gt(-6)
                    water_mask = smoothed.select('Smooth').lt(-17)

                    # Calculate soil moisture
                    Mv = smoothed.select("Smooth").subtract(dry_index).divide(sensitivity)
                    Mv = Mv.updateMask(water_mask.Not()).updateMask(urban_mask.Not())

                    # Upscale and clamp soil moisture
                    Mv_upscaled = Mv.reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True).reproject(crs=Mv.projection(), scale=250)
                    Mv_upscaled_clamped = Mv_upscaled.clamp(0, 0.6)

                    # Calculate mean soil moisture in the region
                    median_ssm = Mv_upscaled_clamped.reduceRegion(
                        reducer=ee.Reducer.median(),
                        geometry=region_geometry,
                        scale=250,
                        maxPixels=1e9
                    ).get('Smooth')

                    # Calculate field capacity
                    median_field_capacity = field_capacity.reduceRegion(
                        reducer=ee.Reducer.median(),
                        geometry=region_geometry,
                        scale=250,
                        maxPixels=1e9
                    ).get('b10')
                    # Return as a feature with aggregated data for the period
                    return ee.Feature(None, {
                        percolation_feature_name: region_id,
                        'MedianSoilMoisture': median_ssm,
                        'MedianFieldCapacity': median_field_capacity
                    })

                # Apply the processing function to each image in the smoothed collection
                soil_moisture_features = s1_collection.map(process_image)

                # Retrieve data for each region and create a DataFrame
                ssm_features = soil_moisture_features.getInfo()
                week_data = [{percolation_feature_name: f['properties'][percolation_feature_name],
                            'MedianSoilMoisture': f['properties']['MedianSoilMoisture'],
                            'MedianFieldCapacity': f['properties']['MedianFieldCapacity']}
                            for f in ssm_features['features']]
                
                # Append the data for the current region to the list
                region_df = pd.DataFrame(week_data)
                region_dataframes.append(region_df)
                logging.info('Estimated Soil Moisture for Region:'+str(region_id)+' For week:'+str(wktime))
            # Combine all regional dataframes into a single DataFrame
            final_df = pd.concat(region_dataframes, ignore_index=True)
            # Group by region and calculate the mean values for each metric
            final_means_df = final_df.groupby(percolation_feature_name).mean().reset_index()
            final_means_df['MedianPercolation'] = (final_means_df['MedianSoilMoisture'] - final_means_df['MedianFieldCapacity']).clip(lower=0)
            os.makedirs(f'{save_data_loc}/percolation', exist_ok = True)
            final_means_df.to_csv(f'{save_data_loc}/percolation/Percolation_{wktime}.csv', index = False)

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            # print(f"Error found: {e}")
            # print('Error on', wktime.capitalize())
            # print(f"Exception in {fname} on line {exc_tb.tb_lineno}")
            logging.error("Week: " + str(wktime) + ", Region: " + str(region_id) + " Showed Error During Processing. Error Line:" + str(exc_tb.tb_lineno)+" See error below.\n")
            logging.error(traceback.format_exc())  # Logging full traceback
            continue    

def ensure_feature_collection(roi):
    if isinstance(roi, ee.geometry.Geometry):
        # If ROI is a Geometry, convert it to a Feature, then to a FeatureCollection
        return ee.FeatureCollection([ee.Feature(roi)])
    elif isinstance(roi, ee.feature.Feature):
        # If ROI is a Feature, convert it to a FeatureCollection
        return ee.FeatureCollection([roi])
    elif isinstance(roi, ee.featurecollection.FeatureCollection):
        # If ROI is already a FeatureCollection, use it directly
        return roi
    else:
        raise ValueError("ROI is neither a Geometry, Feature, nor a FeatureCollection")
    
if __name__ == '__main__':
    makeDirectory()
    if clear_directory_condition:
        clearData()
    if ET_estimation:
        penman_mean_values,sebal_mean_values,irr_mean_values = sebal_eto_Landsat()
        ET_Values(penman_mean_values = penman_mean_values,sebal_mean_values= sebal_mean_values,irr_mean_values= irr_mean_values)
        unziptiffs()
        convertTiffs()
        convertingToETo()
    if precipitation_condition:
        imergprecip()
    if weather_condition:
        gfsdata()
    if percolation_condition:
        percolation_estimation()
    if region_stats_condition:
        union_info()