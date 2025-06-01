# Importing required libraries 
import ee 
import logging
import colorlog
from tqdm import tqdm
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
from ruamel.yaml import YAML
import warnings
warnings.filterwarnings("ignore")

# Initialise Earth Engine
# ee.Authenticate() # If you are running the earth engine first time in your machine, you need to run the ee.Authenticate() command first.
ee.Initialize()

#########*******Module (1) USER INPUT (ROI and Intereset Dates and Planting Date)********########
# Initialize ruamel.yaml YAML object
yaml = YAML()
yaml.preserve_quotes = True  # Optional: preserves quotes around strings

with open('../config_files/Script_Config.yaml', 'r') as f:
    script_config = yaml.load(f)

# Accessing the information from 'Irrigation_cmd_area_shapefile' section
irrigation_cmd_area_path = script_config['Irrigation_cmd_area_shapefile']['path']
feature_name = script_config['Irrigation_cmd_area_shapefile']['feature_name']

bounds_leftlon = float(script_config['Irrigation_cmd_area_shapefile_Bounds']['leftlon'])
bounds_rightlon = float(script_config['Irrigation_cmd_area_shapefile_Bounds']['rightlon'])
bounds_toplat = float(script_config['Irrigation_cmd_area_shapefile_Bounds']['toplat'])
bounds_bottomlat = float(script_config['Irrigation_cmd_area_shapefile_Bounds']['bottomlat'])

save_data_loc = script_config['Save_Data_Location']['save_data_loc']

gee_asset_section = script_config.get('GEE_Asset_ID', {})

if gee_asset_section.get('id'): 
    gee_asset_id = gee_asset_section['id']
elif gee_asset_section.get('shp'):
    gee_asset_id = gee_asset_section['shp']
else:
    raise ValueError(
        "Configuration error: 'GEE_Asset_ID' must contain either 'id' or 'shp'."
        "Both are missing or empty."
    )

# Accessing the information from 'Date_Running' section
start_date = script_config['Date_Running']['start_date']
default_run_week = script_config['Date_Running']['default_run_week']
run_week = (
    ["lastweek", "currentweek"]
    if default_run_week
    else script_config['Date_Running'].get('run_week', [])
)

# Accessing the information from 'OSGEO_Path' section
osgeo_path = script_config['OSGEO_Path']['path']
canal_config_path = script_config['Canal_Config']['path']

# === Conditional Flags ===
clear_directory_condition = script_config['Clean_Directory']['clear_directory_condition']
ET_estimation = script_config['Run_ET_Estimation']['ET_estimation']
glcc_mask = script_config['GLCC_Mask'].get('glcc_mask', False)
precipitation_condition = script_config['Precipitation_Config']['consider_preciptation']
weather_condition = script_config['Weather_Config']['consider_forecasted_weather']
percolation_condition = script_config['Percolation_Config']['consider_percolation']
region_stats_condition = script_config['Region_stats']['estimate_region_stats']
tiff_to_png = script_config['Tiff_to_PNGs'].get('tiff_to_pngs', False)

# === Accessing Secrets Configuration Files ===
secrets_file_path = script_config['Secrets_Path']['path']
with open(rf'{secrets_file_path}', 'r') as f:
    secrets = yaml.load(f)
imerg_username = script_config['IMERG_Account']['username']
imerg_password = script_config['IMERG_Account']['password']

# === SETUP LOGGING ===
os.makedirs(f'{save_data_loc}',exist_ok = True)
os.makedirs(f'{save_data_loc}/logs/',exist_ok = True)
dt_fmt = '%Y%m%d_%H%M%S'
logging.basicConfig(filename=f'{save_data_loc}/logs/{datetime.datetime.today().strftime(dt_fmt)}.log', format='%(levelname)s - %(asctime)s :%(message)s', level=logging.DEBUG)
logging.info("")
logging.info("Location of the Data Folder: "+save_data_loc)
logging.info("")
logging.critical('Start Date:'+str(start_date) )
logging.critical('Run Week:'+ str(run_week) )


irrigation_cmd_area= ee.FeatureCollection(gee_asset_id) 
canal_list = irrigation_cmd_area.reduceColumns(ee.Reducer.toList(1), [feature_name]).get('list').getInfo()
glcc = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")

# """ Reading the canal config file to parse the planting date, 
# crop type and soil coefficient for each command region """

with open(canal_config_path, 'r') as f:
    config = yaml.safe_load(f)

# Extract defaults and check use_default
defaults = config.get('default_settings', {})
use_default = defaults.get('use_default', False)

def read_cmd_area_settings(canal_data, defaults):
    """Read settings for a specific canal using fallback from defaults."""
    if canal_data.get('use_default', False):
        return {
            'planting_date': defaults['planting_date'],
            'crop_type': defaults['crop_type'],
            'soil_coefficient': defaults['soil_coefficient']
        }
    else:
        return {
            'planting_date': canal_data.get('planting_date', defaults['planting_date']),
            'crop_type': canal_data.get('crop_type', defaults['crop_type']),
            'soil_coefficient': canal_data.get('soil_coefficient', defaults['soil_coefficient'])
        }

# Build canal settings dictionary
if use_default:
    canal_settings = {'DEFAULT': {
        'planting_date': defaults['planting_date'],
        'crop_type': defaults['crop_type'],
        'soil_coefficient': defaults['soil_coefficient']
    }}
else:
    canal_settings = {}
    for canal_name, canal_data in config.get('canals', {}).items():
        canal_settings[canal_name] = read_cmd_area_settings(canal_data, defaults)


def read_crop_coefficients(config_path, crop):
    """Read crop coefficient values from YAML config."""
    with open(config_path, 'r') as f:
        crop_config = yaml.safe_load(f)
    
    crop_data = crop_config.get(crop, {})
    coefficients = {}

    for key, value in crop_data.items():
        start_day, end_day = map(int, key.split('-'))
        if isinstance(value, list) and value[0] == 'linear':
            # Format: [linear, start_value, end_value, days]
            coefficients[(start_day, end_day)] = ('linear', float(value[1]), float(value[2]), int(value[3]))
        else:
            coefficients[(start_day, end_day)] = float(value)
    
    return coefficients

def get_growth_kc(coefficients, num_days):
    """Calculate crop coefficient (Kc) based on planting date and satellite overpass date."""
    for day_range, kc in coefficients.items():
        if day_range[0] <= num_days <= day_range[1]:
            if isinstance(kc, tuple) and kc[0] == 'linear':
                # Linear interpolation
                _, start_value, end_value, days = kc
                days_into_range = num_days - day_range[0]
                delta_per_day = (end_value - start_value) / days
                return start_value + delta_per_day * days_into_range
            return kc
    return None  # If num_days is outside defined ranges


penman_mean_values = []
sebal_mean_values = []
irr_mean_values = []

#########*******Module (2) Estimates Penman-Monteith and SEBAL Evapotranspiration********########

def sebal_eto_Landsat(logger = None):
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
      logger.critical('Running Week:'+str(wktime))
      logger.critical("Running Week's Start Date:"+str(startdate))
      logger.critical("Running Week's End Date:"+str(enddate))
      for canal in canal_list:
        regionid = canal[0]
        #logger.info(regionid)
        regionn = canal[0]
        regionid = regionid.replace(" ", "-")
        if os.path.exists(rf"{save_data_loc}/landsat/sebal/" + wktime + r"/sebal_eto_" + regionid + ".zip") is False:
          try:
            table = irrigation_cmd_area.filter(ee.Filter.equals(feature_name, regionn))
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
              logger.info("Planting Date For:"+str(canal_settings[canal[0]])+' is :'+ planting_YY_MM_DD)
              logger.info("Crop Type For:"+str(canal_settings[canal[0]])+' is :'+ crop_type)
            else:
              planting_YY_MM_DD = canal_settings['DEFAULT']['planting_date']
              crop_type = canal_settings['DEFAULT']['crop_type']
              soil_coef = canal_settings['DEFAULT']['soil_coef']
            
            planting_date = datetime.datetime.strptime(planting_YY_MM_DD, '%Y-%m-%d').date()
            num_days = abs((my_datetime - planting_date).days)
            # print(f'Planting Date for command area {canal[0]}:',planting_date)
            # print('Number of days = ',abs((my_datetime - planting_date).days))
            # print((f'Crop grown in command area {canal[0]}:',crop_type))
            logger.info("Planting Date for the command area {area}:{date}".format(area = canal[0],date =planting_date ) )
            logger.info("Number of days = %s",str(abs((my_datetime - planting_date).days)))
            logger.info('Crop grown in command area {canalname}: {crop_type}'.format(canalname = canal[0],crop_type=crop_type))

            crop_config_path = '../Config_files/Crop_Config.ini'
            coefficients = read_crop_coefficients(crop_config_path, crop_type)
            growth_kc = get_growth_kc(coefficients, num_days)
            # print(f"Growth Kc for {crop_type} at {num_days} days: {growth_kc}")  
            logger.info("Growth Kc for {crop_type} at {num_days} days: {growth_kc}".format(crop_type=crop_type,num_days=num_days,growth_kc=growth_kc))
      
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
            seaLevelPressure = 101325 ##// Standard sea level pressure in Pascals
            Rd = 287.05 ##// Gas constant for dry air in J/(kg*K)
            g = 9.80665## // Gravitational acceleration in m/s^2
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
            ##logger.info(centroid)
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
            ##logger.info(tmp1.getInfo())
            tmp2=ee.Number(-1).multiply(ee.Number(cenLAT*math.pi/180).tan())
            ##logger.info(cenLAT)
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
            ##logger.info(Rn)
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
            ##logger.info('tempHot',tempHot)
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
            ##logger.info(tmp)
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
            if glcc_mask:
                ETref_24 = ETref_24.updateMask(glcc_crop)  
            else:
                ETref_24 = ETref_24.clip(ROI) 
            # #logger.info(ETref_24.getInfo())
            penmanET = ee.Image.constant(soil_coef).multiply(ETref_24.multiply(ee.Image.constant(growth_kc))).multiply(ee.Image.constant(dayVal))
            # #logger.info(penmanET_24.getInfo())
            medET_ref= penmanET.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            # #logger.info(medET_ref.getInfo())
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
            # logger.info("Download URL:" + url)
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
            if glcc_mask:
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
            ##logger.info(ETA_24)
            total_etc = totET.get('eta').getInfo()
            params={}
            
            params = {'name': "sebal_eto_" + regionid, 'filePerBand': "false", "scale": ETA_mask.projection().nominalScale(), 'region': ROI.geometry()}
            url = ETA_mask.getDownloadURL(params)
            # #logger.info("Download URL:" + url)
            
            urllib.request.urlretrieve(url, rf"{save_data_loc}/landsat/sebal/" + wktime + r"/sebal_eto_" + regionid + ".zip")
            
            avgET= ETA_mask.reduceRegion(
              reducer= ee.Reducer.mean(),
              geometry= ROI,
              scale= selscale,
              maxPixels= 1e11
              )
            ##logger.info(ETA_24)
            avg_etc = avgET.get('eta').getInfo()
            overirri = ETA_mask.subtract(penmanET)
            if glcc_mask:
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
            logger.info('ET Values For The Processed Command Area - Format: Command Area, Average Penman ET, Average SEBAL ET, Average Overirrigation (SEBAL - Penman)')
            logger.info((canal[0] + ',' + str(median_etref) + "," + str(avg_etc)+ ',' + str(avg_irri)))
            penman_mean_values.append((regionid, median_etref))
            sebal_mean_values.append((regionid, avg_etc))
            irr_mean_values.append((regionid, avg_irri))
            
          except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            # print(f"Error found: {e}")
            logger.error("Week: " + str(wktime) + ", Region: " + str(regionid) + " Showed Error During Processing. Error Line:" + str(exc_tb.tb_lineno)+" See error below.\n")
            logger.error(traceback.format_exc())  # logger full traceback
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
  for canal in tqdm(canal_list, desc="Processing CMD Area TIFF to Readable Format", unit = ' canal'):
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