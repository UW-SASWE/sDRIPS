# Importing required libraries 
import ee 
import logging
from tqdm import tqdm
import geemap
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
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import warnings
warnings.filterwarnings("ignore")

from sdrips.utils.ee_initialize import initialize_earth_engine
from sdrips.utils.utils import (
   load_yaml_config,
    read_cmd_area_settings,
    read_crop_coefficients,
    get_growth_kc
)
from src.sdrips.utils.logging_utils import (
    worker_logger_setup,
    worker_init
)

# Initialise Earth Engine
initialize_earth_engine()

#########*******Module (1) USER INPUT (ROI and Intereset Dates and Planting Date)********########
# Initialize ruamel.yaml YAML object
yaml = YAML()
yaml.preserve_quotes = True  # Optional: preserves quotes around strings

script_config = load_yaml_config('config_files/test_script_config.yaml')

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
cmd_config_path = script_config['Cmd_Area_Config']['path']

# === Conditional Flags ===
clear_directory_condition = script_config['Clean_Directory']['clear_directory_condition']
ET_estimation = script_config['Run_ET_Estimation']['et_estimation']
glcc_mask = script_config['GLCC_Mask'].get('glcc_mask', False)
precipitation_condition = script_config['Precipitation_Config']['consider_preciptation']
weather_condition = script_config['Weather_Config']['consider_forecasted_weather']
percolation_condition = script_config['Percolation_Config']['consider_percolation']
region_stats_condition = script_config['Region_stats']['estimate_region_stats']
tiff_to_png = script_config['Tiff_to_PNGs'].get('tiff_to_pngs', False)

# === multiprocessing ===
max_workers = script_config.get("Multiprocessing", {}).get("max_workers")
worker_count = max_workers if max_workers is not None else multiprocessing.cpu_count() - 1

# === Accessing Secrets Configuration Files ===
secrets_file_path = script_config['Secrets_Path']['path']
secrets = load_yaml_config(rf'{secrets_file_path}')
imerg_username = secrets['IMERG_Account']['username']
imerg_password = secrets['IMERG_Account']['password']


# === SETUP LOGGING ===
# os.makedirs(f'{save_data_loc}',exist_ok = True)
# os.makedirs(f'{save_data_loc}/logs/',exist_ok = True)
# dt_fmt = '%Y%m%d_%H%M%S'
# logging.basicConfig(filename=f'{save_data_loc}/logs/{datetime.datetime.today().strftime(dt_fmt)}.log', format='%(levelname)s - %(asctime)s :%(message)s', level=logging.DEBUG)
# logging.info("")
# logging.info("Location of the Data Folder: "+save_data_loc)
# logging.info("")
# logging.critical('Start Date:'+str(start_date) )
# logging.critical('Run Week:'+ str(run_week) )

if gee_asset_section.get('id'):
    irrigation_cmd_area = ee.FeatureCollection(gee_asset_id) 
    # print(f'GEE Asset ID: {gee_asset_id}')
else:
    irrigation_cmd_area= geemap.shp_to_ee(gee_asset_id)
 
cmd_area_list = irrigation_cmd_area.reduceColumns(ee.Reducer.toList(1), [feature_name]).get('list').getInfo()

def get_cmd_area_list():
  return cmd_area_list
    
glcc = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")

# """ Reading the canal config file to parse the planting date, 
# crop type and soil coefficient for each command region """

config = load_yaml_config(cmd_config_path)

# Extract defaults and check use_default
defaults = config.get('DEFAULT', {})
use_default = defaults.get('use_default', False)

# Build cmd_area settings dictionary
if use_default:
    cmd_area_settings = {'DEFAULT': {
        'planting_date': defaults['planting_date'],
        'crop_type': defaults['crop_type'],
        'soil_coefficient': defaults['soil_coef']
    }}
else:
    cmd_area_settings = {}
    for cmd_area_name, cmd_area_data in config.items():
        if cmd_area_name == 'DEFAULT':
            continue
        cmd_area_settings[cmd_area_name] = read_cmd_area_settings(cmd_area_data, defaults)


penman_mean_values = []
sebal_mean_values = []
irr_mean_values = []

#########*******Module (2) Estimates Penman-Monteith and SEBAL Evapotranspiration********########
def process_single_cmd_area(args):
  logger = logging.getLogger()
  cmd_area, wktime, save_data_loc = args
  try:
    regionid = cmd_area[0].replace(" ", "-")
    regionn = cmd_area[0]
    output_zip = rf"{save_data_loc}/landsat/sebal/{wktime}/sebal_eto_{regionid}.zip"

    if os.path.exists(output_zip):
      return None  # Skip if already processed

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
    if list(cmd_area_settings.keys())[0] != 'DEFAULT':
      planting_YY_MM_DD = cmd_area_settings[cmd_area[0]]['planting_date']
      crop_type = cmd_area_settings[cmd_area[0]]['crop_type']
      soil_coef = cmd_area_settings[cmd_area[0]]['soil_coef']
      logger.info("Planting Date For:"+str(cmd_area_settings[cmd_area[0]])+' is :'+ planting_YY_MM_DD)
      logger.info("Crop Type For:"+str(cmd_area_settings[cmd_area[0]])+' is :'+ crop_type)
    else:
      planting_YY_MM_DD = cmd_area_settings['DEFAULT']['planting_date']
      # print('here is the planting date', planting_YY_MM_DD)
      crop_type = cmd_area_settings['DEFAULT']['crop_type']
      # print('here is the crop_type', crop_type)
      # print('here is the type of data of crop_type:', type(crop_type))
      # print('All canal settings for default:', cmd_area_settings['DEFAULT'])
      soil_coef = cmd_area_settings['DEFAULT']['soil_coefficient']
    
    planting_date = datetime.datetime.strptime(planting_YY_MM_DD, '%Y-%m-%d').date()
    num_days = abs((my_datetime - planting_date).days)
    # print(f'Planting Date for command area {canal[0]}:',planting_date)
    # print('Number of days = ',abs((my_datetime - planting_date).days))
    # print((f'Crop grown in command area {canal[0]}:',crop_type))
    logger.info("Planting Date for the command area {area}:{date}".format(area = cmd_area[0],date =planting_date ) )
    logger.info("Number of days = %s",str(abs((my_datetime - planting_date).days)))
    logger.info('Crop grown in command area {canalname}: {crop_type}'.format(canalname = cmd_area[0],crop_type=crop_type))

    crop_config_path = 'config_files/crop_config.yaml'
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
    os.makedirs(rf"{save_data_loc}/landsat/penman/" + wktime, exist_ok=True)
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
    os.makedirs(rf"{save_data_loc}/landsat/sebal/" + wktime, exist_ok=True)
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
    os.makedirs(rf"{save_data_loc}/landsat/irrigation/" + wktime, exist_ok=True)
    urllib.request.urlretrieve(url, rf"{save_data_loc}/landsat/irrigation/" + wktime + r"/irrigation_" + regionid + ".zip")
    logger.info('ET Values For The Processed Command Area - Format: Command Area, Average Penman ET, Average SEBAL ET, Average Overirrigation (SEBAL - Penman)')
    logger.info((cmd_area[0] + ',' + str(median_etref) + "," + str(avg_etc)+ ',' + str(avg_irri)))
    penman_mean_values.append((regionid, median_etref))
    sebal_mean_values.append((regionid, avg_etc))
    irr_mean_values.append((regionid, avg_irri))
    return (cmd_area[0], penman_mean_values, sebal_mean_values, irr_mean_values)
            
  except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    line_no = exc_tb.tb_lineno

    logger.error(
        f"Exception occurred in {fname} on line {line_no}: {e}"
    )
    logger.error(traceback.format_exc())


def process_cmd_area_parallel(main_logger, log_queue, max_workers=worker_count):
  results = []
  total_iterations = len(run_week) * len(cmd_area_list)
  with tqdm(total=total_iterations, desc="Estimating ET with models", unit=" Command Area") as pbar:
    for wktime in run_week:
      stats_file = rf"{save_data_loc}/landsat/stats_{wktime}.txt"
      os.makedirs(os.path.dirname(stats_file), exist_ok=True)

      with open(stats_file, 'w') as txt:
          txt.write("Region,Penman_ET,Sebal_ET,Irrigation\n")

      main_logger.critical(f"Running Week: {wktime}")
      main_logger.critical(f"Start Date: {start_date}")
      # logger.critical(f"End Date: {(datetime.datetime.strptime(start_date, '%Y-%m-%d') + datetime.timedelta(days=8+2)).strftime('%Y-%m-%d')}")

      args_list = [(cmd_area, wktime, save_data_loc) for cmd_area in cmd_area_list]

      with ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=worker_logger_setup,
                initargs=(log_queue,)
            ) as executor:
        futures = [executor.submit(process_single_cmd_area, args) for args in args_list]
        for future in as_completed(futures):
          try:
              result = future.result()
              if result and not isinstance(result, str):
                  region, penman, sebal, irr = result
                  with open(stats_file, 'a') as txt:
                      txt.write(f"{region},{penman},{sebal},{irr}\n")
                  results.append(result)
              elif isinstance(result, str):
                  main_logger.error(result)
          except Exception:
              main_logger.exception("Exception in worker process")

          pbar.update(1)
    # Group by variable
    penman_mean_values = [r[1] for r in results]
    sebal_mean_values = [r[2] for r in results]
    irr_mean_values = [r[3] for r in results]

    return penman_mean_values, sebal_mean_values, irr_mean_values