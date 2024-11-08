import ee 
import logging
import colorlog
from tqdm import tqdm
# import geemap
import os, datetime, requests, zipfile, time, math
import urllib.request
import datetime,math
import numpy as np
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
import warnings
warnings.filterwarnings("ignore")

## SETUP LOGGING
save_data_loc = "../Data_2023_March_Check2/"
os.makedirs(save_data_loc,exist_ok = True)

dt_fmt = '%Y%m%d_%H%M%S'
logging.basicConfig(filename=f'../logs/{datetime.datetime.today().strftime(dt_fmt)}.log', format='%(levelname)s - %(asctime)s :%(message)s', level=logging.DEBUG)

logging.info("")
logging.info("Location of the Data Folder: "+save_data_loc)
logging.info("")
# def configure_logging():
#     # Create a logger
#     logger = logging.getLogger()
#     logger.setLevel(logging.INFO)

#     # Create a colorlog handler
#     handler = logging.FileHandler(f'../logs/{datetime.datetime.today().strftime(dt_fmt)}.html')
#     handler.setFormatter(colorlog.ColoredFormatter(
#         '%(log_color)s%(levelname)-8s%(reset)s %(message)s',
#         log_colors={
#             'DEBUG': 'cyan',
#             'INFO': 'green',
#             'WARNING': 'yellow',
#             'ERROR': 'red',
#             'CRITICAL': 'bold_green',
#         },
#         secondary_log_colors={},
#         style='%'
#     ))

#     # Add the handler to the logger
#     logger.addHandler(handler)

# Initialise Earth Engine
# ee.Authenticate() # If you are running the earth engine first time in your machine, you need to run the ee.Authenticate() command first.
ee.Initialize()


#########*******Module (1) USER INPUT (ROI and Intereset Dates and Planting Date)********########
script_config = configparser.ConfigParser()

script_config.read('..\Config_files\Script_Config.ini')

# Accessing the information from 'Irrigation_Canals_shapefile' section
irrigation_canals_path = script_config.get('Irrigation_Canals_shapefile', 'path')

# Accessing the information from 'Date_Running' section
start_date = script_config.get('Date_Running', 'start_date')
print(f'Start Date:{start_date}')
default_run_week = script_config.getboolean('Date_Running', 'default_run_week')
if default_run_week:
    run_week = ['lastweek', 'currentweek']
else:
    run_week = eval(script_config.get('Date_Running', 'run_week'))
logging.critical('Start Date:'+str(start_date) )
logging.critical('Run Week:'+ str(run_week) )

# Accessing the information from 'OSGEO_Path' section
osgeo_path = script_config.get('OSGEO_Path', 'path')

# irrigation_canals= ee.FeatureCollection('users/skhan7/PhD/Chapter2/Mapped_Buffer_Canals')
irrigation_canals= ee.FeatureCollection('users/skhan7/PhD/Chapter2/BWDB_Teesta_Command_Area')
canal_list = irrigation_canals.reduceColumns(ee.Reducer.toList(1), ['Canal_ID']).get('list').getInfo()
glcc = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")


# """ Reading the canal config file to parse the planting date, 
# crop type and soil coefficient for each command region """

config = configparser.ConfigParser()

# Read the INI file
config.read('../Config_files/New_Canal_Config.ini')

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



def download_tif_from_ee(image, path, name):
    """Generate download URL from EE and download the TIFF using urllib."""
    
    params = {
        'name':f'{name}',
        'filePerBand': "false",
        'scale': 25000,
        'crs': 'EPSG:4326',
        'fileFormat': 'GeoTIFF',
        'region': ee.Geometry.Rectangle([88.5, 25.5, 89.5, 26.5])
    }
    url = image.getDownloadURL(params)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    urllib.request.urlretrieve(url, path)
    logging.info(f'Downloaded Successfully From GEE: {path}')
    # except Exception as e:
    #     logging.error(f'Failed to download {path}: {str(e)}')

# def gfsdata():
#     startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
#     endDate = startDate + datetime.timedelta(days=1)
#     logging.info(f'GFS Start date: {startDate}')
#     date_difference = (datetime.datetime.today() - startDate).days
#     logging.info(f'Date difference: {date_difference}')

#     if date_difference > 10:
#         logging.info("Switching to Google Earth Engine for data download due to date difference.")
#         dataset = ee.ImageCollection('NOAA/GFS0P25')
#         dataset = dataset.filterDate(startDate.strftime('%Y-%m-%d'), endDate.strftime('%Y-%m-%d')).filterBounds(ee.Geometry.Rectangle([88.5, 25.5, 89.5, 26.5]))

#         # APCP download for next week, 168-hour forecast
#         apcp_image = dataset.select('total_precipitation_surface').filter(ee.Filter.eq('forecast_hours', 168)).first()
#         apcp_path = os.path.join(save_data_loc, "precip", f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.zip")
#         download_tif_from_ee(apcp_image, apcp_path,name =f"precip_gfs_{startDate.strftime('%Y%m%d')}_nextweek" )
#         extract_dir = os.path.join(save_data_loc, "precip")
#         with zipfile.ZipFile(apcp_path, "r") as zip_input:
#             for member in zip_input.namelist():
#                 extracted_path = zip_input.extract(member, extract_dir)
#                 new_path = os.path.join(extract_dir, f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.tif")
#                 os.rename(extracted_path, new_path)

#         # Remove the zip file
#         os.remove(apcp_path)
#         # Variables for weeks and hours
#         weeks = ['currentweek', 'nextweek']
#         params_ee = ['u_component_of_wind_10m_above_ground', 'v_component_of_wind_10m_above_ground', 'temperature_2m_above_ground']
#         total_iterations = len(weeks) * len(params_ee+1) * 14  # 14 comes from the range(12, 169, 12)
#         with tqdm(total=total_iterations, desc="Downloading and Processing GFS Data", unit="file") as pbar:
#             for week in weeks:
#                 for param in params_ee:
#                     week_date = startDate + datetime.timedelta(days=-7) if week == 'currentweek' else startDate
#                     for hours in range(12, 169, 12):
#                         image = dataset.select(param).filter(ee.Filter.eq('forecast_hours', hours)).first()
#                         folder = 'ugrd' if param == 'u_component_of_wind_10m_above_ground' else ('vgrd' if param == 'v_component_of_wind_10m_above_ground' else 'temp')
#                         path = os.path.join(save_data_loc, folder, f"{folder}.{week_date.strftime('%Y%m%d')}.{week}.{str(hours).zfill(3)}.zip")
#                         download_tif_from_ee(image, path,name = f"{folder}_{week_date.strftime('%Y%m%d')}_{week}_{str(hours).zfill(3)}")
#                         if folder != 'temp':
#                             extract_dir = os.path.join(save_data_loc, folder)
#                             with zipfile.ZipFile(path, "r") as zip_input:
#                                 for member in zip_input.namelist():
#                                     extracted_path = zip_input.extract(member, extract_dir)
#                                     gtifpath = os.path.join(extract_dir, f"{folder}.gfs.{startDate.strftime('%Y%m%d')}.{week}{str(hours).zfill(3)}.tif")
#                                     os.rename(extracted_path, gtifpath)
#                             ascpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.asc'
#                             finalpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.tif'
#                             x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
#                             logging.info(x.read())
#                             x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
#                             logging.info(x.read())
#                             logging.info((finalpath))
#                             pbar.update(1)
                        

#                         else:
#                             folders=['tmax','tmin']
#                             for folder in folders:
#                                 extract_dir = os.path.join(save_data_loc, folder)
#                                 with zipfile.ZipFile(path, "r") as zip_input:
#                                     for member in zip_input.namelist():
#                                         extracted_path = zip_input.extract(member, extract_dir)
#                                         gtifpath = os.path.join(extract_dir, f"{folder}.gfs.{startDate.strftime('%Y%m%d')}.{week}{str(hours).zfill(3)}.tif")
#                                         os.rename(extracted_path, gtifpath)
#                                 ascpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.asc'
#                                 finalpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.tif'
#                                 x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
#                                 logging.info(x.read())
#                                 x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
#                                 logging.info(x.read())
#                                 logging.info((finalpath))
#                                 pbar.update(1)
#         logging.info("Using data from NOMADS due to shorter date difference.")






def gfsdata_ee():
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days = 0)
    datestr = startDate.strftime("%Y%m%d")
    endDate = startDate + datetime.timedelta(days=1)
    dataset = ee.ImageCollection('NOAA/GFS0P25')
    dataset = dataset.filterDate(startDate.strftime('%Y-%m-%d'), endDate.strftime('%Y-%m-%d')).filterBounds(ee.Geometry.Rectangle([88.5, 25.5, 89.5, 26.5]))

    # APCP download for next week, 168-hour forecast
    apcp_image = dataset.select('total_precipitation_surface').filter(ee.Filter.eq('forecast_hours', 168)).first()
    apcp_path = os.path.join(save_data_loc, "precip", f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.zip")
    download_tif_from_ee(apcp_image, apcp_path,name =f"precip_gfs_{startDate.strftime('%Y%m%d')}_nextweek" )
    extract_dir = os.path.join(save_data_loc, "precip")
    with zipfile.ZipFile(apcp_path, "r") as zip_input:
        for member in zip_input.namelist():
            extracted_path = zip_input.extract(member, extract_dir)
            new_path = os.path.join(extract_dir, f"precip.gfs.{startDate.strftime('%Y%m%d')}.nextweek.tif")
            os.rename(extracted_path, new_path)
    os.remove(apcp_path)
    weeks = ['currentweek', 'nextweek']
    params_ee = ['u_component_of_wind_10m_above_ground', 'v_component_of_wind_10m_above_ground', 'temperature_2m_above_ground']
    total_iterations = len(weeks) * (len(params_ee)+1) * 14  # 14 comes from the range(12, 169, 12)
    with tqdm(total=total_iterations, desc="Downloading GFS Data From GEE", unit="file") as pbar:
        for week in weeks:
            week_date = startDate + datetime.timedelta(days=-7) if week == 'currentweek' else startDate
            week_date_end = week_date + datetime.timedelta(days=1)
            dataset = ee.ImageCollection('NOAA/GFS0P25')
            dataset = dataset.filterDate(week_date.strftime('%Y-%m-%d'), week_date_end.strftime('%Y-%m-%d')).filterBounds(ee.Geometry.Rectangle([88.5, 25.5, 89.5, 26.5]))
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
                        os.remove(path)
                        ascpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.asc'
                        finalpath = rf"{save_data_loc}" + folder + r'/' + folder + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                        x = os.popen(f'{osgeo_path} gdal_translate -of AAIGrid -b 1 ' + gtifpath + ' ' + ascpath)
                        logging.info(x.read())
                        x = os.popen(f'{osgeo_path} gdal_translate -of GTiff ' + ascpath + ' ' + finalpath)
                        logging.info(x.read())
                        logging.info((finalpath))
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
    path = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_APCP=on&subregion=&leftlon=88.5&rightlon=89.5&toplat=26.5&bottomlat=25.5&dir=%2Fgfs.' + datestr + '%2F00%2Fatmos'
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
    with tqdm(total=total_iterations, desc="Downloading and GFS Data From NOAA Server", unit="file") as pbar:
        for week in weeks:
            for param in params:
                datestr = startDate.strftime("%Y%m%d")
                if week=='currentweek':
                    datestr = (startDate+datetime.timedelta(days=-7)).strftime("%Y%m%d")
                for hours in range(12, 169, 12):
                    paramid = paramids[params.index(param)]
                    serverpath = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_' + paramid + '=on&subregion=&leftlon=88.5&rightlon=89.5&toplat=26.5&bottomlat=25.5&dir=%2Fgfs.' + datestr + '%2F00%2Fatmos'
                    gtifpath = rf"{save_data_loc}" + param + r'/' + param +'.' + datestr + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                    ascpath = rf"{save_data_loc}" + param + r'/' + param + '.' + week + '.' + str(hours).zfill(3) + '.asc'
                    finalpath = rf"{save_data_loc}" + param + r'/' + param + '.' + week + '.' + str(hours).zfill(3) + '.tif'
                    logging.info(f'GFS Link for {param}: {serverpath}')

                    try:
                        # urllib.request.urlretrieve(serverpath, gtifpath)
                        args = f'curl.exe --output {gtifpath} "{serverpath}"'
                        sub.run(args, shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
                    except:
                        serverpath = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.f168&var_' + paramid + '=on&subregion=&leftlon=88.5&rightlon=89.5&toplat=26.5&bottomlat=25.5&dir=%2Fgfs.' + datestr + '%2F00'
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
    with tqdm(total=total_iterations, desc="Processing GFS Data", unit="file") as pbar:
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

# Run the function
gfsdata()