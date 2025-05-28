import geopandas as gpd
import configparser
import argparse
from pathlib import Path

# Function to read shapefile and get canal IDs
def get_canal_ids(shapefile_path):
    gdf = gpd.read_file(shapefile_path)
    return gdf['Canal_ID'].unique()

# Function to create INI file
def create_ini_file(canal_ids, ini_file_path, default_planting_date, default_crop_type,defualt_soil_coef,default_distribution_unif):
    config = configparser.ConfigParser()

    config['DEFAULT'] = {
        'use_default': 'True',
        'default_planting_date': default_planting_date,
        'default_crop_type': default_crop_type,
        'defualt_soil_coef': defualt_soil_coef,
        'default_distribution_unif': default_distribution_unif,
    }

    for canal_id in canal_ids:
        config[canal_id] = {
            'use_default': 'False',
            'planting_date': default_planting_date,  # Replace with specific date
            'crop_type': default_crop_type,  # Replace with specific crop type
            'soil_coef': defualt_soil_coef,
            'distribution_unif': default_distribution_unif,
        }

    with open(ini_file_path, 'w') as configfile:
        config.write(configfile)

def main():
    parser = argparse.ArgumentParser(description="Generate canal config file using irrigation canals shapefile.")
    parser.add_argument("-s","--shp_path", help="Path to the shapefile, should be in double quotes for windows, and single quotes for Mac/Linux.")
    parser.add_argument("-o","--output", help="Path where the canal_config file will be created,should be in double quotes for windows, and single quotes for Mac/Linux.")
    parser.add_argument("-d","--default_planting_date", default="2023-04-01", help="Default planting date in YYYY-MM-DD format.")
    parser.add_argument("-cc","--default_crop_type", default="Rice", help="Default crop type with first letter capital. eg - Rice, Wheat, Corn etc.")
    parser.add_argument("-sc","--defualt_soil_coef", default=0.5, help="Default soil coefficient. Should not be in quotes eg - 0.3,0.5 etc")
    parser.add_argument("-du","--default_distribution_unif", default=1, help="Default distribution uniformity (efficiency). Should not be in quotes eg - 0.3,0.5 etc. Gives information on how well the water will be distributed, depends on method of irrigation such as drip, sprinkler etc." )
    args = parser.parse_args()

    canal_ids = get_canal_ids(args.shp_path)
    create_ini_file(canal_ids, args.output, args.default_planting_date, args.default_crop_type,args.defualt_soil_coef, args.default_distribution_unif)

    print(f"Canal Config file created at {args.output}")

if __name__ == "__main__":
    main()
