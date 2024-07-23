import numpy as np 
import pandas as pd 
import geopandas as gpd 
from shapely.geometry import LineString
import configparser
import argparse
import os

def iterative_buffer_algo(geometry, target_area, initial_guess, tolerance=0.005, max_iterations=1000):
    """
    Adjusts the buffer distance until the area of the buffered feature 
    is close to the target area.
    """
    left, right = 0, initial_guess  
    iteration = 0

    while iteration < max_iterations:
        mid = (left + right) / 2
        area = geometry.buffer(mid).area

        if np.isclose(area, target_area, rtol=tolerance) == True:
            # print(np.isclose(area, target_area, rtol=tolerance))
            return mid
        elif area < target_area:
            left = mid
        else:
            right = mid
            
        
        iteration += 1
    if iteration == max_iterations:
        print('Crossed max iter, mid:',mid)

    return mid

def apply_buffer(shapefile_path,output_path):
    gdf = gpd.read_file(shapefile_path)
    # Apply the function to each row of the dataframe
    gdf['buffer_distance'] = gdf.apply(lambda row: iterative_buffer_algo(geometry = row['geometry'], target_area = row['Command Area (m2)'], initial_guess = gdf['Command Area (m2)'].max()), axis=1)

    # Create the buffered geometries
    gdf['buffered_geometry'] = gdf.apply(lambda row: row['geometry'].buffer(row['buffer_distance']), axis=1)

    # Drop the initial irrigation canal geometries, and save the geometry of the new buffer polygons
    gdf = gdf.drop(columns=['geometry'])
    gdf=gdf.rename(columns = {'buffered_geometry':'geometry'})
    gdf.to_file(output_path+'/Generated_Command_Areas.shp')

def main():
    parser = argparse.ArgumentParser(description="Generates command areas around the network of canals using an iterative approach. By default this scrip looks for the Command Area (m2) attribute name. Users can change it to desired name in this script.")
    parser.add_argument("-s","--shp_path", help="Path to the irrigation canal network shapefile, should be in double quotes for windows, and single quotes for Mac/Linux.")
    parser.add_argument("-o","--output", help="Folder path where the command area shapefile (Generated_Command_Areas.shp) will be created,should be in double quotes for windows, and single quotes for Mac/Linux.")    
    args = parser.parse_args()
    apply_buffer(args.shp_path, args.output)

    print(f"Command Area Shapefile created at {args.output}")

if __name__ == "__main__":
    main()