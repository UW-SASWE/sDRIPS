import numpy as np 
import pandas as pd 
import geopandas as gpd 
from shapely.geometry import LineString
import configparser
import argparse
import os
from pathlib import Path
from typing import Any, Union

def iterative_buffer_algo(geometry: Any,
    target_area: float,
    initial_guess: float,
    tolerance: float = 0.005,
    max_iterations: int = 1000) -> float:
    """Find a buffer distance so buffered area is close to the target.

    Uses a bisection search between 0 and ``initial_guess`` until the area of
    ``geometry.buffer(distance)`` is within the given relative tolerance of
    ``target_area`` or the maximum number of iterations is reached.

    Parameters
    ----------
    geometry : shapely geometry
        Geometry to buffer. Typically a LineString/MultiLineString.
    target_area : float
        Desired area for the buffered geometry (in units of the CRS).
    initial_guess : float
        Initial upper bound for the buffer distance (same units as CRS).
    tolerance : float, optional
        Relative tolerance for area comparison (default is 0.005 → 0.5%).
    max_iterations : int, optional
        Maximum number of bisection iterations (default is 1000).

    Returns
    -------
    float
        The buffer distance found by the search.
    """
    left, right = 0, float(initial_guess)  
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

def create_cmd_area(shape_path: Union[str, Path], 
                    column_name: str ,
                    output_path: Union[str, Path]) -> None: 
    """Apply iterative buffering to each feature and write polygon output.

    This function reads the input dataset, computes a per-row buffer distance
    so that each feature's buffered area approaches the value in the
    ``column_name``, replaces the geometry with the buffered
    polygon, and writes the result to the
    given ``output_path``.
    Note: Use this function only when command area boundaries are not available, as the output
    command areas from this function is an estimation, and may not be accurate. Hence, it should be used with caution.

    Parameters
    ----------
    shape_path : str or pathlib.Path
        Path to the irrigation canal network dataset readable by GeoPandas.
    column_name : str
        Name of the column containing the target area for each feature. Should be in m2.
    output_path : str or pathlib.Path
        Path to the output shapefile (including .shp extension) created by this function.

    Returns
    -------
    pathlib.Path
        Path to the created command area shapefile.
    """
    gdf = gpd.read_file(shape_path)
    if gdf.crs is None:
        raise ValueError("Input data has no CRS. Please set or define a projected CRS before buffering.")
    if column_name not in gdf.columns:
        raise ValueError(f"Column '{column_name}' not found in the shapefile. \n Available columns: {gdf.columns.tolist()}")
    crs = gdf.crs
    # Apply the function to each row of the dataframe
    gdf['buffer_distance'] = gdf.apply(lambda row: iterative_buffer_algo(geometry = row['geometry'], target_area = row[column_name], initial_guess = gdf[column_name].max()), axis=1)

    gdf['buffered_geometry'] = gdf.apply(lambda row: row['geometry'].buffer(row['buffer_distance']), axis=1)
    gdf = gdf.set_geometry("buffered_geometry")
    gdf = gdf.drop(columns=["geometry"])
    gdf = gdf.set_crs(crs, allow_override=True)
    gdf = gdf.rename_geometry("geometry")
    out = Path(os.path.join(output_path))
    os.makedirs(os.path.dirname(output_path), exist_ok=True) 
    gdf.to_file(out)
    return out

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate command-area polygons around a canal/line network using a per-feature "
            "iterative buffer so that each polygon area matches a target area."
        )
    )
    parser.add_argument(
        "-s", "--shp_path", required=True,
        help="Path to the input canal network (e.g., shapefile - shp format)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path for the output shapefile."
    )
    parser.add_argument(
        "-c", "--column_name", default="Command Area (m2)",
        help="Name of the area column in the input (default: 'Command Area (m2)')."
    )
    args = parser.parse_args()
    create_cmd_area(shape_path = args.shp_path, column_name = args.column_name, output_path = args.output)

    print(f"Command Area Shapefile created at {args.output}")

if __name__ == "__main__":
    main()