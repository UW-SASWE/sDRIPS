import geopandas as gpd
import networkx as nx
import pandas as pd
import pandas as pd
import os
import logging
from shapely.wkt import loads
from ruamel.yaml import YAML
from pathlib import Path
from typing import List, Dict, Any


from sdrips.utils.utils import load_yaml_config


yaml = YAML()
yaml.preserve_quotes = True


def calculate_canal_cumulative_discharge(config_path: Path) -> pd.DataFrame:
    """
    Computes cumulative water discharge along a canal network using parent-child relationships
    defined in a shapefile and a stats CSV file.

    Parameters:
        config_path (str): Path to main configuration file.

    Returns:
        pd.DataFrame: Merged DataFrame with cumulative discharge and percentage.
    """
    logger = logging.getLogger()
    logger.critical('Started Canal Water Distribution Calculation')
    config = load_yaml_config(config_path)

    canals_shp_file = config['Canal_water_allotment_config']['canals_shp_file']
    logger.info(f"Loading canal shapefile from {canals_shp_file}")
    canals_feature_name = config['Canal_water_allotment_config']['canals_feature_name']
    numeric_id_name = config['Canal_water_allotment_config']['numeric_id_name']
    parent_cols = config['Canal_water_allotment_config']['parents_canals_nomenclature']
    subdivisions = config['Canal_water_allotment_config']['subdivisions_in_canal_network']
    logger.info(f"Number of sub-parents in the canal network: {subdivisions}")
    supply_data_path = config['Canal_water_allotment_config']['supply_data_path']
    save_data_loc = config['Save_Data_Location']['save_data_loc']
    
    gdf = gpd.read_file(canals_shp_file)
    gdf = gdf.to_crs('EPSG:4326')

    cnlnm_id_dict = {
        row.get(canals_feature_name): row.get(numeric_id_name)
        for _, row in gdf.iterrows()
        if pd.notna(row.get(canals_feature_name)) and pd.notna(row.get(numeric_id_name))
    }

    for col in parent_cols:
        new_col = 'node_' + col
        gdf[new_col] = gdf[col].map(cnlnm_id_dict)


    stats_df = pd.read_csv(supply_data_path)
    merged = gdf.merge(stats_df, on=numeric_id_name)

    def update_nodes(row):
        last_valid = row[numeric_id_name]
        for col in ['node_' + pc for pc in parent_cols]:
            if row.get(col) == last_valid:
                row[col] = None
            else:
                last_valid = row.get(col)
        return row

    merged = merged.apply(update_nodes, axis=1)

    merged['coords'] = merged['geometry'].apply(lambda geom: geom.coords[0] if geom is not None else None)

    G = nx.DiGraph()
    for _, row in merged.iterrows():
        if row['coords']:
            G.add_node(row[numeric_id_name], pos=row['coords'])

    for _, row in merged.iterrows():
        for i in range(subdivisions):
            parent_col = 'node_' + parent_cols[i]
            child_col = numeric_id_name if i == 0 else 'node_' + parent_cols[i - 1]
            parent = row.get(parent_col)
            child = row.get(child_col)
            if pd.notna(parent) and pd.notna(child):
                G.add_edge(parent, child)

    water_dict = {}

    def calculate_water(node):
        if node in water_dict:
            return water_dict[node]

        water_value = merged.loc[merged[numeric_id_name] == node, 'net_water_req_m3'].values[0]

        if G.out_degree(node) == 0:
            water_dict[node] = water_value
            return water_value

        total_water = water_value
        for child in G.successors(node):
            total_water += calculate_water(child)

        water_dict[node] = total_water
        return total_water

    root_nodes = [n for n in G.nodes if G.in_degree(n) == 0]
    for root in root_nodes:
        calculate_water(root)

    dict_df = gpd.GeoDataFrame(list(water_dict.items()), columns=[numeric_id_name, 'Cum_canal_water_req'])
    dict_df['Cum_canal_water_req_m3'] = dict_df['Cum_canal_water_req_m3'].abs()
    merged = merged.merge(dict_df, on = numeric_id_name, how = 'outer')
    merged.to_csv(os.path.join(save_data_loc, 'Landsat_Command_Area_Stats.csv'), index=False)
    columns_to_keep = [numeric_id_name, 'geometry', 'Cum_canal_water_req']
    merged = merged[columns_to_keep]
    merged['Cum_canal_water_req_m3_percentage'] = (merged['Cum_canal_water_req_m3']/merged['Cum_canal_water_req_m3'].max())*100
    merged['geometry_str'] = merged['geometry'].apply(lambda geom: geom.wkt)
    grouped = merged.groupby('geometry_str').agg({'Cum_canal_water_req_percentage': 'max'}).reset_index()
    grouped['geometry'] = grouped['geometry_str'].apply(lambda geom: loads(geom))
    grouped = gpd.GeoDataFrame(grouped, geometry='geometry')

    grouped.set_crs(merged.crs)
    grouped.crs = merged.crs
    grouped = grouped.rename(columns={'Cumulative Discharge percentage': 'Max Discharge Percentage'})

    output_path = os.path.join(save_data_loc, 'canal_distribution.csv')
    grouped.to_csv(output_path, index=False)
    logger.info(f"Canal water distribution saved to: {output_path}")
    logger.critical('Finished Canal Water Distribution Calculation')
    return grouped

    