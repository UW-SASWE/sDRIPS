import os
import argparse
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from shapely import wkt
import numpy as np
from datetime import datetime

# Function to lighten colormap
def lighten_colormap(cmap, lightness_factor=0.5):
    cmap_colors = cmap(np.linspace(0, 1, cmap.N))
    white = np.ones((cmap.N, 4))
    lightened_colors = (1 - lightness_factor) * cmap_colors + lightness_factor * white
    return ListedColormap(lightened_colors)

def generate_water_requirement_plots(input_folder):
    # Set CRS
    btm_proj4 = '+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=-2000000 +ellps=evrst30 +units=m +no_defs +type=crs'
    
    # Define the Plots directory in the output folder
    plots_dir = os.path.join(input_folder, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Find all CSV files in the Distribution_Files folder
    distribution_folder = os.path.join(input_folder, 'distribution_files')
    for file in os.listdir(distribution_folder):
        if file.endswith('.csv'):
            file_path = os.path.join(distribution_folder, file)
            
            # Extract the date from the filename (assuming the format is like 'Distribution_File_Stats_2023_04_07.csv')
            try:
                date_str = file.split('_')[-3:]  # Extract "2023_04_07" part
                date_str = "_".join(date_str).replace('.csv', '')  # Join back and remove ".csv"
                date_obj = datetime.strptime(date_str, "%Y_%m_%d")
                date_formatted = date_obj.strftime("%d %B %Y")
            except ValueError:
                print(f"Filename {file} does not match the expected date format.")
                continue

            # Load data
            data = pd.read_csv(file_path)
            data['geometry'] = data['geometry'].apply(wkt.loads)
            gdf = gpd.GeoDataFrame(data, geometry='geometry')

            # Set and convert CRS
            gdf.set_crs(btm_proj4, inplace=True)
            gdf = gdf.to_crs("EPSG:4326")

            # Convert 'net_water_req (m3)' to numeric and take absolute values
            gdf['net_water_req (m3)'] = pd.to_numeric(gdf['net_water_req (m3)'], errors='coerce').abs()

            # Convert cubic meters to acre-feet
            gdf['net_water_req (acre-feet)'] = gdf['net_water_req (m3)'] / 1233.48

            # Lighten the viridis colormap
            lightened_viridis = lighten_colormap(plt.get_cmap('Oranges'), lightness_factor=0)

            # Create a normalized color scale
            norm = mcolors.Normalize(vmin=gdf['net_water_req (acre-feet)'].min(), vmax=gdf['net_water_req (acre-feet)'].max())

            # Plot the choropleth map
            fig, ax = plt.subplots(1, 1, figsize=(18, 12))
            gdf.boundary.plot(ax=ax, linewidth=1, color='black')  # Plot boundaries for better visualization
            gdf.plot(column='net_water_req (acre-feet)', cmap=lightened_viridis, linewidth=0.8, ax=ax, edgecolor='0.8', norm=norm)

            # Add a color bar
            sm = plt.cm.ScalarMappable(cmap=lightened_viridis, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax)
            cbar.set_label('Net Water Requirement (acre-feet)', fontsize=17)
            cbar.ax.tick_params(axis='both', labelsize=15)

            # Add title and labels
            ax.set_title(f'Net Water Requirement for {date_formatted}', fontsize=20)
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
            ax.tick_params(axis='x', colors='black', labelsize=15)
            ax.tick_params(axis='y', colors='black', labelsize=15)
            ax.set_xlabel("Longitude", fontsize=17, color='black')
            ax.set_ylabel("Latitude", fontsize=17, color='black')

            # Save the plot
            plot_filename = f"Net_Water_Requirement_{date_str}.png"
            plot_path = os.path.join(plots_dir, plot_filename)
            plt.savefig(plot_path, dpi=900, bbox_inches='tight')
            plt.close(fig)
            print(f"Plot saved for {date_formatted} at {plot_path}")

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate net water requirement plots from distribution files.")
    parser.add_argument('-i', '--input_folder', type = str, required=True, help="Path to the input folder containing Distribution_Files.")
    args = parser.parse_args()

    # Run the function with provided arguments
    generate_water_requirement_plots(args.input_folder)
