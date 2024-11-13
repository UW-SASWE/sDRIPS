import os
import argparse
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import numpy as np
from datetime import datetime
from shapely import wkt
import configparser
from matplotlib.colors import Normalize
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# Read configuration
script_config = configparser.ConfigParser()
script_config.read('..\Config_files\Script_Config.ini')
save_data_loc = script_config.get('Save_Data_Location', 'save_data_loc')

# Custom normalization class
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1.]
        return np.ma.masked_array(np.interp(value, x, y))

    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.interp(value, x, y)

def generate_deficit_surplus_plots(input_folder = save_data_loc):
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
                date_str = file.split('_')[-3:]  
                date_str = "_".join(date_str).replace('.csv', '') 
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

            # Convert 'Deficit between requirement and teesta' to acre-feet
            gdf['Deficit between requirement and supply (acre-feet)'] = gdf['Deficit between requirement and supply'] / 1233.48

            # Define column for plotting
            column_name = 'Deficit between requirement and supply (acre-feet)'
            min_value = gdf[column_name].min()
            max_value = gdf[column_name].max()

            # Apply custom normalization
            norm = MidpointNormalize(vmin=min_value, vcenter=0, vmax=max_value)

            # Create the plot
            fig, ax = plt.subplots(1, 1, figsize=(18, 12))
            gdf.boundary.plot(ax=ax, linewidth=1, color='black')
            gdf.plot(column=column_name, ax=ax, legend=False, cmap='RdYlBu', norm=norm, linewidth=0.8, edgecolor='black')

            # Add title and labels
            ax.set_title(f'Deficit(-)/Surplus(+) For {date_formatted}', fontsize=24)
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
            ax.tick_params(axis='x', colors='black', labelsize=17)
            ax.tick_params(axis='y', colors='black', labelsize=17)
            ax.set_xlabel("Longitude", fontsize=19, color='black')
            ax.set_ylabel("Latitude", fontsize=19, color='black')

            # Create and customize colorbar
            cmap = matplotlib.cm.get_cmap('RdYlBu')
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, orientation='vertical')
            cbar.set_label('Deficit(-)/Surplus(+) in acre-feet', fontsize=22)

            # Create rounded ticks for the colorbar
            below_zero = list(np.arange(min_value, 0, 100))
            above_zero = list(np.arange(0, max_value, 100))
            ticks = below_zero + above_zero
            rounded_ticks = [round(tick, -len(str(int(abs(tick)))) + 1) for tick in ticks]
            cbar.set_ticks(ticks)
            cbar.set_ticklabels([f"{int(tick):}" for tick in rounded_ticks])
            cbar.ax.tick_params(axis='both', labelsize=19)

            # Save the plot
            plot_filename = f"Deficit_Surplus_{date_str}.png"
            plot_path = os.path.join(plots_dir, plot_filename)
            plt.tight_layout()
            plt.savefig(plot_path, dpi=900, bbox_inches='tight')
            plt.close(fig)

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate deficit/surplus plots from the distribution files.")
    parser.add_argument('-i', '--input_folder', type=str, required=True, help="Path to the input folder containing distribution_files. By default, it will be the run folder a.k.a save_data_loc from the Script_Config.ini.")
    args = parser.parse_args()

    # Run the function with provided arguments
    generate_deficit_surplus_plots(args.input_folder)
