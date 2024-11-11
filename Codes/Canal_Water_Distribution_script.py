import geopandas as gpd
import pandas as pd
import pandas as pd
import numpy as np
import os
import argparse
import configparser

script_config = configparser.ConfigParser()
script_config.read('..\Config_files\Script_Config.ini')
feature_name = script_config.get('Irrigation_Canals_shapefile', 'feature_name')
save_data_loc = script_config.get('Save_Data_Location', 'save_data_loc')
run_date = script_config.get('Date_Running', 'start_date')
supply_path = '../TBP_Data/TBP_CHR_Cleaned_Pivoted.csv'

def get_discharge(date, df):
    # Get the discharge value for the interested date
    discharge = df[df['Date'] == date]['Discharge'].values[0]
    
    # Check if the discharge is NaN
    if np.isnan(discharge):
        # Find the index of the interested date
        idx = df.index[df['Date'] == date].tolist()[0]
        
        # Get the discharge value just above and below
        if idx > 0 and idx < len(df) - 1:
            above_discharge = df.iloc[idx - 1]['Discharge']
            below_discharge = df.iloc[idx + 1]['Discharge']
            
            # If both above and below values are not NaN, take the average
            if not np.isnan(above_discharge) and not np.isnan(below_discharge):
                discharge = (above_discharge + below_discharge) / 2
            # If only above value is not NaN, use the above value
            elif not np.isnan(above_discharge):
                discharge = above_discharge
            # If only below value is not NaN, use the below value
            elif not np.isnan(below_discharge):
                discharge = below_discharge
            # If both are NaN, keep it as NaN
        # Handle edge cases where the date is at the start or end of the series
        elif idx == 0:
            below_discharge = df.iloc[idx + 1]['Discharge']
            if not np.isnan(below_discharge):
                discharge = below_discharge
        elif idx == len(df) - 1:
            above_discharge = df.iloc[idx - 1]['Discharge']
            if not np.isnan(above_discharge):
                discharge = above_discharge
                
    return discharge

def Irr_status(data_run_folder = save_data_loc, interested_date = run_date,  net_water_req_corr = False, shp_file = '../TBP Shape Files/TCAD_Phase1_DataPartA/Modified_Command_Area_Irrigable.csv' ):
    # cmd = gpd.read_file('../TBP Shape Files/TCAD_Phase1_DataPartA/ShapeFiles/commandarea_teesta_ph1_simplified_modified.shp')
    cmd = pd.read_csv(shp_file)
    cmd['Distribution_Factor'] = cmd['AREA']/cmd['AREA'].sum()
    if os.path.isabs(data_run_folder):
        stats_path = f'{data_run_folder}/Landsat_Command_Area_Stats.csv'
    else:
        stats_path = f'{data_run_folder}/Landsat_Command_Area_Stats.csv'
    stats_file = pd.read_csv(stats_path)
    if net_water_req_corr == False:
        try:
            stats_file['net_water_req'] = stats_file['net_water_req'].apply(lambda x: x if x <= 0 else 0)
        except:
            print('Correction wasnt applied')
    stats_merged = stats_file.merge(cmd,on=[feature_name])
    stats_merged = stats_merged.rename(columns = {'net_water_req':'net_water_req (mm)','Currentweek PPT':'Currentweek PPT (mm)','Nextweek PPT':'Nextweek PPT (mm)','7Day Irrigation':'7Day Irrigation (mm)','AREA_x':'AREA','PERIMETER_x':'PERIMETER'})
    # stats_merged['net_water_req (m3)'] = (stats_merged['net_water_req (mm)']/1000)*(stats_merged['AREA'])
    stats_merged['net_water_req (m3)'] = (stats_merged['net_water_req (mm)']/1000)*(stats_merged['Irrigable_m2'])
    coming_supply = pd.read_csv(supply_path)
    coming_supply['Date'] = pd.to_datetime(coming_supply['Date'], format='%m/%d/%Y')
    # discharge = coming_supply[coming_supply['Date']==interested_date]['Discharge'].values[0]
    discharge = get_discharge(interested_date, coming_supply)
    print(f'Discharge:{discharge}')
    cusecs_to_m3_day = 2446.58
    discharge_m3d = discharge*cusecs_to_m3_day
    print(f'Discharge (m3/day):{discharge_m3d}')
    stats_merged['Water Provided'] = stats_merged['Distribution_Factor']*discharge_m3d
    stats_merged['Deficit between requirement and supply'] = stats_merged['Water Provided']+stats_merged['net_water_req (m3)'] 
    stats_merged['Irrigation Status'] = stats_merged['Deficit between requirement and supply'].apply(lambda x: 'Surplus' if x > 0 else ('Right Amount' if x == 0 else 'Deficit'))
    os.makedirs(f'{data_run_folder}/distribution_files', exist_ok=True)
    stats_merged.to_csv(f'{data_run_folder}/'+'/distribution_files'+f'/Distribution_File_Stats_{interested_date.replace("-", "_")}.csv',index = False)
    return stats_merged

def main():
    parser = argparse.ArgumentParser(description="Generates water distribution csv file based on the latest water supply conditions. Water supply conditions should be in CSV format. The output will be saved in the distribution_files folder.")
    parser.add_argument("-i","--inputfolder", type = str, help="Path to the input folder. By default input folder is the run folder.")
    parser.add_argument("-d","--interested_date", type = str, help="Latest supply date in YYYY-MM-DD format, or the supply date of the run. By default takes the date of the run.")
    parser.add_argument("-s","--shpfile", type = str, help="Path to the shapefile or a csv file with the geometry column in it.")
    args = parser.parse_args()
    Irr_status(args.inputfolder, args.interested_date, args.shpfile)

if __name__ == "__main__":
    main()