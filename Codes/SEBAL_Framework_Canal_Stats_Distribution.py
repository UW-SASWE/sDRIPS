import geopandas as gpd
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


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


def Irr_status(folder_run,interested_date,net_water_req_corr = False,shapefile_path = '../TBP Shape Files/TCAD_Phase1_DataPartA/ShapeFiles/commandarea_teesta_ph1_simplified_modified.shp' ):
    cmd = gpd.read_file(shapefile_path)
    cmd['Distribution_Factor'] = cmd['AREA_HA']/cmd['AREA_HA'].sum()
    stats_file = pd.read_csv(f'../{folder_run}/Landsat_Command_Area_Stats.csv')
    if net_water_req_corr == False:
        try:
            stats_file['net_water_req'] = stats_file['net_water_req'].apply(lambda x: x if x <= 0 else 0)
        except:
            print('Correction wasn\'t  applied')
    stats_merged = stats_file.merge(cmd,on=['CNLNM_ID'])
    stats_merged = stats_merged.drop(['AREA_y','PERIMETER_y','ID_y','CNLSYS_y','AREA_HA_y','SEC_CNLNM_y','CNLNM_y'], axis=1)
    stats_merged = stats_merged.rename(columns = {'net_water_req':'net_water_req (mm)','Currentweek PPT':'Currentweek PPT (mm)','Nextweek PPT':'Nextweek PPT (mm)','7Day Irrigation':'7Day Irrigation (mm)','AREA_x':'AREA','PERIMETER_x':'PERIMETER','ID_x':'ID','CNLSYS_x':'CNLSYS','AREA_HA_x':'AREA_HA','SEC_CNLNM_x':'SEC_CNLNM','CNLNM_x':'CNLNM'})
    stats_merged['net_water_req (m3)'] = (stats_merged['net_water_req (mm)']/1000)*(stats_merged['AREA'])
    teesta_supply = pd.read_csv('../TBP_Data/TBP_CHR_Cleaned_Pivoted.csv')
    teesta_supply['Date'] = pd.to_datetime(teesta_supply['Date'], format='%m/%d/%Y')
    # discharge = teesta_supply[teesta_supply['Date']==interested_date]['Discharge'].values[0]
    discharge = get_discharge(interested_date, teesta_supply)
    print(f'Discharge From Teesta:{discharge}')
    cusecs_to_m3_day = 2446.58
    discharge_m3d = discharge*cusecs_to_m3_day
    print(f'Discharge From Teesta (m3/day):{discharge_m3d}')
    stats_merged['Water Provided via Teesta'] = stats_merged['Distribution_Factor']*discharge_m3d
    stats_merged['Deficit between requirement and teesta'] = stats_merged['Water Provided via Teesta']+stats_merged['net_water_req (m3)'] 
    stats_merged['Irrigation Status'] = stats_merged['Deficit between requirement and teesta'].apply(lambda x: 'Surplus' if x > 0 else ('Right Amount' if x == 0 else 'Deficit'))
    return stats_merged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate Irrigation Distribution Stats File')
    parser.add_argument('-f','--folder_run', type=str, help='Path to the folder where data is stored')
    parser.add_argument('-d','--interested_date', type=str, help='Interested date in the format YYYY-MM-DD, typically the date on which Landsat Image was found')
    parser.add_argument('-nwc','--net_water_req_corr', action='store_true', help='Apply net water requirement correction')

    args = parser.parse_args()
    
    interested_date = pd.to_datetime(args.interested_date)
    result = Irr_status(args.folder_run, interested_date, args.net_water_req_corr)
    print(result)

