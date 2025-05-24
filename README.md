# sD.R.I.P.S - Framework

## Please Note: Major changes are coming to create a Python package from this framework and make it more user friendly and structured. 

### Satellite Data Rendered Irrigation using Penman and SEBAL (sD.R.I.P.S) for Surface Water Irrigation Optimization: 

<div style="text-align: justify;">
This study proposes a water-provider-centric irrigation advisory system designed to manage surface water resources and allocate water efficiently to areas in need, thereby promoting sustainable irrigation practices in the context of a changing climate. The system utilizes satellite remote sensing based SEBAL (Surface Energy Balance Algorithm for Land) and Penman-Monteith evapotranspiration models to estimate crop water use. By integrating the responses from the previous irrigation cycle, current precipitation, forecasted precipitation, and evapotranspiration-based water needs, the framework calculates the net water requirements for command areas within irrigation canal networks. Operating on a weekly basis, the system generates advisories that enable water providers to make informed, science-based decisions about water allocation. These advisories quantify the water needs, giving water providers the flexibility to direct water to areas of higher need based on their judgment. Additionally, the proposed framework can simulate future scenarios by assuming potential policy changes come into effect today and evaluating their impact on water supply to irrigation areas. This allows stakeholders to assess the response of crops to such changes and develop strategies to mitigate adverse effects on water supply conditions.
</div>

### Repository Outline:
<div style="text-align: justify;">
This repository mirrors the directory structure of the sD.R.I.P.S framework prior to setup for any specific region of interest. The repository includes:

- **Codes folder**: Contains the scripts necessary to run the framework.
- **Config_Files folder**: Holds the configuration files required for the framework's operation.
- **Shapefiles folder**: Users can add shapefiles of their regions of interest here.
- **Distribution_Files folder**: This is where the CSV outputs of the sD.R.I.P.S framework will be stored.

</div>

### How to Use:
<div style="text-align: justify;">
The detailed functioning of the framework is outlined in the sD.R.I.P.S Paper. If you encounter any issues, please let us know. Here are the basic steps to set up the sD.R.I.P.S for your region of interest:

1. **Clone the repository**
2. **Storing the shapefiles and creation of `Canal_Config.ini`**. 
    - Store the region of interest shapefiles in the `Shapefile` folder and run the `Canal_Config_Creator.py` script. Users may need to change the attribute name in the `Canal_Config_Creator.py` script to their desired attribute name. Each command area should have a unique name for reference. For example, in the `Canal_Config_Creator.py` script, it is `CNLNM_ID`.
    - For help, users can type:
      ```bash
      python Canal_Config_Creator.py -h
      ```
    Note: If users dont have the command area shapefile, but have the canal network shapefiles, users should run the `Command_Area_Creation.py` to generate an estimate of the command areas boundary.
3. **Set up the configuration files**: In this step, we set up the `Script_Config` and `Secrets.ini` files.
   
    - **`Script_Config.ini` file**: 
        - Users need to provide the following information:
          - The path of the command area shapefile.
          - The bounds of the command area shapefile.
          - Google Earth Engine (GEE) asset ID of the command area shapefile.
          - The date for which the users want to run the framework.
          - The location where the `sDRIPS_Framework.py` will save the raster and CSV outputs.

    - **`Secrets.ini` file**:
        - Users need to fill in the PPS credentials to automatically download the IMERG precipitation dataset.
4. **Run the `sDRIPS_Framework.py` script**:
    - The script will automatically estimate the Penman-Monteith and SEBAL based ET using Landsat 8 and Landsat 9 collections.
    - Additionally, the script will estimate the intensity of precipitation events that occurred in the running week and the forecasted precipitation.
    - The end products of this script are the ET rasters and the CSV files (`Landsat_Command_Area_Stats.csv`) storing the mean ET values and precipitation intensity for each command area.
    - To debug, please use the log file present in the logs folder.
    - After successfully running the script, the directory structure of the folder where `sDRIPS_Framework.py` saved its outputs will be as follows:

      ```
      save_data_loc/
      ├── avgt/
      ├── landsat/
      │   ├── irrigation
      │   │   ├── currentweek
      │   │   └── lasttweek (if ran)   
      │   ├── penman
      │   │   ├── currentweek
      │   │   └── lasttweek (if ran)     
      │   └── sebal
      │       ├── currentweek
      │       └── lasttweek (if ran)
      ├── logs/
      ├── precip/
      ├── temp/
      ├── tmax/
      ├── tmin/
      ├── ugrd/
      ├── uploads/
      ├── vgrd/
      ├── wind/  
      └── Landsat_Command_Area_Stats.csv
      ```
5. **Run the `Canal_Water_Distribution_script.py`**:  
    - This script estimates the deficit and surplus regions after the given water supply.
    - Users might need to change the names of a few columns to match their water supply files.
    - For help, users can type:
      ```bash
      python Canal_Water_Distribution_script.py -h
      ```
</div>

### How to Cite
<div style="text-align: justify;">
If you use the scripts from this repository, please cite the following paper:
Khan, S; et.al (2024). <i>Satellite Data Rendered Irrigation using Penman and SEBAL (sD.R.I.P.S) for Surface Water Irrigation Optimization, Agricultural Water Management.</i> Paper Under Review.
</div>
