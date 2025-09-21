# Command Line Interface (CLI) Guide: Complete sDRIPS Run in Command Area Advisory Mode

In this tutorial, we will set up sDRIPS to generate advisories for command areas. A command area refers to the land irrigated by a particular canal, which may contain multiple farm fields. A similar setup can also be created for any individual farm fields if shapefiles for those fields are available.

!!! lightning_note "Note"
    sDRIPS supports multiple configuration pathways. For clarity, this tutorial demonstrates one representative setup.

Follow the steps below after completing the [***installation***](en/latest/Quick_Start/Installation/#installation).

<span class="preparation_step">Step 1: Initialization</span> <br>
Run the following command in your terminal (make sure your [***conda environment***](en/latest/Quick_Start/Installation/#installation) is activated).
```
sdrips init -d .
```

Expected output
```yaml
Fetching latest configurations:
• config_links.yaml... Configuration Files have been downloaded successfully.
• crop_config.yaml... Configuration Files have been downloaded successfully.
• sdrips_config.yaml... Configuration Files have been downloaded successfully.
• secrets.yaml... Configuration Files have been downloaded successfully.
```

<span class="preparation_step">Step 2: Secrets</span> <br>
Fill in the required details in `secrets.yaml`. A Google Earth Engine (GEE) service account is recommended, though a standard GEE account can also be used.

<span class="preparation_step">Step 3: [Run Tests (Optional)](en/latest/Quick_Start/Installation/#testing) </span> <br>
If you haven’t already, run the test suite. Tests can be run using only satellite and model data configurations, or by including in-situ data along with satellite/model inputs. Check the [***testing guide***](en/latest/Quick_Start/Installation/#testing) for more details.

<span class="preparation_step">Step 4: [Create ca_config.yaml file](en/latest/tutorials/cmd_area_creation/cli/) </span> <br>
Run the following command in your terminal. At this point, you should have a shapefile describing the command areas or farm fields. Here we assume: 

- A uniform crop (Rice) is grown across all fields.
- Planting date is 15 January 2023.
- `CNLNM_ID` is the column containing unique IDs for the command areas. 
```
python -m sdrips.cmd_config -s "Shapefiles/TBP_Command_Area/commandarea_teesta_ph1_simplified_modified_Irrigable_Area_Final.shp" -c CNLNM_ID -d 2023-01-15 -cc Rice
```
!!! tip_note "Tip"
    
    - If the column name contains spaces, wrap it in quotes (e.g., "CNLNM ID").  
    - Entries can be checked at `config_files/ca_config.yaml`. 
    - Crop entries can also be changed manually in `ca_config.yaml`.  
    - For help on inputs, run:
    ```
    python -m sdrips.cmd_config -h
    ```
!!! lightning_note "Note"
    While sDRIPS supports multiple crops, for simplicity we assume rice is grown across all fields in this tutorial.

<span class="preparation_step">Step 5: Update the `sdrips_config.yaml`</span> <br>
Typical sections to update include: 

- Irrigation_cmd_area_shapefile, 
- GEE_Asset_ID, Date_Running,
- Cmd_Area_Config, 
- Irrigation_cmd_area_shapefile_Bounds [Optional]  

For this tutorial, we will generate advisories for farmers without in-situ sensor integration.The high-level controls in `sdrips_config.yaml` should look like this:
```yaml
# -------------------------------
# HIGH-LEVEL MODULE CONTROLS
# -------------------------------
Command_Area_Net_Water_Requirement: true  # Generate advisory on how much water should be allotted to each command area
Canal_water_allotment: false          # Generate advisory on how much water should be allotted to each canal (From primary to tertiary canals)
Insitu_Sensor_integration: false  # Toggle to integrate in-situ sensor data with satellite data for bias correction
Weather_station_integration: false  # Toggle to integrate local weather station data to improve accuracy
```

<span class="preparation_step">Step 6: Run sDRIPS </span> <br>
Finally, run sDRIPS with the updated configuration:
```
sdrips run -c config_files/sdrips_config.yaml 
```
!!! tip_note "Tip"
    Be patient while assets are exported to GEE for processing. Duration depends on the size of the shapefile.

Progress bars will show the current execution step and estimated completion time. At the end, you should see **sDRIPS execution complete. Check logs.**  

Check outputs:  

- Raster outputs → `Data/Landsat/`  
- Percolation Loss → `percolation/Percolation_currentweek.csv`
- Net water requirement results → `Landsat_Command_Area_Stats.csv`


!!! congrats_note "Success"
    Congratulations! You’ve completed a full sDRIPS run in the Command Area Advisory setting. This mode is designed to generate advisories that can be directly used by farmers.