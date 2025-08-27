# sDRIPS Configuration File
The **sDRIPS configuration file** (`sdrips_config.yaml`) serves as the **primary control interface** for running the model. It provides a structured mechanism for tailoring simulations to user needs and computational resources.  

Configuration parameters are divided into two major categories:  

1. **High-level module controls**  

    - Define the overall **functional mode** of sDRIPS (e.g., generating advisories for farmers or canal operators).  
    - Specify whether to integrate **in-situ sensor data**, **local weather station inputs**, or rely exclusively on **satellite and numerical weather model estimates** without bias correction.  
  
2. **Low-level module controls**  
    
    - Configure common parameters such as the **running period**, **boundary polygon location**, **project directory paths**, **data-cleaning options**, etc.  

Together, these controls determine the **scope, inputs, and computational logic** of a given sDRIPS run.  


## High-level module controls
High-level controls govern the *core objectives* of a model run. A high-level control in sDRIPS is defined as a configuration switch that activates a **distinct analytical capability or output dimension** beyond the core water requirement estimation workflow.  
Unlike low-level controls (which specify run parameters, inputs, file structures and output format), high-level controls represent **conceptual modules** that extend the scientific scope of sDRIPS.  

Examples:  
- **New considerations** → Incorporating additional physical processes (e.g., deep groundwater–surface water interaction).  
- **New outputs** → Generating novel advisory products (e.g., crop yield forecasts, fertilizer optimization metrics).  

In essence, a high-level control introduces a **new modeling domain, process pathway, or decision-support output** that meaningfully expands the functionality of sDRIPS. 
  
Users can activate or deactivate modules to suit specific applications such as:  

- **Command area advisories** → estimating water requirements at the irrigation command level.  
- **Canal water allocation** → distributing water resources across canal networks.  
- **Data integration modes** → enabling or disabling in-situ sensor and weather station data to improve accuracy.  

Below is a representative snippet from the configuration file:  

```yaml
# -------------------------------
# HIGH-LEVEL MODULE CONTROLS
# -------------------------------
Command_Area_Net_Water_Requirement: true   # Advisory: water requirement per command area
Canal_water_allotment: false               # Advisory: water allotment per canal (primary → tertiary)
Insitu_Sensor_integration: false           # Integrate in-situ sensor data for bias correction
Weather_station_integration: false         # Integrate local weather station data for improved accuracy
```

### Command Area Water Allotment
This control is designed to generate **irrigation advisories at the command area level** (i.e., farmer-focused outputs).  
When enabled, sDRIPS estimates the **net crop water requirement** within each defined command area and translates it into actionable recommendations for farmers.  

- If the user wishes to incorporate **in-situ sensor data** or **local weather station observations**, those features should be activated in parallel while keeping this control enabled.  
- This functionality is most useful for **field-scale water management and scheduling decisions**.  


### Canal Water Allotment
This control generates **advisories for canal operators and irrigation managers**.  
When enabled, sDRIPS:  

1. Estimates the **net water requirement** for all fields served by a canal command area.  
2. Aggregates demand across the irrigation hierarchy:  
    - **Tertiary canals** (child nodes) →  **Secondary canals** →  **Primary canals** (root canal).  

This hierarchical aggregation ensures that water requirements are properly scaled from field-level needs to the **canal network level**, supporting equitable and efficient distribution.  

- In-situ or weather station integration can be activated alongside this featurefor improved accuracy.  
- This module is particularly relevant for **system-level water allocation and canal scheduling**.  

### In-Situ Sensor Configuration
Activating this control enables the incorporation of **locally deployed sensors** (e.g., air temperature, wind speed, soil moisture) into the analysis.  

- Sensor data are used to **bias-correct satellite and model inputs**, thereby improving the fidelity of water requirement estimates.  
- This feature is intended for applications where **ground validation data are available**.  


### Weather Station Configuration
This control incorporates **local meteorological station data** (e.g., precipitation, air temperature, wind speed, atmospheric pressure) into sDRIPS.  

- Similar to in-situ sensors, weather station inputs enhance **local accuracy** by replacing or correcting global model estimates.  
- This module is particularly valuable in regions with **dense weather station coverage**.  

### Adding a New High-Level Control
To extend the capabilities of sDRIPS, new high-level controls can be defined using the [***Developer version***](/en/latest/Development/Developer_Version/).  

A **high-level control** refers to a feature that introduces a **new analytical capability or output dimension** to the system (e.g., crop yield prediction, fertilizer application effects, groundwater–surface water interactions).  

Adding a new control involves:  
1. Declaring the control in `sdrips_config.yaml` with a **Boolean toggle** (`true` / `false`).  
2. Implementing the corresponding logic in the source code so that the new control activates relevant computations or outputs.  
3. Updating documentation to reflect the new functionality.  
 

## Low-level module controls
Low-level controls govern the **operational details** of an sDRIPS run. These parameters are common across high-level modules and provide essential contextual inputs:  

- **Boundary definition** → Path to the region of interest boundary file (Shapefile or GeoJSON), field name of interest, and unique identifier.  
- **Running period** → Start and end dates for the run.  
- **Project structure** → Directory paths, links to supplementary configuration files, and output storage locations.  
- **Precipitation configuration** → Specification of how precipitation should be incorporated into the analysis (a parameter shared across all high-level modules).  
- **Data handling options** → Directory-cleaning settings, overwrite rules, and log management.  
- **Output options** → Outputs can include PNGs (for practitioners), CSVs (for researchers and practitioners), and rasters (for advanced analysis).



