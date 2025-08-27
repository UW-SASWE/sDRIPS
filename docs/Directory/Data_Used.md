# Data Used

The design of **sDRIPS** emphasizes *minimal user-side data requirements* while leveraging cloud-hosted datasets to ensure scalability and global applicability.  

To execute sDRIPS for a given region of interest (ROI), users must provide only two essential inputs:  

1. A **boundary polygon file** (e.g., Shapefile or GeoJSON) that defines the spatial extent of the analysis.  
2. The **operational mode** in which sDRIPS should run.  

Beyond these minimal requirements, sDRIPS automatically integrates multiple **satellite remote sensing** and **model-based datasets** via [Google Earth Engine (GEE)](/en/latest/Quick_Start/Getting_ready/#setup-gee), including:  

- **Satellite imagery** (e.g., Landsat, Sentinel)  
- **Precipitation estimates** (NASA IMERG)  
- **Meteorological inputs** (NOAA GFS)  

Where available, **in situ measurements** (e.g., local weather station or sensor data) can supplement or replace global datasets. This adaptive design allows sDRIPS to operate seamlessly across varying levels of data availability, from global-scale applications to local field-level analysis.  

- **Table 1** (adapted from *Khan et al., (2025)*) summarizes the minimal data and credentials required to run sDRIPS.  
- **Table 2** (adapted from *Khan and Hossain, (2025)*) outlines the global datasets integrated into sDRIPS for estimation and analysis. 

### Table 1. Minimal data and credentials required to run sDRIPS  

| **Sr. No.** | **Requirement**                               | **Purpose**                                                                 |
|-------------|-----------------------------------------------|-----------------------------------------------------------------------------|
| 1           | Google Earth Engine (GEE) Account             | Access to cloud computing power and hosted satellite datasets               |
| 2           | Precipitation Processing System (PPS) Account | Download NASA IMERG precipitation data                                      |
| 3           | Boundary Polygon (Shapefile / GeoJSON)        | Defines the spatial extent of the analysis                                  |
| 4           | Crop Type and Planting Date                   | Provides crop-specific information needed for accurate ET estimation        |
| 5           | Sensor or Weather Station Data (optional)     | Improves estimation of net water requirement; required only if in-situ data are to be used |


### Table 2. Global datasets integrated into sDRIPS  

| **Dataset**                          | **Dataset ID on GEE**                                    | **Band Name / Derived Products**           | **Description**                                   | **Spatial Resolution** | **Temporal Resolution**         |
|--------------------------------------|----------------------------------------------------------|---------------------------------------------|---------------------------------------------------|------------------------|---------------------------------|
| Global Forecasting System (GFS)      | NOAA/GFS0P25                                             | temperature_2m_above_ground                 | Air Temperature                                  | 25 km                  | 6 hours                         |
|                                      |                                                          | u_component_of_wind_10m_above_ground        | Wind Speed (u component)                         |                        |                                 |
|                                      |                                                          | u_component_of_wind_10m_above_ground        | Wind Speed (v component)                         |                        |                                 |
|                                      |                                                          | specific_humidity_2m_above_ground           | Specific Humidity                                |                        |                                 |
|                                      | ----                                                     | ----                                        | Pressure (Estimated using Hypsometric Equation)   |                        |                                 |
|                                      | GEE and NOAA                                             | total_precipitation_surface                 | 168-hour Precipitation                           |                        |                                 |
| Landsat 8 and Landsat 9 satellite    | LANDSAT/LC08/C02/T1_TOA <br> LANDSAT/LC09/C02/T1_TOA     | B2                                          | Blue                                             | 30 m                   | 16 days individual, 8 days combined |
|                                      |                                                          | B4                                          | Red                                              |                        |                                 |
|                                      |                                                          | B5                                          | NIR                                              |                        |                                 |
|                                      |                                                          | B6                                          | Shortwave IR                                     |                        |                                 |
|                                      |                                                          | B10                                         | Low Gain Thermal                                 |                        |                                 |
|                                      |                                                          | B11                                         | High Gain Thermal                                |                        |                                 |
| Sentinel 1 satellite                 | COPERNICUS/S1_GRD                                        | VV                                          | Single co-polarization, vertical transmit/receive | 10 m                   | 10 days                         |
| Shuttle Radar Terrain Mapping (SRTM) | USGS/SRTMGL1_003                                         | ---                                         | Digital Elevation Map                            | 30 m                   | ---                             |
| Global Land Cover Classification     | COPERNICUS/Landcover/100m/Proba-V-C3/Global                                                        | discrete_classification                      | Land Use Land Cover                              | 100 m                  | ---                             |
| IMERG GPM                            | ---                                                      | ---                                         | Precipitation                                    | 10 km                  | Around 4 hours                  |
