# Command Line Interface (CLI) Guide for Generating Command Area Configuration File

This function reads a vector dataset (e.g., Shapefile, GeoJSON, GeoPackage, GDB file)
using GeoPandas, extracts the unique command‑area identifiers from the specified
attribute column, and writes a YAML config with one section per ID plus a DEFAULT
section containing baseline parameters. Here Shapefile format is highly recommended as it 
can be integrated with Google Earth Engine easily. 

### Basic Usage

```bash
python -m sdrips.cmd_config -s "../../Shapefiles/cmd_area_wgs/Teesta_Command_Areas.shp" -c "CNLNM" -d 2025-07-07 -cc Rice -o '../Data/cmd_config.yaml'
```
### Checking Help
```bash
python -m sdrips.cmd_config -h
```
### Expected Output
```yaml
DEFAULT:
  use_default: true
  planting_date: "2025-07-30"
  crop_type: Wheat
  soil_coef: 0.6
  distribution_unif: 0.8

Teesta Canal A:
  use_default: false
  planting_date: "2025-07-30"
  crop_type: Wheat
  soil_coef: 0.6
  distribution_unif: 0.8

T1T:
  use_default: false
  planting_date: '2025-07-07'
  crop_type: Wheat
  soil_coef: 0.5
  distribution_unif: 1.0
...
```