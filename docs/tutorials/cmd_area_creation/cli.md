# Command Line Interface (CLI) Guide for Generating Command Areas

This function generate command-area polygons around a canal/line network using a per-feature
iterative buffer so that each polygon area matches a target area.

### Basic Usage

```bash
python -m sdrips.cmd_area_creation -s '../../Shapefiles/canal_network/Irrigation_Canal_Networks.shp' -c "CMD_Ar_m2" -o '../Data/created_cmd_area/created_cmd_area.shp'
```
### Checking Help
```bash
python -m sdrips.cmd_area_creation -h
```