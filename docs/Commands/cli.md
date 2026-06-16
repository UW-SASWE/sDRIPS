# Command Line Interface of sDRIPS

This page summarizes the command-line workflows exposed by sDRIPS. Use it as a quick reference after completing the [installation guide](/en/latest/Quick_Start/Installation/).

!!! tip "Recommended workflow"
    Run commands from the root of your sDRIPS project directory. This keeps relative paths such as `config_files/sdrips_config.yaml`, `Shapefiles/...`, and `Data/...` easy to manage.

## Command Map

| Task | Command | Main output |
| --- | --- | --- |
| Initialize a project | `sdrips init` | Project folders and default configuration files |
| Verify the installation | `sdrips test` | Downloaded test data, generated outputs, comparison results |
| Run the model | `sdrips run` | ET, precipitation, percolation, statistics, and advisory outputs |
| Build command-area config | `python -m sdrips.cmd_config` | `ca_config.yaml` |
| Estimate command-area polygons | `python -m sdrips.cmd_area_creation` | Command-area shapefile |

!!! note "Installed command and module commands"
    The installed package exposes `sdrips init`, `sdrips test`, and `sdrips run`.
    Utility workflows such as command-area configuration and command-area polygon creation are currently run with `python -m`.

## Initialize a Project

Use `sdrips init` once for each new project directory.

```bash
sdrips init -d ./sdrips_project
```

The initialization command creates the working folders:

```text
sdrips_project/
  Data/
  Shapefiles/
  config_files/
```

It also downloads the default configuration files into `config_files/`:

| File | Purpose |
| --- | --- |
| `sdrips_config.yaml` | Main workflow controller |
| `secrets.yaml` | Google Earth Engine and IMERG/PPS credentials |
| `crop_config.yaml` | Crop coefficient curves by crop type |
| `config_links.yaml` | External data-source URL patterns |

### Options

| Option | Description |
| --- | --- |
| `-d`, `--dir` | Project directory. Default is the current directory. |
| `-f`, `--force` | Force update of existing configuration files. |

!!! warning "Using `--force`"
    `--force` can replace local configuration files with fresh defaults. Use it only when you intentionally want to refresh the template files.

## Run Verification Tests

Use the test command to check whether sDRIPS can download test data, run the workflow, and compare outputs.

```bash
sdrips test -d ./tests
```

To run the sensor-corrected test case:

```bash
sdrips test -d ./tests -s
```

### Options

| Option | Description |
| --- | --- |
| `-d`, `--dir` | Test directory. Default is `./tests`. |
| `-s`, `--sensor` | Enables the in-situ sensor test mode. |

!!! note "What the test command does"
    The test workflow downloads expected test files, copies project configuration files, updates paths for testing, runs sDRIPS, and compares selected raster and CSV outputs.

## Run sDRIPS

After updating `sdrips_config.yaml`, run the full workflow:

```bash
sdrips run -c config_files/sdrips_config.yaml
```

### Options

| Option | Description |
| --- | --- |
| `-c`, `--config` | Path to the main sDRIPS configuration file. Default is `./config_files/sdrips_config.yaml`. |

The modules that run are controlled from `sdrips_config.yaml`.

| Configuration section | Controls |
| --- | --- |
| `Command_Area_Net_Water_Requirement` | Command-area advisory generation |
| `Canal_water_allotment` | Canal-level water allotment |
| `Insitu_Sensor_integration` | Bias correction with in-situ sensor data |
| `Weather_station_integration` | Bias correction with weather-station data |
| `Run_ET_Estimation` | ET estimation workflow |
| `Precipitation_Config` | IMERG precipitation workflow |
| `Weather_Config` | Forecasted weather data workflow |
| `Percolation_Config` | Percolation estimation |
| `Region_stats` | Command-area statistics |

Typical outputs are written under the configured `Save_Data_Location`.

```text
Data/
  Landsat/
  percolation/
  logs/
  Landsat_Command_Area_Stats.csv
```

!!! tip "Logs"
    At the end of a run, sDRIPS prints the log-file path. Check that log first if a module finishes unexpectedly or a data source fails.

## Generate Command-Area Configuration

Use `sdrips.cmd_config` to generate `ca_config.yaml` from a command-area vector file. This file defines crop and planting information for each command area.

```bash
python -m sdrips.cmd_config \
  -s "Shapefiles/cmd_area/Teesta_Command_Areas.shp" \
  -c "CNLNM" \
  -d 2025-07-07 \
  -cc Rice \
  -o config_files/ca_config.yaml
```

On Windows PowerShell, use backticks for line continuation:

```powershell
python -m sdrips.cmd_config `
  -s "Shapefiles/cmd_area/Teesta_Command_Areas.shp" `
  -c "CNLNM" `
  -d 2025-07-07 `
  -cc Rice `
  -o config_files/ca_config.yaml
```

### Options

| Option | Required | Description |
| --- | --- | --- |
| `-s`, `--shp_path` | Yes | Path to the command-area vector file. |
| `-c`, `--column_name` | Yes | Attribute column containing unique command-area IDs. |
| `-d`, `--default_planting_date` | No | Default planting date in `YYYY-MM-DD` format. |
| `-cc`, `--default_crop_type` | No | Default crop type, such as `Rice`, `Wheat`, or `Corn`. |
| `-sc`, `--default_soil_coef` | No | Default soil coefficient. |
| `-du`, `--default_distribution_unif` | No | Default distribution uniformity. |
| `-o`, `--output_path` | No | Output YAML path. Default is `config_files/ca_config.yaml`. |

Expected structure:

```yaml
DEFAULT:
  use_default: true
  planting_date: 2025-07-07
  crop_type: Rice
  soil_coef: 0.5
  distribution_unif: 1.0

Command_Area_1:
  use_default: false
  planting_date: 2025-07-07
  crop_type: Rice
  soil_coef: 0.5
  distribution_unif: 1.0
```

!!! tip "After generating `ca_config.yaml`"
    Review the generated file before running sDRIPS. You can edit crop type, planting date, soil coefficient, or distribution uniformity for individual command areas.

## Create Command-Area Polygons

Use `sdrips.cmd_area_creation` when command-area boundaries are not available and you need to estimate polygons from a canal or line network.

```bash
python -m sdrips.cmd_area_creation \
  -s "Shapefiles/canal_network/Irrigation_Canal_Networks.shp" \
  -c "CMD_Ar_m2" \
  -o "Shapefiles/created_cmd_area/created_cmd_area.shp"
```

### Options

| Option | Required | Description |
| --- | --- | --- |
| `-s`, `--shp_path` | Yes | Path to the input canal or line-network shapefile. |
| `-o`, `--output` | Yes | Path for the output command-area shapefile. |
| `-c`, `--column_name` | No | Column containing target area values in square meters. Default is `Command Area (m2)`. |

!!! warning "Estimated boundaries"
    This command estimates polygons using an iterative buffer. If measured command-area or field boundaries are available, use those boundaries instead.

## Help Commands

Each command can print its available options:

```bash
sdrips -h
sdrips init -h
sdrips test -h
sdrips run -h
python -m sdrips.cmd_config -h
python -m sdrips.cmd_area_creation -h
```

## Related Tutorials

| Tutorial | Use it when |
| --- | --- |
| [Full CLI guide](/en/latest/tutorials/full_guide/cli/) | You want a complete command-area advisory workflow. |
| [Creating command areas](/en/latest/tutorials/cmd_area_creation/cli/) | You need to estimate command areas from a canal network. |
| [Creating command-area config files](/en/latest/tutorials/config_files/cli/) | You need a focused guide for generating `ca_config.yaml`. |

## Developer and Utility Notes

The current source code does not expose additional terminal commands from `sdrips.utils`. Those modules are internal helpers used by the run pipeline. Useful utilities for notebooks and development are documented in the [notebook interface](/en/latest/Commands/notebook/).
