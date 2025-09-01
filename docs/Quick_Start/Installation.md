# User Guide
To configure sDRIPS for your area of interest, follow the steps in the [***Installation***](#installation) section. sDRIPS is platform-independent and can be installed on **Windows, Linux, or macOS**. Its modular architecture offers flexibility through configuration files.

- Use the Standard Version for most applications. This version contains the complete workflow and is recommended for general use.  
- Use the Developer Version if you plan to extend sDRIPS by adding custom modules or adapting it for research applications by tweaking the source code. Detailed instructions for the developer version can be found in [***Developer Verison***](/en/latest/Developer_Version/Developer_Version).  

## Requirements
To setup the sDRIPS, following things are required:  

- Anaconda/Miniconda - To install the appropriate distribution for your operating system, click [here](https://www.anaconda.com/download/success).  
- Credentials for Google Earth Engine - Required for processing satellite imagery. If you do not have credentials, follow the [***Getting Ready***](/en/latest/Quick_Start/Getting_ready).  
- Credentials to download IMERG-GPM precipitation data - Required for downloading precipitation data via the NASA Precipitation Processing System (PPS). Instructions are provided in the [***Getting Ready***] guide.  

## Installation
<span class="preparation_step">Step 1:</span> <br>
Create an empty project directory.
```
mkdir ./sdrips_project
```
    
<span class="preparation_step">Step 2:</span> <br>
Create a conda environment
```
conda create -n sdrips_env python=3.11 pip
```
  
<span class="preparation_step">Step 3:</span> <br>
Activate your virtual environment using conda.
```
conda activate sdrips_env
```

<span class="preparation_step">Step 4:</span> <br>
Install sDRIPS using either pip (currently recommended) or mamba (under development).  
Using [pip](https://test.pypi.org/project/sdrips/):
```
pip install sdrips
```
Using mamba (under development):
```
conda install mamba -c conda-forge
mamba install sdrips -c conda-forge
```
!!! success_note "Success"
    sDRIPS has been successfully installed. The next step is initialization, which is required before using sDRIPS. 

## Initialization
Every new sDRIPS project must be initialized. This step generates the configuration files required to run workflows.

Run the following command in your project directory (e.g., sdrips_project):
```
sdrips init -d ./sdrips_project
```
Initialization creates four configuration files in the `config_files/` directory. These files form the core of sDRIPS, storing user inputs and controlling workflow execution. At this stage, add your credential details (IMERG-GPM PPS account and Google Earth Engine service account, if applicable) to the `secrets.yaml` file.

!!! warning_note "Warning"
    You only need to initialize the project **once per directory**. Re-running initialization resets configuration files to their default values and may overwrite custom settings. To prevent accidental data loss, sDRIPS will **not overwrite existing configuration files** unless you explicitly use the `--force` option.


## Testing
After installation, you can verify whether sDRIPS has been set up correctly. Running this test suite not only ensures a successful installation but also helps detect any issues in the pipeline when new functionality is introduced.

sDRIPS supports two testing modes for advisory generation:

1. **Satellite + weather model data** – Runs the functional logic using satellite datasets in combination with numerical weather model data.

2. **Satellite + weather model data (bias-corrected with in-situ sensors)** – Extends the first mode by integrating in-situ observations to bias-correct the numerical model data before estimating evapotranspiration and generating advisories.

For more in-depth information on testing for developers, refer to the [***Developer Version***](/en/latest/Development/Developer_Version/#testing).  

!!! tip_note "Tip"
    Make sure to provide your credentials in the `secrets.yaml` file before running tests.

Running Tests
```
sdrips test -d ./tests
```

Expected Output  

A successful installation produces output similar to the following:
```yaml
Tests (3 Tests) for ET based raster outputs completed successfully.

Checks (2 Checks) for precipitation raster outputs completed successfully.

Tests (2 Tests) for CSV outputs completed successfully.

All tests passed!
```

!!! tip_note "Tip"
    Need help with the CLI? Run `sdrips -h` or `sdrips --help` to see all available options. Check [***Commands***](/en/latest/Commands/cli/) section for more details.
