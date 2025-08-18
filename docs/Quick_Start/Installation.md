# User Guide
To configure sDRIPS for your area of interest, follow the steps described in the [***Installation***](#installation) section. sDRIPS is platform-independent and can be installed on **Windows, Linux, or macOS**. It has a modular architecture that provides flexibility through configuration files.

- Use the Standard Version for most applications. This version contains the complete workflow and is recommended for general use.  
- Use the Developer Version if you want to extend sDRIPS by adding custom modules or adapting it for research applications by tweaking the source code. Detailed instructions for the developer version can be found in [***Developer Verison***](/Developer_Version/Developer_Version).  

## Requirements
To setup the sDRIPS, following things are required:  

- Anaconda/Miniconda. To install the appropriate distribution for your operating system, click [here](https://www.anaconda.com/download/success).  
- Credentials for Google Earth Engine. Required for processing satellite imagery. If you do not have credentials, follow the [***Getting Ready***](/Quick_Start/Getting_ready).  
- Credentials to download IMERG-GPM precipitation data, Required for downloading precipitation data via the NASA Precipitation Processing System (PPS). Credentials can be requested by following the [***Getting Ready***](/Quick_Start/Getting_ready).  

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
Install sDRIPS either using pip (recommended at the moment) or mamba (under development).  
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
Initialization is required for every new sDRIPS project. It creates the configuration files necessary to run workflows.

Run the following command in your project directory (e.g., sdrips_project):
```
sdrips init -d ./sdrips_project
```
After initialization, four configuration files will be created in the `config_files/` directory. These files form the core of sDRIPS, as they store user inputs and control the execution of workflows. At this stage, add your credential details (IMERG-GPM PPS account and Google Earth Engine service account, if applicable) to the `secrets.yaml` file.

## Testing
After initialization, you can verify whether sDRIPS has been installed correctly.
This test step is required only for **first-time installations**. 
Under Developement!