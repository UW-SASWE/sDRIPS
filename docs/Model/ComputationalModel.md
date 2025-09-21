# Computational Model of sDRIPS

The schematic flowchart below illustrates the decision logic, where execution pathways are determined by user-defined configuration settings. At a high level, the pipeline is composed of five sequential components:  

1.	Configuration setup (orange patch) – Reading and validating user-defined parameters from the configuration files. Here user also creates the [`ca_config.yaml`](en/latest/tutorials/config_files/cli/). 

    !!! tip_note "Note"
        If sDRIPS is run in a canal-operator setting and command-area boundaries are unavailable, the `cmd_config` module can be used to generate a baseline structure of command areas. This process is iterative, and the resulting command areas may serve as approximations rather than fully accurate representations.  

2.	Evapotranspiration estimation (purple patch) – Computing crop water demand using Penman-Monteith and SEBAL models.  
3.	Percolation loss estimation (pink patch) – Quantifying percolation losses specific to field conditions.  
4.	Precipitation integration (blue patch) – Incorporating both nowcasted and forecasted rainfall from global models or local sources (only nowcast).  
5.	Advisory generation (green patch) – Comparing water balance components to produce net water requirement advisories.  


![sDRIPS Logic](../images/Model/sDRIPS_Logic.png)
<p style="text-align:center"> <sup><i>Flowchart of the sDRIPS pipeline with decision logic defined by user configurations.</i></sup> </p>
