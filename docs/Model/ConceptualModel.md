# Conceptual Model of sDRIPS
**sDRIPS** is a cloud-based, open-source irrigation advisory system designed to optimize surface water use by integrating Earth observation datasets and numerical weather models. It leverages publicly available data from multiple sources, including:

1. **Landsat**, **Sentinel-1**, and **SMAP** satellites for surface and soil moisture monitoring.
2. **Global Forecast System (GFS)** for meteorological forecasts.
3. **Global Precipitation Measurement (GPM)** – IMERG for precipitation estimates. 

At its core, sDRIPS estimates evapotranspiration (ET) using two models:

1. **Penman–Monteith** <a href="https://www.fao.org/4/x0490e/x0490e00.htm" target="_blank">(Allen et al., 1998)</a> method – used as a proxy of potential crop water demand under optimal (non-stressed) conditions.

2. **Surface Energy Balance Algorithm for Land (SEBAL)** <a href="https://doi.org/10.1016/S0022-1694(98)00253-4" target="_blank">(Bastiaanssen et al., 1998)</a> – used to estimate actual ET under prevailing field conditions.

By comparing these two ET estimates, sDRIPS quantifies the gap between water needed by crops and water actually supplied. This information is then combined with nowcast and forecast precipitation data as well as estimated percolation losses to calculate the net irrigation requirement.

sDRIPS applies both energy balance and water balance principles to generate actionable irrigation advisories at multiple spatial scales — from individual farms to large irrigation command areas. 

![sDRIPS Concept Illustration](../images/sDRIPS_Illustration.png)
<p style="text-align:center"> <sup><i>Illustration of how sDRIPS works. Green represents the regions where surplus irrigation has been provided, while red regions are deficit regions needing more irrigation.</i></sup> </p>


## Evapotranspiration 
<span id="penman_et"></span>
### Penman-Monteith ET
The **Penman–Monteith** evapotranspiration (\(ETo\)) represents the **potential water demand** for a well-watered reference crop. In sDRIPS, \(ETo\) is computed for each Landsat overpass date and across all regions of interest, following the methodology of <a href="https://www.fao.org/4/x0490e/x0490e00.htm" target="_blank">(Allen et al., 1998)</a>.  

The **FAO-56 Penman–Monteith** equation is:  

$$ ETo = \frac{0.408\,\Delta\,(R_n - G) + \gamma \,\frac{900}{T + 273} \, u_2 \,(e_s - e_a)}{\Delta + \gamma \,(1 + 0.34\,u_2)}$$
    
Where:  

- $ETo$ — Reference evapotranspiration (mm/day)  
- $R_n$ — Net radiation at the crop surface (MJ m$^{-2}$ day$^{-1}$)  
- $G$ — Soil heat flux density (MJ m$^{-2}$ day$^{-1}$)  
- $T$ — Mean daily air temperature at 2 m height (°C)  
- $u_2$ — Wind speed at 2 m height (m/s)  
- $e_s$ — Saturation vapor pressure (kPa)  
- $e_a$ — Actual vapor pressure (kPa)  
- $(e_s - e_a)$ — Vapor pressure deficit (kPa)  
- $\Delta$ — Slope of the vapor pressure curve (kPa °C$^{-1}$)  
- $\gamma$ — Psychrometric constant (kPa °C$^{-1}$)  

The $ETo$ computed above represents the water use of a reference crop under optimal (non-stressed) conditions. To estimate crop-specific potential evapotranspiration ($ET$), sDRIPS applies crop and water stress coefficients:  

$$
ET = ETo \times K_c \times K_s
$$

Where:  

- $ET$ — Potential evapotranspiration for the specific crop (mm/day)  
- $K_c$ — Crop coefficient (varies by crop type and growth stage)  
- $K_s$ — Soil water stress coefficient    

<span id="sebal_et"></span>
### SEBAL ET
The SEBAL algorithm, developed by <a href="https://doi.org/10.1016/S0022-1694(98)00253-4" target="_blank">Bastiaanssen et al. (1998)</a>, has been successfully implemented in more than 30 countries <a href="https://doi.org/10.1061/(ASCE)0733-9437(2005)131:1(85)" target="_blank">(Bastiaanssen et al., 2005)</a>. It has demonstrated high accuracy in estimating evapotranspiration.

SEBAL solves the surface energy balance equation to estimate evapotranspiration using satellite imagery and meteorological forcing data. It computes an instantaneous ET flux at the time of the Landsat overpass. The model calculates net surface radiation (\(R_n\)), soil heat flux (\(G\)), and sensible heat flux to the air (\(H\)). The residual energy flux, which is the latent heat flux driving evapotranspiration, is obtained by subtracting the soil and sensible heat fluxes from net radiation at the surface. For each pixel, the latent heat flux (\( \lambda E \)) is given by the surface energy balance equation:

\[
\lambda E = R_n - G - H
\]

Where:

- \( \lambda E \) is the latent heat flux (energy used for evapotranspiration),  
- \( R_n \) is net radiation at the surface,  
- \( G \) is soil heat flux, and  
- \( H \) is sensible heat flux.

The SEBAL-based evapotranspiration (SEBAL ET) serves as a proxy for the **actual water consumed by the crop** <a href="https://doi.org/10.1029/2020WR028654" target="_blank">Bose et al. (2021)</a>. In sDRIPS, daily (24-hour) evapotranspiration is estimated assuming that instantaneous ET variations are negligible over the day <a href="https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)" target="_blank">(Allen et al., 2007)</a>. For weekly ET, the evapotranspiration estimated on the Landsat acquisition day is assumed steady for the following seven days until the next image is available.

<span id="percolation"></span>
## Percolation
Percolation represents water that moves beyond the active crop root zone and is therefore no longer available to meet near-term crop water demand. In the updated sDRIPS framework, percolation is estimated by combining fine-scale Sentinel-1 SAR information with coarser soil moisture and soil hydraulic datasets, following the published methodology of <a href="https://doi.org/10.5194/hess-30-3675-2026" target="_blank">Khan et al. (2026)</a>.

Global soil moisture products such as SMAP provide useful volumetric soil moisture information, but their native spatial resolution is too coarse for command-area or field-scale irrigation advisories. Sentinel-1 C-band SAR, available through Google Earth Engine (GEE), provides finer spatial detail at approximately 10 m resolution and is therefore used as the primary source for estimating the spatial pattern of relative soil moisture. Sentinel-1 C-band observations have been widely used for soil moisture estimation and have shown useful sensitivity near the upper soil layer, approximately up to 100 mm depth (<a href="https://doi.org/10.1109/TGRS.2018.2858004" target="_blank">Bauer-Marschallinger et al., 2019</a>; <a href="https://doi.org/10.1016/j.asr.2022.03.019" target="_blank">Bhogapurapu et al., 2022</a>).

The Sentinel-1 VV backscatter signal is converted into a relative soil moisture index by comparing the current backscatter state with dry and wet reference conditions over the analysis window:

\[
SMI = \frac{\sigma^0_{VV} - \sigma^0_{VV,dry}}{\sigma^0_{VV,wet} - \sigma^0_{VV,dry}}
\]

Where:

- \(SMI\) -- Sentinel-1 relative soil moisture index  
- \(\sigma^0_{VV}\) -- Sentinel-1 VV backscatter for the analysis date or composite  
- \(\sigma^0_{VV,dry}\) -- dry reference backscatter condition  
- \(\sigma^0_{VV,wet}\) -- wet reference backscatter condition  

The relative Sentinel-1 soil moisture index is then scaled to volumetric soil moisture using SMAP volumetric soil moisture estimates, and the scaled volumetric estimate is converted to water depth:

\[
SM_{mm} = \theta_{SMAP\text{-}scaled} \times Z
\]

Where:

- \(SM_{mm}\) -- soil water depth in millimeters  
- \(\theta_{SMAP\text{-}scaled}\) -- volumetric soil moisture after scaling the Sentinel-1 relative index with SMAP  
- \(Z\) -- equivalent soil depth used for the percolation estimate  

Soil moisture at field capacity is derived from the ISRIC SoilGrids dataset available on GEE (<a href="https://isric.org/explore/soilgrids" target="_blank">ISRIC SoilGrids</a>; Poggio et al., 2021). Field capacity is converted to the same water-depth basis:

\[
FC_{mm} = \theta_{FC} \times Z
\]

Percolation is calculated for each pixel only when the estimated soil water depth exceeds field capacity:

\[
Percolation =
\begin{cases}
SM_{mm} - FC_{mm}, & SM_{mm} > FC_{mm} \\
0, & SM_{mm} \leq FC_{mm}
\end{cases}
\]

Pixel-level percolation estimates are then summarized over each command area and used as a loss term in the net water requirement calculation.

<span id="precipitation"></span>
## Precipitation
Accurate estimation of net water requirement depends on incorporating precipitation events. Precipitation for the operational week is derived from GPM IMERG early run data.

The cumulative precipitation over the week before the latest Landsat overpass is considered as **nowcast precipitation**. Additionally, the cumulative forecasted precipitation for the next week after the Landsat overpass is included.

## Net Water Requirement
The net water requirement is calculated as:

$$
\text{Net Water Requirement} = \text{SEBAL ET} - \text{PET} - \text{Per} + P_{nc} + P_{fc}
$$

Where:  

- [**SEBAL ET**](#sebal_et): Estimated actual evapotranspiration  
- [**PET**](#penman_et): Potential evapotranspiration (Penman-Monteith method)  
- [**Per**](#percolation): Percolation losses  
- [**P<sub>nc</sub>**](#precipitation): Nowcast precipitation  
- [**P<sub>fc</sub>**](#precipitation): Forecasted precipitation  



**Interpretation:**

1. If net water requirement is **positive**, the water supply (including precipitation and percolation) meets or exceeds crop demand, so no additional irrigation is needed.

2. If net water requirement is **negative**, the water supply is insufficient, indicating additional irrigation is required for the upcoming week.

