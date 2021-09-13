# Project Title

Code supplement for Slater et al (2021) - Increased variability in Greenland Ice Sheet runoff from satellite observations

## Description

MATLAB code used to generate data and figures described in Slater et al (2021) - Increased variability in Greenland Ice Sheet runoff from satellite observations

### Dependencies

* MATLAB
* Various incredibly useful functions available in the [Climate Data Toolbox](https://www.chadagreene.com/CDT/CDT_Getting_Started.html)

### Scripts

#### /main/

##### seasonal_elevation_changes.m

Produce seasonal elevation change time series from CryoSat-2 plane fit data in the Greenland Ice Sheet Ablation Zone

##### compute_runoff.m

Compute runoff from CryoSat-2 plane fit data in the Greenland Ice Sheet within identified ablation and inland runoff zones

##### dz_mean.m

Function to derive average elevation time series and estimated uncertainty within a user defined area

##### dz_to_runoff.m

Function to convert elevation changes from CryoSat-2 elevation data (in x,y,time) to runoff within user defined area, allowing for corrections for the effects of ice advection, dynamics and converting to mass at a defined density

#### /figs/

##### plot_figX.m
Scripts used for generating figures used in the manuscript (X = 1-4)

## Authors

Thomas Slater (t.slater1@leeds.ac.uk)

## Version History

* 1.0
    * Initial Release