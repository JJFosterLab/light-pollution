Light Pollution Forces a Change in Dung Beetle Orientation Behaviour
====================
This repository contains files relevant to the manuscript "_Light Pollution Forces a Change in Dung Beetle Orientation Behaviour_" by James J. Foster, 	Claudia Tocco, Jochen Smolka, Lana Khaldy, Emily Baird, Marcus J. Byrne, Dan-Eric Nilsson and Marie Dacke (submitted).
The source code is [statistical analysis](https://github.com/JJFosterLab/light-pollution/tree/master/Behaviour), written in ```R```, and [image analysis](https://github.com/JJFosterLab/light-pollution/tree/master/Imaging), written in ```MATLAB```, based on the **E**nvironmental **L**ight **F**ield ([ELF](https://github.com/sciencedjinn/elf)) system developed by Smolka & Nilsson (D.-E. Nilsson, J. Smolka, 2021 [Quantifying biologically essential aspects of environmental light.](https://doi.org/10.1098/rsif.2021.0184) _J. R. Soc. Interface_ 18: 20210184).  Image data (```.mat``` files) are available from [Figshare](https://doi.org/10.1098/rsif.2021.0184). 
# File Manifest
## Behaviour
Scripts for organisation, plotting and analysis of exit headings recorded during behavioural experiments. Raw data can be found in ```LP Headings all 20200505.txt``` and in a neater format in ```LPdataOrganised20200505.txt```. Mean vector lengths for each individual can be found in ```LPallRho20200528.txt``` and mean bearings in ```LPallMu20200528.txt```. ```LPHeadings-OrganiseData.R``` is a script to label data according to relevant experimental conditions. Plotting functions can be found in ```LP_PlotFunctions.R``` and rely functions from the [circular](https://cran.r-project.org/web/packages/circular/index.html) package. ```MMRayleigh.test.R``` is an implementation of Moore's modified Rayleigh test ([Moore, 1980](https://doi.org/10.1093/biomet/67.1.175)), for which p-values are estimated numerically rather than using lookup tables.
```RuralUrban-Tests.R``` and ```DirectIndirect-Tests.R``` contain the hypothesis tests for the comparisons of light-polluted and dark skies, and direct and indirect light pollution, respectively.
### Contributions:
All scripts written by JJ Foster. Data collected by JJ Foster, C Tocco, L Khaldy, E Baird, MJ Byrne and M Dacke. All analysis performed by JJ Foster. [ELF](https://doi.org/10.1098/rsif.2021.0184) was developed by J Smolka and DE Nilsson.
## Imaging
Scripts for processing, plotting and analysis of scenes recorded during behavioural experiments. Based the ```night``` and ```polar``` modules of **ELF** developed by Jochen Smolka (https://github.com/sciencedjinn/elf). Polarization images rely on ```polarELF.1.7.1 171108```, originally written for [Foster _et al._ 2019](https://doi.org/10.1242/jeb.188532). Heatmaps of scene radiance rely on [export_fig](https://github.com/altmany/export_fig) by [Yair Altman](https://github.com/altmany). ```RuralUrban-hmax-Tests.R``` and ```DirectIndirect-hmax-Tests.R``` are ```R``` scripts used to calculate the _snapshot index_ from minimum-width at half maximum for rotational image-difference curves.
### Contributions:
All ELF scripts written by J Smolka with conceptual input from D-E Nilsson. ```skypipe.m```, ```rotimdiffsky.m```, ```hmax.m``` and ```colourmapsky.m``` written by JJ Foster. Data collected by JJ Foster, C Tocco, MJ Byrne and M Dacke. All analysis performed by JJ Foster.

# System Requirements
## Hardware requirements
Requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
All analysis was performed on *macOS*: High Sierra (10.13.1) – Catalina (10.15.5)

### R Dependencies

```
circular
CircStats
muStat
beeswarm
```

### MATLAB Dependencies

```
Curve Fitting Toolbox
Image Processing Toolbox
Statistics and Machine Learning Toolbox
```

****
# Installation Guide:

### Install from Github
```
git clone https://github.com/JJFosterLab/light-pollution
```
