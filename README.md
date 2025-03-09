
# Mars IDP Organics Paper

A Study of buried Interplanetary Dust Particle Organics on Mars and Implications on Methane Detection. This repository contains python scripts and data files used to generate the results described in the paper. 


## Documentation

[Documentation](https://linktodocumentation)

Five python scripts to calculate Organic burial at locations: Gale Crater, Arabia Terra, Schiaparelli Crater, Juventae Chasma, and Global Mars.

### Files
1. daytemp.txt - avergage day time temperatures at Gale Crater in Kelvin (Smith et al., 2016)
2. daytemp_arabia.txt - average day time temperatures at Arabia Terra in Kelvin (Source: Mars Climate Database (Forget et al.,
1999; Madeleine et al. 2011; Milllour et al., 2018))

3. daytemp_schiaparelli.txt - average day time temperature at Schiaparelli crater (Mars Climate Database)
4. daytemp_juventae.txt - average day time temperature at Juventae chasma (Mars Climate Database)
5. tempdata.txt - average day temp data at different latitudes on Mars (Mars Climate Database)
6. uvadata.txt, uvbdata.txt, uvcdata.txt - UVA, UVB, and UVC flux on Mars [solar longitude, latitude, flux ] (Source: Smith & Moores 2020 Model)
7. normalizeddust.txt - normalized dust deposition at lat/long in microns (data initially taken from Mars Climate Database https://www-mars.lmd.jussieu.fr/mcd_python/)
8. normalizeddust_kg.txt - normalized dust deposition at lat/long in kg/m2 (Mars Climate Database)

## Equations used for Organic carbon UV photolysis reactions
Refer to Schuerger et al., 2012
Schuerger, A. C., Moores, J. E., Clausen, C. A., Barlow, N. G., & Britt, D. T. (2012). Methane from UV-irradiated carbonaceous chondrites under simulated Martian conditions. Journal of Geophysical Research: Planets, 117(E8), n/a-n/a https://doi.org/10.1029/2011je004023


