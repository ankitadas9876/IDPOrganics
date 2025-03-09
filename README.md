
# Mars IDP Organics Paper

A Study of buried Interplanetary Dust Particle Organics on Mars and Implications on Methane Detection. This repository contains python scripts and data files used to generate the results described in the paper. 


## Documentation

[Documentation](https://linktodocumentation)

Five python scripts to calculate Organic burial at locations: Gale Crater, Arabia Terra, Schiaparelli Crater, Juventae Chasma, and Global Mars.

### Files
1. daytemp.txt - avergage day time temperatures at Gale Crater in Kelvin (Smith et al., 2016)
2. daytemp_arabia.txt - average day time temperatures at Arabia Terra in Kelvin (Source: Mars Climate Database (Forget et al.,
1999; Milllour et al., 2018))

3. daytemp_schiaparelli.txt - average day time temperature at Schiaparelli crater (Mars Climate Database)
4. daytemp_juventae.txt - average day time temperature at Juventae chasma (Mars Climate Database)
5. tempdata.txt - average day temp data at different latitudes on Mars (Mars Climate Database)
6. uvadata.txt, uvbdata.txt, uvcdata.txt - UVA, UVB, and UVC flux on Mars [solar longitude, latitude, flux ] (Source: Smith & Moores 2020 Model)
7. normalizeddust.txt - normalized dust deposition at lat/long in microns (data initially taken from Mars Climate Database https://www-mars.lmd.jussieu.fr/mcd_python/)
8. normalizeddust_kg.txt - normalized dust deposition at lat/long in kg/m2 (Mars Climate Database)

## Equations used for Organic carbon UV photolysis reactions
Refer to Schuerger et al., 2012
Schuerger, A. C., Moores, J. E., Clausen, C. A., Barlow, N. G., & Britt, D. T. (2012). Methane from UV-irradiated carbonaceous chondrites under simulated Martian conditions. Journal of Geophysical Research: Planets, 117(E8), n/a-n/a https://doi.org/10.1029/2011je004023

## References Relevant to the Python scripts
Daytime temperatures at Gale Crater: Smith, M. D., Zorzano, M. P., Lemmon, M., Martín-Torres, J., & Mendaza de Cal, T. (2016). Aerosol optical depth as observed by the Mars Science Laboratory REMS UV photodiodes. Icarus, 280, 234–248. https://doi.org/10.1016/j.icarus.2016.07.012
Mars Climate Database: 
Forget, F., Hourdin, F., Fournier, R., Hourdin, C., Talagrand, O., Collins, M., … Huot, J.-P. (1999). Improved general circulation models of the Martian atmosphere from the surface to above 80 km. Journal of Geophysical Research, 104, 24155–24176. https://doi.org/10.1029/1999JE001025
Millour, E., Forget, F., Spiga, A., Vals, M., Zakharov, V., Montabone, L., … González-Galindo, F. (n.d.). THE MARS CLIMATE DATABASE (VERSION 5.3). Retrieved from https://www.cosmos.esa.int/documents/1499429/1583871/Millour_E.pdf
UV data from Smith and Moores Model:
Smith, C. L., & Moores, J. E. (2020). Modelled small-scale crack orientations in Martian surface clasts caused by differential insolation-mobilized water. Icarus, 338, 113497. https://doi.org/10.1016/j.icarus.2019.113497 
