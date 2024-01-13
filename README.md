# NP_point_sources_model
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL3.0-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Model to estimate Nitrogen (N) and Phosphorus (P) point sources and other N and P emissions from wastewater. 

### How to cite the software

If you would like to use the software, please cite the following publication:

Sarrazin, F. J., Attinger, A., Kumar, R., Gridded dataset of nitrogen and phosphorus point sources from wastewater in Germany (1950-2019), submitted to Earth System Science Data.

DOI for latest sofware version (V1.1): 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10500705.svg)](https://doi.org/10.5281/zenodo.10500705)

### Code description

The model code is contained in two modules:

#### *model_NP_emissions.py*

Module containing functions to estimate domestic and industrial/commercial N and P gross and net emissions from wastewater at the spatial level for which the input data at available (i.e. urban and rural population data, protein data, detergent data, population connection data:

- *domestic_nutrient_gross_emissions*: function to calculate N and P domestic gross human emissions corresponding to human excreta based on protein supply data.

- *domestic_nutrient_net_emissions*:function to calculate domestic N and P net emissions, i.e. it determines the fate of domestic gross emissions (treated and untreated point sources or other fates).

- *industrial_nutrient_gross_emissions*: function to calculates industrial/commercial gross emissions that are collected by the sewer system.

- *industrial_nutrient_net_emissions*: function to calculate industrial/commercial net emissions, i.e. it determines the fate of industrial/commercial gross emissions (treated and untreated point sources or other fates).

- *aggregate_dom_ind_emissions*: function to aggregate domestic and industrial/commercial emissions. 

#### *disaggregation_NP_point_sources_grid.py*

Module containing functions Module to disaggregate N and P emissions from wastewater to grid level:

- *NUTS2grid*: function to estimate the values at grid level of a variable provided at a certain spatial (NUTS) level.

### Example

#### *main_point_sources_N_P_Germany.py*

Main file to estimate the N and P emissions for Germany in the following publication: Sarrazin, F. J., Attinger, A., Kumar, R., Gridded dataset of nitrogen and phosphorus point sources from wastewater in Germany (1950-2019), submitted to Earth System Science Data. 

The necessary input data to run the code are available from: Sarrazin, F. J., Attinger, S., & Kumar, R. (2024). Gridded dataset of nitrogen and phosphorus point sources from wastewater in Germany (1950-2019) (1.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10500535
