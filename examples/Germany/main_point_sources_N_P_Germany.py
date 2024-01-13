
"""

This script provides an application of the N and P point sources model over
Germany.

LICENSE:
    
    Copyright 2024, Fanny J. Sarrazin, Rohini Kumar
    All rights reserved.
    
    This source code is licensed under the GNU General Public License v3.0 
    license  found in the LICENSE file in the root directory of this source tree. 


PLEASE CITE THIS PUBLICATION IF YOU USE THIS CODE:
    
    Sarrazin, F. J., Attinger, A., Kumar, R., Gridded dataset of nitrogen and 
    phosphorus point sources from wastewater in Germany (1950-2019), submitted 
    to Earth System Science Data.


INPUT DATA TO RUN THE CODE:

    The necessary input data to run the code are available from: 
        
    Sarrazin, F. J., Attinger, S., & Kumar, R. (2024). Gridded dataset of nitrogen 
    and phosphorus point sources from wastewater in Germany (1950-2019) (1.1) 
    [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10500535


TABLE OF CONTENTS:

    1) Load data at NUTS-1 level
        1.1) Load parameter sample
        1.2) Load input data
        1.3) Load calibration data
        1.4) Prepare indices for NUTS-1 regions of Berlin and Brandenburg
    
    2) Load data at grid level
        2.1) Load population data
        2.2) Load NUTS-1 map
    
    3) Calculate domestic and industrial/commercial emissions at NUTS-1 level 
       for urban and rural population
       
    4) Aggregate domestic and industrial/commercial emissions at NUTS-1 level  
       for urban and rural areas
       
    5) Process emissions at NUTS-1 level for Berlin and Brandenburg
    
    6) Estimate the model parameters at NUTS-1 level
    
    7) Downscale the N and P point sources emissions from NUTS-1 to grid level 
       using gridded population data

"""

import os
import numpy as np
from netCDF4 import Dataset

# Set directory to emission model 
# (UPDATE WITH YOUR OWN PATH)
os.chdir('PATH_TO_MODEL') 

from model_NP_emissions_nuts import domestic_nutrient_gross_emissions,\
    domestic_nutrient_net_emissions, industrial_nutrient_gross_emissions,\
        industrial_nutrient_net_emissions,aggregate_dom_ind_emissions

from disaggregation_NP_point_sources_grid import NUTS2grid
        

# Download input data from a zenodo repository:
# https://doi.org/10.5281/zenodo.10500535

# Set directory to downloaded input data
# (UPDATE WITH YOUR OWN PATH)
dir_data = 'PATH_TO_DATA'

# Data files at NUTS-1 level
# parameter sample
parameter_file = dir_data + r'parameter_sample_v1.0.csv'  
# input data            
input_data_file = dir_data + r'data_input_nuts1_v1.0.csv'  
# calibration data          
calibration_data_file = dir_data + r'data_calibration_nuts1_v1.0.csv' 

# Data files at grid level
 # population data
population_data_file = dir_data + r'population_count_v1.0.nc'  
#  NUTS-1 map      
nuts1_map_file = dir_data + r'nuts1_map_v1.0.nc' 

# no data value
no_val = -9999

#%% ---------------------------------------------------------------------------
# 1) Load data at NUTS-1 level
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1.1) Load parameter sample
# -----------------------------------------------------------------------------

Xpar = np.loadtxt(parameter_file, delimiter=',', skiprows=2)

# Fraction of protein supply wasted at the distribution and consumption level (-) 
f_pro_waste = Xpar[:, 1]
# N content in protein (kg N/kg) 
f_N_pro = Xpar[:, 2]
# N:P ratio for human intake (kg N/kg P) 
f_N_P = Xpar[:, 3]
# Fraction of nutrient intake lost via sweat, hair and blood (-) 
f_lossHum = Xpar[:, 4]
# Fraction of N and P emissions lost during wastewater collection and transport 
# (losses in the sewer system, and in cesspits and during transportation of 
# wastewater collected in cesspits to wwtps) (-)
f_lossTransport = Xpar[:, 5]
# Efficiency of N removal for primary treatment (-) 
eff_1_N = Xpar[:, 6]
# Efficiency of N removal for secondary treatment and tertiary treatment without 
# targeted N removal (-) 
eff_23noN_N = Xpar[:, 7]
# Efficiency of N removal for tertiary treatment with targeted N removal (-)
eff_3_N = Xpar[:, 8]
# Efficiency of P removal for primary treatment (-)
eff_1_P = Xpar[:, 9]
# Efficiency of P removal for secondary treatment and tertiary treatment without targeted P removal (-)
eff_23noP_P = Xpar[:, 10]
# Efficiency of P removal for tertiary treatment with targeted P removal (-) 
eff_3_P = Xpar[:, 11]
# Fraction of industrial/commercial (indirect release) to domestic human 
# (physiological) N gross emissions in 1950 (-)
f_ind_dom_1950_N = Xpar[:, 12]
# Fraction of industrial/commercial (indirect release) to domestic human 
# (physiological) N gross emissions for the period 2000–2019 (-)
f_ind_dom_2000_N = Xpar[:, 13]
# Fraction of industrial/commercial (indirect release) to domestic human 
# (physiological) P gross emissions in 1950 (-)  
f_ind_dom_1950_P = Xpar[:, 14]
# Fraction of industrial/commercial (indirect release) to domestic human 
# (physiological) P gross emissions for the period 2000–2019 (-)
f_ind_dom_2000_P = Xpar[:, 15]
# Fraction of industrial/commercial (indirect release) to domestic dishwasher 
# detergent phosphate use for the period 1950-2012 (-)
f_DD_phosphate_ind_dom = Xpar[:, 16]
# Fraction of industrial/commercial (indirect release) to domestic detergent 
# phosphonate use for the period 1950-2015(-)
f_D_phosphonate_ind_dom = Xpar[:, 17]
# Industrial/commercial (indirect release) laundry booster (kg P/capita/yr)  
LD_phosphate_ind_cap = Xpar[:, 18]

# P intake as a fraction of protein supply intake
f_P_pro = f_N_pro/f_N_P

# Number of parameter samples
Npar = len(f_pro_waste)
# Number of parameters
M = Xpar.shape[1]
# years for f_ind_dom
years_f_ind = np.array([1950, 2000])

# -----------------------------------------------------------------------------
# 1.2) Load input data at NUTS-1 level
# -----------------------------------------------------------------------------

# Initialise population connection variable (for tot, urban and rural)
data = np.loadtxt(input_data_file, delimiter=',', skiprows=2)

# id of the NUTS-1 regions
# 1:   Schleswig-Holstein
# 2:   Hamburg
# 3:   Lower Saxony 
# 4:   Bremen
# 5:   North Rhine-Westphalia
# 6:   Hesse
# 7:   Rhineland-Palatinate
# 8:   Baden-Württemberg
# 9:   Bavaria
# 10:  Saarland
# 11:  Berlin
# 12:  Brandenburg
# 13:  Mecklenburg-Vorpommern
# 14:  Saxony
# 15:  Saxony-Anhalt
# 16:  Thuringia
# 111: Part of the population of Berlin whose wastewater is handled in Berlin
# 112: Part of the population of Berlin whose wastewater is handled in Brandenburg

id_nuts1 = data[:, 0]
# years
years_nuts1 = data[:, 1]
# urban and rural population count
pop_nuts1 = data[:, 4:6]
# protein supply per capita at the distribution level (kg/capita/yr)
protein_nuts1 = data[:, 6]
# laundry detergent phosphate P use per capita (kg P/capita/yr) corresponding to: 
# (1) total (both domestic and industrial/commercial) use until 1990 for the
#  NUTS-1 regions of West Germany and until 1991 for the NUTS-1 regions of East 
# Germany 
# (2) domestic use only from 1991 for West Germany and from 1992 for East Germany.
LD_phosphate_tot_P_nuts1 = data[:, 7]
# domestic dishwasher detergent phosphate P use per capita (kg P/capita/yr)
DD_phosphate_dom_P_nuts1 = data[:, 8]
# domestic detergent phosphonate P use per capita (kg P/capita/yr)
D_phosphonate_dom_P_nuts1 = data[:,9]
# percentage of urban and rural population connected to sewage farms (%)
TFarm_nuts1 = data[:, 11:13]/100
# percentage of urban and rural population connected to the sewer system but not 
# to wastewater treatment plants (WWTPs) (%)
T0_nuts1 = data[:, 14:16]/100
# percentage of total population connected to primary (mechanical) treatment 
T1_nuts1 = data[:, 17:19]/100
# percentage of urban and rural population connected to secondary (biological) treatment 
# and tertiary (advanced) treatment without targeted N removal (%)
T23_noN_nuts1 = data[:, 20:22]/100
# percentage of urban and rural population connected to secondary (biological) treatment 
# and tertiary (advanced) treatment without targeted P removal (%)
T23_noP_nuts1 = data[:, 23:25]/100
# percentage of urban and rural population connected to tertiary treatment with 
# targeted N removal (%)
T3_N_nuts1 = data[:, 26:28]/100
# percentage of urban and rural population connected to tertiary treatment with 
# targeted P removal (%)
T3_P_nuts1= data[:, 29:31]/100

# Number of data records
N_data = len(id_nuts1)

# time period
years = np.unique(years_nuts1)
Ny = len(years)

#id of NUTS-1 regions 
id_nuts1_unique = np.unique(id_nuts1)
# number of NUTS-1 regions
Nnuts1 = len(id_nuts1_unique)

# id of states of Western Germany and Berlin
id_west = np.concatenate((np.arange(1, 12), np.array([111, 112])))

# indices of records for Western Germany and Berlin (this will be used to 
# to prepare the industrial/commercial detergent data)
idx_west = np.full((N_data, ), False)
for ii in range(N_data):
    idx_west[ii] = np.any(id_nuts1[ii]==id_west) 

# Urban population connection to sewer or wwtp for Berlin
# (This is used to disaggregate Berlin gross emissions into the part of Berlin 
# handled in Berlin and in Brandenburg)
T_sewerWwtp_urb = T0_nuts1[:, 0] + TFarm_nuts1[:, 0] + T1_nuts1[:, 0] +\
                  T23_noN_nuts1[:, 0] + T3_N_nuts1[:, 0]   
                  
# -----------------------------------------------------------------------------
# 1.3) Load calibration data at NUTS-1 level
# -----------------------------------------------------------------------------

data_calib = np.loadtxt(calibration_data_file, delimiter=',', skiprows=2)
# Set no data value to nan
data_calib[data_calib==no_val] = np.nan

# WWTP outgoing N load
J_N_psWwtpOut_nuts1_obs_tmp = data_calib[:, 2:4]
# WWTP outgoingP load
J_P_psWwtpOut_nuts1_obs_tmp = data_calib[:, 4:6]
# Year variable
years_obs = np.unique(data_calib[:, 1])
Ny_obs = len(years_obs)

# Change the years so that they are consistent with the simulations (1950-2019)
# Indentify first year with observations
idx_year_obs_ini = np.where(years_obs[0] == years)[0][0]
# Initialise variables with nan
J_N_psWwtpOut_nuts1_obs = np.nan*np.zeros((Nnuts1-2, Ny, 2))
J_P_psWwtpOut_nuts1_obs = np.nan*np.zeros((Nnuts1-2, Ny, 2))
# Set values
J_N_psWwtpOut_nuts1_obs[:, idx_year_obs_ini:Ny, 0:2] =\
    np.reshape(J_N_psWwtpOut_nuts1_obs_tmp, (Nnuts1-2, Ny_obs, 2))[:, 0:Ny_obs-1, :]
J_P_psWwtpOut_nuts1_obs[:, idx_year_obs_ini:Ny, 0:2] =\
    np.reshape(J_P_psWwtpOut_nuts1_obs_tmp, (Nnuts1-2, Ny_obs, 2))[:, 0:Ny_obs-1, :]
    
# -----------------------------------------------------------------------------
# 1.4) Prepare indices for NTS-1 regions of Berlin and Brandenburg
# -----------------------------------------------------------------------------

# Indices of Berlin
idx_BE = np.where(id_nuts1==11)[0]
# Indices of part of Berlin handled in Berlin
idx_BEinBE = np.where(id_nuts1==111)[0]
# Indices of part of Berlin handled in Brandenburg
idx_BEinBB = np.where(id_nuts1==112)[0]
# Indices of Brandenburg
idx_BB = np.where(id_nuts1==12)[0]

# Indices of all states but part of Berlin handled in Berlin and Brandenburg
idx_nuts1_16 = np.arange(0, Ny*(Nnuts1-2))

#%% -----------------------------------------------------------------------------
# 2) Load data at grid level
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2.1) Load population data at grid level
# -----------------------------------------------------------------------------

# Load netcdf file
nc_pop_tot = Dataset(population_data_file)
# urban population count
pop_urb = np.ma.filled(nc_pop_tot.variables['pop_urban_count'][:], np.nan)
# rural population count
pop_rur = np.ma.filled(nc_pop_tot.variables['pop_rural_count'][:], np.nan)
# time variable
time_pop = np.ma.filled(nc_pop_tot.variables['time'][:], np.nan)
# latitude
lat_pop = np.ma.filled(nc_pop_tot.variables['lat'][:], np.nan)
# logitude
lon_pop = np.ma.filled(nc_pop_tot.variables['lon'][:], np.nan)
# Close netcdf file
nc_pop_tot.close()


# -----------------------------------------------------------------------------
# 2.1) Load NUTS-1 map at grid level
# -----------------------------------------------------------------------------

# Load netcdf file
nc_nuts1 = Dataset(nuts1_map_file)
# NUTS-1 map
nuts1_grid = np.ma.filled(nc_nuts1.variables['nuts1_map'][:], np.nan)
# Close netcdf file
nc_nuts1.close()


#%% -----------------------------------------------------------------------------
# 3) Calculate domestic and industrial/commercial emissions at NUTS-1 level 
#  (for urban and rural population)
# -----------------------------------------------------------------------------

# Number of population groups (urban, rural)
Npop = 2

# Initialize variables

# N domestic gross emissions (human excreta) (kg/yr)
J_N_gross_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop)) 
# N WWTP outgoing domestic emissions (point sources) (kg/yr)
J_N_psWwtpOut_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions not treated (kg/yr)
J_N_psNoTreat_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N WWTP incoming domestic emissions (kg/yr) 
J_N_wwtpIn_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions treated in sewage farms (kg/yr) (incoming emissions to 
# sewage farms without accounting for treatment removal)
J_N_farm_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions not collected by the sewer system nor treated in wwtp (kg/yr)
# (these emissions are domestic only, as all industrial/commercial emissions
# are collected by the sewer system or treated in wwtps)
J_N_noSewerWwtp_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions lost during collection and transport  (kg/yr)
J_N_transportLoss_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions removed in WWTPs (kg/yr)     
J_N_wwtpRem_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))

# P domestic gross emissions (human excreta, excluding detergents) (kg/yr)
J_P_gross_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop)) 
# P domestic gross detergent emissions
J_P_gross_det_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# P WWTP outgoing domestic emissions (point sources) (kg/yr)
J_P_psWwtpOut_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# P domestic emissions not treated (kg/yr)
J_P_psNoTreat_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# P WWTP incoming domestic emissions (kg/yr) 
J_P_wwtpIn_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# P domestic emissions treated in sewage farms (kg/yr) (incoming emissions to 
# sewage farms without accounting for treatment removal)
J_P_farm_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# P domestic emissions not collected by the sewer system nor treated in wwtps (kg/yr)
# (these emissions are domestic only, as all industrial/commercial emissions
# are collected by the sewer system or treated in wwtps)
J_P_noSewerWwtp_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions lost during collection and transport  (kg/yr)
J_P_transportLoss_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))
# N domestic emissions removed in WWTPs (kg/yr)     
J_P_wwtpRem_dom_nuts1 = np.nan * np.zeros((N_data, Npar, Npop))

# N industrial/commercial gross emissions (kg/yr)
J_N_gross_ind_nuts1 = np.nan * np.zeros((N_data, Npar)) 
# N WWTP outgoing industrial/commercial emissions (point sources) (kg/yr)
J_N_psWwtpOut_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercial emissions not treated (kg/yr)
J_N_psNoTreat_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N WWTP incoming industrial/commercial emissions (kg/yr) 
J_N_wwtpIn_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercialemissions treated in sewage farms (kg/yr) (incoming 
#  emissions to sewage farms without accounting for treatment removal)
J_N_farm_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercial  emissions lost during collection and transport  (kg/yr)
J_N_transportLoss_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercial emissions removed in WWTPs (kg/yr)     
J_N_wwtpRem_ind_nuts1 = np.nan * np.zeros((N_data, Npar))

# P industrial/commercial gross emissions (excluding detergents) (kg/yr)
J_P_gross_ind_nuts1 = np.nan * np.zeros((N_data, Npar)) 
# P industrial/commercial gross detergent emissions
J_P_gross_det_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# P WWTP outgoing industrial/commercial emissions (point sources) (kg/yr)
J_P_psWwtpOut_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercial emissions not treated (kg/yr)
J_P_psNoTreat_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# P WWTP incoming industrial/commercial emissions (kg/yr) 
J_P_wwtpIn_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercialemissions treated in sewage farms (kg/yr) (incoming 
#  emissions to sewage farms without accounting for treatment removal)
J_P_farm_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# N industrial/commercial  emissions lost during collection and transport  (kg/yr)
J_P_transportLoss_ind_nuts1 = np.nan * np.zeros((N_data, Npar))
# P industrial/commercial emissions removed in WWTPs (kg/yr)     
J_P_wwtpRem_ind_nuts1 = np.nan * np.zeros((N_data, Npar))

for nn in range(Npar): # loop on the parameter sets
            
    for jj in range(Npop): # loop on the population groups (urban, rural)
    
        #------------------------------ domestic N emissions----------------------
        # Calculate domestic gross N emissions
        J_N_gross_dom_nuts1[:, nn, jj] =\
            domestic_nutrient_gross_emissions(np.array([f_pro_waste[nn], f_N_pro[nn], f_lossHum[nn]]),\
                                              pop_nuts1[:, jj], protein_nuts1)    
        # Calculate domestic net N emissions
        J_N_psWwtpOut_dom_nuts1[:, nn, jj],\
        J_N_psNoTreat_dom_nuts1[:, nn, jj],\
        J_N_wwtpIn_dom_nuts1[:, nn, jj],\
        J_N_farm_dom_nuts1[:, nn, jj],\
        J_N_noSewerWwtp_nuts1[:, nn, jj],\
        J_N_transportLoss_dom_nuts1[:, nn, jj],\
        J_N_wwtpRem_dom_nuts1[:, nn, jj] =\
            domestic_nutrient_net_emissions(np.array([f_lossTransport[nn], eff_1_N[nn],\
                                                      eff_23noN_N[nn], eff_3_N[nn]]),\
                                            J_N_gross_dom_nuts1[:, nn, jj], T0_nuts1[:, jj],\
                                            TFarm_nuts1[:, jj], T1_nuts1[:, jj],\
                                            T23_noN_nuts1[:, jj], T3_N_nuts1[:, jj])
                
        #------------------------------ domestic P emissions----------------------
        # Calculate domestic gross P emissions
        J_P_gross_dom_nuts1[:, nn, jj] =\
            domestic_nutrient_gross_emissions(np.array([f_pro_waste[nn], f_P_pro[nn], f_lossHum[nn]]),\
                                              pop_nuts1[:, jj], protein_nuts1)
        
        # Calculate domestic gross detergent P emissions
        J_P_gross_det_dom_nuts1[:, nn, jj] =\
            (LD_phosphate_tot_P_nuts1 + DD_phosphate_dom_P_nuts1 +\
             D_phosphonate_dom_P_nuts1)*pop_nuts1[:, jj]
         
        # Calculate domestic net P emissions
        J_P_psWwtpOut_dom_nuts1[:, nn, jj],\
        J_P_psNoTreat_dom_nuts1[:, nn, jj],\
        J_P_wwtpIn_dom_nuts1[:, nn, jj],\
        J_P_farm_dom_nuts1[:, nn, jj],\
        J_P_noSewerWwtp_nuts1[:, nn, jj],\
        J_P_transportLoss_dom_nuts1[:, nn, jj],\
        J_P_wwtpRem_dom_nuts1[:, nn, jj] =\
            domestic_nutrient_net_emissions(np.array([f_lossTransport[nn], eff_1_P[nn],\
                                                      eff_23noP_P[nn], eff_3_P[nn]]),\
                                            J_P_gross_dom_nuts1[:, nn, jj] + J_P_gross_det_dom_nuts1[:, nn, jj] ,\
                                            T0_nuts1[:, jj],  TFarm_nuts1[:, jj], T1_nuts1[:, jj],\
                                            T23_noP_nuts1[:, jj], T3_P_nuts1[:, jj])
           
    #---------------------Industrial/commercial detergents---------------------
    # Industrial/commercial laundry boosters
    LD_boost_phosphate_ind_P_nuts1 = np.zeros((N_data, ))
    # Consider laundry boosters from 1992 onwards
    LD_boost_phosphate_ind_P_nuts1[years_nuts1>=1992] = LD_phosphate_ind_cap[nn]
    # for states of West Germany also consider laundry boosters in 1991
    LD_boost_phosphate_ind_P_nuts1[np.logical_and(years_nuts1==1991, idx_west)] =\
                                                      LD_phosphate_ind_cap[nn]
    
    # Industrial/commercial dishwasher phosphate detergents
    DD_phosphate_ind_P_nuts1 = np.nan * np.zeros((N_data, ))
    # Fraction of domestic emissions up to 2012
    DD_phosphate_ind_P_nuts1[years_nuts1<=2012] =\
        DD_phosphate_dom_P_nuts1[years_nuts1<=2012]*f_DD_phosphate_ind_dom[nn]
    # Constant value from 2012
    DD_phosphate_ind_P_nuts1[years_nuts1>2012] =\
        np.tile(np.expand_dims(DD_phosphate_ind_P_nuts1[years_nuts1==2012], axis=1),\
                                                           (2019-2012)).flatten()
    
    # Industrial/commercial phosphonate detergents
    D_phosphonate_ind_P_nuts1 = np.nan * np.zeros((N_data, ))
    # Fraction of domestic emissions up to 2012
    D_phosphonate_ind_P_nuts1[years_nuts1<=2015] =\
        D_phosphonate_dom_P_nuts1[years_nuts1<=2015] *f_D_phosphonate_ind_dom[nn]
    # Constant value from 2015
    D_phosphonate_ind_P_nuts1[years_nuts1>2015] =\
        np.tile(np.expand_dims(D_phosphonate_ind_P_nuts1[years_nuts1==2015], axis=1),\
                                                           (2019-2015)).flatten()
     
    #---------------------------- industrial N emissions-----------------------    
    # Assess fraction of industrial/commercial to domestic emissions for each record
    # (the value is kept constant after 2000)
    f_ind_dom_interp_N = np.interp(years, years_f_ind, np.array([f_ind_dom_1950_N[nn],\
                                                                 f_ind_dom_2000_N[nn]]))   
    f_ind_dom_N = np.tile(f_ind_dom_interp_N, (Nnuts1))
    
    # Calculate industrial gross N emissions (use total population)
    J_N_gross_ind_nuts1[:, nn] = \
        industrial_nutrient_gross_emissions(np.array([f_pro_waste[nn], f_N_pro[nn], f_lossHum[nn]]),\
                                        pop_nuts1[:, 0]+pop_nuts1[:, 1], protein_nuts1, f_ind_dom_N)
                    
    # For Berlin we distribute the industrial gross emissions to the urban population 
    # connected to sewer or wwtp for the population of Berlin whose wastewater 
    # is handled in Berlin and Brandenburg
    # Population of Berlin whose wastewater is handled in Berlin
    J_N_gross_ind_nuts1[idx_BEinBE, nn] = J_N_gross_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBE]*pop_nuts1[idx_BEinBE, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])
    # Population of Berlin whose wastewater is handled in Brandenburg
    J_N_gross_ind_nuts1[idx_BEinBB, nn] = J_N_gross_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBB]*pop_nuts1[idx_BEinBB, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])

    # Calculate industrial net N emissions (use connection for urban population)
    J_N_psWwtpOut_ind_nuts1[:, nn],\
    J_N_psNoTreat_ind_nuts1[:, nn],\
    J_N_wwtpIn_ind_nuts1[:, nn],\
    J_N_farm_ind_nuts1[:, nn],\
    J_N_transportLoss_ind_nuts1[:, nn],\
    J_N_wwtpRem_ind_nuts1[:, nn] =\
    industrial_nutrient_net_emissions(np.array([f_lossTransport[nn], eff_1_N[nn],\
                                              eff_23noN_N[nn], eff_3_N[nn]]),\
                                      J_N_gross_ind_nuts1[:, nn], T0_nuts1[:, 0],\
                                      TFarm_nuts1[:, 0], T1_nuts1[:, 0],\
                                      T23_noN_nuts1[:, 0], T3_N_nuts1[:, 0])
                
    #---------------------------- industrial P emissions-----------------------
    # Assess fraction of industrial/commercial to domestic emissions for each record
    f_ind_dom_interp_P = np.interp(years, years_f_ind, np.array([f_ind_dom_1950_P[nn],\
                                                                 f_ind_dom_2000_P[nn]]))   
    f_ind_dom_P = np.tile(f_ind_dom_interp_P, (Nnuts1))
    
    # Calculate industrial gross NP emissions (use total population)
    J_P_gross_ind_nuts1[:, nn] = \
       industrial_nutrient_gross_emissions(np.array([f_pro_waste[nn], f_P_pro[nn], f_lossHum[nn]]),\
                                       pop_nuts1[:, 0]+pop_nuts1[:, 1], protein_nuts1, f_ind_dom_P)
    J_P_gross_det_ind_nuts1[:, nn] =\
       (LD_boost_phosphate_ind_P_nuts1 +\
         DD_phosphate_ind_P_nuts1 + D_phosphonate_ind_P_nuts1)*(pop_nuts1[:, 0]+pop_nuts1[:, 1])
                   
    # For Berlin we distribute the industrial gross emissions over the part treated
    # in Berlin and the part treated in Brandenburg using the urban population
    # Population of Berlin whose wastewater is handled in Berlin 
    J_P_gross_ind_nuts1[idx_BEinBE, nn] = J_P_gross_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBE]*pop_nuts1[idx_BEinBE, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])
    J_P_gross_ind_nuts1[idx_BEinBB, nn] = J_P_gross_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBB]*pop_nuts1[idx_BEinBB, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])
    # Population of Berlin whose wastewater is handled in Brandenburg
    J_P_gross_det_ind_nuts1[idx_BEinBE, nn] = J_P_gross_det_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBE]*pop_nuts1[idx_BEinBE, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])
    J_P_gross_det_ind_nuts1[idx_BEinBB, nn] = J_P_gross_det_ind_nuts1[idx_BE, nn] *\
        (T_sewerWwtp_urb[idx_BEinBB]*pop_nuts1[idx_BEinBB, 0])/\
            (T_sewerWwtp_urb[idx_BE]*pop_nuts1[idx_BE, 0])
    
    # Calculate industrial net P emissions (use connection for urban population)
    J_P_psWwtpOut_ind_nuts1[:, nn],\
    J_P_psNoTreat_ind_nuts1[:, nn],\
    J_P_wwtpIn_ind_nuts1[:, nn],\
    J_P_farm_ind_nuts1[:, nn],\
    J_P_transportLoss_ind_nuts1[:, nn],\
    J_P_wwtpRem_ind_nuts1[:, nn] =\
    industrial_nutrient_net_emissions(np.array([f_lossTransport[nn], eff_1_P[nn],\
                                              eff_23noP_P[nn], eff_3_P[nn]]),\
                                      J_P_gross_ind_nuts1[:, nn] + J_P_gross_det_ind_nuts1[:, nn],\
                                      T0_nuts1[:, 0], TFarm_nuts1[:, 0], T1_nuts1[:, 0],\
                                      T23_noP_nuts1[:, 0], T3_P_nuts1[:, 0])
        
#%% ---------------------------------------------------------------------------       
# 4) Aggregate domestic and industrial/commercial emissions at NUTS-1 level for 
#  urban and rural areas
# -----------------------------------------------------------------------------

J_N_gross_nuts1 =\
    aggregate_dom_ind_emissions(J_N_gross_dom_nuts1, J_N_gross_ind_nuts1,\
                                'J_N_gross', axis_pop=2, idx_urb=0, idx_rur=1)
J_P_gross_nuts1 = \
    aggregate_dom_ind_emissions(J_P_gross_dom_nuts1 + J_P_gross_det_dom_nuts1,\
                                          J_P_gross_ind_nuts1 + J_P_gross_det_ind_nuts1,\
                                              'J_P_gross', axis_pop=2, idx_urb=0, idx_rur=1)
J_N_psNoTreat_nuts1 =\
    aggregate_dom_ind_emissions(J_N_psNoTreat_dom_nuts1,J_N_psNoTreat_ind_nuts1,\
                                'J_N_psNoTreat', axis_pop=2, idx_urb=0, idx_rur=1) 
J_P_psNoTreat_nuts1 =\
    aggregate_dom_ind_emissions(J_P_psNoTreat_dom_nuts1, J_P_psNoTreat_ind_nuts1,\
                                'J_P_psNoTreat', axis_pop=2, idx_urb=0, idx_rur=1)

J_N_psWwtpOut_nuts1 =\
    aggregate_dom_ind_emissions(J_N_psWwtpOut_dom_nuts1, J_N_psWwtpOut_ind_nuts1,\
                                'J_N_psWwtpOut', axis_pop=2, idx_urb=0, idx_rur=1)
J_P_psWwtpOut_nuts1 =\
    aggregate_dom_ind_emissions(J_P_psWwtpOut_dom_nuts1, J_P_psWwtpOut_ind_nuts1,\
                                'J_P_psWwtpOut', axis_pop=2, idx_urb=0, idx_rur=1)

J_N_wwtpIn_nuts1 =\
    aggregate_dom_ind_emissions(J_N_wwtpIn_dom_nuts1, J_N_wwtpIn_ind_nuts1,\
                                'J_N_wwtpIn', axis_pop=2, idx_urb=0, idx_rur=1) 
J_P_wwtpIn_nuts1 =\
    aggregate_dom_ind_emissions(J_P_wwtpIn_dom_nuts1, J_P_wwtpIn_ind_nuts1,\
                                'J_P_wwtpIn', axis_pop=2, idx_urb=0, idx_rur=1)

J_N_transportLoss_nuts1 =\
    aggregate_dom_ind_emissions(J_N_transportLoss_dom_nuts1, J_N_transportLoss_ind_nuts1,\
                                'J_N_transportLoss', axis_pop=2, idx_urb=0, idx_rur=1)
J_P_transportLoss_nuts1 =\
    aggregate_dom_ind_emissions(J_P_transportLoss_dom_nuts1, J_P_transportLoss_ind_nuts1,\
                                'J_P_transportLoss', axis_pop=2, idx_urb=0, idx_rur=1)

J_N_wwtpRem_nuts1 =\
    aggregate_dom_ind_emissions(J_N_wwtpRem_dom_nuts1, J_N_wwtpRem_ind_nuts1,\
                                'J_N_wwtpRem', axis_pop=2, idx_urb=0, idx_rur=1)
J_P_wwtpRem_nuts1 =\
    aggregate_dom_ind_emissions(J_P_wwtpRem_dom_nuts1, J_P_wwtpRem_ind_nuts1,\
                                'J_P_wwtpRem', axis_pop=2, idx_urb=0, idx_rur=1)

J_N_farm_nuts1 =\
    aggregate_dom_ind_emissions(J_N_farm_dom_nuts1, J_N_farm_ind_nuts1,\
                                             'J_N_farm', axis_pop=2, idx_urb=0, idx_rur=1) 
J_P_farm_nuts1 =\
    aggregate_dom_ind_emissions(J_P_farm_dom_nuts1, J_P_farm_ind_nuts1,\
                                             'J_P_farm', axis_pop=2, idx_urb=0, idx_rur=1)

# Delete variables to free memory
del J_N_gross_dom_nuts1
del J_P_gross_dom_nuts1
del J_P_gross_det_dom_nuts1
del J_N_gross_ind_nuts1
del J_P_gross_ind_nuts1
del J_P_gross_det_ind_nuts1

del J_N_psNoTreat_dom_nuts1
del J_N_psNoTreat_ind_nuts1
del J_P_psNoTreat_dom_nuts1   
del J_P_psNoTreat_ind_nuts1  

del J_N_wwtpIn_dom_nuts1
del J_N_wwtpIn_ind_nuts1
del J_P_wwtpIn_dom_nuts1
del J_P_wwtpIn_ind_nuts1

del J_N_psWwtpOut_dom_nuts1
del J_N_psWwtpOut_ind_nuts1
del J_P_psWwtpOut_dom_nuts1
del J_P_psWwtpOut_ind_nuts1

del J_N_farm_dom_nuts1
del J_N_farm_ind_nuts1
del J_P_farm_dom_nuts1
del J_P_farm_ind_nuts1

del J_N_transportLoss_dom_nuts1
del J_N_transportLoss_ind_nuts1
del J_P_transportLoss_dom_nuts1
del J_P_transportLoss_ind_nuts1

#%% ---------------------------------------------------------------------------       
# 5) Process emissions at NUTS-1 level for Berlin and Brandenburg
# -----------------------------------------------------------------------------

# Add part of Berlin handled in Brandenburg to Brandenburg and set Berlin to
# the part of Berlin handled in Berlin. 
# Then remove the part of Berlin handled in Berlin and Brandenburg from the data

J_N_gross_nuts1[idx_BB, :] = J_N_gross_nuts1[idx_BB, :] +\
                             J_N_gross_nuts1[idx_BEinBB, :]
J_N_gross_nuts1[idx_BE, :] = J_N_gross_nuts1[idx_BEinBE, :]
J_N_gross_nuts1 = J_N_gross_nuts1[idx_nuts1_16, :]

J_P_gross_nuts1[idx_BB, :] = J_P_gross_nuts1[idx_BB, :] +\
                             J_P_gross_nuts1[idx_BEinBB, :]
J_P_gross_nuts1[idx_BE, :] = J_P_gross_nuts1[idx_BEinBE, :]
J_P_gross_nuts1 = J_P_gross_nuts1[idx_nuts1_16, :]

J_N_psNoTreat_nuts1[idx_BB, :] = J_N_psNoTreat_nuts1[idx_BB, :] +\
                                 J_N_psNoTreat_nuts1[idx_BEinBB, :]
J_N_psNoTreat_nuts1[idx_BE, :] = J_N_psNoTreat_nuts1[idx_BEinBE, :]
J_N_psNoTreat_nuts1 = J_N_psNoTreat_nuts1[idx_nuts1_16, :]

J_P_psNoTreat_nuts1[idx_BB, :] = J_P_psNoTreat_nuts1[idx_BB, :] +\
                                 J_P_psNoTreat_nuts1[idx_BEinBB, :]
J_P_psNoTreat_nuts1[idx_BE, :] = J_P_psNoTreat_nuts1[idx_BEinBE, :]
J_P_psNoTreat_nuts1 = J_P_psNoTreat_nuts1[idx_nuts1_16, :]

J_N_psWwtpOut_nuts1[idx_BB, :] = J_N_psWwtpOut_nuts1[idx_BB, :] +\
                                 J_N_psWwtpOut_nuts1[idx_BEinBB, :]
J_N_psWwtpOut_nuts1[idx_BE, :] = J_N_psWwtpOut_nuts1[idx_BEinBE, :]
J_N_psWwtpOut_nuts1 = J_N_psWwtpOut_nuts1[idx_nuts1_16, :]

J_P_psWwtpOut_nuts1[idx_BB, :] = J_P_psWwtpOut_nuts1[idx_BB, :] +\
                                 J_P_psWwtpOut_nuts1[idx_BEinBB, :]
J_P_psWwtpOut_nuts1[idx_BE, :] = J_P_psWwtpOut_nuts1[idx_BEinBE, :]
J_P_psWwtpOut_nuts1 = J_P_psWwtpOut_nuts1[idx_nuts1_16, :]

J_N_wwtpIn_nuts1[idx_BB, :] = J_N_wwtpIn_nuts1[idx_BB, :] +\
                              J_N_wwtpIn_nuts1[idx_BEinBB, :]
J_N_wwtpIn_nuts1[idx_BE, :] = J_N_wwtpIn_nuts1[idx_BEinBE, :]
J_N_wwtpIn_nuts1 = J_N_wwtpIn_nuts1[idx_nuts1_16, :]

J_P_wwtpIn_nuts1[idx_BB, :] = J_P_wwtpIn_nuts1[idx_BB, :] +\
                              J_P_wwtpIn_nuts1[idx_BEinBB, :]
J_P_wwtpIn_nuts1[idx_BE, :] = J_P_wwtpIn_nuts1[idx_BEinBE, :]
J_P_wwtpIn_nuts1 = J_P_wwtpIn_nuts1[idx_nuts1_16, :]

J_N_transportLoss_nuts1[idx_BB, :] = J_N_transportLoss_nuts1[idx_BB, :] +\
                                     J_N_transportLoss_nuts1[idx_BEinBB, :]
J_N_transportLoss_nuts1[idx_BE, :] = J_N_transportLoss_nuts1[idx_BEinBE, :]
J_N_transportLoss_nuts1 = J_N_transportLoss_nuts1[idx_nuts1_16, :]

J_P_transportLoss_nuts1[idx_BB, :] = J_P_transportLoss_nuts1[idx_BB, :] +\
                                     J_P_transportLoss_nuts1[idx_BEinBB, :]
J_P_transportLoss_nuts1[idx_BE, :] = J_P_transportLoss_nuts1[idx_BEinBE, :]
J_P_transportLoss_nuts1 = J_P_transportLoss_nuts1[idx_nuts1_16, :]

J_N_wwtpRem_nuts1[idx_BB, :] = J_N_wwtpRem_nuts1[idx_BB, :] +\
                               J_N_wwtpRem_nuts1[idx_BEinBB, :]
J_N_wwtpRem_nuts1[idx_BE, :] = J_N_wwtpRem_nuts1[idx_BEinBE, :]
J_N_wwtpRem_nuts1 = J_N_wwtpRem_nuts1[idx_nuts1_16, :]

J_P_wwtpRem_nuts1[idx_BB, :] = J_P_wwtpRem_nuts1[idx_BB, :] +\
                               J_P_wwtpRem_nuts1[idx_BEinBB, :]
J_P_wwtpRem_nuts1[idx_BE, :] = J_P_wwtpRem_nuts1[idx_BEinBE, :]
J_P_wwtpRem_nuts1 = J_P_wwtpRem_nuts1[idx_nuts1_16, :]

J_N_farm_nuts1[idx_BB, :] = J_N_farm_nuts1[idx_BB, :] +\
                            J_N_farm_nuts1[idx_BEinBB, :]
J_N_farm_nuts1[idx_BE, :] = J_N_farm_nuts1[idx_BEinBE, :]
J_N_farm_nuts1 = J_N_farm_nuts1[idx_nuts1_16, :]

J_P_farm_nuts1[idx_BB, :] = J_P_farm_nuts1[idx_BB, :] +\
                            J_P_farm_nuts1[idx_BEinBB, :]
J_P_farm_nuts1[idx_BE, :] = J_P_farm_nuts1[idx_BEinBE, :]
J_P_farm_nuts1 = J_P_farm_nuts1[idx_nuts1_16, :]

J_N_noSewerWwtp_nuts1[idx_BB, :] = J_N_noSewerWwtp_nuts1[idx_BB, :] +\
                                   J_N_noSewerWwtp_nuts1[idx_BEinBB, :]
J_N_noSewerWwtp_nuts1[idx_BE, :] = J_N_noSewerWwtp_nuts1[idx_BEinBE, :]
J_N_noSewerWwtp_nuts1 = J_N_noSewerWwtp_nuts1[idx_nuts1_16, :]

J_P_noSewerWwtp_nuts1[idx_BB, :] = J_P_noSewerWwtp_nuts1[idx_BB, :] +\
                                   J_P_noSewerWwtp_nuts1[idx_BEinBB, :]
J_P_noSewerWwtp_nuts1[idx_BE, :] = J_P_noSewerWwtp_nuts1[idx_BEinBE, :]
J_P_noSewerWwtp_nuts1 = J_P_noSewerWwtp_nuts1[idx_nuts1_16, :]

#%% -----------------------------------------------------------------------------
# 6) Estimate the model parameters at NUTS-1 level
# -----------------------------------------------------------------------------

# Set number of posterior (behavioural) simulations to be selected
Nb_behav = 100 

# Reshape simulated target variable (wwtp output load for sum of urban and rural)
J_N_psWwtpOut_nuts1_rshp = np.reshape(J_N_psWwtpOut_nuts1[:, :, 0] +\
                                      J_N_psWwtpOut_nuts1[:, :, 1],\
                                          (Nnuts1-2, Ny, Npar))
J_P_psWwtpOut_nuts1_rshp = np.reshape(J_P_psWwtpOut_nuts1[:, :, 0] +\
                                      J_P_psWwtpOut_nuts1[:, :, 1],\
                                          (Nnuts1-2, Ny, Npar))

# Prepare observation data
J_N_psWwtpOut_nuts1_obs_all =\
    np.tile(np.expand_dims(J_N_psWwtpOut_nuts1_obs, axis=2), (1, 1, Npar, 1))
J_P_psWwtpOut_nuts1_obs_all =\
    np.tile(np.expand_dims(J_P_psWwtpOut_nuts1_obs, axis=2), (1, 1, Npar, 1))
    
# Calculate the distance between simulations and observations consider the lower
# and upper bounds of the observations
d_N = np.maximum(J_N_psWwtpOut_nuts1_rshp-J_N_psWwtpOut_nuts1_obs_all[:, :, :, 1], 0) +\
      np.maximum(J_N_psWwtpOut_nuts1_obs_all[:, :, :, 0]-J_N_psWwtpOut_nuts1_rshp, 0)

d_P = np.maximum(J_P_psWwtpOut_nuts1_rshp-J_P_psWwtpOut_nuts1_obs_all[:, :, :, 1], 0) +\
      np.maximum(J_P_psWwtpOut_nuts1_obs_all[:, :, :, 0]-J_P_psWwtpOut_nuts1_rshp, 0)

# Calculate objective function (Root Mean Square Error -RMSE for N and P output 
# load from wwtp
fun_obj_N = np.sqrt(np.nanmean((d_N)**2, axis=1))/\
    np.nanmean(J_N_psWwtpOut_nuts1_obs_all, axis=(1,3))
fun_obj_P = np.sqrt(np.nanmean((d_P)**2, axis=1))/\
    np.nanmean(J_P_psWwtpOut_nuts1_obs_all, axis=(1,3))
fun_obj_all = (fun_obj_N + fun_obj_P)/2 # mean value between N and P

# Indices of posterior (behavioural) simulations (1st 100 best objective functions)
# Sort objective function
idx_sort = np.argsort(fun_obj_all, axis=1)
# Select 100 best realizations
idxb = idx_sort[:, 0:Nb_behav]

# Extract the posterior (behavioural) emissions
# Initialize variables
J_N_gross_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_gross_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_psNoTreat_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_psNoTreat_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_psWwtpOut_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_psWwtpOut_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_wwtpIn_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_wwtpIn_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_transportLoss_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_transportLoss_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_wwtpRem_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_wwtpRem_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_farm_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_farm_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_N_noSewerWwtp_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))
J_P_noSewerWwtp_nuts1_b = np.nan * np.zeros(((Nnuts1-2)*Ny, Nb_behav, Npop))

for ii in range(Nnuts1-2):
    J_N_gross_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_gross_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_gross_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_gross_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_N_psNoTreat_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_psNoTreat_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_psNoTreat_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_psNoTreat_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :] 
        
    J_N_psWwtpOut_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_psWwtpOut_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_psWwtpOut_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_psWwtpOut_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_N_wwtpIn_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_wwtpIn_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_wwtpIn_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_wwtpIn_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :] 
        
    J_N_transportLoss_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_transportLoss_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_transportLoss_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_transportLoss_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_N_wwtpRem_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_wwtpRem_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_wwtpRem_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_wwtpRem_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_N_farm_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_farm_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_farm_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_farm_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_N_noSewerWwtp_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_N_noSewerWwtp_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]
        
    J_P_noSewerWwtp_nuts1_b[ii*Ny:(ii+1)*Ny, :] =\
        J_P_noSewerWwtp_nuts1[ii*Ny:(ii+1)*Ny, idxb[ii, :], :]

# Calculate total emissions (aggregated over urban and rural areas)  
J_N_gross_nuts1_tot_b = np.sum(J_N_gross_nuts1_b, axis=2)
J_P_gross_nuts1_tot_b = np.sum(J_P_gross_nuts1_b, axis=2)
J_N_psNoTreat_nuts1_tot_b = np.sum(J_N_gross_nuts1_b, axis=2)
J_P_psNoTreat_nuts1_tot_b = np.sum(J_P_gross_nuts1_b, axis=2)
J_N_psWwtpOut_nuts1_tot_b_tot = np.sum(J_N_psWwtpOut_nuts1_b, axis=2)
J_P_psWwtpOut_nuts1_tot_b = np.sum(J_P_gross_nuts1_b, axis=2)
J_N_wwtpIn_nuts1_tot_b = np.sum(J_N_gross_nuts1_b, axis=2)
J_P_wwtpIn_nuts1_tot_b = np.sum(J_P_gross_nuts1_b, axis=2)
J_N_transportLoss_nuts1_tot_b = np.sum(J_N_gross_nuts1_b, axis=2)
J_P_transportLoss_nuts1_tot_b = np.sum(J_P_gross_nuts1_b, axis=2)
J_N_wwtpRem_nuts1_tot_b = np.sum(J_N_wwtpRem_nuts1_b, axis=2)
J_P_wwtpRem_nuts1_tot_b = np.sum(J_P_wwtpRem_nuts1_b, axis=2)
J_N_farm_nuts1_tot_b = np.sum(J_N_farm_nuts1_b, axis=2)
J_P_farm_nuts1_tot_b = np.sum(J_P_farm_nuts1_b, axis=2)
J_N_noSewerWwtp_nuts1_tot_b = np.sum(J_N_noSewerWwtp_nuts1_b, axis=2)
J_P_noSewerWwtp_nuts1_tot_b = np.sum(J_P_noSewerWwtp_nuts1_b, axis=2)

#%% -----------------------------------------------------------------------------
# 7) Downscale the N and P point sources emissions (treated and untreated) from 
#    NUTS-1 to grid level using gridded population data 
# -----------------------------------------------------------------------------

# Reshape NUTS-1 emissions
J_N_psNoTreat_nuts1_b_rshp = np.reshape(J_N_psNoTreat_nuts1_b,\
                                        (Nnuts1-2, Ny, Nb_behav, Npop))
J_P_psNoTreat_nuts1_b_rshp = np.reshape(J_P_psNoTreat_nuts1_b,\
                                        (Nnuts1-2, Ny, Nb_behav, Npop))
J_N_psWwtpOut_nuts1_b_rshp = np.reshape(J_N_psWwtpOut_nuts1_b,\
                                        (Nnuts1-2, Ny, Nb_behav, Npop))
J_P_psWwtpOut_nuts1_b_rshp = np.reshape(J_P_psWwtpOut_nuts1_b,\
                                        (Nnuts1-2, Ny, Nb_behav, Npop))
    
for bb in range(Nb_behav): # loop on the posterior realisations

    # Disaggregate N urban untreated point sources based on urban population
    J_N_psNoTreat_disagPop_urb =\
        NUTS2grid(J_N_psNoTreat_nuts1_b_rshp[:, :, bb, 0],\
                  id_nuts1_unique[0:Nnuts1-2], pop_urb, nuts1_grid)
    # Disaggregate P urban untreated point sources based on urban population    
    J_P_psNoTreat_disagPop_urb =\
    NUTS2grid(J_P_psNoTreat_nuts1_b_rshp[:, :, bb, 0],\
              id_nuts1_unique[0:Nnuts1-2], pop_urb, nuts1_grid)
    # Disaggregate N rural untreated point sources based on rural population
    J_N_psNoTreat_disagPop_rur =\
        NUTS2grid(J_N_psNoTreat_nuts1_b_rshp[:, :, bb, 1],\
                  id_nuts1_unique[0:Nnuts1-2], pop_rur, nuts1_grid)
    # Disaggregate P rural untreated point sources based on rural population    
    J_P_psNoTreat_disagPop_rur = \
        NUTS2grid(J_P_psNoTreat_nuts1_b_rshp[:, :, bb, 1],\
                  id_nuts1_unique[0:Nnuts1-2], pop_rur, nuts1_grid)    
        
    # Disaggregate N urban treated point sources based on urban population
    J_N_psWwtpOut_disagPop_urb =\
        NUTS2grid(J_N_psWwtpOut_nuts1_b_rshp[:, :, bb, 0],\
                  id_nuts1_unique[0:Nnuts1-2], pop_urb, nuts1_grid)
    # Disaggregate P urban treated point sources based on urban population    
    J_P_psWwtpOut_disagPop_urb =\
        NUTS2grid(J_P_psWwtpOut_nuts1_b_rshp[:, :, bb, 0],\
                  id_nuts1_unique[0:Nnuts1-2], pop_urb, nuts1_grid)
    # Disaggregate N rural treated point sources based on rural population
    J_N_psWwtpOut_disagPop_rur =\
        NUTS2grid(J_N_psWwtpOut_nuts1_b_rshp[:, :, bb, 1],\
                  id_nuts1_unique[0:Nnuts1-2], pop_rur, nuts1_grid)
    # Disaggregate P rural treated point sources based on rural population    
    J_P_psWwtpOut_disagPop_rur =\
        NUTS2grid(J_P_psWwtpOut_nuts1_b_rshp[:, :, bb, 1],\
                  id_nuts1_unique[0:Nnuts1-2], pop_rur, nuts1_grid)  
        
    # Calculate total emissions over urban and rural areas
    # N untreated point sources
    J_N_psNoTreat_disagPop = J_N_psNoTreat_disagPop_urb +\
                             J_N_psNoTreat_disagPop_rur
    # P untreated point sources                         
    J_P_psNoTreat_disagPop = J_P_psNoTreat_disagPop_urb +\
                             J_P_psNoTreat_disagPop_rur                         
    # N treated point sources
    J_N_psWwtpOut_disagPop = J_N_psWwtpOut_disagPop_urb +\
                             J_N_psWwtpOut_disagPop_rur
    # P treated point sources                         
    J_P_psWwtpOut_disagPop = J_P_psWwtpOut_disagPop_urb +\
                             J_P_psWwtpOut_disagPop_rur                            