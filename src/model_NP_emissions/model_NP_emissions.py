"""

Module to estimate domestic and industrial/commercial Nitrogen (N) and 
Phosphorus (P) gross and net emissions from wastewater

LICENSE:
    
    Copyright 2024, Fanny J. Sarrazin, Rohini Kumar
    All rights reserved.
    
    This source code is licensed under the GNU General Public License v3.0 
    found in the LICENSE file in the root directory of this source tree. 


PLEASE CITE THIS PUBLICATION IF YOU USE THIS CODE:
    
    Sarrazin, F. J., Attinger, A., Kumar, R., Gridded dataset of nitrogen and 
    phosphorus point sources from wastewater in Germany (1950-2019), submitted 
    to Earth System Science Data.
    
    
CONTENTS:
    
    - domestic_nutrient_gross_emissions: 
      function to calculate N and P domestic gross human emissions corresponding 
      to human excreta based on protein supply data.

    - domestic_nutrient_net_emissions:
      function to calculate domestic N and P net emissions, i.e. it determines the 
      fate of domestic gross emissions (treated and untreated point sources or 
      other fates).
    
    - industrial_nutrient_gross_emissions: 
      function to calculates industrial/commercial gross emissions that are 
      collected by the sewer system.
    
    - industrial_nutrient_net_emissions: 
      function to calculate industrial/commercial net emissions, i.e. it determines 
      the fate of industrial/commercial gross emissions (treated and untreated 
      point sources or other fates).
    
    - aggregate_dom_ind_emissions: 
      function to aggregate domestic and industrial/commercial emissions.

"""
    
import numpy as np
from warnings import warn


def domestic_nutrient_gross_emissions(Xpar, pop, protein):
    '''
    This function calculates nutriment (N or P) domestic gross human emissions
    corresponding to human excreta (without the detergent component) based on
    on protein supply data.
      
      inputs:
                     Xpar = coefficients
                  Xpar[0] = Fraction of protein supply wasted at the distribution
                            and consumption level f_protein_waste (-)
                  Xpar[1] = Nutriment intake as a fraction of protein intake 
                            f_Nutri_pro (-)
                  Xpar[2] = Fraction of nutrient intake lost via sweat, hair and 
                            blood f_lossHum (-)
                      pop = population count (capita)
                  protein = protein supply at the retail level (kg/capita/yr)
       outputs:
        J_Nutri_gross_dom = Nutrient domestic gross emissions (human excreta) (kg/yr)
    '''
    
    #%% Check inputs
    # Protein, population are >0
    if np.any(protein<0):
        raise ValueError('protein<0')
    if np.any(pop<0):
        raise ValueError('pop<0')
        
    # Check that parameter values are between 0 and 1
    if np.any(np.logical_or(Xpar<0, Xpar>1)):
        raise ValueError('Xpar<0 or Xpar>1')

    f_pro_waste = Xpar[0]
    f_Nutri_pro = Xpar[1]
    f_lossHum = Xpar[2]

    #%% Calculate domestic emissions 
    
    # Gross emissions
    J_Nutri_gross_dom = pop*protein*f_Nutri_pro*(1-f_pro_waste)*(1-f_lossHum)
    
     #%% Check results
    
    # Check that all variables are >0
    if np.any(J_Nutri_gross_dom<0):
        raise ValueError('J_Nutri_gross_dom<0')
           
    # check that there are no nan 
    if np.any(np.isnan(J_Nutri_gross_dom)):
        raise ValueError('J_Nutri_gross_dom is nan')
       
    return J_Nutri_gross_dom

    
def domestic_nutrient_net_emissions(Xpar, J_Nutri_gross_dom, T0, TFarm, T1,\
                                    T23_noNutri, T3_Nutri, err_mb=10**(-8)):
    '''
    This function calculates domestic nutriment (N or P) net emissions, i.e. 
    it determines the fate of domestic gross emissions (treated and untreated point 
    sources or other fates).
      
      inputs:
                  Xpar = coefficients    
                  Xpar[0] = Fraction of N and P emissions lost during wastewater 
                            collection and transport (losses in the sewer system, 
                            and in cesspits and duringtransportation of wastewater 
                            collected in cesspits to wwtps) f_lossTransport (-)
                  Xpar[1] = Efficiency of nutrient removal for primary treatment 
                            eff_1_Nutri (-)
                  Xpar[2] = Efficiency of nutrient removal for secondary treatment  
                            and tertiary treatment without targeted nutrient
                            removal eff_23_noNutri (-)
                  Xpar[3] = Efficiency of nutrient removal for tertiary treatment 
                            with targeted nutrient removal eff_3_Nutri (-)        
        J_Nutri_gross_dom = Nutrient domestic gross emissions (kg/yr)
                       T0 = fraction of population connected to the sewer but 
                            not to public wastewater treatment (-)
                    TFarm = fraction of population connected to sewer farms (-)     
                       T1 = fraction of population connected to primary wastewater 
                            treatment (-)
              T23_noNutri = fraction of population connected to secondary and
                            tertiary wastewater treatment without targeted 
                            nutrient removal (-)
                 T3_Nutri = fraction of population connected to tertiary wastewater 
                            treatment with targeted nutrient removal (-)   
      optional input:                      
                  err_mb = tolerance on mass balance error
                
      outputs:
   J_Nutri_psWwtpOut_dom = Nutrient WWTP outgoing point sources domestic emissions  (kg/yr)
   J_Nutri_psNoTreat_dom = Nutrient domestic point sources emissions not treated (kg/yr)
      J_Nutri_wwtpIn_dom = Nutrient WWTP incoming domestic emissions (kg/yr) 
        J_Nutri_farm_dom = Nutrient domestic emissions treated in sewage farms
                           (incoming emissions to sewage farms without accounting
                            for treatment removal)     
 J_Nutri_noSewerWwtp_dom = Nutrient domestic emissions not collected by the sewer 
                            system nor treated in wwtp (kg/yr)
J_Nutri_transportLoss_dom = Nutrient domestic emissions lost during collection
                            and transport  (kg/yr)
      J_Nutri_wwtpRem_dom = Nutrient domestic emissions removed in WWTPs (kg/yr)     
                 
     '''
    
    #%% Check inputs
    
    # T0, T1, T23_noNutri, T3_Nutri between 0 and 1
    if np.any(np.logical_or(T0<0, T0>1)):
        raise ValueError('T0<0 or T0>1')
    if np.any(np.logical_or(T1<0, T1>1)):
        raise ValueError('T1<0 or T1>1')
    if np.any(np.logical_or(T23_noNutri<0, T23_noNutri>1)):
        raise ValueError('T23_noNutri<0 or T23_noNutri>1')
    if np.any(np.logical_or(T3_Nutri<0, T3_Nutri>1)):
        raise ValueError('T3_Nutri<0 or T3_Nutri>1') 
    if np.any(np.logical_or(TFarm<0, TFarm>1)):
        raise ValueError('TFarm<0 or TFarm>1')         
    # TFarm+T0+T1+T2+T3<1
    if np.any(TFarm+T0+T1+T23_noNutri+T3_Nutri>1.000001):
        raise ValueError('T0+T1+T23_noNutri+T3_Nutri>1')
    # J_Nutri_dom_gross is >0
    if np.any(J_Nutri_gross_dom<0):
        raise ValueError('J_Nutri_gross_dom<0')
        
    # Check that parameter values are between 0 and 1
    if np.any(np.logical_or(Xpar<0, Xpar>1)):
        raise ValueError('Xpar<0 or Xpar>1')

    f_lossTransport = Xpar[0]
    eff_1_Nutri = Xpar[1]
    eff_23_noNutri = Xpar[2]
    eff_3_Nutri = Xpar[3]

    #%% Calculate domestic emissions 
    
    # Fraction of population connected to sewer and/or wastewater treatment
    T_sewer_wwtp = np.minimum(TFarm+T0+T1+T23_noNutri+T3_Nutri, 1) # ensure that it is <=1
    # (could be slightly >1 due to roundoff errors)
        
    # WWTP incoming load
    J_Nutri_wwtpIn_dom = J_Nutri_gross_dom*(1-f_lossTransport)*(T1 + T23_noNutri+ T3_Nutri)
        
    # WWTP outgoing load
    J_Nutri_psWwtpOut_dom = J_Nutri_gross_dom*(1-f_lossTransport)*(T1*(1-eff_1_Nutri) +\
                           T23_noNutri*(1-eff_23_noNutri) + T3_Nutri*(1-eff_3_Nutri))
    
    # Treated in sewage farms
    J_Nutri_farm_dom = J_Nutri_gross_dom*(1-f_lossTransport)*TFarm
        
    # Untreated point sources from the sewer system
    J_Nutri_psNoTreat_dom = J_Nutri_gross_dom*(1-f_lossTransport)*T0
    
    # Emissions not collected by the sewer or treated in wwtps
    J_Nutri_noSewerWwtp_dom = J_Nutri_gross_dom*(1-T_sewer_wwtp)
    
    # Emissions lost during wastewater collection and transport
    J_Nutri_transportLoss_dom = J_Nutri_gross_dom*T_sewer_wwtp*f_lossTransport
    
    # Removal in WWTP
    J_Nutri_wwtpRem_dom = J_Nutri_wwtpIn_dom - J_Nutri_psWwtpOut_dom
        
    #%% Check results
    
    # Check that the sum of all path is equal to the gross emissions   
    idx_pos = J_Nutri_gross_dom!=0 # index of non-zero values of gross emissions
    if np.any(abs(J_Nutri_psWwtpOut_dom[idx_pos] + J_Nutri_psNoTreat_dom[idx_pos] +\
              J_Nutri_farm_dom[idx_pos] + J_Nutri_noSewerWwtp_dom[idx_pos] +\
              J_Nutri_transportLoss_dom[idx_pos] + J_Nutri_wwtpRem_dom[idx_pos] -\
              J_Nutri_gross_dom[idx_pos])/J_Nutri_gross_dom[idx_pos]>err_mb):
        raise ValueError('Mismatch between gross and net emissions')
    idx_zero = J_Nutri_gross_dom==0 # index of zero values of gross emissions
    if np.any(abs(J_Nutri_psWwtpOut_dom[idx_zero] + J_Nutri_psNoTreat_dom[idx_zero] +\
              J_Nutri_farm_dom[idx_zero] + J_Nutri_noSewerWwtp_dom[idx_zero] +\
              J_Nutri_transportLoss_dom[idx_zero] + J_Nutri_wwtpRem_dom[idx_zero] -\
              J_Nutri_gross_dom[idx_zero])>err_mb):
        raise ValueError('Mismatch between gross and net emissions')

    # Check that all variables are >0
    if np.any(J_Nutri_psNoTreat_dom<0):
        raise ValueError('J_Nutri_psNoTreat_dom<0')
    if np.any(J_Nutri_psWwtpOut_dom<0):
        raise ValueError('J_Nutri_psWwtpOut_dom<0') 
    if np.any(J_Nutri_wwtpIn_dom<0):
        raise ValueError('J_Nutri_wwtpIn_dom<0') 
    if np.any(J_Nutri_farm_dom<0):
        raise ValueError('J_Nutri_farm_dom<0')    
    if np.any(J_Nutri_noSewerWwtp_dom<0):
       raise ValueError('J_Nutri_noSewerWwtp_dom<0')     
    if np.any(J_Nutri_transportLoss_dom<0):
        raise ValueError('J_Nutri_transportLoss_dom<0')    
    if np.any(J_Nutri_wwtpRem_dom<0):
        raise ValueError('J_Nutri_wwtpRem_dom<0')    
        
    # check that there are no nan 
    if np.any(np.isnan(J_Nutri_psNoTreat_dom)):
       raise ValueError('J_Nutri_psNoTreat_dom is nan')
    if np.any(np.isnan(J_Nutri_psWwtpOut_dom)):
        raise ValueError('J_Nutri_psWwtpOut_dom is nan')   
    if np.any(np.isnan(J_Nutri_wwtpIn_dom)):
        raise ValueError('J_Nutri_wwtpIn_dom is nan')  
    if np.any(np.isnan(J_Nutri_farm_dom)):
        raise ValueError('J_Nutri_farm_dom is nan')    
    if np.any(np.isnan(J_Nutri_noSewerWwtp_dom)):
        raise ValueError('J_Nutri_noSewerWwtp_dom is nan') 
    if np.any(np.isnan(J_Nutri_noSewerWwtp_dom)):
        raise ValueError('J_Nutri_noSewerWwtp_dom is nan')     
    if np.any(np.isnan(J_Nutri_wwtpRem_dom)):
        raise ValueError('J_Nutri_wwtpRem_dom is nan') 
        
        
    return J_Nutri_psWwtpOut_dom, J_Nutri_psNoTreat_dom, J_Nutri_wwtpIn_dom,\
           J_Nutri_farm_dom, J_Nutri_noSewerWwtp_dom, J_Nutri_transportLoss_dom,\
           J_Nutri_wwtpRem_dom 
           
           
           
def industrial_nutrient_gross_emissions(Xpar, pop, protein, f_ind_dom):
    '''
    This function calculates industrial/commercial gross emissions that are 
    collected by the sewer system (i.e. that do not end up in stabilization ponds 
    or are lost via volatilization). We assume that gross emissions are equal to 
    a time varying fraction of domestic gross emissions (f_indus_dom). This 
    function does not account for the detergent part.
  
      inputs:
                    Xpar = coefficients
                  Xpar[0] = Fraction of protein supply wasted at the distribution
                            and consumption level f_protein_waste (-)
                  Xpar[1] = Nutriment intake as a fraction of protein intake 
                            f_Nutri_pro (-)
                  Xpar[2] = Fraction of nutrient intake lost via sweat, hair and 
                            blood f_lossHum (-)
                  protein = protein supply at the retail level (kg/capita/yr)
                f_ind_dom = fraction of industrial/commercial to domestic 
                            nutrient gross emissions (-)
          
      outputs:

        J_Nutri_ind_gross = Nutrient industrial gross emissions (kg/yr)
    '''
    
    #%% Check inputs
    
    # Protein, population, f_indus_dom are >0
    if np.any(protein<0):
        raise ValueError('protein<0')
    if np.any(pop<0):
        raise ValueError('pop<0')
    if np.any(f_ind_dom<0):
        raise ValueError('f_ind_dom<0')
        
    # Check that parameter values are between 0 and 1
    if np.any(np.logical_or(Xpar<0, Xpar>1)):
        raise ValueError('Xpar<0 or Xpar>1')

    f_pro_waste = Xpar[0]
    f_Nutri_pro = Xpar[1]
    f_lossHum = Xpar[2]

    #%% Calculate industrial emissions
    
    # Gross Industrial emissions  
    J_Nutri_gross_ind = f_ind_dom*pop*protein*f_Nutri_pro*(1-f_pro_waste)*(1-f_lossHum)
    
    #%% Check results
    
    # Check that all variables are >0
    if np.any(J_Nutri_gross_ind<0):
        raise ValueError('J_Nutri_gross_ind<0')
        
    # check that there are no nan    
    if np.any(np.isnan(J_Nutri_gross_ind)):
        raise ValueError('J_Nutri_gross_ind is nan')
    
        
    return J_Nutri_gross_ind


def industrial_nutrient_net_emissions(Xpar, J_Nutri_gross_ind, T0, TFarm, T1,\
                                      T23_noNutri, T3_Nutri, err_mb=10**(-8)):
    '''
    This function calculates industrial/commercial net emissions assuming 
    that all of these industrial emissions are collected by the sewer system
    are treated with the same treatment efficiencies as domestic emissions.
    This function determines the fate of industrial/commercial gross emissions 
    (treated and untreated point sources or other fates).
      
      inputs:
                     Xpar = coefficients
                  Xpar[0] = Fraction of N and P emissions lost during wastewater 
                            collection and transport (losses in the sewer system, 
                            and in cesspits and duringtransportation of wastewater 
                            collected in cesspits to wwtps) f_lossTransport (-)
                  Xpar[1] = Efficiency of nutrient removal for primary treatment 
                            eff_1_Nutri (-)
                  Xpar[2] = Efficiency of nutrient removal for secondary treatment  
                            and tertiary treatment without targeted nutrient
                            removal eff_23_noNutri (-)
                  Xpar[3] = Efficiency of nutrient removal for tertiary treatment 
                            with targeted nutrient removal eff_3_Nutri (-)
        J_Nutri_gross_ind = Nutrient industrial gross emissions (kg/yr)
                f_ind_dom = ratio of industrial to domestic nutrient gross 
                            emissions (-)
                       T0 = fraction of population connected to the sewer but 
                            not to public wastewater treatment (-)
                    TFarm = fraction of population connected to sewer farms (-)  
                       T1 = fraction of population connected to primary wastewater 
                            treatment (-)
              T23_noNutri = fraction of population connected to secondary and
                            tertiary wastewater treatment without targeted 
                            nutrient removal (-)
                 T3_Nutri = fraction of population connected to tertiary wastewater 
                            treatment with targeted nutrient removal (-)    
      optional input                        
                   err_mb = tolerance on mass balance error
          
      outputs:
   J_Nutri_psWwtpOut_ind = Nutrient WWTP outgoing point sources industrial/commercial  
                           emissions (kg/yr)
   J_Nutri_psNoTreat_ind = Nutrient industrial/commercial point sources emissions
                           not treated (kg/yr)
      J_Nutri_wwtpIn_ind = Nutrient WWTP incoming industrial/commercial emissions 
                           (kg/yr) 
J_Nutri_transportLoss_ind = Nutrient industrial/commercial emissions lost during 
                            collection and transport  (kg/yr)
      J_Nutri_wwtpRem_ind = Nutrient domestic emissions removed in WWTPs (kg/yr)     
    
    '''
    
     
    #%% Check inputs
    # Total connection to sewer and wwtp
    T_sewer_wwtp = TFarm + T0 + T1 + T23_noNutri + T3_Nutri
     
    # T0, T1, T23_noNutri, T3_Nutri between 0 and 1
    if np.any(np.logical_or(T0<0, T0>1)):
        raise ValueError('T0<0 or T0>1')
    if np.any(np.logical_or(TFarm<0, TFarm>1)):
        raise ValueError('TFarm<0 or TFarm>1')
    if np.any(np.logical_or(T1<0, T1>1)):
        raise ValueError('T1<0 or T1>1')
    if np.any(np.logical_or(T23_noNutri<0, T23_noNutri>1)):
        raise ValueError('T23_noNutri<0 or T23_noNutri>1')
    if np.any(np.logical_or(T3_Nutri<0, T3_Nutri>1)):
        raise ValueError('T3_Nutri<0 or T3_Nutri>1')    
    if np.any(np.logical_or(T_sewer_wwtp<0, T_sewer_wwtp>1.000001)):
        raise ValueError('Ttot<0 or Ttot>1')    
    # J_Nutri_indus_gross is >0
    if np.any(J_Nutri_gross_ind<0):
        raise ValueError('J_Nutri_indus_gross<0')

    # Check that parameter values are between 0 and 1
    if np.any(np.logical_or(Xpar<0, Xpar>1)):
        raise ValueError('Xpar<0 or Xpar>1')

    f_lossTransport = Xpar[0]
    eff_1_Nutri = Xpar[1]
    eff_23_noNutri = Xpar[2]
    eff_3_Nutri = Xpar[3]
    
    Nval = len(J_Nutri_gross_ind) # number of records
    
    # Check whether the T_sewer_wwtp contains <0 values
    if np.any(T_sewer_wwtp==0):
        warn('T_sewer_wwtp has zero values. When it is zero the entire'+\
             ' ''J_Nutri_gross_ind'' ends up as untreated point sources (''J_Nutri_psNoTreat_ind'')'+\
             ' while accounting for transport losses (''J_Nutri_transportLoss_ind'')')
        
    #%% Calculate industrial emissions
    
    # Initialise net industrial emissions
    J_Nutri_wwtpIn_ind = np.nan * np.zeros((Nval, ))
    J_Nutri_psWwtpOut_ind = np.nan * np.zeros((Nval, ))
    J_Nutri_psNoTreat_ind = np.nan * np.zeros((Nval, ))
    J_Nutri_transportLoss_ind = np.nan * np.zeros((Nval, ))
    J_Nutri_farm_ind = np.nan * np.zeros((Nval, ))
    
    # Case 1: T_sewer_wwtp is not equal to 0    
    # We normalise T0, T1, T23_noNutri and T3_Nutri so that the sum of the 
    # different treatment types is equal to 1 (unlike domestic emissions, 
    # industrial emissions are assumed to be all collected by the sewer)

    # Untreated point sources from the sewer system
    J_Nutri_psNoTreat_ind[T_sewer_wwtp!=0] =\
        J_Nutri_gross_ind[T_sewer_wwtp!=0]*\
        (1-f_lossTransport)*T0[T_sewer_wwtp!=0]/T_sewer_wwtp[T_sewer_wwtp!=0]
    
    # WWTP incoming load
    J_Nutri_wwtpIn_ind[T_sewer_wwtp!=0] =\
        J_Nutri_gross_ind[T_sewer_wwtp!=0]*(1-f_lossTransport)*\
            (T1[T_sewer_wwtp!=0] +\
             T23_noNutri[T_sewer_wwtp!=0]+\
             T3_Nutri[T_sewer_wwtp!=0])/T_sewer_wwtp[T_sewer_wwtp!=0]
    
    # WWTP outgoing load
    J_Nutri_psWwtpOut_ind[T_sewer_wwtp!=0] =\
        J_Nutri_gross_ind[T_sewer_wwtp!=0]*(1-f_lossTransport)*\
            (T1[T_sewer_wwtp!=0]*(1-eff_1_Nutri) +\
             T23_noNutri[T_sewer_wwtp!=0]*(1-eff_23_noNutri) +\
             T3_Nutri[T_sewer_wwtp!=0]*(1-eff_3_Nutri))/T_sewer_wwtp[T_sewer_wwtp!=0]
                
    J_Nutri_farm_ind[T_sewer_wwtp!=0] =\
    J_Nutri_gross_ind[T_sewer_wwtp!=0]*\
    (1-f_lossTransport)*TFarm[T_sewer_wwtp!=0]/T_sewer_wwtp[T_sewer_wwtp!=0]
    
    # Emissions lost during wastewater collection and transport
    J_Nutri_transportLoss_ind[T_sewer_wwtp!=0] = J_Nutri_gross_ind[T_sewer_wwtp!=0]*\
             f_lossTransport*\
            (T0[T_sewer_wwtp!=0] +\
             TFarm[T_sewer_wwtp!=0]+\
             T1[T_sewer_wwtp!=0] +\
             T23_noNutri[T_sewer_wwtp!=0]+\
             T3_Nutri[T_sewer_wwtp!=0])/T_sewer_wwtp[T_sewer_wwtp!=0]
    
    # Case 2: if TFarm+T0+T1+T2+T3 = 0 assume no treatment, but that all emissions are
    # transported via pipes to surface waters
    J_Nutri_wwtpIn_ind[T_sewer_wwtp==0] = 0
    J_Nutri_psWwtpOut_ind[T_sewer_wwtp==0] = 0
    J_Nutri_farm_ind[T_sewer_wwtp==0] = 0
    J_Nutri_psNoTreat_ind[T_sewer_wwtp==0] =\
        J_Nutri_gross_ind[T_sewer_wwtp==0]*(1-f_lossTransport)
    
    J_Nutri_transportLoss_ind[T_sewer_wwtp==0] = J_Nutri_gross_ind[T_sewer_wwtp==0]*f_lossTransport
    # Variable that are calculated the same way for case 1 and 2 above 
    # Removal in WWTP
    J_Nutri_wwtpRem_ind = J_Nutri_wwtpIn_ind - J_Nutri_psWwtpOut_ind
        
    #%% Check results
                       
    # Check that the sum of all path is equal to the gross emissions
    idx_T_sewer_wwtp_pos = np.where(T_sewer_wwtp!=0)[0]
    J_Nutri_in = J_Nutri_gross_ind[idx_T_sewer_wwtp_pos]*\
        (TFarm[idx_T_sewer_wwtp_pos]+T0[idx_T_sewer_wwtp_pos]+T1[idx_T_sewer_wwtp_pos]+\
         T23_noNutri[idx_T_sewer_wwtp_pos]+T3_Nutri[idx_T_sewer_wwtp_pos])/T_sewer_wwtp[idx_T_sewer_wwtp_pos]
    
    idx_pos = J_Nutri_in!=0
    if np.any(abs(J_Nutri_psWwtpOut_ind[idx_T_sewer_wwtp_pos[idx_pos]] +\
                  J_Nutri_psNoTreat_ind[idx_T_sewer_wwtp_pos[idx_pos]] +\
                  J_Nutri_transportLoss_ind[idx_T_sewer_wwtp_pos[idx_pos]] +\
                  J_Nutri_wwtpRem_ind[idx_T_sewer_wwtp_pos[idx_pos]] +\
                  J_Nutri_farm_ind[idx_T_sewer_wwtp_pos[idx_pos]]-\
                  J_Nutri_in[idx_pos])/J_Nutri_in[idx_pos]>err_mb):
        raise ValueError('Mismatch between gross and net emissions')
    
    idx_zero = J_Nutri_in==0
    if np.any(abs(J_Nutri_psWwtpOut_ind[idx_T_sewer_wwtp_pos[idx_zero]] +\
                  J_Nutri_psNoTreat_ind[idx_T_sewer_wwtp_pos[idx_zero]] +\
                  J_Nutri_transportLoss_ind[idx_T_sewer_wwtp_pos[idx_zero]] +\
                  J_Nutri_wwtpRem_ind[idx_T_sewer_wwtp_pos[idx_zero]] +\
                  J_Nutri_farm_ind[idx_T_sewer_wwtp_pos[idx_zero]]-\
                  J_Nutri_in[idx_zero])>err_mb):
        raise ValueError('Mismatch between gross and net emissions')
        
    # Check that all variables are >0
    if np.any(J_Nutri_psWwtpOut_ind<0):
        raise ValueError('J_Nutri_psWwtpOut_ind<0')
    if np.any(J_Nutri_psNoTreat_ind<0):
        raise ValueError('J_Nutri_psNoTreat_ind<0')
    if np.any(J_Nutri_transportLoss_ind<0):
        raise ValueError('J_Nutri_transportLoss_ind<0') 
    if np.any(J_Nutri_wwtpRem_ind<0):
        raise ValueError('J_Nutri_wwtpRem_ind<0')
    if np.any(J_Nutri_wwtpIn_ind<0):
        raise ValueError('J_Nutri_wwtpIn_ind<0') 
    if np.any(J_Nutri_farm_ind<0):
        raise ValueError('J_Nutri_farm_ind<0') 
        
    # check that there are no nan
    if np.any(np.isnan(J_Nutri_psWwtpOut_ind)):
        raise ValueError('J_Nutri_psWwtpOut_ind is nan')  
    if np.any(np.isnan(J_Nutri_psNoTreat_ind)):
        raise ValueError('J_Nutri_psNoTreat_ind is nan')  
    if np.any(np.isnan(J_Nutri_transportLoss_ind)):
        raise ValueError('J_Nutri_transportLoss_ind is nan')      
    if np.any(np.isnan(J_Nutri_wwtpRem_ind)):
        raise ValueError('J_Nutri_wwtpRem_ind is nan')  
    if np.any(np.isnan(J_Nutri_wwtpIn_ind)):
        raise ValueError('J_Nutri_wwtpIn_ind is nan') 
    if np.any(np.isnan(J_Nutri_farm_ind)):
        raise ValueError('J_Nutri_farm_ind is nan') 
        
    return J_Nutri_psWwtpOut_ind, J_Nutri_psNoTreat_ind, J_Nutri_wwtpIn_ind,\
        J_Nutri_farm_ind, J_Nutri_transportLoss_ind, J_Nutri_wwtpRem_ind 


def aggregate_dom_ind_emissions(X_dom, X_ind, VarName, axis_pop, idx_urb, idx_rur):
    
     '''
     This function aggregates domestic and industrial/commercial emissions.
     Domestic emissions should be provided for urban and rural areas, while
     all industrial/commercial emissions are attributed to the urban fraction.
     
     inputs:
        X_dom = domestic emissions (contains emissions for total, urban and rural
                population along the axis axis_pop)
        X_ind = industrial/commercial emissions
     axis_pop = axis for urban and rural emissions   
      idx_urb = index for urban emissions along the axis axis_pop
      idx_rur = index for rural emissions along the axis axis_pop
      
     output 
        X_agg = aggregated domestic and industrial/commercial emissions 
                (contains emissions for urban and rural
                        population along the axis axis_pop)
     '''
    
     # Move axes of X so that axis_pop is the first dimension
     X_dom_perm = np.moveaxis(X_dom, axis_pop, 0)
    
      # Initialize variable
     X_agg  = np.nan * np.zeros((X_dom_perm.shape))
    
     # Urban values (urban part of domestic emissions + industry emissions)
     X_agg[idx_urb, :] = X_dom_perm[idx_urb, :] + X_ind
     # Rural values (rural part of domestic emissions)
     X_agg[idx_rur, :] = X_dom_perm[idx_rur, :]

     # Move axis so that axis_pop goes back to its original position
     X_agg = np.moveaxis(X_agg, 0, axis_pop)  
    
     return X_agg