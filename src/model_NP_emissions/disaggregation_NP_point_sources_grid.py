# -*- coding: utf-8 -*-
"""
Module to disaggregate Nitrogen (N) and Phosphorus (P) emissions from wastewater
to grid level

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
    
    NUTS2grid:
    function to estimate the values at grid level of a variable provided at a 
    certain spatial (NUTS) level.    
    
"""

import numpy as np


def NUTS2grid(X_nuts, id_nuts, w_grid, nuts_grid):
    '''
    This function estimate the values at grid level of a variable provided at 
    NUTS level (X_nuts) based on gridded (weight) variable (w_grid)-
    
    inputs:
        X_nuts = variable to disaggregate to grid level. 
                  Array of shape (N_nuts, Nyear) when N_nuts is the number of NUTS 
                  unit and Nyear is the number of years
        id_nuts = id of the NUTS units along axis 0 of X_nuts. 
                  Array of shape (N_nuts, )
        w_grid = gridded variable use to disaggregate X_nuts at grid level 
                  Array of shape (Nyear, Ny, Nx) when Ny is the number of cells 
                  along the y direction nad Nx if the number of cells along the 
                  x direction.
                  or Array of shape (Ny, Nx)
      nuts_grid = gridded map of NUTS regions. The values of the cells are the 
                  id of the NUTS regions.
                  Array of shape (Ny, Nx)
        
    output:
        X_grid = disaggregated variable to grid level.
                  Array of shape (Nyear, Ny, Nx).
                    
    '''
    #%% Check variables
    
    # Check that the number of elements along axis 0 is the same for X_nuts
    # and id_nuts
    shp_X = X_nuts.shape
    shp_id_nuts = id_nuts.shape
    if shp_X[0]!=shp_id_nuts[0]:
        raise ValueError('X_nuts and id_nuts should have the same number of elements '+\
                          'along axis 0')
    
    shp_nuts_grid = nuts_grid.shape
    shp_w_grid = w_grid.shape
    
    if len(shp_w_grid)==2:# if w_grid has 2 axes
       # Check the w_grid and nuts_grid have the same dimensions
       if shp_w_grid[0]!=shp_nuts_grid[0] or shp_w_grid[1]!=shp_nuts_grid[1]:
           raise ValueError('The dimensions of w_grid and nuts_grid should be consistent')
    
    elif len(shp_w_grid)==3:# if w_grid has 3 axes
        if shp_w_grid[1]!=shp_nuts_grid[0] or shp_w_grid[2]!=shp_nuts_grid[1]:
            raise ValueError(' the dimensions of w_grid and nuts_grid should be consistent')
        if shp_w_grid[0]!=shp_X[1]:
            raise ValueError('The number of elements in w_grid along axis 0 should equal to the number of elements in X_nuts along axis 1')    
    else:
        raise ValueError('shp_w_grid should have 2 or 3 axes')
 
    # Check that the values of id_nuts and the unique values of nuts_grid are the same
    nuts_grid_u = np.unique(nuts_grid[np.logical_not(np.isnan(nuts_grid))])
    id_nuts_sort = np.sort(id_nuts)
    if len(nuts_grid_u)!=len(id_nuts_sort):
        raise ValueError('The unique values in nuts_grid should be the values in id_nuts')       
    if np.any(nuts_grid_u!=id_nuts_sort):
        raise ValueError('The unique values in nuts_grid should be the values in id_nuts')
    
    # Number of NUTS regions
    Nnuts = shp_X[0]
    # Number of years
    Nyear = shp_X[1]
    # Number of cells along the y direction
    Ny = shp_nuts_grid[0]
    # Number of cells along the x direction
    Nx = shp_nuts_grid[1]
    
    #%% Disaggregate X_nuts

    # Initialize gridded variable
    X_grid = np.nan*np.zeros((Nyear, Ny, Nx))  
    
    # Loop over the NUTS regions
    for ii in range(Nnuts):
        # Indices of the ii-th NUTS region in the NUTS map
        idx_ii = id_nuts[ii]==nuts_grid
        # Number of cells for the ii-th NUTS region
        Ncell_ii = np.sum(idx_ii)
        
        # X_nuts is disaggregated using w_grid as a weight, i.e. we multiply 
        # the X_nuts by w_grid normalized by its sum over the NUTS region       
        X_nuts_ii = np.tile(np.expand_dims(X_nuts[ii, :], axis=1), (1, Ncell_ii))
        
        # If w_grid does not have a time dimension
        if len(shp_w_grid)==2: # no time dimension in w_grid
            # Normalize w_grid_ii by the total value and repeat the same value 
            # for each year
            w_grid_ii = np.tile(np.expand_dims(w_grid[idx_ii]/np.sum(w_grid[idx_ii]), axis=0), (Nyear, 1))
            # Disaggregation
            X_grid[:, idx_ii] = X_nuts_ii*w_grid_ii
        
        if len(shp_w_grid)==3: # time dimension in w_grid
            # Normalize w_grid_ii by the total value and repeat the same value 
            # for each year
            w_grid_ii = w_grid[:, idx_ii]/\
                np.tile(np.expand_dims(np.sum(w_grid[:, idx_ii], axis=1), axis=1), (1, Ncell_ii))
            # Disaggregation
            X_grid[:, idx_ii] = X_nuts_ii*w_grid_ii
            
        # Check that the sum of the gridded variable is equal to the NUTS total
        idx_pos = X_nuts[ii, :]!=0
        if np.any(idx_pos):
            if np.max(abs(np.sum(X_grid[:, idx_ii], axis=1)[idx_pos] -\
                                 X_nuts[ii, idx_pos])/X_nuts[ii, idx_pos])>10**(-4):
                print(np.max(abs(np.sum(X_grid[:, idx_ii], axis=1)[idx_pos] -\
                                     X_nuts[ii, idx_pos])/X_nuts[ii, idx_pos]))
                raise ValueError('The sum of disaggregated variable is not equal to the NUTS value')
        idx_zero = X_nuts[ii, :]==0
        if np.any(idx_zero):
            if np.max(abs(np.sum(X_grid[:, idx_ii], axis=1)[idx_zero] -\
                                 X_nuts[ii, idx_zero]))>10**(-4):
                print(np.max(abs(np.sum(X_grid[:, idx_ii], axis=1)[idx_zero] -\
                                     X_nuts[ii, idx_zero])))
                raise ValueError('The sum of disaggregated variable is not equal to the NUTS value')
        
    # Check that result is positive
    if np.any(X_grid<0):
        raise ValueError('X_grid<0')
    # Check that the NAN of X_grid are the same as nuts_grid
    for yy in range(Nyear):
        if np.any(np.isnan(nuts_grid)!=np.isnan(X_grid[yy, :, :])):
            raise ValueError('NANs should be the same in nuts_grid and X_grid')
        
    
    return X_grid

