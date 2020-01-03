# -*- coding: utf-8 -*-

'''
******************************************************************************

     =========================================================================
     Subroutine to calculate crater dimensions at each step and for the transient
     and final crater diameters. 
    ==========================================================================
    

    
    ==========================================================================
    
    TODO: 
        I should give the possibility to define the transient crater diameter
        at the maximum depth during crater evolution. So we can choose for 
        example this mode in the case of low velocity impact and problem related
        to the somehow sinking of the impact structure.

******************************************************************************
'''

# Basic modules are loaded
import numpy as np
import os
import sys
import copy
import re

# path to in-house scripts and pySALEPlot
path_pySALEPlot = '/work/nilscp/iSALE/Dellen/lib'

# path to pySALEPlot is loaded
sys.path.append(path_pySALEPlot)
import pySALEPlot as psp

'''
***********************************************************************
'''

def tokenize(filename):
    '''
    Function to list filenames in correct order (see
    http://stackoverflow.com/questions/5997006/sort-a-list-of-files-using-python)

    :param filename:
    :return:
    '''
    digits = re.compile(r'(\d+)')

    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

'''
***********************************************************************
'''

def movingaverage(array, window):
    '''
    description:
    average values along a moving window

    param:
    : array : 1-D numpy array
    : window : window size (in n number of values)
        
    returns:
    : ma : moving average
    '''

    weights = np.repeat(1.0, window)/window
    ma = np.convolve(array, weights, 'valid')
    return ma


'''
***********************************************************************
'''


def ixdiam(model, step, id_mat = 0):
    
    '''
    description:
    return the first non-empty cell along the pre-impact surface (x-axis)

    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)      
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : ix (int): moving average   
    '''
    
    # get the y-index that correponds to the pre-impact surface
    oridy = np.where(model.yy == 0.)[0][0]
    
    # get the concentration of material in a cell along the pre-impact surface
    conc = step.cmc[id_mat][:, oridy]

    # get a moving average over 20 cells (to avoid detection of ejected materials)
    ma = movingaverage(conc, 20)

    # if the amount of material fill more than 99.9% of the cell, get the index
    if np.max(ma > 0.999):
        ix = np.where(ma > 0.999)[0][0]
        
    # otherwise get the index where the first maximum is detected
    else:
        ix = np.where(ma == np.max(ma))[0][-1]

    return ix


'''
***********************************************************************
'''


def appvolume(model, step, id_mat = 0):
    
    '''
    description:
    calculate the volume of the apparent crater diameter (from the pre-impact
    surface, Va)
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : Vapp (float): Apparent crater volume   

    '''
    # if this is the first time step, the volume is then equal to 0
    if step.time == 0:
        Vapp = 0
            
    else:        
        try:
            # get the y-index that correponds to the pre-impact surface
            oridy = np.where(model.yy == 0.)[0][0]
            
            # get the first non-empty cell along the pre-impact surface (x-axis)
            ix = ixdiam(model, step, id_mat)
                        
            # get the concentration in all cells up to the pre-impact surface (y)
            # and up to the apparent crater diameter (x)
            conc = step.cmc[id_mat][:ix, :oridy]
    
            # select only cells that are not completelly filled with materials
            empty_cells = np.where(conc < 1)
    
            # select the x-index, x-index + 1 and y-index
            i = empty_cells[0]
            ii = empty_cells[0] + 1
            j = empty_cells[1]
    
            # r is the distance from the axis of symmetry
            r = model.x[i, j]
    
            # h is the height of the cells, which is the same for all the cells
            h = model.dy
    
            # the distance along the x-axis between two cells
            # takes into account the concentration of materials
            R = model.x[i, j] + ((1 - step.cmc[id_mat][i, j])
                                 * (model.x[ii, j] - model.x[i, j]))
    
            # calculation of the volume (volume of a cylinder)
            vol = np.pi * h * ((R**2) - (r**2))
            Vapp = np.sum(vol)
            
        except:
            Vapp = np.nan

    return Vapp

'''
***********************************************************************
'''


def rimrim(model, step, id_mat = 0):
    
    '''
    description:
    calculate crater dimensions up to the rim crest
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : alt (float): Height of the rim over the pre-impact surface
    : drim (float): Rim-to-rim diameter
    : ix_rim (int): location of the crater rim (index along x-axis) 
    : iy_rim (int): location of the crater rim (index along y-axis)   
    ''' 

    # if time equal to 0, then all crater dimensions and indexes are equal to 0
    if step.time == 0:
        alt = 0
        drim = 0
        ix_rim = 0
        iy_rim = 0

    else:        
        try:                        
            # we neglect the 50 closest cells to the axis of symmetry
            conc = step.cmc[id_mat][:-50, :]
    
            # and we take only the partially filled cells (conc = concentration of
            # matters)
            empty_cells = np.where((conc > 0) & (conc < 1))
    
            # x- and y-indexes
            i = empty_cells[0]
            j = empty_cells[1]
    
            # creation of new indexes
            ik = np.array([])
            jk = np.array([])
    
            for idx, values in np.ndenumerate(j):
    
                # we loop through every index, and take the mean of 50 cells (i,j-50:j)
                # along the y-axis
                thres_d = np.nanmean(step.cmc[0][i[idx[0]], values-50:values])
    
                # if the concentration mean is larger than 99% in average and that
                # it is not close the axis symmetry, append index
                if (thres_d > 0.99) and (i[idx[0]] > 50):
    
                    jk = np.append(jk, values)
                    ik = np.append(ik, i[idx[0]])
    
            # if empty, the height and rim-to-rim crater diameter is equal to 0
            if (ik.size == 0) or (jk.size == 0):
                alt = np.nan
                drim = np.nan
                ix_rim = np.nan
                iy_rim = np.nan
    
            # otherwise, we take the maximum y-index (with the biggest height)
            else:
                idx_depth = np.where(jk == np.max(jk))[0]
                n = np.shape(idx_depth)[0]
                ix_depth = idx_depth[0]
    
                # if the maximum height is reached in one cell, we just take the value
                if n < 2:
                    ix_rim = ik[ix_depth].astype(int)
                    iy_rim = jk[ix_depth].astype(int)
                # if the maximum height is reached in several cells, we take the mean
                else:
                    ix_rim = np.int(np.mean(ik[idx_depth]))
                    iy_rim = jk[ix_depth].astype(int)
    
                # calculate the rim-to-rim crater diameter and rim height
                alt = np.abs((model.y[ix_rim, iy_rim+1]) -
                             ((1 - step.cmc[0][ix_rim, iy_rim]) * model.dy))
                drim = (np.abs(model.xc[ix_rim, iy_rim])) * 2.
                
        except:
            alt = np.nan
            drim = np.nan
            ix_rim = np.nan
            iy_rim = np.nan
            

    # return also the indexes
    return alt, drim, ix_rim, iy_rim


'''
***********************************************************************
'''


def rimvolume(model, step, id_mat=0):
    
    '''
    description:
    calculate the volume of the crater diameter up to the rim crest, Vrim
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : Vrim (float): Rim-to-rim crater volume
    ''' 
    
    # if this is the first time step, the volume is then equal to 0
    if step.time == 0:
        Vrim = 0
    else:
        try:           
            # get the indexes (ix, iy) for the rim location 
            __, __, ix, iy = rimrim(model, step, id_mat)
            
            # get the concentration in all cells up to the rim height (y)
            # and up to the rim-rim crater diameter (x)
            conc = step.cmc[id_mat][:ix, :iy]
    
            # select only cells that are not completelly filled with materials
            empty_cells = np.where(conc < 1)  # here I have the index
    
            # select the x-index, x-index + 1 and y-index
            i = empty_cells[0]
            ii = empty_cells[0] + 1
            j = empty_cells[1]
    
            # r is the distance from the axis of symmetry
            r = model.x[i, j]
    
            # h is the height of the cells, which is the same for all the cells
            h = model.dy
    
            # the distance along the x-axis between two cells
            # takes into account the concentration of materials
            R = model.x[i, j] + ((1 - step.cmc[id_mat][i, j])
                                 * (model.x[ii, j] - model.x[i, j]))
    
            # calculation of the volume (volume of a cylinder)
            vol = np.pi * h * ((R**2) - (r**2))
            Vrim = np.sum(vol)
            
        except:
            Vrim = np.nan

    return Vrim


'''
***********************************************************************
'''


def maxappdepth(model, step, id_mat=0):
    
    '''
    description:
    calculate the maximum apparent (from pre-impact surface) crater depth (da)
    across the whole cross section
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : da (float): maximum apparent crater depth
    '''     

    # if the simulation did not start (t = 0) or the apparent crater diameter
    # is not detected ix = 0, then the max app da is equal to 0.
    if step.time == 0:
        da = 0

    else:
        try:
            # get the y-index that correponds to the pre-impact surface
            oridy = np.where(model.yy == 0.)[0][0]
            
            # get the first non-empty cell along the pre-impact surface (x-axis)
            ix = ixdiam(model, step, id_mat)
                        
            # store how much the cells are filled with materials
            conc = step.cmc[id_mat][:ix, :oridy]
    
            # select cells, which are partially filled
            empty_cells = np.where((conc > 0) & (
                conc < 1))  # here I have the index
    
            # if there are no partially filled cells then da is equal to 0
            # I guess this condition is not likely to be True
            if empty_cells[0].size == 0:
                da = np.nan
    
            else:
                # get x- and y-index
                i = empty_cells[0]
                j = empty_cells[1]
    
                # create empty new x- and y-index
                ik = np.array([])
                jk = np.array([])
    
                for idx, values in np.ndenumerate(j):
    
                    # a threshold is required
                    # we use a 20-cell window to exclude possible ejected materials
                    thres_d = np.nanmean(step.cmc[0][i[idx[0]], values-20:values])
    
                    # if the 20 cells below the selected cells have more than 99%
                    # of materials in average, then i and j are added to ik and jk
                    if (thres_d > 0.99):
    
                        jk = np.append(jk, values)
                        ik = np.append(ik, i[idx[0]])
    
                # the index where the maximum da is reached can now be selected
                idx_da = np.where(jk == np.min(jk))[0][0]
    
                # conversion to integer
                iif = ik[idx_da].astype(int)
                ijf = jk[idx_da].astype(int)
    
                # the maximum apparent da is calculated, taking into account
                # the concentration of materials in the cell.
                da = np.abs((model.y[iif, ijf+1]) -
                               ((1 - step.cmc[0][iif, ijf]) * model.dy))
        except:
            da = np.nan

    return da


'''
***********************************************************************
'''


def midappdepth(model, step, id_mat=0):
    
    '''
    description:
    calculate the maximum apparent (from pre-impact surface) crater depth (da)
    along the y-axis of symmetry (first column). The result derived from midappdepth
    might be affected by modelling artificats along the axis of symmetry..
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : das (float): maximum apparent crater depth along the y-axis of symmetry
    '''
    
    # if first time step, then the depth is obviously equal to 0
    if step.time == 0:
        das = 0
    else:
        try:
            # the concentration of materials is only check along the axis of symmetry
            conc = step.cmc[id_mat][0, :]
    
            # a 20-cell window average is used
            ma = movingaverage(conc, 20)
    
            # as low indexes represent the deep material, we here use a threshold
            # where the boundary is detected when less than 1% of material is filled
            # in average over the 20-cell average (not the best way)
            ijf = np.where(ma < 0.01)[0][0]
    
            # the depth is then taken
            das = np.abs((model.y[0, ijf]) -
                           ((1 - step.cmc[0][0, ijf-1]) * model.dy))
        except:
            das = np.nan

    return das


'''
***********************************************************************
'''


def appdiameter(model, step, id_mat = 0):
    '''
    description:
    calculate the apparent crater diameter along the pre-impact surface (Da)
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : diam (float): apparent crater diameter (along the pre-impact surface)
    '''
    
    # if first time step, then the diameter is obviously equal to 0
    if step.time == 0:
        diam = 0
    else:
        try:
            # get the y-index that correponds to the pre-impact surface
            oridy = np.where(model.yy == 0.)[0][0]
            
            # get the first non-empty cell along the pre-impact surface (x-axis)
            ix = ixdiam(model, step, id_mat)
            
            # the boundary of the crater is previously detected in function ixdiam
            # the concentration of material(s) in the cell is here calculated
            cmc_diam = (1. - step.cmc[id_mat][ix-1, oridy]) * model.dx
    
            # the apparent crater diameter is then calculated
            diam = (np.abs(model.x[ix-1, oridy]) + cmc_diam) * 2.
    
            # if the apparent crater diameter calculated is larger than the
            # hi resolution than set value equal to 0
            if diam > (model.xhires[1] * 2.):
                diam = np.nan  # overwrite previously calculated value
        except:
            diam = np.nan

    return diam


'''
***********************************************************************
'''


def emptycells(model, step, id_mat = 0):
    
    '''
    description:
    calculate the number of:
        - empty cells (empty of materials, = void)
        - cells with concentration lower than 1 (partially filled)
        - completely filled (full of materials)
        
    This is only calculated for a distance up to the pre-impact surface
    
    param:
    : model: iSALE model (jdata.dat file)
    : step : step in iSALE (read with model.readStep)       
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.  
        
    returns:
    : empty (int): number of empty cells
    : partially_filled (int): number of partially filled cells
    : completely_filled (int): number of completely filled cells
    
    '''
    #
    if step.time == 0:
        empty = 0  # empty cells
        partially_filled = 0  # partially filled
        completely_filled = 0  # less than 1

    else:
        try:
            # get the y-index that correponds to the pre-impact surface
            oridy = np.where(model.yy == 0.)[0][0]
                       
            # get all cells under the surface
            cells_under_surface = step.cmc[id_mat][:, :oridy]
            
            # get number of cells that are ...
            empty = np.sum((cells_under_surface == 0))
            partially_filled = np.sum(np.logical_and(cells_under_surface > 0, 
                                                     cells_under_surface < 1))
    
            completely_filled = np.sum(cells_under_surface == 1)
        except:
            empty = np.nan 
            partially_filled = np.nan 
            completely_filled = np.nan 

    return (empty, partially_filled, completely_filled)



'''
***********************************************************************
'''


def dimensions(model, id_mat = 0):
    '''
    description:
    calculate all final crater dimensions for 1-layer target

    param:
    : model: iSALE model (jdata.dat file)     
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.
      
    returns:
    : xxx (float): all crater dimensions
    '''

    # define all crater variables
    vol = np.zeros(model.nsteps)
    depth = np.zeros(model.nsteps)
    mdepth = np.zeros(model.nsteps)
    diam = np.zeros(model.nsteps)
    time = np.zeros(model.nsteps)
    empty = np.zeros(model.nsteps)
    partially_filled = np.zeros(model.nsteps)
    completely_filled = np.zeros(model.nsteps)
    alt = np.zeros(model.nsteps)
    drim = np.zeros(model.nsteps)
    Vrim = np.zeros(model.nsteps)

    # loop through all saved timesteps
    for u in range(model.nsteps):

        # each timestep are saved in the variable step
        den = model.readStep('Den', u)
        step = copy.deepcopy(den)

        # time in crater evolution
        time[u] = step.time

        # if t = 0
        if time[u] == 0:
            vol[u] = 0
            depth[u] = 0
            mdepth[u] = 0
            diam[u] = 0
            empty[u] = 0
            partially_filled[u] = 0
            completely_filled[u] = 0
            alt[u] = 0
            drim[u] = 0
            Vrim[u] = 0
        else:   
            # Apparent crater dimensions
            diam[u] = appdiameter(model, step, id_mat)
            vol[u] = appvolume(model, step, id_mat)

            depth[u] = midappdepth(model, step, id_mat)
            mdepth[u] = maxappdepth(model, step, id_mat)

            # Rim-to-rim dimensions
            alt[u], drim[u], __, __ = rimrim(model, step, id_mat)
            Vrim[u] = rimvolume(model, step, id_mat)
            empty[u], partially_filled[u], completely_filled[u] = emptycells(model, step, id_mat)

    return (time, depth, mdepth, diam, vol, empty, partially_filled, 
            completely_filled, alt, drim, Vrim)


'''
***********************************************************************
'''

def to_txt(path_txt, output, header_txt):
    
    '''
    description:
        save output to text file
    param:
    : path_txt : absolute path to txt file
    : output : numpy array containing m x t values (m: parameters, t: timesteps)
    : header_txt : header txt to be saved in txt file
    '''
    
    fmt_list = ['%1.6e'] * len(header_txt.split('\t'))
    
    np.savetxt(path_txt, output, header=header_txt, delimiter='\t', fmt=fmt_list)

'''
***********************************************************************
'''


def main(path_jdata, path_tosave, id_mat = 0, transient = 0):
    
    
    '''
    description:
    calculate all final crater dimensions for 1-layer targets

    inputs:
    path: path containing folders with models (jdata in each modelled case)
    folders: all models
    paths: saving directory
    mode: mostly equal to 0 (if artifacts along the boundary, can be set to 1)
    id_mat: material id (0: all, 1: projectile, 2: first layer, 3: second ...)
    transient: 0 (by volume)
             : 1 (by the maximum depth along the crater cavity)
    outputs:
    1 .txt file with all crater dimensions every saved timesteps
    1. txt file with all transient crater dimensions
    1. txt file with all final crater dimensions

    each of these text files are saved in different folders
    so 3 txt files are created per models
    '''
    
    # In case only a single model is specified, transform string to list of strings
    if type(path_jdata) == str:
        path_jdata = [path_jdata]        
    else:
        None
                
    # header txts
    header_txt_evo = "time\ttime_norm\tdepth\tmdepth\tdiameter\tvolume\tempty\tpf\tcf\talt\tDrim\tVrim"
    header_txt_tr = "tr_time\ttr_depth\ttr_mdepth\ttr_diameter\ttr_vol\ttr_alt\ttr_Drim\ttr_Vrim\ttr_idx"
    
    # loop through models in folders
    for jdata_folder in path_jdata:

        # print modelname
        modelname = jdata_folder.split('/')[-2]
        print (modelname)

        # open model with pySALEplot
        jdata = os.path.join(jdata_folder, 'jdata.dat')
        model = psp.opendatfile(jdata)

        # calculate all final crater dimensions
        (time, depth, mdepth, diam, vol, empty, partially_filled, 
         completely_filled, alt, drim, Vrim) = dimensions(
            model, id_mat)
        
        # create plots directory if not existing
        for f in ["evolution", "transient"]:
            if not os.path.exists(os.path.join(path_tosave, modelname, f)):
                os.makedirs(os.path.join(path_tosave, modelname, f))
                
        # the transient crater is defined at the maximum volume
        if transient == 0:
            idx_tr = np.nanargmax(vol)
        # the transient crater is defined at the maximum depth    
        else:
            idx_tr = np.nanargmax(mdepth)  # if maximum depth is used
            
        # normalized time
        tnorm = time / time[idx_tr]
                
        #######################################################################
        # CALCULATE CRATER DIMENSIONS DURING CRATER EVOLUTION (EACH STEP)
        #######################################################################

        # we need to define the name of the txt file that will be saved
        fname = os.path.join(path_tosave, modelname, 'evolution', modelname + '_data.txt')

        # output file
        output = np.column_stack(
            (time, tnorm, depth, mdepth, diam, vol, empty, partially_filled, 
             completely_filled, alt, drim, Vrim))
        
        to_txt(fname, output, header_txt_evo)
        
        #######################################################################
        # CALCULATE CRATER DIMENSIONS AT THE TRANSIENT CRATER
        #######################################################################

        fname = os.path.join(path_tosave, modelname, 'transient', modelname + '_tr.txt')
            
        # output file
        output = np.column_stack(
            (time[idx_tr], depth[idx_tr], mdepth[idx_tr], diam[idx_tr], 
             vol[idx_tr], alt[idx_tr], 
             drim[idx_tr], Vrim[idx_tr], idx_tr))
        
        to_txt(fname, output, header_txt_tr)            
            
        # close model file
        model.closeFile()


'''
***********************************************************************
'''


def find_nearest(array, value):
    '''
    description:
    find the nearest value in an array

    output:
    value and index where the nearest value is found
    '''

    idx = np.nanargmin((np.abs(array-value)))
    return array[idx], idx


 
'''
***********************************************************************
'''

def crossSection(model, step, id_mat = 0):
    
    '''
    description:
        get cross section at the step "step" up to the end of high resolution zone
    '''
        
    # get the cross section across the whole high resolution zone
    ix = np.int(model.xhires[1] / model.dx)

    # get the concentration
    conc = step.cmc[0][:ix, :]  # + 1 before

    # define the x- and y-indexes
    ii = np.arange(ix)
    jj = np.array([])

    # values are append to y-index if.....
    for values in ii:

        ma = movingaverage(conc[values, :], 20)
        lll = np.where(ma < 0.01)[0]
        if lll.size == 0:
            tijf = np.where(ma == np.min(ma))[0][0]
        else:
            tijf = lll[0]
        jj = np.append(jj, tijf)

    # create empty arrays
    jj = jj.astype(int)
    altr = np.array([])
    rrim = np.array([])

    # calculate the height and distance
    for values in np.arange(len(jj)):

        a1 = (model.y[ii[values], jj[values]+1]) - \
            ((1 - step.cmc[0][ii[values], jj[values]]) * model.dy)
        a2 = (np.abs(model.xc[ii[values], jj[values]]))
        altr = np.append(altr, a1)
        rrim = np.append(rrim, a2)
            
            
    return (rrim, altr)
'''
***********************************************************************
'''


def crossSection_atstep(path_jdata, path_tosave, idx, id_mat = 0):
    '''
    description:
        extract elevation and distance points for a step in iSALE or for a
        specific time (normalized to the transient crater time) and save it to 
        the crossSection folder of this modelrun.
        
    param:
    : path_jdata (list) : list of absolute paths to jdata.dat files (i.e., models)
    : path_tosave (str) : absolute path to the main plot folder
    : idx : iSALE time step to save (if >= 0) or normalized time (norm to transient)
      if equal to -1
    : id_mat (int): material id in iSALE, 0 corresponds to all materials, if 
      id_mat is different than 0 (e.g., id_mat = 1), only this material is
      selected within the target. set to 0 as default.

    '''
    
    # In case only a single model is specified, transform string to list of strings
    if type(path_jdata) == str:
        path_jdata = [path_jdata]        
    else:
        None
        
    # header_txt
    header_txt = "distance\tdepth"
                    
    # loop through models in folders
    for jdata_folder in path_jdata:

        # print modelname
        modelname = jdata_folder.split('/')[-2]
        print (modelname)
        
        
        # check if the parameters for the evolution have been calculated
        evo_file = os.path.join(path_tosave, 
                                       modelname, 
                                       'evolution', 
                                       modelname + '_data.txt')
        
        path_cS = os.path.join(path_tosave, modelname, 'crossSections')
        if not os.path.exists(path_cS):
                os.makedirs(path_cS)
        else:
            None
                    
        
        if os.path.isfile(evo_file):
            
            jdata = os.path.join(jdata_folder, 'jdata.dat')
            model = psp.opendatfile(jdata)
            
            # load evolution data
            data = np.loadtxt(evo_file)
            tnorm = data[:,1]
            
            # name of the cross section               
            if idx == -1:
                ix = []
                crossSection_fname = []
                for i in range(1,np.int(np.nanmax(tnorm)) + 1):
                    ix.append(find_nearest(tnorm,i)[1])
                    crossSection_fname.append(os.path.join(path_cS, 'crossSection_tnorm_' + str(i).zfill(1) + '.txt'))                    
            else:
                if type(idx) == str: 
                    ix = [idx]
                else:
                    ix = idx
                crossSection_fname = []
                for i in ix:    
                    crossSection_fname.append(os.path.join(path_cS, 'crossSection_idx_' + str(i).zfill(3) + '.txt'))  
                    
            # open model with pySALEplot
            jdata = os.path.join(jdata_folder, 'jdata.dat')
            model = psp.opendatfile(jdata)
            
            for k, i in enumerate(ix): 
                # step
                step = model.readStep('Den',i)
                (rrim, altr) = crossSection(model, step, id_mat)
                
                output = np.column_stack((rrim, altr))
                
                to_txt(crossSection_fname[k], output, header_txt)
    
            # close model file
            model.closeFile()
                
        else:
            print ("crater.py has to be run first")



    print ("crater profiles have been saved")
    
'''
***********************************************************************
'''

def step(jdata, nDT_SAVE, norm=False):
    
    '''
    returns the step corresponding to the wished saved timestep. 
    This function is quite pratical if you want to work with a specific timestep
    in for example jupyter
    
    params:
    : jdata (str):
    : nDT_SAVE (int): 
        
    Should have a normalization by the transient crater time option 
    (if path_tosave is specified and the transient crater has been calculated)
    then it is pretty straightforward ....
    '''
    
    model = psp.opendatfile(jdata)
    
    # if we want to get normalized data
    if norm:
        None
              
    # load all the data available (except V_ and tracers)
    fieldiSALE = []
    
    for f, f2 in model.fieldlist:
        if not ((f.startswith('V_')) or (f.startswith('Tr'))):
            fieldiSALE.append(f)
            
    step = model.readStep(fieldiSALE, nDT_SAVE)
    
    return (model, step)
    
    
    
    
    
    
