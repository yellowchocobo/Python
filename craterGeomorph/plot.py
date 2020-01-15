# -*- coding: utf-8 -*-

# loading of basic Python's module
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import glob


pathm = ['/uio/kant/geo-ceed-u1/nilscp/Nils/Python/craterGeomorph', 
         '/work/nilscp/iSALE/Dellen/lib']

for pat in pathm:
    sys.path.append(pat)

import pySALEPlot as psp
import subprocess
'''
***********************************************************************
'''

# Switch on the use of Latex text
rc('text', usetex=True)
rc('axes', linewidth=2)
rc('font', weight='bold')

'''
***********************************************************************
'''


def make_video(path, delay_number, loop_number, videolabel):
    '''
    description:

    path = '/mypath/' # note the slash at the end
    delay_number: in hundreds of a second
    loop: = 0 infinite looping, 1 = once, 2 = twice ....
    videoname = 'myvideo.gif'  

    example:
    path = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILWPO_L100/evolution/plots/Distension/all/'
    delay_number = 10
    loop_number = 1
    videoname = 'distension'
    make_video(path, delay_number, loop_number, videoname)
    '''
    os.chdir(path)

    videoname = videolabel + '_delay' + \
        str(int(delay_number)) + '_loop' + str(loop_number) + '.gif'

    command = ("convert -delay " + str(delay_number) + " -loop " + str(loop_number) + " " + str(path) +
               "*.png " + videoname)

    subprocess.Popen(command.split(), cwd=path)

    print ('Steven Spielberg')


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

def zoom(model, zoom_level, scaling_factor):
    
    '''
    description:
        
    params:
        
    returns:
    extentx
    extenty
    
    
    '''
    # in case for 
    if type(zoom_level) == tuple:
        if zoom_level[0][1] > 0:
            extentx = zoom_level[0]
            extenty = zoom_level[1]
        else:
            extentx = (0.0, np.abs(zoom_level[0][1]) * scaling_factor)
            
            # first is negative
            extenty = (zoom_level[1][0] * scaling_factor, 
                       np.abs(zoom_level[1][1]) * scaling_factor)
            
    elif type(zoom_level) == str:
        
        if zoom_level == 'complete_grid':
            extentx = (0.0, np.max(model.x) / scaling_factor)
            
        elif zoom_level == 'hires_grid':
            extentx = (0.0, model.xhires[1] / scaling_factor)
            extenty = (model.yhires[0]/ scaling_factor, model.yhires[1] / scaling_factor)
            
        elif zoom_level == 'medium':
            extentx = (0.0, (model.xhires[1] / scaling_factor) / 2.0)
            extenty = ((model.yhires[0]/2.0) / scaling_factor, model.yhires[1] / scaling_factor)
            
        elif zoom_level == 'close_up':
            extentx = (0.0, (model.xhires[1] / scaling_factor) / 4.0)
            extenty = ((model.yhires[0]/4.0) / scaling_factor, model.yhires[1]/ scaling_factor)
        
        else:
            print ("zoom level not recognized. Please choose from complete grid"
                   ", hires_grid, medium, close_up or specify values as tuple")
            
    return ((extentx, extenty))

'''
***********************************************************************
'''

def truncate(n, decimals=0):
    
    '''
    see website: https://realpython.com/python-rounding/
    
    >>> truncate(-5.963, 1)
    -5.9
    
    
    >>> truncate(-1374.25, -3)
    -1000.0
    '''
    
    multiplier = 10 ** float(decimals)
    return int(n * multiplier) / multiplier

   
'''
***********************************************************************
'''

def getTransientvalues(path_tosave, modelname):
    
    '''
    description:
    
    param:
    :path_tosave (str):
    : modelname (str):
        
    returns:
    : tr_time : transient crater time
    : tr_depth : transient crater depth (along the axis of symmetry)
    : tr_mdepth : transient crater maximum depth (across the whole cross section)
    : tr_diameter : transient crater diameter
    : tr_vol : transient crater volume
    : tr_alt : altitude of the transient rim
    : tr_Drim: rim-to-rim transient diameter
    : tr_Vrim : rim-to-rim transient volume
    : tr_idx : timestep at which the transient crater is reached
    '''
    
    try:
        data = np.loadtxt(path_tosave + modelname + '/transient/' + 
                          modelname + '_tr.txt',
                          delimiter='\t', comments='#')
                    
        return (data)
    
    except:
        
        if os.path.isdir(path_tosave + modelname):
            print ("transient crater dimensions have not yet been calculated")
        else:
            print ("modelname does not exists")
                          
'''
***********************************************************************
'''            
    
def getvalues(path_tosave, modelname, step):
    
    '''
    description:
    
    param:
    :path_tosave (str):
    : modelname (str):
    : step (int): if negative, normalized to the transient crater time
    
    returns:
    : time : absolute time (since impact)
    : time_norm : normalized time (by transient crater time)
    : depth : crater depth (along the axis of symmetry)
    : mdepth : crater maximum depth (across the whole cross section)
    : diameter : crater diameter at the pre-impact surface
    : volume : crater volume at the pre-impact surface
    : empty : Number of void cells within the apparent crater
    : pf : Number of partially filled cells within the apparent crater
    : cf : Number of completely filled cells within the apparent crater
    : alt : Altitude of the crater rim
    : DRim : Diameter of the crater rim
    : Vrim : Volume of the crater rim

    '''
    
    try:
        data = np.loadtxt(path_tosave + modelname + '/evolution/' + 
                          modelname + '_data.txt',
                          delimiter='\t', comments='#')
        
        if step < 0:
            
            # load normalized time
            tnorm = data[:,1]
            
            # should have something that check whether it is a close value or not
            # get the closest value to the specified step
            __, tstep = find_nearest(tnorm, np.abs(step))
            
        else:
            tstep = np.int(step)
            
        return (data[tstep, :])
            
    except:
        
        if os.path.isdir(path_tosave + modelname):
            print ("crater dimensions have not yet been calculated")
        else:
            print ("modelname does not exists")
            
'''
***********************************************************************
''' 

def normalize(jdata, model, path_tosave, norm, norm_tend):
    
    '''
    '''
    
    
    # modelname
    modelname = jdata.split('/')[-2]
    
    # loading of the scaling factor   
    if norm == 0:
        scaling_factor = 1.0
                
    # normalized by the radius of the projectile
    elif norm == 1:
        scaling_factor = model.cppr[0] * model.dx
    
    # normalized by the apparent transient crater diameter     
    elif norm == 2:
        scaling_factor = getTransientvalues(path_tosave, modelname)[3]
        
    # normalized by the final rim-to-rim crater diameter     
    elif norm == 3:
        scaling_factor = getTransientvalues(path_tosave, modelname, norm_tend)[-2]
        
    # scale by a certain value (if for example if you want to change m to km)
    elif norm == 4:
        scaling_factor = 1000.0
        
    # normalized by the final rim-to-rim crater diameter of another crater
    else:
        try:
            scaling_factor = getTransientvalues(path_tosave, norm.split('/')[-2], norm_tend)[-2]
        except:
            print ("modelname does not exist")
            
            
    return (scaling_factor)
    

'''
***********************************************************************
''' 

def getTransientcrossSection(path_tosave, modelname):
    
    '''
    get transient crater cross sections and return the x and y coordinates for
    the surface. Only values smaller than the apparent crater diameter are shown
    '''
    
    try:
        path = os.path.join(path_tosave, modelname, 'crossSections')
        fname = os.path.join(path, 'crossSection_tnorm_1.txt')
        data = np.loadtxt(fname,
                          delimiter='\t', comments='#')
        
        transient_app_diameter = getTransientvalues(path_tosave, modelname)[3]
        
        X = data[:, 0]
        Y = data[:, 1]

        return (X[X <= transient_app_diameter /2.0], 
                Y[X <= transient_app_diameter /2.0])          
                                              
    except:
        if os.path.isdir(path):
            print ('the cross section for the transient crater '
                   'has not yet been generated')


'''
***********************************************************************
'''

def time_lbl(path_tosave, modelname, time, norm, shownormTime):
    
    '''
    Need to modify this part if you don't want to plot the time and norm time
    
    what if we want to have the time related to the penetration time?
    '''
    
    tr_time = getTransientvalues(path_tosave, modelname)[0]
    tr_idx = np.int(getTransientvalues(path_tosave, modelname)[-1])
    
    if shownormTime:
    
        lbl_time = (r"$\mathit{t} = " + "{:.2f}".format(time) + " " + " s, $" + " " +
                                  "$\zeta$ = " + "{:.2f}".format(time/tr_time))
        
    else:
        lbl_time = (r"$\mathit{t} = " + "{:.2f}".format(time) + " " + " s$")
    
    
    return (tr_idx, lbl_time)
    
    
    
'''
***********************************************************************
'''

def figurename(jdata, path_tosave, step, param, norm, zoom_level, dual):
    
    '''
    description:
        give a figurename based on modelname, parameter and type of normalization
        
        
    ppp_nnnnn_Zzzz_xxx.jpg 

    ppp : parameter (DEN, TMP .....)
    nnnnn : type of normalization (TRANS, RPROJ, NORMA, FINAL, FINCO)
    zzzzz : type of zoom (COMPL, HIRES, MIDDL, CLOSE, MANUA)
    xxx : timestep
    
    DEN_TRANS_HIRES_000.png
    
    In evolution/plots/parameter
    
    params:
    : jdata :
    : path_tosave :
    : param :
    : norm :
        
    dual 
    '''
    if dual:
        try:
            P =(param[0] + param[1]).upper()
        except:
            P =(param[0] + param[0]).upper()
    else:
        # parameter name
        P = param.upper()
        
    # modelname
    modelname = jdata.split('/')[-2]
        
    # step
    S = str(step).zfill(3)
    
    # normalization name
    if ((norm == 0) or (norm == 4)):
        N = "NORMA"                
    elif norm == 1:
        N = "RPROJ" 
    elif norm == 2:
        N = "TRANS" 
    elif norm == 3:
        N = "FINAL"        
    else:
        N = "FINCO"
        
    # zoom name
    if type(zoom_level) == tuple:
        Z = "MANUA"
           
    elif type(zoom_level) == str:
        
        if zoom_level == 'complete_grid':
            Z = "COMPL"            
        elif zoom_level == 'hires_grid':
            Z = "HIRES"           
        elif zoom_level == 'medium':
            Z = "MIDDL"                  
        elif zoom_level == 'close_up':
            Z = "CLOSE"           
   

    dirname = os.path.join(path_tosave,
                           modelname,
                           'evolution',
                           'plots', P)
    
    # create if non-existing
    if os.path.isdir(dirname):
        None
    else:
        os.makedirs(dirname)
                           
    figname = os.path.join(dirname, 
                           modelname + "_" + P + "_" + N + "_Z" + Z + "_" + S + ".png") 
    
    return (figname)
           
'''
***********************************************************************
''' 
def adjust_ticks_and_labels(axes, x_lbl, y_lbl):
    
    '''
    detect automatically if one or two axes
    
    '''
    
    try:
        n = len(axes)       
    except:
        n = 1
        axes = [axes]
  
    
    # remove division in the middle for dual plots
    if n == 2:
        ax1, ax2 = axes
        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax1.yaxis.set_label_position("left")
        ax2.yaxis.set_label_position("right")
        #ax1.set_ylabel(y_lbl, fontsize=16, labelpad=-26) #r"$y / R_{r}$"
        #ax2.set_ylabel(y_lbl, rotation=270, fontsize=16, labelpad=-10) #r"$y / R_{r}$"
        #ax2.set_xlabel(x_lbl, fontsize=16)

    else:
        None
                
    # set labelsize to 14
    for ax in axes:
        ax.minorticks_off()
        ax.tick_params('both', labelsize=14, length=5)
        
        # set labels
        ax.set_xlabel(x_lbl, fontsize=16)
        ax.set_ylabel(y_lbl, fontsize=16)

'''
***********************************************************************
'''     
def adjust_time(axes, t_lbl):
    
    '''
    detect automatically if one ax or two axes

    

    Parameters
    ----------
    axes : TYPE
        DESCRIPTION.
    t_lbl : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    try:
        n = len(axes)       
    except:
        n = 1
        ax = axes
    
    # if dual, it would be better to put the time in the middle
    if n == 2:
        ax1, ax2 = axes
        ax2.set_title(t_lbl, position=(0.0, 0.9), fontsize=18)
    else:
        ax.set_title(t_lbl, position=(0.5, 0.9), fontsize=18)

'''
***********************************************************************
''' 

def adjust_zoom(axes, extentx, extenty):
    
    '''
    '''
    
    
    try:
        n = len(axes)       
    except:
        n = 1
        ax = axes
    
    if n == 2:
        ax1, ax2 = axes
        
        ax1.set_xlim((-extentx[1], extentx[0]))
        ax2.set_xlim(extentx)
        
        ax1.set_ylim(extenty)
        ax2.set_ylim(extenty) 
        
    else:
        ax.set_xlim(extentx)
        ax.set_ylim(extenty)        
        
'''
***********************************************************************
'''     
def adjust_colorbar(fig, ax, pco, sign, colormap_lbl, dual):
    
    '''
    Adjust colorbar depending on type of plot    
    '''
    
    # colorbar, colorbar positioning 
    ax1_divider = make_axes_locatable(ax)
    if sign > 0:
        if dual:
            ax.yaxis.tick_right() # maybe should go somewhere else
            cax = ax1_divider.append_axes("right", size="5%", pad="15%")
            cbar = fig.colorbar(pco, cax=cax)
            cax.yaxis.set_ticks_position("right")
            cbar.ax.set_title('$' + colormap_lbl + '$', fontsize=14)
        else:
            cbar = fig.colorbar(pco)
            cbar.ax.set_title('$' + colormap_lbl + '$', fontsize=14)
    else:
        
        cax = ax1_divider.append_axes("left", size="5%", pad="15%")
        cbar = fig.colorbar(pco, cax=cax)
        cax.yaxis.set_ticks_position("left")
        cbar.ax.set_title('$' + colormap_lbl + '$', fontsize=14)
        #cbar.set_label('$' + colormap_lbl + '$', labelpad=-6.0, 
        #       y=1.07, fontsize=14, rotation=0)
        
    cbar.ax.tick_params(labelsize=14)
    
'''
***********************************************************************
''' 

def p(path_tosave, modelname, model, step, param, fig, ax, scaling_factor, 
      vmin, vmax, colormap, colormap_lbl, param_unit_factor, sign, 
      stepiSALE, stepTransient, dual, 
      showSurface, showTransient, showTracers):
    
    # plot data
    pco = ax.pcolormesh(sign * model.x/scaling_factor, 
                         model.y/scaling_factor, 
                         param/param_unit_factor,
                             vmin=vmin,
                             vmax=vmax,
                             cmap=colormap, zorder=2)
    
    # showing surface
    if showSurface:
        ax.contour(sign * model.xc/scaling_factor, 
                    model.yc/scaling_factor,
                    step.cmc[0], 1, colors='k', linewidths=2, 
                    zorder=4)
    
    # showing transient only if below transient
    if showTransient:
        if stepiSALE < stepTransient:
            None
        else:
            # this values need to be scaled
            xtr_surf, ytr_surf = getTransientcrossSection(path_tosave, modelname)
            
            ax.plot(sign * xtr_surf/scaling_factor, ytr_surf/scaling_factor,
                'r', linewidth=2, zorder=5)
    
    # showing tracers        
    if showTracers:                    
        for u in range(model.tracer_numu):
            tru = model.tru[u]
            # Plot the tracers in horizontal lines, every 20 lines
            for l in np.arange(0, len(tru.xlines), showTracers):
                ax.plot(sign *step.xmark[tru.xlines[l]]/scaling_factor,
                         step.ymark[tru.xlines[l]]/scaling_factor,
                         c='k', marker='.', linestyle='None', 
                         markersize=1.0, zorder=3)

            # Plot the tracers as vertical lines, every 20 lines
            for l in np.arange(0, len(tru.ylines), showTracers):
                ax.plot(sign *step.xmark[tru.ylines[l]]/scaling_factor,
                         step.ymark[tru.ylines[l]]/scaling_factor,
                         c='k', marker='.', linestyle='None', 
                         markersize=1.0, zorder=3)
                
    
    # colorbar, colorbar positioning
    adjust_colorbar(fig, ax, pco, sign, colormap_lbl, dual)
    
    return (pco)
    
        
'''
***********************************************************************
''' 
    

def field(jdata, path_tosave, steps, params, norm, zoom_level, vmin, vmax, title_lbl,
           x_lbl, y_lbl, 
           colormap, colormap_lbl, param_unit_factor, norm_tend,
           showhline, showTransient, showSurface, 
           showTracers, shownormTime, dual):
    
    '''
    Important!
    vmin has to take several inputs
    vmax as well
    title_lbl as well
    colormap
    colormap_lbl
    param_unit_factor
    
    param:
    : jdata (str):
    : path_tosave (str):
    : steps (int):
    : param (list of str): ['Den', 'Tmp']
    : norm (str or float): if str specify jdata path
    : zoom (str or float):
    : vmin (float):
    : vmax (float):
    : lbl (str):
    : showhline (flag):
    : showTransient (flag):
    : showSurface (flag):
    : showTracers (flag) to add
    : dual (flag) to add
        
    returns:
    plot(s) in path_tosave
    
    
    example:
    jdata = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/structuralUplift/SHROEDINGER_CRUST_20km/jdata.dat'
    path_tosave = '/uio/kant/geo-ceed-u1/nilscp/Desktop/astra/iSALE/structuralUplift/'
    steps = -1
    params = ['Den', 'Tmp']
    param_unit_factor = [1.0, 1.0]
    dual = True
    norm = 1
    zoom_level = 'medium'
    vmin = [0, 1000]
    vmax = [3000, 2000]
    title_lbl = 'Shroedinger (T = 20 km), density'
    x_lbl = 'x'
    y_lbl = 'y'
    colormap = ['viridis', 'viridis']
    showhline = True
    showTransient = True
    showSurface = True
    showTracers = False # or int
    colormap_lbl = ['kg/m^{3}', 'K']
    
    field(jdata, path_tosave, steps, params, norm, zoom_level, vmin, vmax, title_lbl,
           x_lbl, y_lbl, 
           colormap, colormap_lbl, param_unit_factor, norm_tend = -5.0,
           showhline = False, showTransient = False, showSurface = False, 
           showTracers = False, dual = True)
    
    params = ['Den']
    param_unit_factor = 1.0
    dual = True
    norm = 1
    zoom_level = 'medium'
    vmin = 0
    vmax = 3000
    title_lbl = 'Shroedinger (T = 20 km), density'
    x_lbl = 'x'
    y_lbl = 'y'
    colormap = 'viridis'
    showhline = True
    showTransient = False
    showSurface = True
    showTracers = 20 # or int
    colormap_lbl = 'kg/m^{3}'
    dual = False
    
    field(jdata, path_tosave, steps, params, norm, zoom_level, vmin, vmax, title_lbl,
           x_lbl, y_lbl, 
           colormap, colormap_lbl, param_unit_factor, norm_tend = -5.0,
           showhline, showTransient, showSurface, 
           showTracers, dual)
    
    What if we want to plot only the cmc, it is a bit tricky:
        could have a list of parameters Cm1, Cm2 and so forth..
        automatically detect all the Cm fieldvar and merge it to a single array
        
    cmc = []
    nmat = []
    values = []
    for i, j in model.fieldlist:
        if i.startswith('Cm'):
            cmc.append(i)
            nmat.append(i[-1])
    
    step.readStep
    # not sure because it can be percentage of materials
    '''
    
    # figures does not pop up
    plt.ioff()
    
    # open jdata
    model = psp.opendatfile(jdata)
    
    # modelname
    modelname = jdata.split('/')[-2]
    
    # get scaling factor based on if it will be normalized or not
    scaling_factor = normalize(jdata, model, path_tosave, norm, norm_tend)
       
    # get extent
    extentx, extenty = zoom(model, zoom_level, scaling_factor)
        
    # managing values for steps
    if steps == -1:
        steps = range(model.nsteps)
    elif type(steps) == int:
        steps = [steps]
    else:
        None
    
    for i in steps:
            step = model.readStep(params, i)  
            t1 = step.time
            
            # set time label
            stepTransient, t_lbl = time_lbl(path_tosave, modelname, t1, norm, 
                                            shownormTime)

            
            # dual 
            if dual:
                
                # get figurename based on the different selected parameters (have to include somehow dual)
                figname = figurename(jdata, path_tosave, i, params, norm, zoom_level, dual) # should only be two
                    
                # fix plot and other related plot layout
                fig, ax = plt.subplots(figsize=(12, 6), nrows = 1, ncols = 2)
                fig.subplots_adjust(wspace=0.0)
                
                # plots based on preferences
                cpo1 = p(path_tosave, modelname, model, step, step.data[0], 
                        fig, ax[0], scaling_factor, vmin[0], vmax[0], colormap[0], 
                        colormap_lbl[0], param_unit_factor[0], -1.0, 
                        i, stepTransient, dual, 
                        showSurface, showTransient, showTracers)
                
                # dual case where it is the same param in left and right panels
                if len(params) == 1:
                    cpo2 =p(path_tosave, modelname, model, step, step.data[0], 
                            fig, ax[1], scaling_factor, vmin[0], vmax[0], colormap[0], 
                            colormap_lbl[0], param_unit_factor[0], 1.0, 
                            i, stepTransient, dual, 
                            showSurface, showTransient, showTracers)
                    
                # left and right panels have a different parameter
                else:
                    cpo2 =p(path_tosave, modelname, model, step, step.data[1], 
                            fig, ax[1], scaling_factor, vmin[1], vmax[1], colormap[1], 
                            colormap_lbl[1], param_unit_factor[1], 1.0, 
                            i, stepTransient, dual, 
                            showSurface, showTransient, showTracers)
                
                # adjust zoom (for dual case)
                adjust_zoom(ax, extentx, extenty)
                
                # adjust time
                adjust_time(ax, t_lbl)

                # adjust ticks and labels
                adjust_ticks_and_labels(ax, x_lbl, y_lbl)
                
                # tight_layout               
                fig.tight_layout()
                fig.subplots_adjust(wspace=0.0)
                
                # main title
                st = fig.suptitle(title_lbl, fontsize=20) # main title
                st.set_y(0.97)
                st.set_x(0.5)
                fig.subplots_adjust(top=0.9)
                
                # saving figure
                fig.savefig(figname, dpi=300)   
                plt.close()
                
            else:
                
                for iparam, param in enumerate(params):
                    
                    # one figure per param
                    fig, ax = plt.subplots(figsize=(6, 6), nrows = 1, ncols = 1)
                    

                    # get figurename based on the different selected parameters
                    figname = figurename(jdata, path_tosave, i, param, norm, zoom_level, dual)
                    
                    p(path_tosave, modelname, model, step, step.data[iparam], 
                      fig, ax, scaling_factor, vmin[iparam], vmax[iparam], colormap[iparam], 
                      colormap_lbl[iparam], param_unit_factor[iparam], 1.0, 
                      i, stepTransient, dual, 
                      showSurface, showTransient, showTracers)
                    
                    # adjust zoom
                    adjust_zoom(ax, extentx, extenty)
                    
                    # adjust time
                    adjust_time(ax, t_lbl)
    
                    # adjust ticks and labels
                    adjust_ticks_and_labels(ax, x_lbl, y_lbl)
                    
                    # tight_layout               
                    fig.tight_layout()
                    fig.subplots_adjust(wspace=0.0)
                    
                    # main title
                    st = fig.suptitle(title_lbl[iparam], fontsize=20) # main title
                    st.set_y(0.97)
                    st.set_x(0.5)
                    fig.subplots_adjust(top=0.9)
                    
                    # saving figure
                    fig.savefig(figname, dpi=300)  
                    plt.close()
                    
    model.closeFile()


'''
***********************************************************************
'''

def peakShockPressureZone(jdata):
        
    # Open the datafile
    model=psp.opendatfile(jdata)
        
    step0 = model.readStep(['Den','TrP'],0)
    endstep = model.readStep(['Den','TrP'],model.laststep)
    
    # I could save these two variables in a text file (quite lot of data though)
        
    return model, step0, endstep

'''
***********************************************************************
'''

def plotPressureZone(model,step0,endstep,maxP):

    # [model.tru[u].truInfo() for u in range(model.tracer_numu)]
    
    
    xli = ((np.min(step0.xmark),np.max(step0.xmark)))
    yli = ((np.min(step0.ymark),np.max(step0.ymark)))
    
    plt.scatter(step0.xmark,step0.ymark,c=endstep.data[1]*1e-9,
                cmap='jet',vmin=0,vmax=maxP,s=4,linewidths=0)
    
    #plt.scatter(endstep.xmark,endstep.ymark,c=endstep.data[1]*1e-9,
    #            cmap='jet',vmin=0,vmax=maxP,s=4,linewidths=0)
                
    plt.xlim(xli)
    plt.ylim(yli)
    plt.colorbar()
    
'''
***********************************************************************
'''

def t(model, step, param, idxTracers, fig, ax, scaling_factor, 
      vmin, vmax, colormap, colormap_lbl, param_unit_factor, sign, 
      dual):
    
    """
    fig, ax = plt.subplots(figsize=(6, 6), nrows = 1, ncols = 1)
    
    
    """
    
        # plot data
    if idxTracers:
        pco = ax.scatter(sign * step.xmark[idxTracers]/scaling_factor, 
                             step.ymark[idxTracers]/scaling_factor, 
                             c = param/param_unit_factor,
                                 vmin=vmin,
                                 vmax=vmax,
                                 cmap=colormap, zorder=2)
    else:
        # showing tracers every ..(start from here again, does not work ....)
        if showTracers:                    
            for u in range(model.tracer_numu):
                tru = model.tru[u]
                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.xlines), showTracers):
                        ax.scatter(sign *step.xmark[tru.xlines[l]]/scaling_factor,
                                 step.ymark[tru.xlines[l]]/scaling_factor,
                                 c = param[l]/param_unit_factor,
                                     vmin=vmin,
                                     vmax=vmax,
                                     cmap='viridis', zorder=2) 
    
                # Plot the tracers as vertical lines, every 20 lines
                for l in np.arange(0, len(tru.ylines), showTracers):
                    ax.scatter(sign *step.xmark[tru.ylines[l]]/scaling_factor,
                             step.ymark[tru.ylines[l]]/scaling_factor,
                             c = param[l]/param_unit_factor,
                                 vmin=vmin,
                                 vmax=vmax,
                                 cmap=colormap, zorder=2)
        
        pco = ax.scatter(sign * step.xmark[:]/scaling_factor, 
                             step.ymark[:]/scaling_factor, 
                             c = param/param_unit_factor,
                                 vmin=vmin,
                                 vmax=vmax,
                                 cmap=colormap, zorder=2)        
    
                
    
    # colorbar, colorbar positioning
    adjust_colorbar(fig, ax, pco, sign, colormap_lbl, dual)
    

def tracers():
    
    """
    norm
    norm_tend
    shownormTime
    zoom_level
    param_unit_facor
    vmin
    vmax
    title_lbl
    colormap
    colormap_lbl
    dual (but only showing the same)
    figsize (also in field)
    fontsize (also in field)
    """
    None

'''
***********************************************************************
'''

# TO DO: improve

def excavated(pathdata, modelname, pathplots):
    '''
    description:
    routine that plots the excavated crater dimensions

    inputs:
    pathdata: path where the text files have been stored (text files that are obtained
    from running the main.py script)

    folders: all modelcases that need to be plotted

    pathplots: saving plot directory

    example:
    pathdata = '/media/nilscp/Zell/Collapse/data/'
    modelname = 'C00P20F08_L250'
    pathplots = '/media/nilscp/Zell/Collapse/rplots/'
    excavated(pathdata,modelname,pathplots)
    '''

    # should be outside of the loop
    if not os.path.exists(pathplots):
        os.makedirs(pathplots)

    os.chdir(pathdata + modelname + '/excavated/')

    # load excavated data
    # possible miss something at the end
    Ve, de, De = np.loadtxt(modelname + '_excavated.txt',
                            delimiter='\t', comments="#")
    Z = np.loadtxt(modelname + '_tracersXY_excavated.txt', delimiter='\t',
                   comments="#")  # possible miss something at the end
    X = Z[:, 0]
    Y = Z[:, 1]

    # load transient profile
    os.chdir(pathdata + modelname + '/transient/')
    dataXY = np.loadtxt(modelname + '_XYtransientprofile.txt',
                        delimiter='\t', comments='#')
    X2 = dataXY[:, 0]
    Y2 = dataXY[:, 1]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.plot(X, Y, "bo")
    ax.plot(X2, Y2, "ro")
    ax.hlines(de, np.min(X2), np.max(X2), 'b', linewidth=3,
              label='de: apparent excavated depth')
    ax.plot(De/2., 0, 'yo', ms=10, label='Re: apparent excavated radius')
    ax.legend(loc='lower right')
    ax.set_xlabel('X (m)', fontsize=18)
    ax.set_ylabel('Y (m)', fontsize=18)
    ax.tick_params('both', labelsize=16, length=10, width=1., which='major')
    fig.tight_layout()
    st = fig.suptitle(
        modelname + ' excavated and transient profile', fontsize=18)
    st.set_y(0.95)
    fig.subplots_adjust(top=0.9)
    figtitle = modelname + ' excavated_profile.png'
    fig.savefig(pathplots+figtitle, dpi=300)
    plt.close()

'''
***********************************************************************
'''


'''
example:
param = 'Por'
paths = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILWPO_L100/'
zoom_id = ['mid', 'hires']

fld_param, fld_name, fld_cmap, fld_unit, fld_factor, npath = field_definition(model, zoom_id, paths, param = 'all')

#Updated field_definition the 24th June 2019



fld = np.array(['Cm1', 'Cm2', 'Den', 'Tmp', 'Pre', 'Sie', 'Yld',
                'Dam', 'Alp', 'TPS', 'YAc'])

fld_param = np.array(['Cm1', 'Cm2', 'Den', 'Tmp', 'Pre', 'Sie', 'Yld',
                      'Dam', 'Alp', 'TPS', 'YAc'])

fld_name = np.array(['material 1', 'material 2', 'density', 'temperature', 'pressure', 'sie',
                     'YS', 'damage', 'distension', 'TPS',
                     'YAc'])

fld_folder = np.array(['Cm1', 'Cm2', 'Density', 'Temperature', 'Pressure', 'Specific_internal_energy',
                       'Yield_strength', 'Damage', 'Distension', 'Total_plastic_strain',
                       'Acoustic_fluidisation_strength'])

fld_cmap = np.array(['BrBG_r', 'BrBG_r', 'viridis', 'autumn', 'bwr', 'bwr',
                     'summer', 'winter', 'gray_r', 'Reds_r',
                     'summer'])

fld_unit = np.array(['Conc. material','Conc. material','kg/m^{3}', 'Kelvin', 'GPa', 'MJ',
                     'GPa', 'Damage', 'Distention', 'Strain',
                     'MPa'])

fld_factor = np.array([1., 1., 1., 1., 1e9, 1e6, 1e9, 1., 1., 1., 1e6])

'''


