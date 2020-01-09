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


def morphology(path1, paths, norm, normpath, timei, showtransient, showhline, lbl):
    '''    
    path1 = '/work/nilscp/iSALE/isaleruns/testruns/ejecta/L20km_part1_HYDRO_CPPR1000/results/L20km_HYDRO_part1g/'
    paths = '/work/nilscp/iSALE/isaleruns/data/ejecta/L20km_HYDRO_part1g/'
    normpath = '/work/nilscp/iSALE/isaleruns/data/ejecta/L20km_HYDRO_part1g/'
    norm = 0.5 # 0.5 by the radius of the projectile?
    showtransient = False
    showhline = False
    timei = 'pen' #'transient' or 'norm' or 'pen' (for penetration time)

    lbl = r"$\mathit{U} = 12.7 km/s, \mathit{L} = 20 km$"

    morphology(path1, paths, norm, normpath, timei, showtransient, showhline, lbl)



    #loop
    L = [100,200,300,400,500,600,700,800,900,1000]
    name = 'C00P10F06LONG_U04_L'

    for l in L:
        tmp_name = name + str(l)
         path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/velocity/results/' + tmp_name
         paths = '/run/media/nilscp/Zell/velocity/data/'  + tmp_name
         lbl = r'$\mathit{U} = 400 m/s, \mathit{L} = ' + str(l) + ' m$"
         evoplot(path1,paths, path1tr, pathstr, norm, extentx, extenty, lbl)


    ############################3

    What if it is nan?

    What if we want to plot the transient crater in addition

    Maybe I could also have a function without having the need to have extentx and extenty
    taking 20% of the final drim and so much of the maximum depth.... or a fraction of the final drim
    we known that d ~ 0.2D

    norm = 0.5 can only be used for showtransient == False
    '''

    plt.ioff()  # figures does not pop up

    os.chdir(path1)
    mod1 = psp.opendatfile('jdata.dat')

    modelname = path1.split('/')[-2]

    if not os.path.exists(paths + 'evolution/plots/morphology'):
        os.makedirs(paths + 'evolution/plots/morphology')

    if (norm == 0):
        if not os.path.exists(paths + 'evolution/plots/morphology/absolute/'):
            os.makedirs(paths + 'evolution/plots/morphology/absolute')
    elif (norm == 0.5):
        if not os.path.exists(paths + 'evolution/plots/morphology/norm_proj/'):
            os.makedirs(paths + 'evolution/plots/morphology/norm_proj')
    elif (norm == 1):
        if not os.path.exists(paths + 'evolution/plots/morphology/norm_its/'):
            os.makedirs(paths + 'evolution/plots/morphology/norm_its')
    else:
        if not os.path.exists(paths + 'evolution/plots/morphology/norm/'):
            os.makedirs(paths + 'evolution/plots/morphology/norm')

    if norm == 1:
        os.chdir(paths + 'transient/')

        # load data for the transient crater 1
        t_trg, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        os.chdir(paths + 'final/')
        t_f, __, __, __, __, drim_f, __ = np.loadtxt(
            modelname + '_final.txt', delimiter='\t', comments='#')
        rinner1 = drim_f/2.

    elif norm == 2:

        os.chdir(paths + 'transient/')

        # load data for the transient crater
        t_trg, __, __, __, __, __, __, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        modelnametr = normpath.split('/')[-2]
        os.chdir(normpath + 'transient/')

        # load data for the transient crater modtr
        t_tr1, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelnametr + '_tr.txt', delimiter='\t', comments='#')

        os.chdir(normpath + 'final/')
        t_fnorm, __, __, __, __, drim_fnorm, __ = np.loadtxt(
            modelnametr + '_final.txt', delimiter='\t', comments='#')
        rinner1 = drim_fnorm/2.

    else:
        os.chdir(paths + 'transient/')

        # load data for the transient crater 1
        t_trg, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        os.chdir(paths + 'final/')
        t_f, __, __, __, __, drim_f, __ = np.loadtxt(
            modelname + '_final.txt', delimiter='\t', comments='#')

    # if we want to show the transient crater
    if showtransient:
        os.chdir(paths + 'transient/')
        dataXY = np.loadtxt(
            modelname + '_XYtransientprofile.txt', delimiter='\t', comments='#')
        X = dataXY[:, 0]
        Y = dataXY[:, 1]

    os.chdir(path1)
    for i in range(mod1.nsteps):
        step = mod1.readStep('Den', i)

        t1 = step.time
        # plotting of the data
        fig = plt.figure(figsize=(6, 6))
        ax1 = fig.add_subplot(111)

        if norm >= 1:
            ax1.pcolormesh(mod1.x/rinner1, mod1.y/rinner1, step.mat,
                           cmap='BrBG_r', vmin=-5., vmax=5., zorder=2)
            ax1.contour(mod1.xc/rinner1, mod1.yc/rinner1,
                        step.cmc[0], 1, colors='b', linewidths=4, zorder=4)
            if t1 >= t_trg:
                if showtransient:
                    # I need to show it only for a range
                    ax1.plot(X/rinner1, Y/rinner1, "ro")
                else:
                    None
            #ij1 = np.where((data_tr1[:,1]<=500.) & (data_tr1[:,0]<=Dr_tr1/2.))
            #ij2 = np.where((data_tr2[:,1]<=500.) & (data_tr2[:,0]<=Dr_tr2/2.))

            if showhline:
                ax1.hlines(0, 0, 1.2, 'r', linewidths=4, zorder=1)
            else:
                None

            ax1.set_xlim([0, 1.2])  # -2000,0 ## -1500,1500

            for ax in [ax1]:
                ax.minorticks_off()
                ax.tick_params('both', labelsize=16, length=10,
                               width=2., which='major')
                ax.set_ylim([-0.8, 0.8])

            for u in range(mod1.tracer_numu):
                tru = mod1.tru[u]
                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.xlines), 20):
                    ax1.plot(step.xmark[tru.xlines[l]]/rinner1,
                             step.ymark[tru.xlines[l]]/rinner1,
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.ylines), 20):
                    ax1.plot(step.xmark[tru.ylines[l]]/rinner1,
                             step.ymark[tru.ylines[l]]/rinner1,
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

        elif norm == 0.5:
            ax1.pcolormesh(mod1.x/(mod1.cppr[0]*mod1.dx), mod1.y/(
                mod1.cppr[0]*mod1.dx), step.mat, cmap='BrBG_r', vmin=-5., vmax=5., zorder=2)
            ax1.contour(mod1.xc/(mod1.cppr[0]*mod1.dx), mod1.yc/(
                mod1.cppr[0]*mod1.dx), step.cmc[0], 1, colors='b', linewidths=4, zorder=4)
            if t1 >= t_trg:
                if showtransient:
                    # I need to show it only for a range
                    ax1.plot(X/(mod1.cppr[0]*mod1.dx),
                             Y/(mod1.cppr[0]*mod1.dx), "ro")
                else:
                    None

            if showhline:
                ax1.hlines(
                    0, 0, mod1.xhires[1]/(mod1.cppr[0]*mod1.dx), 'r', linewidths=4, zorder=1)
            else:
                None
            # -2000,0 ## -1500,1500
            ax1.set_xlim([0, mod1.xhires[1]/(mod1.cppr[0]*mod1.dx)])

            for ax in [ax1]:
                ax.set_ylim([mod1.yhires[0]/(mod1.cppr[0]*mod1.dy),
                             mod1.yhires[1]/(mod1.cppr[0]*mod1.dy)])
                ax.minorticks_off()
                ax.tick_params('both', labelsize=16, length=10,
                               width=2., which='major')

            for u in range(mod1.tracer_numu):
                tru = mod1.tru[u]
                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.xlines), 20):
                    ax1.plot(step.xmark[tru.xlines[l]]/(mod1.cppr[0]*mod1.dx),
                             step.ymark[tru.xlines[l]]/(mod1.cppr[0]*mod1.dx),
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.ylines), 20):
                    ax1.plot(step.xmark[tru.ylines[l]]/(mod1.cppr[0]*mod1.dx),
                             step.ymark[tru.ylines[l]]/(mod1.cppr[0]*mod1.dx),
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

        else:
            ax1.pcolormesh(mod1.x, mod1.y, step.mat,
                           cmap='BrBG_r', vmin=-5., vmax=5., zorder=2)
            ax1.contour(mod1.xc, mod1.yc,
                        step.cmc[0], 1, colors='b', linewidths=4, zorder=4)
            if t1 >= t_trg:
                if showtransient:
                    ax1.plot(X, Y, "ro")  # I need to show it only for a range
                else:
                    None

            if showhline:
                ax1.hlines(0, 0, (drim_f/2.)*1.2, 'r', linewidths=4, zorder=1)
            else:
                None
            ax1.set_xlim([0, (drim_f/2.)*1.2])  # -2000,0 ## -1500,1500

            for ax in [ax1]:
                ax.set_ylim([-(drim_f/2.)*0.8, (drim_f/2.)*0.8])
                ax.minorticks_off()
                ax.tick_params('both', labelsize=16, length=10,
                               width=2., which='major')

            for u in range(mod1.tracer_numu):
                tru = mod1.tru[u]
                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.xlines), 20):
                    ax1.plot(step.xmark[tru.xlines[l]],
                             step.ymark[tru.xlines[l]],
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

                # Plot the tracers in horizontal lines, every 20 lines
                for l in np.arange(0, len(tru.ylines), 20):
                    ax1.plot(step.xmark[tru.ylines[l]],
                             step.ymark[tru.ylines[l]],
                             c='k', marker='.', linestyle='None', markersize=1.0, zorder=3)

        if norm == 0:
            ax1.set_xlabel("x (m)", fontsize=20)
            ax1.set_ylabel("y (m)", fontsize=20)
            fig.tight_layout()

        elif norm == 0.5:
            ax1.set_xlabel(r"$x / a$", fontsize=20)
            ax1.set_ylabel(r"$y / a$", fontsize=20)
            fig.tight_layout()

        else:
            ax1.set_xlabel(r"$x / R_{r}$", fontsize=20)
            ax1.set_ylabel(r"$y / R_{r}$", fontsize=20)
            fig.tight_layout()

        if timei == 'normal':
            ax1.set_title(r"$\mathit{t} = " + str(np.around(t1, decimals=2)
                                                  ) + " " + " s$", position=(0.5, 0.9), fontsize=18)
        elif timei == 'pen':
            ax1.set_title(r"$\mathit{t} = " + str(np.around(t1, decimals=2)) + " " + " s, $" + " " + "$t / t_s = $" + str(
                np.around(t1/((mod1.cppr[0]*mod1.dx*2.) / (mod1.objvel)), decimals=2)), position=(0.5, 0.9), fontsize=18)
        else:
            ax1.set_title(r"$\mathit{t} = " + str(np.around(t1, decimals=2)) + " " + " s, $" + " " +
                          "$\zeta$ = $" + str(np.around(t1/t_trg, decimals=2)), position=(0.5, 0.9), fontsize=18)

        st = fig.suptitle(lbl, fontsize=22)
        st.set_y(0.97)
        st.set_x(0.55)
        fig.subplots_adjust(top=0.9)

        if norm == 0:
            figtitle = modelname + '_' + str(int(i)).zfill(3) + ".png"
            fig.savefig(
                paths + 'evolution/plots/morphology/absolute/' + figtitle, dpi=300)
        elif norm == 0.5:
            figtitle = modelname + '_' + str(int(i)).zfill(3) + ".png"
            fig.savefig(
                paths + 'evolution/plots/morphology/norm_proj/' + figtitle, dpi=300)
        elif norm == 1:
            figtitle = 'itself_norm_' + modelname + \
                '_' + str(int(i)).zfill(3) + ".png"
            fig.savefig(
                paths + 'evolution/plots/morphology/norm_its/' + figtitle, dpi=300)
        else:
            figtitle = 'norm_' + modelname + \
                '_' + str(int(i)).zfill(3) + ".png"
            fig.savefig(
                paths + 'evolution/plots/morphology/norm/' + figtitle, dpi=300)
        plt.close()

    mod1.closeFile()


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




def field_definition(model, zoom_id, paths, param = 'all'):
    '''
    example:
    param = 'Por'
    paths = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILWPO_L100/'
    zoom_id = ['mid', 'hires']

    fld_param, fld_name, fld_cmap, fld_unit, fld_factor, npath = field_definition(model, zoom_id, paths, param = 'all')
    
    #Updated field_definition the 24th June 2019
    
    '''

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
    
    
    # In case for all the parameters 
    if param == 'all':
        
        # get all the field available for the model
        iSALEparam0 = model.fieldlist
        iSALEparam = []
    
        for param in iSALEparam0:
            if ((param[0].startswith('V_')) | (param[0].startswith('Tr'))):
                None
            else:
                iSALEparam.append(param[0])
                
        # get the lists
        param1 = []
        name = []
        cmap = []
        unit = []
        factor = []
        npathl = []
        
        for parami in iSALEparam:
            
            idx = np.where(fld == parami)[0][0]
        
            param1.append(fld_param[idx])
            name.append(fld_name[idx])
            cmap.append(fld_cmap[idx])
            unit.append(fld_unit[idx])
            factor.append(fld_factor[idx])
            
            
            for z in zoom_id:
                npath = paths + 'evolution/plots/' + fld_folder[idx] + '/' + z + '/'
                npathl.append(npath)
                
                if not os.path.exists(npath):
                    os.makedirs(npath)
                    
    
    else:
        idx = np.where(fld == param)[0][0]
    
        param1 = fld_param[idx]
        name = fld_name[idx]
        cmap = fld_cmap[idx]
        unit = fld_unit[idx]
        factor = fld_factor[idx]
        npathl = []
        
        for z in zoom_id:
            npath = paths + 'evolution/plots/' + fld_folder[idx] + '/' + z + '/'
            npathl.append(npath)
            
            if not os.path.exists(npath):
                os.makedirs(npath)

    return param1, name, cmap, unit, factor, npathl


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

def time_lbl(path_tosave, modelname, time, norm):
    
    '''
    Need to modify this part if you don't want to plot the time and norm time
    
    what if we want to have the time related to the penetration time?
    '''
    
    tr_time = getTransientvalues(path_tosave, modelname)[0]
    tr_idx = np.int(getTransientvalues(path_tosave, modelname)[-1])
    
    lbl_time = (r"$\mathit{t} = " + "{:.2f}".format(time) + " " + " s, $" + " " +
                              "$\zeta$ = " + "{:.2f}".format(time/tr_time))
    
    
    return (tr_idx, lbl_time)
    
    
    
'''
***********************************************************************
'''

def figurename(jdata, path_tosave, step, param, norm, zoom_level):
    
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
    
    # modelname
    modelname = jdata.split('/')[-2]
    
    # parameter name
    P = param.upper()
    
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
        ax1.set_ylabel(y_lbl, fontsize=16, labelpad=-26) #r"$y / R_{r}$"
        ax2.set_ylabel(y_lbl, rotation=270, fontsize=16, labelpad=-10) #r"$y / R_{r}$"
        ax2.set_xlabel(x_lbl, fontsize=16)

    else:
        None
                
    # set labelsize to 14
    for ax in axes:
        ax.minorticks_off()
        ax.tick_params('both', labelsize=14, length=5)
        
        # set labels
        ax.set_xlabel(x_lbl, fontsize=16)

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
        ax2.set_title('blablabla', position=(0.0, 0.9), fontsize=18)
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

def p(path_tosave, modelname, model, step, iparam, fig, ax, scaling_factor, 
      vmin, vmax, colormap, colormap_lbl, param_unit_factor, sign, 
      stepiSALE, stepTransient, dual, 
      showSurface, showTransient, showTracers):
    
    # plot data
    pco = ax.pcolormesh(sign * model.x/scaling_factor, 
                         model.y/scaling_factor, 
                         step.data[iparam]/param_unit_factor,
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
    
        

    

def field(jdata, path_tosave, steps, params, norm, zoom_level, vmin, vmax, title_lbl,
           x_lbl, y_lbl,
           colormap, colormap_lbl, param_unit_factor = 1.0,
           showhline = False, showTransient = False, showSurface = False, 
           showTracers = False, dual = False,
           norm_tend = -5.0):
    
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
    params = ['Den', ]
    norm = 1
    zoom_level = 'medium'
    vmin = 0
    vmax = 3000
    title_lbl = 'Shroedinger (T = 40 km), density'
    x_lbl = 'x'
    y_lbl = 'y'
    colormap = 'viridis'
    showhline = True
    showTransient = True
    showSurface = True
    showTracers = 'int' (number of tracers space between) or False
    colormap_lbl = 'kg/m^{3}'
    crater(jdata, path_tosave, steps, params, norm, zoom_level, vmin, vmax, title_lbl,
           x_lbl, y_lbl, colormap, colormap_lbl, 1.0, 1.0, showhline, showTransient, showSurface, showTracers)
    
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
    
    # should load the data to 
    
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
            
            # dual 
            if dual:
                    
                # fix plot and other related plot layout
                fig, ax = plt.subplots(figsize=(12, 6), nrows = 1, ncols = 2)
                fig.subplots_adjust(wspace=0.0)
                
                # plots based on preferences
                cpo1 = p() # step0
                cpo2 = p() # step1
                
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
                st = fig.suptitle(title_lbl, fontsize=20) # main title
                st.set_y(0.97)
                st.set_x(0.5)
                fig.subplots_adjust(top=0.9)
                
                # saving figure
                fig.savefig(figname, dpi=300)   
                
            else:
                
                fig, ax = plt.subplots(figsize=(6, 6), nrows = 1, ncols = 1)
                
                for iparam, param in enumerate(params):
                    
                    # get figurename based on the different selected parameters
                    figname = figurename(jdata, path_tosave, i, param, norm, zoom_level)
                    
                    # set time label
                    idx_tr, t_lbl = time_lbl(path_tosave, modelname, t1, norm)
                    
                    #                     
                    cpo = p()

                    
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
                    st = fig.suptitle(title_lbl, fontsize=20) # main title
                    st.set_y(0.97)
                    st.set_x(0.5)
                    fig.subplots_adjust(top=0.9)
                    
                    # saving figure
                    fig.savefig(figname, dpi=300)  
                    
                    
                    
                    
                    # distance between the two plots
                    # 
                else:
                    fig = plt.figure(figsize=(6, 6))
                    ax1 = fig.add_subplot(111)
                
                # plot data
                cax = ax1.pcolormesh(x_direction * model.x/scaling_factor, 
                                     model.y/scaling_factor, 
                                     step.data[iparam]/param_unit_factor,
                                         vmin=vmin,
                                         vmax=vmax,
                                         cmap=colormap, zorder=2)
                
                # showing surface
                if showSurface:
                    ax1.contour(x_direction * model.xc/scaling_factor, 
                                model.yc/scaling_factor,
                                step.cmc[0], 1, colors='k', linewidths=2, 
                                zorder=4)
                
                # showing transient only if below transient
                if showTransient:
                    if i < idx_tr:
                        None
                    else:
                        # this values need to be scaled
                        xtr_surf, ytr_surf = getTransientcrossSection(path_tosave, modelname)
                        
                        ax1.plot(x_direction * xtr_surf/scaling_factor, ytr_surf/scaling_factor,
                            'r', linewidth=2, zorder=5)
                
                # showing tracers        
                if showTracers:                    
                    for u in range(model.tracer_numu):
                        tru = model.tru[u]
                        # Plot the tracers in horizontal lines, every 20 lines
                        for l in np.arange(0, len(tru.xlines), showTracers):
                            ax1.plot(x_direction *step.xmark[tru.xlines[l]]/scaling_factor,
                                     step.ymark[tru.xlines[l]]/scaling_factor,
                                     c='k', marker='.', linestyle='None', 
                                     markersize=1.0, zorder=3)
        
                        # Plot the tracers as vertical lines, every 20 lines
                        for l in np.arange(0, len(tru.ylines), showTracers):
                            ax1.plot(x_direction *step.xmark[tru.ylines[l]]/scaling_factor,
                                     step.ymark[tru.ylines[l]]/scaling_factor,
                                     c='k', marker='.', linestyle='None', 
                                     markersize=1.0, zorder=3)

                    
                for ax in [ax1]:
                    ax.minorticks_off()
                    ax.tick_params('both', labelsize=16, length=10,
                                   width=2., which='major')
                    
                cbar = fig.colorbar(cax)
                cbar.ax.tick_params(labelsize=15)
                cbar.set_label('$' + colormap_lbl + '$', labelpad=-
                               42.5, y=1.07, fontsize=15, rotation=0)
                
                # set labels
                ax1.set_xlabel(x_lbl, fontsize=20)
                ax1.set_ylabel(y_lbl, fontsize=20) #r"$y / R_{r}$"
                
                # set title            
                ax1.set_title(t_lbl, position=(0.5, 0.9), fontsize=18)
                
                # related to the zoom
                ax1.set_xlim(extentx)
                for ax in [ax1]:
                    ax.set_ylim(extenty)
                fig.tight_layout()
                
                st = fig.suptitle(title_lbl, fontsize=20) # main title
                st.set_y(0.97)
                st.set_x(0.5)
                fig.subplots_adjust(top=0.9)
                
                fig.savefig(figname, dpi=300)                                       
                plt.close()

# TODO: delete

def field(path1, paths, norm, normpath, zoom_id, lbl, vmiin, vmaax, manualx, manualy, param = 'all'):
    '''
    param = 'Den'  
    path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/CDILWPO/CDILWPO_L100_dev_larger_cells/'
    paths = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILWPO_L100_dev_larger_cells/'
    normpath = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILWPO_L100/'
    norm = 0
    zoom_id = ['close', 'mid', 'hires', 'all', 'manual'] # 'mid'; 'hires'; 'all' #depending on the zoom
    param = 'all'
    vmiin = [] # or vmiin = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e-5, 0.0] # need to take into account the field factor
    vmaax = [] # or vmaax = [1.0, 1.0, 2000.0, 400.0, 1.0, 1.0, 1.0, 1.0, 1.25,1.0, 1.0]

    lbl = r"$\mathit{L} = 100 m"

    field(path1, paths, norm, normpath, zoom_id, lbl, vmiin, vmaax, param)
    
    I should replace vmiin and vmaax by percentile values that are rounded to the hundreds


    Field list:
    [['Cm1', 'Concentration, mat. 1'],
     ['Cm2', 'Concentration, mat. 2'],
     ['Den', 'Density'],
     ['Tmp', 'Temperature'],
     ['Pre', 'Pressure'],
     ['Sie', 'Specific internal energy'],
     ['Yld', 'Yield strength'],
     ['Dam', 'Damage'],
     ['Alp', 'Distension'],
     ['TPS', 'Total plastic strain (shear)'],
     ['YAc', 'Acoustic fluidisation strength'],
     ['V_x', 'Velocity (x)'],
     ['V_y', 'Velocity (y)'],
     ['Trx', 'Tracer x position'],
     ['Try', 'Tracer y position'],
     ['TrP', 'Tracer peak pressure'],
     ['TrT', 'Tracer peak temperature']]
    '''

    plt.ioff()  # figures does not pop up

    os.chdir(path1)
    mod1 = psp.opendatfile('jdata.dat')

    # modelname
    modelname = path1.split('/')[-2]

    # loading ot the data
    ###########################################################################
    if norm == 1:
        os.chdir(paths + 'transient/')

        # load data for the transient crater 1
        t_trg, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        os.chdir(paths + 'final/')
        t_f, __, __, __, __, drim_f, __ = np.loadtxt(
            modelname + '_final040.txt', delimiter='\t', comments='#') # put to 040 here
        
        rinner1 = drim_f/2.

    elif norm == 2:

        os.chdir(paths + 'transient/')

        # load data for the transient crater
        t_trg, __, __, __, __, __, __, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        modelnametr = normpath.split('/')[-2]
        os.chdir(normpath + 'transient/')

        # load data for the transient crater modtr
        t_tr1, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelnametr + '_tr.txt', delimiter='\t', comments='#')
        
        # is it normal that it only load the data for the norm path?
        # maybe I guess so
        os.chdir(normpath + 'final/')
        t_fnorm, __, __, __, __, drim_fnorm, __ = np.loadtxt(
            modelnametr + '_final040.txt', delimiter='\t', comments='#') # put to 040 here
        
        rinner1 = drim_fnorm/2.

    else:
        os.chdir(paths + 'transient/')

        # load data for the transient crater 1
        t_trg, da_tr1, Da_tr1, V_tr1, h_tr1, Dr_tr1, Vr_tr1, __ = np.loadtxt(
            modelname + '_tr.txt', delimiter='\t', comments='#')

        os.chdir(paths + 'final/')
        t_f, __, __, __, __, drim_f, __ = np.loadtxt(
            modelname + '_final040.txt', delimiter='\t', comments='#') # put to 040 here
            
        rinner1 = drim_f/2.

    ###########################################################################
    
    try:
        # choosing the degree of the zooming on the figure
        # get the extent
        extentx, extenty, zoom_id_up = zoom(zoom_id, norm, mod1, rinner1, manualx, manualy)
        
        # field list (should I make a function out of it?)
        fld_param, fld_name, fld_cmap, fld_unit, fld_factor, npath = field_definition(mod1, zoom_id_up, paths, param)

        # update label
        nlbl = []
        
        # small change (for case with only one selected parameter)
        if param != 'all':
            nlbl.append(lbl + ", " + fld_name) #+ "$" # the $ introduces some errors ...
            fld_param = [fld_param]
            fld_factor = [fld_factor]
            fld_name = [fld_name]
            fld_cmap = [fld_cmap]
            fld_unit = [fld_unit]
        else:        
            for fname in fld_name:
                nlbl.append(lbl + ", " + fname) #+ "$" # the $ introduces some errors ...
        
        # print fld_param
        print (fld_param)


        

        
        #vmin and vmax (update here)
        if len(vmiin) == 0:
            vmiin, vmaax = vmin_vmax(path1, mod1, fld_param, fld_factor, param)
        else:
            None
    
        # looping of
        os.chdir(path1)
        for i in range(mod1.nsteps):
            step = mod1.readStep(fld_param, i)
    
            t1 = step.time
            # plotting of the data
    
            # changing the vmin and vmax every step will change the value everytime
            # maybe I should define manually
            
            if norm >= 1:
                
                '''
                After the change including all parameters this statement will never
                be fullfill expect when param = 'Por'
                '''
                
                # don't like the 'Por' thing
                #if param == 'Por':
                #    cax = ax1.pcolormesh(mod1.x/rinner1, mod1.y/rinner1, (1. - (1./step.data[0]))*100.,
                #                         vmin=vmiin,
                #                         vmax=vmaax,
                #                         cmap=fld_cmap, zorder=2)
                
                count = 0
                
                for ifld, fld in enumerate(fld_param):
                    
                    for iz, z in enumerate(zoom_id_up):
                        
                        fig = plt.figure(figsize=(6, 6))
                        ax1 = fig.add_subplot(111)
                        
                        # plot data
                        cax = ax1.pcolormesh(mod1.x/rinner1, mod1.y/rinner1, step.data[ifld]/fld_factor[ifld],
                                                 vmin=vmiin[ifld],
                                                 vmax=vmaax[ifld],
                                                 cmap=fld_cmap[ifld], zorder=2)
        
                        ax1.contour(mod1.xc/rinner1, mod1.yc/rinner1,
                                step.cmc[0], 1, colors='k', linewidths=2, zorder=4)
                    
                    
        
                        for ax in [ax1]:
                            ax.minorticks_off()
                            ax.tick_params('both', labelsize=16, length=10,
                                           width=2., which='major')
                            
                        cbar = fig.colorbar(cax)
                        cbar.ax.tick_params(labelsize=15)
                        cbar.set_label('$' + fld_unit[ifld] + '$', labelpad=-
                                       42.5, y=1.07, fontsize=15, rotation=0)
                        
                        # set labels
                        ax1.set_xlabel(r"$x / R_{r}$", fontsize=20)
                        ax1.set_ylabel(r"$y / R_{r}$", fontsize=20)
                        
                        # set title            
                        ax1.set_title(r"$\mathit{t} = " + "{:.2f}".format(t1) + " " + " s, $" + " " +
                                      "$\zeta$ = " + "{:.2f}".format(t1/t_trg), position=(0.5, 0.9), fontsize=18)
                        
                        # related to the zoom
                        ax1.set_xlim(extentx[iz])
                        for ax in [ax1]:
                            ax.set_ylim(extenty[iz])
                        fig.tight_layout()
                        
                        st = fig.suptitle(nlbl[ifld], fontsize=20)
                        st.set_y(0.97)
                        st.set_x(0.5)
                        fig.subplots_adjust(top=0.9)
                            
                                                    
                        if norm == 1:
                            figtitle = 'itself_norm_' + modelname + \
                                '_' + str(int(i)).zfill(3) + ".png"
                            fig.savefig(npath[count] + figtitle, dpi=300)
                            
                        elif norm == 2:
                            figtitle = 'norm_' + modelname + \
                                '_' + str(int(i)).zfill(3) + ".png"
                            fig.savefig(npath[count] + figtitle, dpi=300)
                            
                        plt.close()
                        
                        count = count + 1
    
    
            else:
                count = 0
                
                for ifld, fld in enumerate(fld_param):
                    
                    for iz, z in enumerate(zoom_id_up):
                        
                        fig = plt.figure(figsize=(6, 6))
                        ax1 = fig.add_subplot(111)
                        
                        # plot data
                        cax = ax1.pcolormesh(mod1.x, mod1.y, step.data[ifld]/fld_factor[ifld],
                                                 vmin=vmiin[ifld],
                                                 vmax=vmaax[ifld],
                                                 cmap=fld_cmap[ifld], zorder=2)
        
                        ax1.contour(mod1.xc, mod1.yc,
                                step.cmc[0], 1, colors='k', linewidths=2, zorder=4)
                    
                    
        
                        for ax in [ax1]:
                            ax.minorticks_off()
                            ax.tick_params('both', labelsize=16, length=10,
                                           width=2., which='major')
                            
                        cbar = fig.colorbar(cax)
                        cbar.ax.tick_params(labelsize=15)
                        cbar.set_label('$' + fld_unit[ifld] + '$', labelpad=-
                                       42.5, y=1.07, fontsize=15, rotation=0)
                        
                        # set labels
                        ax1.set_xlabel(r"$x$", fontsize=20)
                        ax1.set_ylabel(r"$y$", fontsize=20)
                        
                        # set title            
                        ax1.set_title(r"$\mathit{t} = " + "{:.2f}".format(t1) + " " + " s$", position=(0.5, 0.9), fontsize=18)
                        
                        # related to the zoom
                        ax1.set_xlim(extentx[iz])
                        for ax in [ax1]:
                            ax.set_ylim(extenty[iz])
                        fig.tight_layout()
                        
                        st = fig.suptitle(nlbl[ifld], fontsize=20)
                        st.set_y(0.97)
                        st.set_x(0.5)
                        fig.subplots_adjust(top=0.9)
                            
                        figtitle = modelname + '_' + str(int(i)).zfill(3) + ".png"
                        fig.savefig(npath[count] + figtitle, dpi=300)
                        
                            
                        plt.close()
                        
                        count = count + 1
    except:
        None

    mod1.closeFile()


'''
***********************************************************************
'''

# TO DO: delete

def mainevo(folders, path1, paths, norm, normpath, zoom_id, param):
    
    '''
    path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/'
    paths = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/'
    normpath = ''
    norm = 0
    zoom_id = ['close', 'all', 'mid', 'hires'] # 'mid'; 'hires'; 'all' #depending on the zoom
    param = 'all'
    lbl = r"$\mathit{L} = 100 m"
    L = [100]
    
    # maybe detect automatically from what stands behind L 
    
    Need to be tested
    
    
    '''
    os.chdir(path1)
    
    folders = glob.glob('*')
    
    for f in folders:
        os.chdir(path1 + '/' + f)
        modelnames = glob.glob('*')
        print (modelnames)
     
        for im, modelname in enumerate(modelnames):
            
            pathsa = paths + modelname + '/'
            path1a = path1 + f + '/' + modelname + '/'
            
            L = modelname.split('_L')[1]
            Ln = L.split('_')[0]
            
            #normpath should be the same
            lbla = r"$\mathit{L} = " + Ln + " m"
                        
            field(path1a, pathsa, norm, normpath, zoom_id, lbla, param = 'all')
        

# TO DO: delete        

def main(pathdata, folders, pathplots):
    '''
    description:
    routine that plots both the transient, excavated and final cavities.
    This routine is useful to see whether the different crater dimensions by
    the scripts in crater.py and ejecta.py managed to perform well (i.e., this is
    visual quality assesment of the results)

    inputs:
    pathdata: path where the text files have been stored (text files that are obtained
    from running the main.py script)

    folders: all modelcases that need to be plotted

    pathplots: saving plot directory
    '''

    for modelname in folders:
        transient(pathdata, modelname, pathplots)
        final(pathdata, modelname, pathplots)
        excavated(pathdata, modelname, pathplots)

    print ("ROCK'N ROLL")


'''
***********************************************************************
'''

# TO DO: improve

def transient(pathdata, modelname, pathplots):
    '''
    description:
    routine that plots transient crater dimensions

    inputs:
    pathdata: path where the text files have been stored (text files that are obtained
    from running the main.py script)

    folders: all modelcases that need to be plotted

    pathplots: saving plot directory

    example:
    pathdata = '/media/nilscp/Zell/Collapse/data/'
    modelname = 'C00P20F08_L250'
    pathplots = '/media/nilscp/Zell/Collapse/rplots/'
    transient(pathdata,modelname,pathplots)
    '''

    # should be outside of the loop
    if not os.path.exists(pathplots):
        os.makedirs(pathplots)

    os.chdir(pathdata + modelname + '/transient/')

    t, da, Da, V, h, Dr, Vr, __ = np.loadtxt(
        modelname + '_tr.txt', delimiter='\t', comments='#')
    dataXY = np.loadtxt(modelname + '_XYtransientprofile.txt',
                        delimiter='\t', comments='#')
    X = dataXY[:, 0]
    Y = dataXY[:, 1]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.plot(X, Y, "ro", zorder=2)
    ax.hlines(0, np.min(X), np.max(X), 'k', linewidth=3)
    ax.plot(Da/2., 0, 'yo', ms=10, label='Ra: apparent radius')
    ax.plot(Dr/2., h, 'co', ms=10, label='Drim and hrim')
    ax.hlines(-da, np.min(X), np.max(X), 'b',
              linewidth=3, label='da: apparent depth')
    ax.legend(loc='center right')
    ax.set_xlabel('X (m)', fontsize=18)
    ax.set_ylabel('Y (m)', fontsize=18)
    ax.tick_params('both', labelsize=16, length=10, width=1., which='major')
    fig.tight_layout()
    st = fig.suptitle(modelname + ' transient profile', fontsize=18)
    st.set_y(0.95)
    fig.subplots_adjust(top=0.9)
    figtitle = modelname + ' transient_profile.png'
    fig.savefig(pathplots+figtitle, dpi=300)
    plt.close()


'''
***********************************************************************
'''

# TO DO: improve

def final(pathdata, modelname, pathplots):
    '''
    description:
    routine that plots final crater dimensions

    inputs:
    pathdata: path where the text files have been stored (text files that are obtained
    from running the main.py script)

    folders: all modelcases that need to be plotted

    pathplots: saving plot directory

    example:
    pathdata = '/media/nilscp/Zell/Collapse/data/'
    modelname = 'C00P20F08_L250'
    pathplots = '/media/nilscp/Zell/Collapse/rplots/'
    final(pathdata,modelname,pathplots)
    '''

    # should be outside of the loop
    if not os.path.exists(pathplots):
        os.makedirs(pathplots)

    os.chdir(pathdata + modelname + '/final/')

    t, da, Da, V, h, Dr, Vr = np.loadtxt(
        modelname + '_final.txt', delimiter='\t', comments='#')
    dataXY = np.loadtxt(modelname + '_XYfinalprofile.txt',
                        delimiter='\t', comments='#')
    X = dataXY[:, 0]
    Y = dataXY[:, 1]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.plot(X, Y, "ro", zorder=2)
    ax.hlines(0, np.min(X), np.max(X), 'k', linewidth=3)
    ax.plot(Da/2., 0, 'yo', ms=10, label='Ra: apparent radius')
    ax.plot(Dr/2., h, 'co', ms=10, label='Drim and hrim')
    ax.hlines(-da, np.min(X), np.max(X), 'b',
              linewidth=3, label='da: apparent depth')
    ax.legend(loc='center right')
    ax.set_xlabel('X (m)', fontsize=18)
    ax.set_ylabel('Y (m)', fontsize=18)
    ax.tick_params('both', labelsize=16, length=10, width=1., which='major')
    fig.tight_layout()
    st = fig.suptitle(modelname + ' final profile', fontsize=18)
    st.set_y(0.95)
    fig.subplots_adjust(top=0.9)
    figtitle = modelname + ' final_profile.png'
    fig.savefig(pathplots+figtitle, dpi=300)
    plt.close()


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


