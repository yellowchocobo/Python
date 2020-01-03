# -*- coding: utf-8 -*-

# Import basic Python modules
import numpy as np
import os
import sys
from scipy.optimize import curve_fit
from datetime import datetime

# path pySALEPlot
path_pySALEPlot = '/work/nilscp/iSALE/Dellen/lib'

# append path to pySALEPlot to system
sys.path.append(path_pySALEPlot)

# Import pySALEPlot
import pySALEPlot as psp

'''
***********************************************************************
'''


def loadej(path_tosave, modelname):
    '''
    description:
    load txt file containing ejected parameters information based on path_tosave
    and the modelname
    
    params:
    : path_tosave (str): main plotting folder
    : modelname (str) : modelname
    
    returns:
    : tracer_idx (int): index of ejected lagrangian tracers
    : ts (float): time of ejection (in s)
    : v (float): ejection velocity (in m/s)
    : angle (float): ejection angle (in degrees)
    : xpos (float): ejection position (in meters)
    : tair (float): time in the air (s)
    : dland (float): landing position of ejected tracers (in meters)
    : n (int): number of detected ejected tracers per timestep
    '''

    data = np.loadtxt(path_tosave + modelname + '/ejtracers/' + 
                      modelname + '_ejtracers.txt',
                      delimiter=';', comments='#')

    # extract the columns from the data
    tracer_ix = data[:, 0]
    ts = data[:, 1]
    v = data[:, 2]
    angle = data[:, 3]
    xpos = data[:, 4]
    tair = data[:, 5]
    dland = data[:, 6]

    n = np.loadtxt(path_tosave + modelname + '/ejtracers/' + 
                   modelname + '_ntracers.txt', delimiter='\t', comments='#')

    return (tracer_ix, ts, v, angle, xpos, tair, dland, n)


'''
***********************************************************************
'''


def posTracer(model, method, threshold):
    '''
    description:
    calculate the positions of ejected materials based on the threshold
    (defined in terms of height)

    params:
    model (iSALE model): model output from psp
    method: either use the 2-dot (method = 2) or 3-dot (method = 3) methods.
    Give slightly different results if two or three saved timesteps are used to 
    determine the impact velocity, and ejection angle.
    threshold (float): in meters above the pre-impact surface

    returns:
    xx (float) : positions of detected tracers (x) - (either 2 or 3 positions)
    yy (float): positions of detected tracers (y) - (either 2 or 3 positions)
    tracer_idx : tracer_idx (int) : index of ejected materials (lagrangian tracer indexs)
    timestep : saved timestep when tracers are first detected
    dt : time increment (saved timestep)
    n : number detected tracers each saved timestep
    '''

    # start time to calculate the time it takes to run the script
    startTime = datetime.now()

    # create empty arrays
    x0, x1, x2, y0, y1, y2, tracer_idx, timestep, n = (np.array([]),)*9

    # loop through saved timesteps
    for i in range(model.nsteps-2):

        # if equal to 0, we store DTSAVE (i.e., the time between saved timesteps)
        if i == 0:
            step1 = model.readStep(['Den', 'TrP'], i+1)
            dt = step1.time
            n = np.append(n, 0.)

        # else we read the step before i, i, i+1 and i+2 (i being a saved timestep)
        else:
            step0 = model.readStep(['Den', 'TrP'], i-1)
            step1 = model.readStep(['Den', 'TrP'], i)
            step2 = model.readStep(['Den', 'TrP'], i+1)
            step3 = model.readStep(['Den', 'TrP'], i+2)

            ''' find in step1 and step2 values above pre-impact surface 
            (tracers from the projectile are also selected)'''

            # select where tracer heights are higher than the diameter of the
            # projectile (needs to fill this criteria for the three saved timesteps)
            ix = np.where((step1.ymark > threshold) & (step2.ymark > threshold) & (step3.ymark > threshold)
                          & (step3.ymark > step2.ymark) & (step2.ymark > step1.ymark))[0]

            # as we do this step every saved timestep, we want to only select
            # the same tracer only one time. To do so, we compare the detected
            # tracers at a timestep with tracers that
            # previously met this this criteria. Only new tracers are appended
            # to the array tracer_idx (or calculated)
            mask = np.setdiff1d(ix, tracer_idx)
            mask = mask.astype(int)

            # For the new detected tracers, positions are calculated
            # for step0, step1, step2 and step3
            x0 = np.append(x0, step0.xmark[mask])
            x1 = np.append(x1, step1.xmark[mask])
            x2 = np.append(x2, step2.xmark[mask])
            y0 = np.append(y0, step0.ymark[mask])
            y1 = np.append(y1, step1.ymark[mask])
            y2 = np.append(y2, step2.ymark[mask])

            # append only newly detected tracers
            tracer_idx = np.append(tracer_idx, mask)

            # the timestep when they were detected is also stored
            timestep = np.append(timestep, (i,) * len(mask))

            # the numbers of tracers per timestep are also stored
            n = np.append(n, len(mask))
    '''
    The flag (method) allows to use either two or three positions to calculate
    the ejection velocity and angles.
    '''
    # 3 dots method
    if method == 3:
        x = np.array((x0, x1, x2))
        y = np.array((y0, y1, y2))
        xx = np.swapaxes(x, 0, 1)
        yy = np.swapaxes(y, 0, 1)

    # 2 dots method
    if method == 2:
        x = np.array((x0, x1))
        y = np.array((y0, y1))
        xx = np.swapaxes(x, 0, 1)
        yy = np.swapaxes(y, 0, 1)

    # ellapsed time since the script has been started
    print datetime.now() - startTime

    # return an array with the positions (for a 3-dot or 2-dot methods)
    # the detected tracers, when they were detected, the time between saved
    # timesteps and the number of detected tracers per saved timestep
    return xx, yy, tracer_idx, timestep, dt, n


'''
***********************************************************************
'''


def linearfit(x, a, b):
    '''
    description:
    linear fit
    
    param:
    : x (1D-array): values
    : a (float) : Slope of linear equation
    : b (float): Intercept of linear equation
    '''

    return a*x + b


'''
***********************************************************************
'''


def parameters(x, y, tracer_idx , method, dt):
    '''
    description:
    calculate the ejection velocity, angle and position based on either the two
    or three positions (one before and one or two positions above the treshold)

    param:
    : x (float): positions of detected tracers (x) - (either 2 or 3 positions)
    : y (float): positions of detected tracers (y) - (either 2 or 3 positions)
    : tracer_idx (int) : index of ejected materials (lagrangian tracer indexs)
    : method (int) : 2 or 3 (two or three-dot methods)
    : dt : timestep (in absolute time in seconds)
        
    returns:
    : v (float): ejection velocity (in m/s)
    : angle (float): ejection angle (in degrees)
    : xpos (float): ejection position (in meters)
    '''

    # start the time at which the script has been started
    startTime = datetime.now()

    # creaty empty arrays
    xpos, v, angle = (np.array([]),)*3

    # get the number of tracers?
    n = np.shape(x)[0]

    # see
    # https://stackoverflow.com/questions/12473406/scipy-optimize-leastsq-returns-best-guess-parameters-not-new-best-fit
    x = x.astype(np.float64)
    y = y.astype(np.float64)

    # loop through all tracers (all materials that has been ejected during crater evolution)
    for idx in np.arange(n):

        # for each tracer (and the two or three positions)
        # we fit a linear equation
        popt, pcov = curve_fit(linearfit, x[idx], y[idx])

        # generate a fit through points
        xs = np.linspace(np.min(x[idx]), np.max(x[idx]), 10)
        ys = linearfit(xs, *popt)

        # launch position
        pos = - popt[1]/popt[0]

        # get the dx and dy
        xmin, xmax, ymin, ymax = np.min(xs), np.max(xs), np.min(ys), np.max(ys)

        # calculate the distance with Pythagore
        L = np.sqrt(((xmax - xmin)**2.) + ((ymax - ymin)**2.))

        # calculate the velocity of the ejecta! That should depend on the method
        # if it is three dots, than it should be 3.*dt
        v = np.append(v, (L / (method*dt)))  # should be model.dt

        # calculate the angle
        teta = np.arctan((ymax - ymin) / (xmax - xmin)) * (180. / np.pi)

        angle = np.append(angle, teta)

        # calculate the original position
        xpos = np.append(xpos, pos)

    # print the time ellapsed since we have started the script
    print datetime.now() - startTime

    # return data
    return v, angle, xpos


'''
***********************************************************************
'''


def wrap(model, method, thresholdf, g):
    '''
    description:
    get the extent of how much of the materials get excavated during crater evolution
    The volume,  depth and diameter of excavated materials are calculated
    
    params:
    : model (iSALE model) :
    : method (int): selection of two- (method = 2) or three-dot (3) methods
    : thresholdf (float): threshold value (which is multiplied by the CPPR x cellsize)
    : g (float): surface gravity of the planetary body (m2/s)
        

    returns:
    : Ve (float) : volume of excavated material
    : de (float) : maximum depth of excavated material
    : De (float) : diameter of excavated material
    : X (float) : original position of ejected materials (x)
    : Y (float) : original position of ejected materials (y)
    : Xc (float) : boundary hinge streamline (x)
    : Yc (float) : boundary hinge streamline (y)
    : timesteps (float) : timesteps when the materials are ejected
    '''
    threshold = (model.cppr[0]*model.dx) * thresholdf
    x, y, tracer_idx, timesteps, dt, n = posTracer(model, method, threshold)
    v, angle, xpos = parameters(x, y, tracer_idx, method, dt)
    tracer_idx = tracer_idx.astype(int)
    step = model.readStep(['Den', 'TrP'], 0)
    stepend = model.readStep(['Den', 'TrP'], model.nsteps-1)
    Ve, de, De = excavation(model, step, stepend, tracer_idx)
    X, Y, Xc, Yc = excavationProfile(model, step, tracer_idx)
    tair, xland = ballistic(v, angle, xpos, tracer_idx, g)

    return (Ve, de, De, X, Y, Xc, Yc, timesteps, tracer_idx, n, v, angle, xpos,
            tair, xland)


'''
***********************************************************************
'''


def main(path_jdata, path_tosave, method, thresholdf, g):
    '''
    description:
    loading of the data + all calculation in the function wrap
    
    params:
    : path_jdata (list) :
    : path_tosave (str) :
    : method (int) :
    : thresholdf (float) :
    : g (float) :
        
    returns:
    Three .txt files are saved:

    one text file that contains the depth, diameter and volume of the excavated
    cavity, name : modelname + '_excavated.txt'

    one text file that contains the tracers and their original position
    name : modelname + '_tracersXY_excavated.txt'

    one text file that contains the boundary between excavated and non excavated
    materials (i.e, the hinge streamline)
    name : modelname + '_tracersXY_contour_excavated.txt'
    '''

    # In case only a single model is specified, transform string to list of strings
    if type(path_jdata) == str:
        path_jdata = [path_jdata]        
    else:
        None


    for jdata_folder in path_jdata:

        # print modelname
        modelname = jdata_folder.split('/')[-2]
        print (modelname)
        
        # open model with pySALEplot
        jdata = os.path.join(jdata_folder, 'jdata.dat')
        model = psp.opendatfile(jdata)

        (Ve, de, Re, X, Y, Xc, Yc, timesteps, tracer_idx, n, v, angle, xpos,
         tair, xland) = wrap(model, method, thresholdf, g)

        # create plots directory
        path_data = path_tosave + modelname + "/excavated/"

        if not os.path.exists(path_data):
            os.makedirs(path_data)

        os.chdir(path_data)

        # DATA
        header_txt = "Volume_exc;depth_exc;diameter_exc"
        # we need to define the name of the txt file that will be saved
        fname = modelname + '_excavated.txt'
        # time, depth, diameter, volume
        output = np.column_stack((Ve, de, Re*2.))
        np.savetxt(fname, output, header=header_txt,
                   delimiter=';', fmt=['%1.6e', '%1.6e', '%1.6e'])

        # DATA TRACERS XY ALL
        fname = modelname + '_tracersXY_excavated.txt'
        header_txt = "X;Y"
        output = np.column_stack((X, Y))
        np.savetxt(fname, output, header=header_txt,
                   delimiter=';', fmt=['%1.6e', '%1.6e'])

        # DATA TRACERS XY contours
        fname = modelname + '_tracersXY_contour_excavated.txt'
        header_txt = "XC;YC"
        output = np.column_stack((Xc, Yc))
        np.savetxt(fname, output, header=header_txt,
                   delimiter=';', fmt=['%1.6e', '%1.6e'])

        # tracers
        path_ej = path_tosave + modelname + "/ejtracers/"

        if not os.path.exists(path_ej):
            os.makedirs(path_ej)

        os.chdir(path_ej)

        header_txt = "tracer_idx;timesteps;v;angle;xpos;tair;xland"
        # we need to define the name of the txt file that will be saved
        fname_ej = modelname + '_ejtracers.txt'
        # time, depth, diameter, volume
        output = np.column_stack((tracer_idx, timesteps, v, angle, xpos,
                                  tair, xland))
        np.savetxt(fname_ej, output, header=header_txt,
                   delimiter=';', fmt=['%1.6e', '%1.6e', '%1.6e',
                                       '%1.6e', '%1.6e', '%1.6e', '%1.6e'])

        header_txt = "ntracers"
        # we need to define the name of the txt file that will be saved
        fname_ej = modelname + '_ntracers.txt'
        # time, depth, diameter, volume
        output = n
        np.savetxt(fname_ej, output, header=header_txt,
                   delimiter=';', fmt='%1.6e')

        # close model file
        model.closeFile()


'''
***********************************************************************
'''


def parabolafit(x, a, b, c):
    '''
    description:
    parabola fit
    '''

    y = (a*(x**2)) + (b*x) + c

    return y


'''
***********************************************************************
'''


def ballistic(v, angle, xpos, tracer_idx, g):
    '''
    description:
    preliminary script to calculate where the material ejected will land
    (in order to calculate the thickness of the ejecta)
    More work is required here!!!
    '''

    # flight time
    tes = (2. * v * np.sin(angle * (np.pi/180))) / g

    # position where it landed
    land = xpos + (v*tes*np.cos(angle * (np.pi/180.)))

    # only take ejected materials that landed outside of the crater
    #mask = np.where(land > radius)
    #idx = mask[0]
    # I dunno, this things does not work so well
    return tes, land


'''
***********************************************************************
'''


def excavation(model, step, stepend, tracer_idx):
    '''
    description:
    calculate the volume, depth, and diameter of the exavated materials
    This task is done using position of tracers

    inputs:
    model: modelcase
    step: initial step
    stepend: final step
    tracer_idx: tracers detected as excavated materials

    outputs:
    Ve: excavated volume
    de: maximum excavated depth
    De: diameter of excavated cavity

    example:
    Ve, de, De = excavation(model,step,stepend,tracer_idx)
    '''

    mat1_start = model.tru[1].start
    mat1_end = model.tru[1].end

    # index for tracers of material 1
    ix_mat1 = np.where((tracer_idx >= mat1_start) & (tracer_idx <= mat1_end))

    # tracer number for material 1 and 2
    tr_idx1 = tracer_idx[ix_mat1]

    # calculate volume
    r1 = step.xmark[tr_idx1]
    rii1 = step.xmark[tr_idx1]+model.dx
    h = model.dy  # same values as in vimod (take that)
    R1 = r1 + ((1.) * (rii1 - r1))
    vol1 = np.pi * h * ((R1**2) - (r1**2))
    VT1 = np.sum(vol1)

    # calculate maximum depth of excavation
    de = np.min(step.ymark[tr_idx1])

    # calculate diameter of excavation
    De = np.max(step.xmark[tr_idx1])

    # return data
    return VT1, de, De


'''
***********************************************************************
'''


def excavationProfile(model, step, tracer_idx):
    '''
    description:
    calculate the volume, depth, and diameter of the exavated materials
    This task is done using position of tracers

    inputs:
    model: modelcase
    step: initial step
    stepend: final step
    tracer_idx: tracers detected as excavated materials

    outputs:
    X: return the original x-position of excavated materials
    Y: return the original y-position of excavated materials
    Xcontour: return the boundary of the excavated cavity (x)
    Ycontour: return the boundary of the excavated cavity (y)

    example:
    X, Y, Xcontour, Ycontour = excavationProfile(model,step,tracer_idx)
    '''

    mat1_start = model.tru[1].start
    mat1_end = model.tru[1].end

    # index for tracers of material 1
    ix_mat1 = np.where((tracer_idx >= mat1_start) & (tracer_idx <= mat1_end))

    # tracer number for material 1 and 2
    tr_idx1 = tracer_idx[ix_mat1]
    X = step.xmark[tr_idx1]
    Y = step.ymark[tr_idx1]

    # only take unique distance
    Xcontour = np.unique(X)

    # get the x-index
    xidx = Xcontour/(model.dx/2.)
    xidx = xidx.astype(int)

    # create empty matrix y-index
    Ycontour = np.ones(len(Xcontour))

    # the minimum origin emplacement of ejected materials are calculated
    # in order to get the boundary
    for i, var in np.ndenumerate(xidx):
        ii = i[0]
        idx_test = np.where((X == Xcontour[ii]))
        Ycontour[ii] = np.min(Y[idx_test])

    return X, Y, Xcontour, Ycontour


'''
***********************************************************************
PLAYGROUND


# load matplotlib
import matplotlib.pyplot as plt 

# make few plots
plt.plot(step.xmark[tracer_idx],step.ymark[tracer_idx],"ro",zorder=2)
plt.plot(step.xmark,step.ymark,"ko",zorder=1)




# get rid of the tracers from the projectile and plot it
mat0_start = model.tru[0].start
mat0_end = model.tru[0].end

ifff = np.where(tracer_idx > mat0_end)
tr_idx = tracer_idx[ifff]

plt.plot(step.xmark[tr_idx],step.ymark[tr_idx],"ro",zorder=2)
plt.plot(step.xmark,step.ymark,"ko",zorder=1)

# should use ix_mat2 if if you have more than one target
## calculate the volume and the density (I should get the tracer with the density) I can have the density based on the tru_u

# start and end tracer number for material 1 (lower layer) and material 2 (upper layer)
mat1_start = model.tru[1].start
mat1_end = model.tru[1].end

#mat2_start = model.tru[2].start
#mat2_end = model.tru[2].end

# index for tracers of material 1 and 2
ix_mat1 = np.where((tracer_idx >= mat1_start) & (tracer_idx <= mat1_end))
#ix_mat2 = np.where((tracer_idx >= mat2_start) & (tracer_idx <= mat2_end))

# tracer number for material 1 and 2
tr_idx1 = tracer_idx[ix_mat1]
#tr_idx2 = tracer_idx[ix_mat2]

#Let's check if we got the right things
plt.plot(step.xmark[tr_idx1],step.ymark[tr_idx1],"ro",zorder=2)
plt.plot(step.xmark,step.ymark,"ko",zorder=1)

# get the original density of material 1 and 2 (much simpler if it is only one material)


if tr_idx1.size == 0:
    dens1 = np.nan
else:
    x_mat1 , y_mat1 = step.xmark[tr_idx1[0]], step.ymark[tr_idx1[0]]
    # the closest to
    ix_mat1 = np.where(((np.round(model.xc,decimals=2)) == np.round(x_mat1,decimals=2)) &
                       ((np.round(model.yc,decimals=2)) == np.round(y_mat1,decimals=2)))

    dens1 = step.data[0][ix_mat1].data[0]
    

if tr_idx2.size == 0:
    dens2 = np.nan
else:
    x_mat2 , y_mat2 = step.xmark[tr_idx2[0]], step.ymark[tr_idx2[0]]
    # the closest to (use the argmin)
    ix_mat2 = np.where(((np.round(model.xc,decimals=2)) == np.round(x_mat2,decimals=2)) &
                       ((np.round(model.yc,decimals=2)) == np.round(y_mat2,decimals=2)))
    dens2 = step.data[0][ix_mat2].data[0]

# need to calculate density
if tr_idx1.size == 0:
    overall_dens = dens2
elif tr_idx2.size == 0:
    overall_dens = dens1
else:
    # calculate the volume of mat1 and mat2 (i guess I can calculate from the volume the density)
    r1 = step.xmark[tr_idx1]
    rii1 = step.xmark[tr_idx1]+model.dx
    h = model.dy # same values as in vimod (take that)
    R1 = r1 + ((1.) * (rii1 - r1))
    vol1 = np.pi * h * ((R1**2) - (r1**2))
    VT1 = np.sum(vol1)
    
    r2 = step.xmark[tr_idx2]
    rii2 = step.xmark[tr_idx2]+model.dx
    R2 = r2 + ((1.) * (rii2 - r2))
    vol2 = np.pi * h * ((R2**2) - (r2**2))
    VT2 = np.sum(vol2)
    
    #volumetric density
    ovearall_dens = ((VT1/(VT2+VT1)) * dens1) + ((VT2/(VT2+VT1)) * dens2)
    


# we use the initial density to calculate the excavated volume

dens = step.data[0][0,0] # we do not actually need the density for calculating the volume

r1 = step.xmark[tr_idx1]
rii1 = step.xmark[tr_idx1]+model.dx
h = model.dy # same values as in vimod (take that)
R1 = r1 + ((1.) * (rii1 - r1))
vol1 = np.pi * h * ((R1**2) - (r1**2))
VT1 = np.sum(vol1)

print VT1/1e9 #in km3

# Volume excavated
Ve = VT1/1e9

# maximum depth excavated
de = np.min(step.ymark[tr_idx1])
De = np.max(step.xmark[tr_idx1])

'''
