# -*- coding: utf-8 -*-
"""
Created on Sun Sep 09 15:39:47 2018

@author: Chocobo
"""
import glob
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import copy

sys.path.append('M:/Nils/Python/arcpy/') 
import geomorphCraters as wk

path = 'X:/Moon/downloading/STEPMED/ASCIIreproj/'
pathplot = 'X:/Moon/downloading/STEPMED/plots_detection_ndata2_05112018/'
pathdata = 'X:/Moon/downloading/STEPMED/ndata2_05112018/'


'''
**************************************************************************
'''

def ifnot_mkdir(path):

    if os.path.isdir(path):
        pass
    else:
        os.mkdir(path)
        print ("A new directory has been created")
'''
**************************************************************************
'''

def loadgis(path,filename):
    
    os.chdir(path)
    data = np.loadtxt(filename, delimiter=";", skiprows=1)
    
    return data
       
'''
**************************************************************************
'''
def run1(path, filenameXY, filenamecrater, pathplot, pathdata):
    
    '''
    filenameXY = 'data.txt'
    filenamecrater = 'crater_id.txt'
    '''

    # the ASCII file (converted from the raster in ArcGIS) is loaded
    
    
    os.chdir(path)
        
    filenames = glob.glob("crater*.asc")
    filenames.sort(key=wk.tokenize)
    dataXYD = np.loadtxt(filenameXY,delimiter=";",comments="#")
    xcrater = dataXYD[:,0]
    ycrater = dataXYD[:,1]
    diam0 = dataXYD[:,2]
    diam0 = diam0 * 1000.
    r0 = diam0 /2.
    
    #it takes way too much time this way
    #xarcout = np.array([])
    #yarcout = np.array([])
    
    #if it is the last step
    
    for indf, filename in enumerate(filenames):
        
        print (indf)
        
        # should be good this way
        data = np.loadtxt(path + filename, skiprows=6)
        data2 = np.rot90(data)
        data3 = np.rot90(data2)
        data4 = np.rot90(data3)
        
        del data
        del data2
        del data3
        data = copy.deepcopy(data4*0.50) #multiplying factor of the elevation
        del data4
        
        #correspond to orgin in the lower-left corner and data[x,y] 
        # with x = ncols, y = nrows !!! important
                
        # the number of columns, rows and extent of the raster is read
        (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path,filename)
        

        # if not the same number of rows or columns
        if ncols > nrows:
            data = data[:nrows,:nrows]
            ncols = nrows
        elif nrows > ncols:
            data = data[:ncols,:ncols]
            nrows = ncols
        else:
            pass
        
                
        #print diam
        # x and y axes are created (x and y have the same size so it is not a big deal
        # if we mix up here)
        x = np.arange(0,((ncols)*cellsize),cellsize) 
        y = np.arange(0,((nrows)*cellsize),cellsize)
        
        # x and y at the center of the cell
        xe = np.arange(cellsize/2.,((ncols)*cellsize),cellsize)
        ye = np.arange(cellsize/2.,((nrows)*cellsize),cellsize)
        
        if len(x) > len(xe):
            x = x[:len(xe)]
            y = y[:len(xe)]
        
        # create my own matrices for the plotting with pcolor or pcolormesh
        xc = np.zeros_like(data)
        yc = np.zeros_like(data)
        xce = np.zeros_like(data)
        yce = np.zeros_like(data)
        
        for i in range(ncols):
            xc[:,i] = x
            xce[:,i] = xe
            
        for i in range(nrows):
            yc[i,:] = y
            yce[i,:] = ye
        
        # get the cell the closest to the center of the crater (this part has to be changed as the proj is changed)
        #__, ncenterx = wk.find_nearest(xllcorner+xe, xcrater[indf]) #correspond to columns and x
        
        # because yllcorner is from the top left now
        #__, ncentery = wk.find_nearest(yllcorner+ye,ycrater[indf]) #correspond to rows and y
        
        # middle of the dtm
        ncentery = int(nrows/2) #Python 3 770/2 is not directly integer anymore
        ncenterx = int(nrows/2)
        
        # Average height at 1R, 2R, 3R, 3.5R and 4R from the centre of the crater
        # (need to be replace with the right centery and x)
        
        #other places where ncenterx and ncentery is required 
        x1, y1 = wk.xy_circle(1.0*r0[indf], xe[ncenterx], ye[ncentery]) # xe and ye are equals
        #x2, y2 = wk.xy_circle(2.0*r0[indf], xe[ncenterx], ye[ncentery])
        #x3, y3 = wk.xy_circle(3.0*r0[indf], xe[ncenterx], ye[ncentery])
        #x35, y35 = wk.xy_circle(3.5*r0[indf], xe[ncenterx], ye[ncentery])
        #x4, y4 = wk.xy_circle(4.0*r0[indf], xe[ncenterx], ye[ncentery])
        
        # First detrending (a plane is fitted through the elevations taken at circles
        # 2.0 and 3.0R from the center of the crater )
        
        stduse = True #use standard deviation removal
        ndata = wk.detrending(xc, yc, r0[indf], cellsize, ncenterx, ncentery, data, stduse)
                              
        # txt name
        name_crater_txt = filename.split('.asc')[0] + 'XY.txt'
                        
        #and if XY is not found in a certain path
        if os.path.isfile(pathdata + name_crater_txt):
            pass
        
        else:
    
            # the Maximum elevation are selected a first time
            first_run = True
            mingrade = 0.05
            minclust = 0.05
            slen = 0.1
            
            (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME) = (
             wk.rim(xc, yc, ncenterx, ncentery, ndata, r0[indf], 
                    cellsize, slen, minclust, mingrade, first_run))
            
            # we get rid of the nan-values for the detrending
            ixnan = np.where(np.isnan(col_coord_ME) == False)
            
            # Maximum elevations without NaNs
            coldetrend = col_coord_ME[ixnan]
            rowdetrend = row_coord_ME[ixnan]
            elev4detrend = elev_ME[ixnan]
            
            # the second detrending through the selected maximum elevation points are done
            stduse = False
            Z2 = wk.linear3Ddetrending(coldetrend,rowdetrend,elev4detrend, xc, yc, stduse)
            
            # the detrended plane is substracted to the DEM
            ndata2 = ndata - Z2
            
            # the first detrending has ben run so we change now the first_run flag to false
            first_run = False
            
            # second run of the function rim
            (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME,
             col_coord_LE, row_coord_LE, col_cells_LE, row_cells_LE, elev_LE, prof_LE,
             col_coord_BS, row_coord_BS, col_cells_BS, row_cells_BS, elev_BS, prof_BS) = (
             wk.rim(xc, yc, ncenterx, ncentery, ndata2, r0[indf], 
                    cellsize, slen, minclust, mingrade, first_run))
            
            
            # takes only nonnan values from *BS data and merge it to local elevation data
            ixfinite = np.where(np.isfinite(col_coord_BS) == True)
            col_coord_BS = col_coord_BS[ixfinite] 
            row_coord_BS = row_coord_BS[ixfinite] 
            col_cells_BS = col_cells_BS[ixfinite] 
            row_cells_BS = row_cells_BS[ixfinite] 
            elev_BS = elev_BS[ixfinite] 
            prof_BS = prof_BS[ixfinite] 
            
            use_break_in_slope = False
            
            if use_break_in_slope:
                # concatenate local maxima and slope change!
                colint = np.concatenate((col_coord_LE, col_coord_BS))
                rowint = np.concatenate((row_coord_LE, row_coord_BS))
                colmap = np.concatenate((col_cells_LE, col_cells_BS)) 
                rowmap = np.concatenate((row_cells_LE, row_cells_BS))
                elevint = np.concatenate((elev_LE, elev_BS))
                profint = np.concatenate((prof_LE, prof_BS))
            
            else:
                # or using only local maxima
                colint = col_coord_LE
                rowint = row_coord_LE
                colmap = col_cells_LE
                rowmap = row_cells_LE
                elevint = elev_LE
                profint = prof_LE
                
            # Maximum allowed radial discontinuity Drad (I should convert these values in cells)
            Drad = 0.1 * r0[indf]
            
            # Distance of interest (searching distance)
            Dint = 0.05 * r0[indf]
            
            # Maximum angular discontinuity (avoid unnecessary large gap angle in the 
            # data)
            angle = 2.0 #(in degrees)
            
            stangle = [0,45,90,135,180,225,270,315]
            
            contloop = True
            siftRedundant = True
            kpstitch = False
            
            
            OptRims, Omegas, gap, maxradf = wk.rim_composite(col_coord_ME, row_coord_ME, col_cells_ME,
                                                          row_cells_ME, elev_ME, prof_ME,
                                                          colint, rowint, colmap, rowmap, xc, yc, elevint, 
                                                          profint, ncenterx, ncentery,
                                                          angle, stangle, Drad, Dint, 
                                                          contloop, siftRedundant, kpstitch)
            
            # we want to export the data as X, Y, Z
            # find the best fit
            for i in range(np.shape(OptRims)[0]):
                a = np.mean(maxradf[i]) + (0.5*Omegas[i])
                if i == 0:
                    b = a
                    c = i
                else:
                    if a < b:
                        b = a
                        c = i
            
            
            # I should get rid of the 0,0 counts and should not append it takes way too much time
            
            # temporary name
            name_crater_png = filename.split('.asc')[0] + 'XY.png'
            
            #xarcout = np.append(xarcout, np.array(OptRims[c][0,:]) + xllcorner)
            #yarcout = np.append(yarcout, np.array(OptRims[c][1,:]) + yllcorner)
            
            plt.ioff()
            plt.figure()
            plt.pcolormesh(xc, yc, ndata2)
            plt.title(filename.split('.asc')[0],fontsize=16)
            plt.colorbar()
            
            xarcgis = np.array(OptRims[c][0,:])
            xarcgis = xarcgis[np.nonzero(xarcgis)] 
            
            yarcgis = np.array(OptRims[c][1,:])
            yarcgis = yarcgis[np.nonzero(yarcgis)]
            
            zarcgis = np.array(OptRims[c][2,:])
            zarcgis = zarcgis[np.nonzero(yarcgis)]
            
            profgis = np.array(OptRims[c][3,:])
            profgis = profgis[np.nonzero(yarcgis)]
            
            flaggis = np.array(OptRims[c][4,:])
            flaggis = flaggis[np.nonzero(yarcgis)]
            
            #xarcgis = np.array(OptRims[c][0,:]) + xllcorner
            #xarcgis = xarcgis[np.nonzero(xarcgis - xllcorner)] 
            
            #yarcgis = np.array(OptRims[c][1,:]) + yllcorner
            #yarcgis = yarcgis[np.nonzero(yarcgis - yllcorner)]
            
            #zarcgis = np.array(OptRims[c][2,:])
            #zarcgis = zarcgis[np.nonzero(yarcgis - yllcorner)]
            
            #profgis = np.array(OptRims[c][3,:])
            #profgis = profgis[np.nonzero(yarcgis - yllcorner)]
            
            #flaggis = np.array(OptRims[c][4,:])
            #flaggis = flaggis[np.nonzero(yarcgis - yllcorner)]
            
            #print np.shape(OptRims)
            plt.plot(xarcgis,yarcgis,"ko")
            plt.plot(x1,y1,"b")
            
            
            # save plot
            plt.savefig(pathplot + name_crater_png,dpi=300)
            plt.close()

            #export the data to a .csv file with ";"
            #xarcgis = np.array(OptRims[c][0,:]) + xllcorner
            #yarcgis = np.array(OptRims[c][1,:]) + yllcorner
            #zarcgis = np.array(OptRims[c][2,:])
            header_txt = 'X;Y;Z;prof;flag'
            
            arcout = np.column_stack((xarcgis, yarcgis, zarcgis, profgis, flaggis))
            np.savetxt(pathdata + name_crater_txt, arcout, delimiter = ";", header=header_txt,fmt='%10.5f', comments='')
            
            
            
                
'''
**************************************************************************
'''

def run2(path, filenameXY, filenamecrater, pathdata):
    
    '''
    filenameXY = 'data.txt'
    filenamecrater = 'crater_id.txt'
    '''        
    os.chdir(path)
    
    filenames = glob.glob("crater*.asc")
    filenames.sort(key=wk.tokenize)    
    crater_id = np.genfromtxt(filenamecrater,skip_header=1,dtype=None)
    dataXYD = np.loadtxt(filenameXY,delimiter=";",comments="#")
    xcrater = dataXYD[:,0]
    ycrater = dataXYD[:,1]
    diam0 = dataXYD[:,2]
    diam0 = diam0 * 1000.
    r0 = diam0 /2.    
  
    unc_cse = np.ones(len(crater_id))
    med_cse  = np.ones(len(crater_id))
    cse_25 = np.ones(len(crater_id))
    cse_75 = np.ones(len(crater_id))
    cse_min = np.ones(len(crater_id))
    cse_max = np.ones(len(crater_id))
    
    
    unc_mcw  = np.ones(len(crater_id))
    med_mcw = np.ones(len(crater_id))
    mcw_25  = np.ones(len(crater_id))
    mcw_75 = np.ones(len(crater_id))
    mcw_min  = np.ones(len(crater_id))
    mcw_max  = np.ones(len(crater_id))

    
    unc_ucw = np.ones(len(crater_id))
    med_ucw = np.ones(len(crater_id))
    ucw_25  = np.ones(len(crater_id))
    ucw_75 = np.ones(len(crater_id))
    ucw_min  = np.ones(len(crater_id))
    ucw_max  = np.ones(len(crater_id))
    
    unc_h = np.ones(len(crater_id))
    med_h = np.ones(len(crater_id))
    h_25  = np.ones(len(crater_id))
    h_75 = np.ones(len(crater_id))
    h_min  = np.ones(len(crater_id))
    h_max  = np.ones(len(crater_id))    
    
    depthf = np.ones(len(crater_id))
    
    diamf = np.ones(len(crater_id))
    med_diam = np.ones(len(crater_id))
    diam_25  = np.ones(len(crater_id))
    diam_75 = np.ones(len(crater_id))
    diam_min  = np.ones(len(crater_id))
    diam_max  = np.ones(len(crater_id))     
    
    for indf, filename in enumerate(filenames):
    
        print (indf)
        
        # load data
        data = np.loadtxt(path + filename, skiprows=6)
        data2 = np.rot90(data)
        data3 = np.rot90(data2)
        data4 = np.rot90(data3)
        
        del data
        del data2
        del data3
        data = copy.deepcopy(data4*0.50) # multiplying by scaling factor
        del data4
        
        # txt name
        name_crater_txt = filename.split('.asc')[0] + 'XY.txt'
        

        # the number of columns, rows and extent of the raster is read
        (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path,filename)
        
        # if not the same number of rows or columns
        if ncols > nrows:
            data = data[:nrows,:nrows]
            ncols = nrows
        elif nrows > ncols:
            data = data[:ncols,:ncols]
            nrows = ncols
        else:
            pass
        
        # define a few variables
        x = np.arange(0,((ncols)*cellsize),cellsize) 
        y = np.arange(0,((nrows)*cellsize),cellsize)
        
        # x and y at the center of the cell
        xe = np.arange(cellsize/2.,((ncols)*cellsize),cellsize)
        ye = np.arange(cellsize/2.,((nrows)*cellsize),cellsize)
        
        if len(x) > len(xe):
            x = x[:len(xe)]
            y = y[:len(xe)]
        
        # create my own matrices for the plotting with pcolor or pcolormesh
        xc = np.zeros_like(data)
        yc = np.zeros_like(data)
        xce = np.zeros_like(data)
        yce = np.zeros_like(data)
        
        for i in range(ncols):
            xc[:,i] = x
            xce[:,i] = xe
            
        for i in range(nrows):
            yc[i,:] = y
            yce[i,:] = ye
        
        #load the data
        datagis = loadgis(pathdata, name_crater_txt)
        xarcgis = datagis[:,0]
        yarcgis = datagis[:,1]
        zarcgis = datagis[:,2]
        profgis = datagis[:,3]
        flaggis = datagis[:,4]
               
        #get rid of outliers 
        #Q1 = np.percentile(zarcgis,25.)
        #Q3 = np.percentile(zarcgis,75.)
        #IQR = Q3 - Q1
        #outliers = (zarcgis < (Q1 - 1.5 * IQR)) | (zarcgis > (Q3 + 1.5 * IQR))
        #not_outliers = np.where(outliers == False)
        
        
        x_not_outliers = xarcgis 
        y_not_outliers = yarcgis
        z_not_outliers = zarcgis
        prof_not_outliers = profgis
        
        # old center
        # get the cell the closest to the center of the crater
        #__, ncenterx_old = wk.find_nearest(xllcorner+xe, xcrater[indf]) #correspond to columns and x
        
        # because yllcorner is from the top left now
        #__, ncentery_old = wk.find_nearest(yllcorner+ye,ycrater[indf]) #correspond to rows and y
        ncenterx_old = nrows/2
        ncentery_old = nrows/2
        
        #fit a circle
        xnewcenter, ynewcenter, rnew, residu = wk.leastsq_circle(x_not_outliers,y_not_outliers)
                    
        # get the cell the closest to the center of the crater (new calculation)
        __, ncenterx = wk.find_nearest(xe,xnewcenter)
        __, ncentery = wk.find_nearest(ye,ynewcenter)

                   
        #other places where ncenterx and ncentery is required (new calculation)
        #x1, y1 = wk.xy_circle(1.0*rnew, xe[ncentery], ye[ncenterx])

        
        # First detrending (a plane is fitted through the elevations taken at circles
        # with new center of the circle
        try:
            stduse = True
            ndata = wk.detrending(xc, yc, rnew, cellsize, ncenterx_old, ncentery_old, data, stduse) #detrending
            # with old data otherwise sometimes the dem area is not large enough
            
            # the second detrending through the selected maximum elevation points are done
            stduse = False
            Z2 = wk.linear3Ddetrending(x_not_outliers,y_not_outliers,z_not_outliers, xc, yc, stduse)
                
            # the detrended plane is substracted to the DEM (with the newly detected)
            ndata2 = ndata - Z2
            
            
            # I should actually rerun the rim detection things one  more time
            # as the centre of the crater has been moved (not the same profiles anymore)
            
            
            #need to put some indices there (maybe, they will have different sizes)
            (R_upcw, R_ufrc, cse, slope_mcw, slope_ucw, slope_fsa, slope_lrs, slope_urs,
            h, depth, diam, nnn, prf, crossSections) = (wk.calculation(xc, yc, x_not_outliers, 
                              y_not_outliers, z_not_outliers, prof_not_outliers,
                rnew, ndata2, cellsize, ncenterx, ncentery, xllcorner, yllcorner))
            
			# crossSections has been added
            
            # the problem is that there are several values for profile = 0? for some reasons?
            # I don't understand that, I guess it must be profiles where....
            # np.where(prf != 0) to get rid of it but it would be nice if the problem would be detected
            
            #just need to make the calculation, we have the final circle
            # need to find the closest values that does not overlap
            
            unc_cse[indf] = np.nanstd(cse)/np.sqrt(nnn)
            med_cse[indf] = np.nanmedian(cse)
            cse_25[indf] = np.nanpercentile(cse,25)
            cse_75[indf] = np.nanpercentile(cse,75)
            cse_min[indf] = np.nanmin(cse)
            cse_max[indf] = np.nanmax(cse)
        
            unc_mcw[indf] = np.nanstd(slope_mcw)/np.sqrt(nnn)
            med_mcw[indf] = np.nanmedian(slope_mcw)
            mcw_25[indf] = np.nanpercentile(slope_mcw,25)
            mcw_75[indf] = np.nanpercentile(slope_mcw,75)
            mcw_min[indf] = np.nanmin(slope_mcw)
            mcw_max[indf] = np.nanmax(slope_mcw)
            
            unc_ucw[indf] = np.nanstd(slope_ucw)/np.sqrt(nnn)
            med_ucw[indf] = np.nanmedian(slope_ucw)
            ucw_25[indf] = np.nanpercentile(slope_ucw,25)
            ucw_75[indf] = np.nanpercentile(slope_ucw,75)
            ucw_min[indf] = np.nanmin(slope_ucw)
            ucw_max[indf] = np.nanmax(slope_ucw)
            
            unc_h[indf] = np.nanstd(h)/np.sqrt(nnn)
            med_h[indf] = np.nanmedian(h)
            h_25[indf] = np.nanpercentile(h,25)
            h_75[indf] = np.nanpercentile(h,75)
            h_min[indf] = np.nanmin(h)
            h_max[indf] = np.nanmax(h)
        
            depthf[indf] = depth
            
            diamf[indf] = np.nanstd(diam)/np.sqrt(nnn)
            med_diam[indf] = np.nanmedian(diam)
            diam_25[indf] = np.nanpercentile(diam,25)
            diam_75[indf] = np.nanpercentile(diam,75)
            diam_min[indf] = np.nanmin(diam)
            diam_max[indf] = np.nanmax(diam)
                    
            header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
                          'depth;' + 
                          'mh;uh;h25;h75;hmin;hmax;' +
                          'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
                          'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
                          'mcse;ucse;cse_25;cse_75;cse_min;cse_max;')
            
        except:
            unc_cse[indf] = np.nan
            med_cse[indf] = np.nan
            cse_25[indf] = np.nan
            cse_75[indf] = np.nan
            cse_min[indf] = np.nan
            cse_max[indf] = np.nan
            
            unc_mcw[indf] = np.nan
            med_mcw[indf] = np.nan
            mcw_25[indf] = np.nan
            mcw_75[indf] = np.nan
            mcw_min[indf] = np.nan
            mcw_max[indf] = np.nan
            
            unc_ucw[indf] = np.nan
            med_ucw[indf] = np.nan
            ucw_25[indf] = np.nan
            ucw_75[indf] = np.nan
            ucw_min[indf] = np.nan
            ucw_max[indf] = np.nan
            
            unc_h[indf] = np.nan
            med_h[indf] = np.nan
            h_25[indf] = np.nan
            h_75[indf] = np.nan
            h_min[indf] = np.nan
            h_max[indf] = np.nan
            
            depthf[indf] = np.nan
            
            diamf[indf] = np.nan
            med_diam[indf] = np.nan
            diam_25[indf] = np.nan
            diam_75[indf] = np.nan
            diam_min[indf] = np.nan
            diam_max[indf] = np.nan
            
            header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
              'depth;' + 
              'mh;uh;h25;h75;hmin;hmax;' +
              'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
              'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
              'mcse;ucse;cse_25;cse_75;cse_min;cse_max;')
            
        # txt name
        name_crater_f = filename.split('.asc')[0] + '_res.txt'
            
        arc = np.column_stack((med_diam[indf], diamf[indf], diam_25[indf], diam_75[indf], diam_min[indf], diam_max[indf],
                               depthf[indf], 
                               med_h[indf], unc_h[indf], h_25[indf], h_75[indf], h_min[indf], h_max[indf], 
                               med_mcw[indf], unc_mcw[indf], mcw_25[indf], mcw_75[indf], mcw_min[indf], mcw_max[indf],
                               med_ucw[indf], unc_ucw[indf], ucw_25[indf], ucw_75[indf], ucw_min[indf], ucw_max[indf],
                               med_cse[indf], unc_cse[indf], cse_25[indf], cse_75[indf], cse_min[indf], cse_max[indf]))
        
        np.savetxt(pathdata + name_crater_f, arc, delimiter = ";", header=header_txt,fmt='%10.5f', comments='#')
		
        # added part about saving cross sections. I still need to test if this good or not
        name_crater_cross = filename.split('.asc')[0] + '_cross_sections.txt'
        
        # in order to save it I have to transform the dictionnary into an array
        # I need to get the maximum length
        maxlength_crossSections = 0
        for xx_tmp, x_tmp in crossSections.items():
            if len(x_tmp) > maxlength_crossSections:
                maxlength_crossSections = len(x_tmp)
            else:
                None
                
        number_crossSections = len(crossSections)
        
        # create array
        array_crossSections = np.empty((maxlength_crossSections, number_crossSections))
        array_crossSections[:] = np.nan
                
        # fill array with values
        for xx_tmp, x_tmp in crossSections.items():
            maxlength_crossSections = len(x_tmp)
            array_crossSections[:maxlength_crossSections,xx_tmp] = x_tmp
        
        np.savetxt(pathdata + name_crater_cross, array_crossSections, delimiter = ";",fmt='%10.5f', comments='#')
		
    return (med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max)

'''
**************************************************************************
'''

def mergedata(path,pathdata,filenamecrater):
    
    '''
    path = "X:/Moon/downloading/STEPMED/ASCII/"
    pathdata = "X:/Moon/downloading/STEPMED/firstrun2/"
    filenamecrater = 'crater_id.txt'
    (med_diam, diamf, depthf, med_h, unc_h, med_mcw, unc_mcw, med_ucw, unc_ucw, med_cse, unc_cse) = mergedata(path,pathdata,filenamecrater)
    '''
    
    
    os.chdir(path)
    
    filenames = glob.glob("crater*.asc")
    filenames.sort(key=wk.tokenize)    
    crater_id = np.genfromtxt(filenamecrater,skip_header=1,dtype=None)
    
    unc_cse = np.ones(len(crater_id))
    med_cse  = np.ones(len(crater_id))
    cse_25 = np.ones(len(crater_id))
    cse_75 = np.ones(len(crater_id))
    cse_min = np.ones(len(crater_id))
    cse_max = np.ones(len(crater_id))
    
    
    unc_mcw  = np.ones(len(crater_id))
    med_mcw = np.ones(len(crater_id))
    mcw_25  = np.ones(len(crater_id))
    mcw_75 = np.ones(len(crater_id))
    mcw_min  = np.ones(len(crater_id))
    mcw_max  = np.ones(len(crater_id))

    
    unc_ucw = np.ones(len(crater_id))
    med_ucw = np.ones(len(crater_id))
    ucw_25  = np.ones(len(crater_id))
    ucw_75 = np.ones(len(crater_id))
    ucw_min  = np.ones(len(crater_id))
    ucw_max  = np.ones(len(crater_id))
    
    unc_h = np.ones(len(crater_id))
    med_h = np.ones(len(crater_id))
    h_25  = np.ones(len(crater_id))
    h_75 = np.ones(len(crater_id))
    h_min  = np.ones(len(crater_id))
    h_max  = np.ones(len(crater_id))    
    
    depthf = np.ones(len(crater_id))
    
    diamf = np.ones(len(crater_id))
    med_diam = np.ones(len(crater_id))
    diam_25  = np.ones(len(crater_id))
    diam_75 = np.ones(len(crater_id))
    diam_min  = np.ones(len(crater_id))
    diam_max  = np.ones(len(crater_id))   
    
    for indf, filename in enumerate(filenames):
    
        name_crater_f = filename.split('.asc')[0] + '_res.txt'
        
        (med_diam[indf], diamf[indf], diam_25[indf], diam_75[indf], diam_min[indf], diam_max[indf], 
         depthf[indf], 
         med_h[indf], unc_h[indf], h_25[indf], h_75[indf], h_min[indf], h_max[indf],
         med_mcw[indf], unc_mcw[indf], mcw_25[indf], mcw_75[indf], mcw_min[indf], mcw_max[indf],
         med_ucw[indf], unc_ucw[indf], ucw_25[indf], ucw_75[indf], ucw_min[indf], ucw_max[indf],
         med_cse[indf], unc_cse[indf], cse_25[indf], cse_75[indf], cse_min[indf], cse_max[indf]) = np.loadtxt(pathdata + name_crater_f,delimiter=";", comments="#")
        
    return (med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max)
      
'''
**************************************************************************
'''

# define the folder containing the DTMs, and the folders where the plot
# and data routines will save data
path = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/ASCII_28022019/'
pathplot = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/plots_28022019/'
pathdata = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/data_28022019/'

# Location + cratername
filenameXY = 'data.txt'
filenamecrater = 'crater_id.txt'

#make directory if not existing
ifnot_mkdir(pathplot)
ifnot_mkdir(pathdata)

# run step 1
run1(path, filenameXY, filenamecrater, pathplot, pathdata)

# run step 2
(med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max) = run2(path, filenameXY, filenamecrater, pathdata)


## load data 

(med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max) = mergedata(path,pathdata,filenamecrater)


header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
                          'depth;' + 
                          'mh;uh;h25;h75;hmin;hmax;' +
                          'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
                          'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
                          'mcse;ucse;cse_25;cse_75;cse_min;cse_max;dD')

dD = (med_h-depthf) / med_diam

name_crater_f = 'final_res.txt'
            
arc = np.column_stack((med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max, dD))
        
np.savetxt(pathdata + name_crater_f, arc, delimiter = ";", header=header_txt,fmt='%10.5f', comments='#') 
        