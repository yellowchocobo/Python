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
import pandas as pd

#sys.path.append('M:/Nils/Python/arcpy/')

sys.path.append('/uio/kant/geo-ceed-u1/nilscp/Nils/Python/arcpy/') 
import geomorphCraters as wk
import calculate_volume as cv

#path = 'X:/Moon/downloading/STEPMED/ASCIIreproj/'
#pathplot = 'X:/Moon/downloading/STEPMED/plots_detection_ndata2_05112018/'
#pathdata = 'X:/Moon/downloading/STEPMED/ndata2_05112018/'


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
def run1(path, filenameXY, filenamecrater, pathplot, pathdata, scaling_factor, centre_craters = ''):
    
    '''
    filenameXY = 'data.txt'
    filenamecrater = 'crater_id.txt'
    '''

    # the ASCII file (converted from the raster in ArcGIS) is loaded
    
    
    os.chdir(path)
    
    dataXYD = np.loadtxt(filenameXY,delimiter=";",comments="#")
    crater_id = np.genfromtxt(filenamecrater,skip_header=1,dtype=str)
    
    # centre of craters
    if centre_craters.endswith('.csv'):
        ccrater = np.loadtxt(centre_craters,delimiter=",",comments="#")
        xccrater = ccrater[:,0]
        yccrater = ccrater[:,1]
        flagccrater = ccrater[:,2]
    else:
        None                      
                         
    diam0 = dataXYD[:,2]
    #diam0 = diam0 * 1000. # comment for coldspots
    r0 = diam0 /2.
    
    # sort the files based on the size so the smallest files will be taken first    
    #size_of_asciis = np.zeros(len(crater_id))
    #
    #for indf, filename in enumerate(crater_id):
    #    
    #    size_of_asciis[indf] = os.path.getsize(path + filename + ".asc")
        
    # sort from the smallest to the largest
    #idx_sort = np.argsort(size_of_asciis)
    #crater_id_new = crater_id[idx_sort]
    #r0_new = r0[idx_sort]
    
    crater_id_new = crater_id
    r0_new = r0
    
    for indf, filename in enumerate(crater_id_new):
        
        print (indf)
        
        name_crater_txt = filename + 'XY.txt'
        
        
        # if the flag is equal to 0 then let's do nothing
        if filename == 'cpcrater0143':
            None

        elif filename == 'cpcrater0199':
            None
            
        elif filename == 'cpcrater0288':
            None
            
        elif filename == 'cpcrater0998':
            None
            
        elif filename == 'cpcrater1315':
            None    
        elif filename == 'cpcrater1387':
            None   
        elif filename == 'cpcrater1562':
            None
        elif filename == 'cpcrater1875':
            None
        elif filename == 'cpcrater2005':
            None  
        elif filename == 'cpcrater2106':
            None  
        # if the file exists from a previous run, let's do nothing    
        elif os.path.isfile(pathdata + name_crater_txt):
            None
        
        elif (centre_craters.endswith('.csv')):
            
            if flagccrater[indf] == 0:
                None
                
            else:
                # should be good this way
                data = pd.read_csv(path + filename + ".asc", delimiter=' ', skiprows=6).values
                data2 = np.rot90(data)
                data3 = np.rot90(data2)
                data4 = (np.rot90(data3)) * scaling_factor #scaling factor # only for sldem
                
                del data
                del data2
                del data3
                data = copy.deepcopy(data4) #multiplying factor of the elevation (not for kaguya, 0.5 for Kaguya + LOLA)
                del data4
                
                #correspond to orgin in the lower-left corner and data[x,y] 
                # with x = ncols, y = nrows !!! important
                        
                # the number of columns, rows and extent of the raster is read
                (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path,filename + ".asc")
        
                # take the minimum of ncols, nrows and ncols and nrows of data
                nrowsd, ncolsd = np.shape(data)
                ncolsnrows = np.min([ncols, nrows, ncolsd, nrowsd])
                data = data[:ncolsnrows, :ncolsnrows]
                
                print (ncols, nrows, ncolsd, nrowsd)
                
                x = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
                y = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
                
                xe = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
                ye = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
                
                # create my own matrices for the plotting with pcolor or pcolormesh
                xc = np.zeros_like(data)
                yc = np.zeros_like(data)
                xce = np.zeros_like(data)
                yce = np.zeros_like(data)
                
                for i in range(ncolsnrows):
                    xc[:,i] = x
                    xce[:,i] = xe
                    
                for i in range(ncolsnrows):
                    yc[i,:] = y
                    yce[i,:] = ye
                
                # This need to be changed!
                
                ncenterx = np.int(np.round(xccrater[indf]/cellsize))
                ncentery = np.int(np.round(yccrater[indf]/cellsize))
                
                #other places where ncenterx and ncentery is required 
                x1, y1 = wk.xy_circle(1.0*r0_new[indf], xe[ncenterx], ye[ncentery]) # xe and ye are equals
        
                
                # First detrending (a plane is fitted through the elevations taken at circles
                # 2.0 and 3.0R from the center of the crater )
                
                stduse = True #use standard deviation removal
                ndata = wk.detrending(xc, yc, r0_new[indf], cellsize, ncenterx, ncentery, data, stduse) # I think it's okay to use xc, yc here
                                      
                
                # the Maximum elevation are selected a first time
                first_run = True
                mingrade = 0.05
                minclust = 0.05
                slen = 0.1
                
                (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME) = (
                 wk.rim(xc, yc, ncenterx, ncentery, ndata, r0_new[indf], 
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
                
                #wk.detrending_rim
                
                # the detrended plane is substracted to the DEM
                ndata2 = ndata - Z2
                
                # the first detrending has ben run so we change now the first_run flag to false
                first_run = False
                
                # second run of the function rim
                (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME,
                 col_coord_LE, row_coord_LE, col_cells_LE, row_cells_LE, elev_LE, prof_LE,
                 col_coord_BS, row_coord_BS, col_cells_BS, row_cells_BS, elev_BS, prof_BS) = (
                 wk.rim(xc, yc, ncenterx, ncentery, ndata2, r0_new[indf], 
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
                Drad = 0.1 * r0_new[indf]
                
                # Distance of interest (searching distance)
                Dint = 0.05 * r0_new[indf]
                
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
                name_crater_png = filename + 'XY.png'
                
                #xarcout = np.append(xarcout, np.array(OptRims[c][0,:]) + xllcorner)
                #yarcout = np.append(yarcout, np.array(OptRims[c][1,:]) + yllcorner)
                
                plt.ioff()
                plt.figure()
                plt.pcolormesh(xc, yc, ndata2)
                plt.title(filename,fontsize=16)
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
                
                #fit a circle
                xnewcenter, ynewcenter, rnew, residu = wk.leastsq_circle(xarcgis,yarcgis)
                                           
                #other places where ncenterx and ncentery is required (new calculation)
                x1n, y1n = wk.xy_circle(1.0*rnew, xnewcenter, ynewcenter)
                plt.plot(x1n,y1n,"r")
                
                
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
                
        
        # in the case where centre_craters does not exist or flagcrater == 1 or 2
        # run the whole thing
        else: 
        
            # should be good this way
            data = pd.read_csv(path + filename + ".asc", delimiter=' ', skiprows=6).values
            data2 = np.rot90(data)
            data3 = np.rot90(data2)
            data4 = (np.rot90(data3)) * scaling_factor #scaling factor # only for sldem
            
            del data
            del data2
            del data3
            data = copy.deepcopy(data4) #multiplying factor of the elevation (not for kaguya, 0.5 for Kaguya + LOLA)
            del data4
            
            #correspond to orgin in the lower-left corner and data[x,y] 
            # with x = ncols, y = nrows !!! important
                    
            # the number of columns, rows and extent of the raster is read
            (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path,filename + ".asc")
    
            # take the minimum of ncols, nrows and ncols and nrows of data
            nrowsd, ncolsd = np.shape(data)
            ncolsnrows = np.min([ncols, nrows, ncolsd, nrowsd])
            data = data[:ncolsnrows, :ncolsnrows]
            
            print (ncols, nrows, ncolsd, nrowsd)
            
            x = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
            y = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
            
            xe = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
            ye = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
            
            # create my own matrices for the plotting with pcolor or pcolormesh
            xc = np.zeros_like(data)
            yc = np.zeros_like(data)
            xce = np.zeros_like(data)
            yce = np.zeros_like(data)
            
            for i in range(ncolsnrows):
                xc[:,i] = x
                xce[:,i] = xe
                
            for i in range(ncolsnrows):
                yc[i,:] = y
                yce[i,:] = ye
            
            # This need to be changed!
            
            # middle of the dtm (this has to change!)
            ncentery = int(nrows/2) # changed
            ncenterx = int(ncols/2) # changed
                    
            #other places where ncenterx and ncentery is required 
            x1, y1 = wk.xy_circle(1.0*r0_new[indf], xe[ncenterx], ye[ncentery]) # xe and ye are equals
    
            
            # First detrending (a plane is fitted through the elevations taken at circles
            # 2.0 and 3.0R from the center of the crater )
            
            stduse = True #use standard deviation removal
            ndata = wk.detrending(xc, yc, r0_new[indf], cellsize, ncenterx, ncentery, data, stduse) # I think it's okay to use xc, yc here
                                  
            
            # the Maximum elevation are selected a first time
            first_run = True
            mingrade = 0.05
            minclust = 0.05
            slen = 0.1
            
            (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME) = (
             wk.rim(xc, yc, ncenterx, ncentery, ndata, r0_new[indf], 
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
            
            #wk.detrending_rim
            
            # the detrended plane is substracted to the DEM
            ndata2 = ndata - Z2
            
            # the first detrending has ben run so we change now the first_run flag to false
            first_run = False
            
            # second run of the function rim
            (col_coord_ME, row_coord_ME, col_cells_ME, row_cells_ME, elev_ME, prof_ME,
             col_coord_LE, row_coord_LE, col_cells_LE, row_cells_LE, elev_LE, prof_LE,
             col_coord_BS, row_coord_BS, col_cells_BS, row_cells_BS, elev_BS, prof_BS) = (
             wk.rim(xc, yc, ncenterx, ncentery, ndata2, r0_new[indf], 
                    cellsize, slen, minclust, mingrade, first_run))
            
            
            # takes only nonnan values from *BS data and merge it to local elevation data
            ixfinite = np.where(np.isfinite(col_coord_BS) == True)
            col_coord_BS = col_coord_BS[ixfinite] 
            row_coord_BS = row_coord_BS[ixfinite] 
            col_cells_BS = col_cells_BS[ixfinite] 
            row_cells_BS = row_cells_BS[ixfinite] 
            elev_BS = elev_BS[ixfinite] 
            prof_BS = prof_BS[ixfinite] 
            
            use_break_in_slope = True # changed
            
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
            Drad = 0.1 * r0_new[indf]
            
            # Distance of interest (searching distance)
            Dint = 0.05 * r0_new[indf]
            
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
            name_crater_png = filename + 'XY.png'
            
            #xarcout = np.append(xarcout, np.array(OptRims[c][0,:]) + xllcorner)
            #yarcout = np.append(yarcout, np.array(OptRims[c][1,:]) + yllcorner)
            
            plt.ioff()
            plt.figure()
            plt.pcolormesh(xc, yc, ndata2)
            plt.title(filename,fontsize=16)
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
            
            #fit a circle
            xnewcenter, ynewcenter, rnew, residu = wk.leastsq_circle(xarcgis,yarcgis)
                                       
            #other places where ncenterx and ncentery is required (new calculation)
            x1n, y1n = wk.xy_circle(1.0*rnew, xnewcenter, ynewcenter)
            plt.plot(x1n,y1n,"r")
            
            
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

def run2(path, path8R, filenameXY, filenamecrater, pathdata, scaling_factor):
    
    '''
    filenameXY = 'data.txt'
    filenamecrater = 'crater_id.txt'
    '''        
    os.chdir(path)
    
  
    crater_id = np.genfromtxt(filenamecrater,skip_header=1,dtype=str)
 
  
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
    
    unc_hr = np.ones(len(crater_id))
    med_hr = np.ones(len(crater_id))
    hr_25  = np.ones(len(crater_id))
    hr_75 = np.ones(len(crater_id))
    hr_min  = np.ones(len(crater_id))
    hr_max  = np.ones(len(crater_id))   
    
    med_depth = np.ones(len(crater_id))
    depthf = np.ones(len(crater_id))
    depth_25  = np.ones(len(crater_id))
    depth_75 = np.ones(len(crater_id))
    depth_min  = np.ones(len(crater_id))
    depth_max  = np.ones(len(crater_id)) 
    
    diamf = np.ones(len(crater_id))
    med_diam = np.ones(len(crater_id))
    diam_25  = np.ones(len(crater_id))
    diam_75 = np.ones(len(crater_id))
    diam_min  = np.ones(len(crater_id))
    diam_max  = np.ones(len(crater_id))

    unc_crdl = np.ones(len(crater_id))
    med_crdl = np.ones(len(crater_id))
    crdl_25 = np.ones(len(crater_id))
    crdl_75 = np.ones(len(crater_id))
    crdl_min = np.ones(len(crater_id))
    crdl_max = np.ones(len(crater_id))
    
    unc_frdl = np.ones(len(crater_id))
    med_frdl = np.ones(len(crater_id))
    frdl_25 = np.ones(len(crater_id))
    frdl_75 = np.ones(len(crater_id))
    frdl_min = np.ones(len(crater_id))
    frdl_max = np.ones(len(crater_id))
    
    unc_rupcw = np.ones(len(crater_id))
    med_rupcw = np.ones(len(crater_id))
    rupcw_25 = np.ones(len(crater_id))
    rupcw_75 = np.ones(len(crater_id))
    rupcw_min = np.ones(len(crater_id))
    rupcw_max = np.ones(len(crater_id))
    
    unc_rufrc = np.ones(len(crater_id))
    med_rufrc = np.ones(len(crater_id))
    rufrc_25 = np.ones(len(crater_id))
    rufrc_75 = np.ones(len(crater_id))
    rufrc_min = np.ones(len(crater_id))
    rufrc_max = np.ones(len(crater_id))
    
    unc_lrs = np.ones(len(crater_id))
    med_lrs = np.ones(len(crater_id))
    lrs_25 = np.ones(len(crater_id))
    lrs_75 = np.ones(len(crater_id))
    lrs_min = np.ones(len(crater_id))
    lrs_max = np.ones(len(crater_id))
    
    unc_urs = np.ones(len(crater_id))
    med_urs = np.ones(len(crater_id))
    urs_25 = np.ones(len(crater_id))
    urs_75 = np.ones(len(crater_id))
    urs_min = np.ones(len(crater_id))
    urs_max = np.ones(len(crater_id))
    
    unc_fsa = np.ones(len(crater_id))
    med_fsa = np.ones(len(crater_id))
    fsa_25 = np.ones(len(crater_id))
    fsa_75 = np.ones(len(crater_id))
    fsa_min = np.ones(len(crater_id))
    fsa_max = np.ones(len(crater_id))
    
    
    # sort the files based on the size so the smallest files will be taken first    
    #size_of_asciis = np.zeros(len(crater_id))
    
    #for indf, filename in enumerate(crater_id):
    #    
    #    size_of_asciis[indf] = os.path.getsize(path + filename + ".asc")
        
    # sort from the smallest to the largest
    #idx_sort = np.argsort(size_of_asciis)
    #crater_id_new = crater_id[idx_sort]
    
    crater_id_new = crater_id
       
    for indf, filename in enumerate(crater_id_new):
        
        print (crater_id_new[indf], indf)
        
        flag8R = 0
        
        if os.path.isfile(pathdata + filename + '_res.txt'):
            None
        
        # added this part in case, we did not get any data for a crater
        elif not os.path.isfile(pathdata + filename + 'XY.txt'):
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
            
            unc_hr[indf] = np.nan
            med_hr[indf] = np.nan
            hr_25[indf] = np.nan
            hr_75[indf] = np.nan
            hr_min[indf] = np.nan
            hr_max[indf] = np.nan
            
            depthf[indf] = np.nan
            med_depth[indf] = np.nan
            depth_25[indf] = np.nan
            depth_75[indf] = np.nan
            depth_min[indf] = np.nan
            depth_max[indf] = np.nan
            
            diamf[indf] = np.nan
            med_diam[indf] = np.nan
            diam_25[indf] = np.nan
            diam_75[indf] = np.nan
            diam_min[indf] = np.nan
            diam_max[indf] = np.nan
            
            ##
            unc_rupcw[indf] = np.nan
            med_rupcw[indf] = np.nan
            rupcw_25[indf] = np.nan
            rupcw_75[indf] = np.nan
            rupcw_min[indf] = np.nan
            rupcw_max[indf] = np.nan
            
            unc_rufrc[indf] = np.nan
            med_rufrc[indf] = np.nan
            rufrc_25[indf] = np.nan
            rufrc_75[indf] = np.nan
            rufrc_min[indf] = np.nan
            rufrc_max[indf] = np.nan
            
            unc_lrs[indf] = np.nan
            med_lrs[indf] = np.nan
            lrs_25[indf] = np.nan
            lrs_75[indf] = np.nan
            lrs_min[indf] = np.nan
            lrs_max[indf] = np.nan
            
            unc_urs[indf] = np.nan
            med_urs[indf] = np.nan
            urs_25[indf] = np.nan
            urs_75[indf] = np.nan
            urs_min[indf] = np.nan
            urs_max[indf] = np.nan
            
            unc_fsa[indf] = np.nan
            med_fsa[indf] = np.nan
            fsa_25[indf] = np.nan
            fsa_75[indf] = np.nan
            fsa_min[indf] = np.nan
            fsa_max[indf] = np.nan
            
            unc_crdl[indf] = np.nan
            med_crdl[indf] = np.nan
            crdl_25[indf] = np.nan
            crdl_75[indf] = np.nan
            crdl_min[indf] = np.nan
            crdl_max[indf] = np.nan
            
            unc_frdl[indf] = np.nan
            med_frdl[indf] = np.nan
            frdl_25[indf] = np.nan
            frdl_75[indf] = np.nan
            frdl_min[indf] = np.nan
            frdl_max[indf] = np.nan  
            
        else:
                  
            # txt name
            name_crater_txt = filename + 'XY.txt'
            
            #load the data
            datagis = loadgis(pathdata, name_crater_txt)
            xarcgis = datagis[:,0]
            yarcgis = datagis[:,1]
            zarcgis = datagis[:,2]
            profgis = datagis[:,3]
            flaggis = datagis[:,4]
            
            # new get rid of nan values
            xarcgis =  xarcgis[~np.isnan(xarcgis)]
            yarcgis =  yarcgis[~np.isnan(yarcgis)]
            zarcgis =  zarcgis[~np.isnan(zarcgis)]
            profgis =  profgis[~np.isnan(profgis)]
            flaggis =  flaggis[~np.isnan(flaggis)]
                   
    
            
            
            x_not_outliers = xarcgis 
            y_not_outliers = yarcgis
            z_not_outliers = zarcgis
            prof_not_outliers = profgis
            
            # old center
            # get the cell the closest to the center of the crater
            #__, ncenterx_old = wk.find_nearest(xllcorner+xe, xcrater[indf]) #correspond to columns and x
            
            # because yllcorner is from the top left now
            #__, ncentery_old = wk.find_nearest(yllcorner+ye,ycrater[indf]) #correspond to rows and y
            #ncenterx_old = nrows/2
            #ncentery_old = nrows/2
            
            #fit a circle
            xnewcenter, ynewcenter, rnew, residu = wk.leastsq_circle(x_not_outliers,y_not_outliers)
            
            
            '''
            Does the ascii cover 3 times the radius of the circle (from the center of the crater)?
            '''
            
            # the number of columns, rows and extent of the raster is read
            (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path, filename + ".asc")        
            
            # check if the ascii fullfill this criteria along the x- and y-axis
            alongx  = xnewcenter + 3.0*(rnew)
            alongy = ynewcenter + 3.0*(rnew)
            
            
            # load the data for 4R anyway
            data = pd.read_csv(path + filename + ".asc", delimiter=' ', skiprows=6).values
            nrowsd, ncolsd = np.shape(data)
            ncolsnrows = np.min([ncols, nrows, ncolsd, nrowsd])
            print (ncols, nrows, ncolsd, nrowsd)
            data = data[:ncolsnrows, :ncolsnrows]
            
            if ((alongx >  (ncols * cellsize)) | (alongy >  (nrows * cellsize))):
                
                # load 8R
                
                ncolsnrows_old = ncolsnrows
                
                data = pd.read_csv(path8R + filename + ".asc", delimiter=' ', skiprows=6).values
            
                # the number of columns, rows and extent of the raster is read
                (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(path8R, filename + ".asc")
                
                nrowsd, ncolsd = np.shape(data)
                ncolsnrows = np.min([ncols, nrows, ncolsd, nrowsd])
                print (ncols, nrows, ncolsd, nrowsd)
                data = data[:ncolsnrows, :ncolsnrows]
                flag8R = 1
                
            else:
                # if not load 4R (all okay)
                None
    
            
            
            # load data
            data2 = np.rot90(data)
            data3 = np.rot90(data2)
            data4 = np.rot90(data3) * scaling_factor
            
            del data
            del data2
            del data3
            data = copy.deepcopy(data4) # multiplying by scaling factor (*0.50 for Kaguya+LOLA)
            del data4
            
            # here we have to do some trigonometry
            if flag8R == 1:
                
                print ("8R ascii is used")
                
                offsetorigin = np.int((ncolsnrows / 8) * 2)
                offsetorigin2 = ncolsnrows - (offsetorigin + ncolsnrows_old)
                
                x = np.linspace(-offsetorigin*cellsize,(ncolsnrows_old + offsetorigin2 -1)*cellsize,ncolsnrows)
                y = np.linspace(-offsetorigin*cellsize,(ncolsnrows_old + offsetorigin2 -1)*cellsize,ncolsnrows)
                
                xe = np.linspace((-offsetorigin*cellsize) + cellsize/2.0,((ncolsnrows_old + offsetorigin2 -1)*cellsize) + cellsize/2.0,ncolsnrows)
                ye = np.linspace((-offsetorigin*cellsize) + cellsize/2.0,((ncolsnrows_old + offsetorigin2 -1)*cellsize) + cellsize/2.0,ncolsnrows)
            
            else:
                x = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
                y = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
            
                xe = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
                ye = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
            
            # create my own matrices for the plotting with pcolor or pcolormesh
            xc = np.zeros_like(data)
            yc = np.zeros_like(data)
            xce = np.zeros_like(data)
            yce = np.zeros_like(data)
            
            for i in range(ncolsnrows):
                xc[:,i] = x
                xce[:,i] = xe
                
            for i in range(ncolsnrows):
                yc[i,:] = y
                yce[i,:] = ye        
            
            # get the cell the closest to the center of the crater (new calculation)
            # (not good, use xnewcenter and ynewcenter)
            __, ncenterx = wk.find_nearest(xe,xnewcenter)
            __, ncentery = wk.find_nearest(ye,ynewcenter)
            
            # new
            #ncenterx = np.int(np.round(xnewcenter,decimals=0))
            #ncentery = np.int(np.round(ynewcenter,decimals=0))
                       
            #other places where ncenterx and ncentery is required (new calculation)
            #x1, y1 = wk.xy_circle(1.0*rnew, xe[ncentery], ye[ncenterx])
    
            
            # First detrending (a plane is fitted through the elevations taken at circles
            # with new center of the circle
            
            try:
                stduse = True
                
                # if we have data for 2R and 3R: no problem!
                ndata = wk.detrending(xc, yc, rnew, cellsize, ncenterx, ncentery, data, stduse) #detrending
                    
                # the second detrending through the selected maximum elevation points are done
                stduse = False
                
                ndata2 = wk.detrending_rim(xc, yc, rnew, cellsize, ncenterx, ncentery, ndata, stduse)
            
            
                                
                #need to put some indices there (maybe, they will have different sizes)
                (R_upcw, R_ufrc, cse, slope_mcw, slope_ucw, slope_fsa, slope_lrs, slope_urs,  crdl, frdl,
                 h, hr, depth, diam, nnn, prf, crossSections, YSections, XSections) = (wk.calculation(xc, yc, x_not_outliers, 
                 y_not_outliers, z_not_outliers, prof_not_outliers,
                 rnew, ndata2, cellsize, ncenterx, ncentery, xllcorner, yllcorner))
                
    			# crossSections has been added
                
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
                
                unc_hr[indf] = np.nanstd(hr)/np.sqrt(nnn)
                med_hr[indf] = np.nanmedian(hr)
                hr_25[indf] = np.nanpercentile(hr,25)
                hr_75[indf] = np.nanpercentile(hr,75)
                hr_min[indf] = np.nanmin(hr)
                hr_max[indf] = np.nanmax(hr)
            
                depthf[indf] = np.nanstd(depth)/np.sqrt(nnn)
                med_depth[indf] = np.nanmedian(depth)
                depth_25[indf] = np.nanpercentile(depth,25)
                depth_75[indf] = np.nanpercentile(depth,75)
                depth_min[indf] = np.nanmin(depth)
                depth_max[indf] = np.nanmax(depth)
                
                diamf[indf] = np.nanstd(diam)/np.sqrt(nnn)
                med_diam[indf] = np.nanmedian(diam)
                diam_25[indf] = np.nanpercentile(diam,25)
                diam_75[indf] = np.nanpercentile(diam,75)
                diam_min[indf] = np.nanmin(diam)
                diam_max[indf] = np.nanmax(diam)
                
                ##
                unc_rupcw[indf] = np.nanstd(R_upcw)/np.sqrt(nnn)
                med_rupcw[indf] = np.nanmedian(R_upcw)
                rupcw_25[indf] = np.nanpercentile(R_upcw,25)
                rupcw_75[indf] = np.nanpercentile(R_upcw,75)
                rupcw_min[indf] = np.nanmin(R_upcw)
                rupcw_max[indf] = np.nanmax(R_upcw)
                
                unc_rufrc[indf] = np.nanstd(R_ufrc)/np.sqrt(nnn)
                med_rufrc[indf] = np.nanmedian(R_ufrc)
                rufrc_25[indf] = np.nanpercentile(R_ufrc,25)
                rufrc_75[indf] = np.nanpercentile(R_ufrc,75)
                rufrc_min[indf] = np.nanmin(R_ufrc)
                rufrc_max[indf] = np.nanmax(R_ufrc)
                
                unc_lrs[indf] = np.nanstd(slope_lrs)/np.sqrt(nnn)
                med_lrs[indf] = np.nanmedian(slope_lrs)
                lrs_25[indf] = np.nanpercentile(slope_lrs,25)
                lrs_75[indf] = np.nanpercentile(slope_lrs,75)
                lrs_min[indf] = np.nanmin(slope_lrs)
                lrs_max[indf] = np.nanmax(slope_lrs)
                
                unc_urs[indf] = np.nanstd(slope_urs)/np.sqrt(nnn)
                med_urs[indf] = np.nanmedian(slope_urs)
                urs_25[indf] = np.nanpercentile(slope_urs,25)
                urs_75[indf] = np.nanpercentile(slope_urs,75)
                urs_min[indf] = np.nanmin(slope_urs)
                urs_max[indf] = np.nanmax(slope_urs)
                
                unc_fsa[indf] = np.nanstd(slope_fsa)/np.sqrt(nnn)
                med_fsa[indf] = np.nanmedian(slope_fsa)
                fsa_25[indf] = np.nanpercentile(slope_fsa,25)
                fsa_75[indf] = np.nanpercentile(slope_fsa,75)
                fsa_min[indf] = np.nanmin(slope_fsa)
                fsa_max[indf] = np.nanmax(slope_fsa)
                
                unc_crdl[indf] = np.nanstd(crdl)/np.sqrt(nnn)
                med_crdl[indf] = np.nanmedian(crdl)
                crdl_25[indf] = np.nanpercentile(crdl,25)
                crdl_75[indf] = np.nanpercentile(crdl,75)
                crdl_min[indf] = np.nanmin(crdl)
                crdl_max[indf] = np.nanmax(crdl)
                
                unc_frdl[indf] = np.nanstd(frdl)/np.sqrt(nnn)
                med_frdl[indf] = np.nanmedian(frdl)
                frdl_25[indf] = np.nanpercentile(frdl,25)
                frdl_75[indf] = np.nanpercentile(frdl,75)
                frdl_min[indf] = np.nanmin(frdl)
                frdl_max[indf] = np.nanmax(frdl)
                
                        
                header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
                              'mdepth;udepth;depth25;depth75;depth_min;depth_max;' + 
                              'mh;uh;h25;h75;hmin;hmax;' +
                              'mhr;uhr;hr25;hr75;hrmin;hrmax;' +
                              'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
                              'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
                              'mcse;ucse;cse_25;cse_75;cse_min;cse_max;'
                              'mrupcw;urupcw;rupcw_25;rupcw_75;rupcw_min;rupcw_max;'
                              'mrufrc;urufrc;rufrc_25;rufrc_75;rufrc_min;rufrc_max;'
                              'mlrs;ulrs;lrs_25;lrs_75;lrs_min;lrs_max;'
                              'murs;uurs;urs_25;urs_75;urs_min;urs_max;'
                              'mfsa;ufsa;fsa_25;fsa_75;fsa_min;fsa_max;'
                              'mcrdl;ucrdl;crdl_25;crdl_75;crdl_min;crdl_max;'
                              'mfrdl;ufrdl;frdl_25;frdl_75;frdl_min;frdl_max')
                
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
                
                unc_hr[indf] = np.nan
                med_hr[indf] = np.nan
                hr_25[indf] = np.nan
                hr_75[indf] = np.nan
                hr_min[indf] = np.nan
                hr_max[indf] = np.nan
                
                depthf[indf] = np.nan
                med_depth[indf] = np.nan
                depth_25[indf] = np.nan
                depth_75[indf] = np.nan
                depth_min[indf] = np.nan
                depth_max[indf] = np.nan
                
                diamf[indf] = np.nan
                med_diam[indf] = np.nan
                diam_25[indf] = np.nan
                diam_75[indf] = np.nan
                diam_min[indf] = np.nan
                diam_max[indf] = np.nan
                
                ##
                unc_rupcw[indf] = np.nan
                med_rupcw[indf] = np.nan
                rupcw_25[indf] = np.nan
                rupcw_75[indf] = np.nan
                rupcw_min[indf] = np.nan
                rupcw_max[indf] = np.nan
                
                unc_rufrc[indf] = np.nan
                med_rufrc[indf] = np.nan
                rufrc_25[indf] = np.nan
                rufrc_75[indf] = np.nan
                rufrc_min[indf] = np.nan
                rufrc_max[indf] = np.nan
                
                unc_lrs[indf] = np.nan
                med_lrs[indf] = np.nan
                lrs_25[indf] = np.nan
                lrs_75[indf] = np.nan
                lrs_min[indf] = np.nan
                lrs_max[indf] = np.nan
                
                unc_urs[indf] = np.nan
                med_urs[indf] = np.nan
                urs_25[indf] = np.nan
                urs_75[indf] = np.nan
                urs_min[indf] = np.nan
                urs_max[indf] = np.nan
                
                unc_fsa[indf] = np.nan
                med_fsa[indf] = np.nan
                fsa_25[indf] = np.nan
                fsa_75[indf] = np.nan
                fsa_min[indf] = np.nan
                fsa_max[indf] = np.nan
                
                unc_crdl[indf] = np.nan
                med_crdl[indf] = np.nan
                crdl_25[indf] = np.nan
                crdl_75[indf] = np.nan
                crdl_min[indf] = np.nan
                crdl_max[indf] = np.nan
                
                unc_frdl[indf] = np.nan
                med_frdl[indf] = np.nan
                frdl_25[indf] = np.nan
                frdl_75[indf] = np.nan
                frdl_min[indf] = np.nan
                frdl_max[indf] = np.nan            
                
                header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
                              'mdepth;udepth;depth25;depth75;depth_min;depth_max;' + 
                              'mh;uh;h25;h75;hmin;hmax;' +
                              'mhr;uhr;hr25;hr75;hrmin;hrmax;' +
                              'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
                              'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
                              'mcse;ucse;cse_25;cse_75;cse_min;cse_max;'
                              'mrupcw;urupcw;rupcw_25;rupcw_75;rupcw_min;rupcw_max;'
                              'mrufrc;urufrc;rufrc_25;rufrc_75;rufrc_min;rufrc_max;'
                              'mlrs;ulrs;lrs_25;lrs_75;lrs_min;lrs_max;'
                              'murs;uurs;urs_25;urs_75;urs_min;urs_max;'
                              'mfsa;ufsa;fsa_25;fsa_75;fsa_min;fsa_max;'
                              'mcrdl;ucrdl;crdl_25;crdl_75;crdl_min;crdl_max;'
                              'mfrdl;ufrdl;frdl_25;frdl_75;frdl_min;frdl_max')
                
                #if it fails cross sections does not work
                
            # txt name
            name_crater_f = filename + '_res.txt'
                
            arc = np.column_stack((med_diam[indf], diamf[indf], diam_25[indf], diam_75[indf], diam_min[indf], diam_max[indf],
                                   med_depth[indf], depthf[indf], depth_25[indf], depth_75[indf], depth_min[indf], depth_max[indf],
                                   med_h[indf], unc_h[indf], h_25[indf], h_75[indf], h_min[indf], h_max[indf],
                                   med_hr[indf], unc_hr[indf], hr_25[indf], hr_75[indf], hr_min[indf], hr_max[indf],
                                   med_mcw[indf], unc_mcw[indf], mcw_25[indf], mcw_75[indf], mcw_min[indf], mcw_max[indf],
                                   med_ucw[indf], unc_ucw[indf], ucw_25[indf], ucw_75[indf], ucw_min[indf], ucw_max[indf],
                                   med_cse[indf], unc_cse[indf], cse_25[indf], cse_75[indf], cse_min[indf], cse_max[indf],
                                   med_rupcw[indf], unc_rupcw[indf], rupcw_25[indf], rupcw_75[indf], rupcw_min[indf], rupcw_max[indf],
                                   med_rufrc[indf], unc_rufrc[indf], rufrc_25[indf], rufrc_75[indf], rufrc_min[indf], rufrc_max[indf],
                                   med_lrs[indf], unc_lrs[indf], lrs_25[indf], lrs_75[indf], lrs_min[indf], lrs_max[indf],
                                   med_urs[indf], unc_urs[indf], urs_25[indf], urs_75[indf], urs_min[indf], urs_max[indf],
                                   med_fsa[indf], unc_fsa[indf], fsa_25[indf], fsa_75[indf], fsa_min[indf], fsa_max[indf],
                                   med_crdl[indf], unc_crdl[indf], crdl_25[indf], crdl_75[indf], crdl_min[indf], crdl_max[indf],
                                   med_frdl[indf], unc_frdl[indf], frdl_25[indf], frdl_75[indf], frdl_min[indf], frdl_max[indf]))
            
            np.savetxt(pathdata + name_crater_f, arc, delimiter = ";", header=header_txt,fmt='%10.5f', comments='#')
    		
            # added part about saving cross sections. I still need to test if this good or not
            name_crater_crossX = filename + '_cross_sectionsX.txt'
            name_crater_crossY = filename + '_cross_sectionsY.txt'
            name_crater_crossZ = filename + '_cross_sectionsZ.txt'
            
            
            
            try:
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
                
                arrayY_crossSections = np.empty((maxlength_crossSections, number_crossSections))
                arrayY_crossSections[:] = np.nan  
                
                arrayX_crossSections = np.empty((maxlength_crossSections, number_crossSections))
                arrayX_crossSections[:] = np.nan           
                
                # fill array with values
                for xx_tmp, x_tmp in crossSections.items():
                    maxlength_crossSections = len(x_tmp)
                    array_crossSections[:maxlength_crossSections,xx_tmp] = x_tmp
                    arrayY_crossSections[:maxlength_crossSections,xx_tmp] = YSections[xx_tmp]
                    arrayX_crossSections[:maxlength_crossSections,xx_tmp] = XSections[xx_tmp]
                
                np.savetxt(pathdata + name_crater_crossX, arrayX_crossSections, delimiter = ";",fmt='%10.5f', comments='#')
                np.savetxt(pathdata + name_crater_crossY, arrayY_crossSections, delimiter = ";",fmt='%10.5f', comments='#')
                np.savetxt(pathdata + name_crater_crossZ, array_crossSections, delimiter = ";",fmt='%10.5f', comments='#')

            except:
                None                              
		
		
    return (med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            med_depth, depthf, depth_25, depth_75, depth_min, depth_max, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_hr, unc_hr, hr_25, hr_75, hr_min, hr_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max,
            med_rupcw, unc_rupcw, rupcw_25, rupcw_75, rupcw_min, rupcw_max,
            med_rufrc, unc_rufrc, rufrc_25, rufrc_75, rufrc_min, rufrc_max,
            med_lrs, unc_lrs, lrs_25, lrs_75, lrs_min, lrs_max,
            med_urs, unc_urs, urs_25, urs_75, urs_min, urs_max,
            med_fsa, unc_fsa, fsa_25, fsa_75, fsa_min, fsa_max,
            med_crdl, unc_crdl, crdl_25, crdl_75, crdl_min, crdl_max,
            med_frdl, unc_frdl, frdl_25, frdl_75, frdl_min, frdl_max)
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
    
    # be careful here for the first line
    crater_id = np.genfromtxt(filenamecrater,comments="#",dtype=str)
    
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
    
    unc_hr = np.ones(len(crater_id))
    med_hr = np.ones(len(crater_id))
    hr_25  = np.ones(len(crater_id))
    hr_75 = np.ones(len(crater_id))
    hr_min  = np.ones(len(crater_id))
    hr_max  = np.ones(len(crater_id))
    
    med_depth = np.ones(len(crater_id))
    depthf = np.ones(len(crater_id))
    depth_25  = np.ones(len(crater_id))
    depth_75 = np.ones(len(crater_id))
    depth_min  = np.ones(len(crater_id))
    depth_max  = np.ones(len(crater_id))  
    
    diamf = np.ones(len(crater_id))
    med_diam = np.ones(len(crater_id))
    diam_25  = np.ones(len(crater_id))
    diam_75 = np.ones(len(crater_id))
    diam_min  = np.ones(len(crater_id))
    diam_max  = np.ones(len(crater_id))
    
    med_rupcw  = np.ones(len(crater_id))
    unc_rupcw  = np.ones(len(crater_id))
    rupcw_25  = np.ones(len(crater_id))
    rupcw_75  = np.ones(len(crater_id))
    rupcw_min  = np.ones(len(crater_id))
    rupcw_max  = np.ones(len(crater_id))
    
    med_rufrc  = np.ones(len(crater_id))
    unc_rufrc  = np.ones(len(crater_id))
    rufrc_25  = np.ones(len(crater_id))
    rufrc_75  = np.ones(len(crater_id))
    rufrc_min  = np.ones(len(crater_id))
    rufrc_max  = np.ones(len(crater_id))
    
    med_lrs  = np.ones(len(crater_id))
    unc_lrs  = np.ones(len(crater_id))
    lrs_25  = np.ones(len(crater_id))
    lrs_75  = np.ones(len(crater_id))
    lrs_min  = np.ones(len(crater_id))
    lrs_max  = np.ones(len(crater_id))
    
    med_urs  = np.ones(len(crater_id))
    unc_urs  = np.ones(len(crater_id))
    urs_25  = np.ones(len(crater_id))
    urs_75  = np.ones(len(crater_id))
    urs_min  = np.ones(len(crater_id))
    urs_max  = np.ones(len(crater_id))
    
    med_fsa  = np.ones(len(crater_id))
    unc_fsa  = np.ones(len(crater_id))
    fsa_25  = np.ones(len(crater_id))
    fsa_75  = np.ones(len(crater_id))
    fsa_min  = np.ones(len(crater_id))
    fsa_max  = np.ones(len(crater_id))
    
    med_crdl  = np.ones(len(crater_id))
    unc_crdl  = np.ones(len(crater_id))
    crdl_25  = np.ones(len(crater_id))
    crdl_75  = np.ones(len(crater_id))
    crdl_min  = np.ones(len(crater_id))
    crdl_max  = np.ones(len(crater_id))
    
    med_frdl  = np.ones(len(crater_id))
    unc_frdl  = np.ones(len(crater_id))
    frdl_25  = np.ones(len(crater_id))
    frdl_75  = np.ones(len(crater_id))
    frdl_min  = np.ones(len(crater_id))
    frdl_max  = np.ones(len(crater_id))
    
    for indf, filename in enumerate(crater_id):
    
        name_crater_f = filename + '_res.txt'
        
        (med_diam[indf], diamf[indf], diam_25[indf], diam_75[indf], diam_min[indf], diam_max[indf], 
         med_depth[indf], depthf[indf], depth_25[indf], depth_75[indf], depth_min[indf], depth_max[indf], 
         med_h[indf], unc_h[indf], h_25[indf], h_75[indf], h_min[indf], h_max[indf],
         med_hr[indf], unc_hr[indf], hr_25[indf], hr_75[indf], hr_min[indf], hr_max[indf],
         med_mcw[indf], unc_mcw[indf], mcw_25[indf], mcw_75[indf], mcw_min[indf], mcw_max[indf],
         med_ucw[indf], unc_ucw[indf], ucw_25[indf], ucw_75[indf], ucw_min[indf], ucw_max[indf],
         med_cse[indf], unc_cse[indf], cse_25[indf], cse_75[indf], cse_min[indf], cse_max[indf],
         med_rupcw[indf], unc_rupcw[indf], rupcw_25[indf], rupcw_75[indf], rupcw_min[indf], rupcw_max[indf],
         med_rufrc[indf], unc_rufrc[indf], rufrc_25[indf], rufrc_75[indf], rufrc_min[indf], rufrc_max[indf],
         med_lrs[indf], unc_lrs[indf], lrs_25[indf], lrs_75[indf], lrs_min[indf], lrs_max[indf],
         med_urs[indf], unc_urs[indf], urs_25[indf], urs_75[indf], urs_min[indf], urs_max[indf],
         med_fsa[indf], unc_fsa[indf], fsa_25[indf], fsa_75[indf], fsa_min[indf], fsa_max[indf],
         med_crdl[indf], unc_crdl[indf], crdl_25[indf], crdl_75[indf], crdl_min[indf], crdl_max[indf],
         med_frdl[indf], unc_frdl[indf], frdl_25[indf], frdl_75[indf], frdl_min[indf], frdl_max[indf]) = np.loadtxt(pathdata + name_crater_f,delimiter=";", comments="#")
        
    return (med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            med_depth, depthf, depth_25, depth_75, depth_min, depth_max, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_hr, unc_hr, hr_25, hr_75, hr_min, hr_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max,
            med_rupcw, unc_rupcw, rupcw_25, rupcw_75, rupcw_min, rupcw_max,
            med_rufrc, unc_rufrc, rufrc_25, rufrc_75, rufrc_min, rufrc_max,
            med_lrs, unc_lrs, lrs_25, lrs_75, lrs_min, lrs_max,
            med_urs, unc_urs, urs_25, urs_75, urs_min, urs_max,
            med_fsa, unc_fsa, fsa_25, fsa_75, fsa_min, fsa_max,
            med_crdl, unc_crdl, crdl_25, crdl_75, crdl_min, crdl_max,
            med_frdl, unc_frdl, frdl_25, frdl_75, frdl_min, frdl_max)
      
'''
**************************************************************************
'''

# define the folder containing the DTMs, and the folders where the plot
# and data routines will save data
#path = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_ASCII/'
#pathplot = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_PLOT2D/'
#pathdata = 'X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_DATA2D/'

#path = 'D:/ANALYSIS/SIMPLECRATERS_MOON/may_2019/ascii/'
#pathplot = 'D:/ANALYSIS/SIMPLECRATERS_MOON/may_2019/double_detrending/plots/'
#pathdata = 'D:/ANALYSIS/SIMPLECRATERS_MOON/may_2019/double_detrending/data/'

#path = 'D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/ASCII/'
#path8R = 'D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/ASCII/'
#pathplot = 'D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/double_detrending_few_changes/plots/'
#pathdata = 'D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/double_detrending_few_changes/data/'

#path = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii4R/'
#path8R = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii8R/'
#pathplot = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/double_detrending_4R8Rlast/plots/'
#pathdata = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/double_detrending_4R8Rlast/data/'
#scaling_factor = 1.0

path = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_16R/'
path8R = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_16R/'
pathplot = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending_slope/plots/'
pathdata = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending_slope/data/'
scaling_factor = 1.0
centre_craters = 'centre_crater.csv'

#path = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_16R/'
#path8R = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_16R/'
#pathplot = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/plots/'
#pathdata = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/data/'
#scaling_factor = 1.0
#centre_craters = 'centre_crater.csv'

#path = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii4R/'
#path8R = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii/'

#pathplot = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/double_detrending/plots/'
#pathdata = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/double_detrending/data/'

#path = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii/'
#path8R = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii/'

#pathplot = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/plots/'
#pathdata = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/data/'

# Location + cratername
filenameXY = 'data.txt'
filenamecrater = 'crater_id.txt'

#make directory if not existing
ifnot_mkdir(pathplot)
ifnot_mkdir(pathdata)

# run step 1
#run1(path, filenameXY, filenamecrater, pathplot, pathdata, scaling_factor, centre_craters)
run1(path8R, filenameXY, filenamecrater, pathplot, pathdata, scaling_factor, centre_craters)


# run step 2
(med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            med_depth, depthf, depth_25, depth_75, depth_min, depth_max, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_hr, unc_hr, hr_25, hr_75, hr_min, hr_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max,
            med_rupcw, unc_rupcw, rupcw_25, rupcw_75, rupcw_min, rupcw_max,
            med_rufrc, unc_rufrc, rufrc_25, rufrc_75, rufrc_min, rufrc_max,
            med_lrs, unc_lrs, lrs_25, lrs_75, lrs_min, lrs_max,
            med_urs, unc_urs, urs_25, urs_75, urs_min, urs_max,
            med_fsa, unc_fsa, fsa_25, fsa_75, fsa_min, fsa_max,
            med_crdl, unc_crdl, crdl_25, crdl_75, crdl_min, crdl_max,
            med_frdl, unc_frdl, frdl_25, frdl_75, frdl_min, frdl_max) = run2(path, path8R, filenameXY, filenamecrater, pathdata, scaling_factor)


## load data 
os.chdir(pathdata)
filenamec = glob.glob("cra*_res.txt")

filenamecup = []

for f in filenamec:
    filenamecup.append(f.split('_res.txt')[0])
    
filenamecrater = path + 'crater_id_current.txt'

with open(filenamecrater, "w") as ff:
    for f in filenamecup:
        ff.writelines(f + "\n")

(med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            med_depth, depthf, depth_25, depth_75, depth_min, depth_max, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_hr, unc_hr, hr_25, hr_75, hr_min, hr_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max,
            med_rupcw, unc_rupcw, rupcw_25, rupcw_75, rupcw_min, rupcw_max,
            med_rufrc, unc_rufrc, rufrc_25, rufrc_75, rufrc_min, rufrc_max,
            med_lrs, unc_lrs, lrs_25, lrs_75, lrs_min, lrs_max,
            med_urs, unc_urs, urs_25, urs_75, urs_min, urs_max,
            med_fsa, unc_fsa, fsa_25, fsa_75, fsa_min, fsa_max,
            med_crdl, unc_crdl, crdl_25, crdl_75, crdl_min, crdl_max,
            med_frdl, unc_frdl, frdl_25, frdl_75, frdl_min, frdl_max) = mergedata(path,pathdata,filenamecrater)




## volume



header_txt = ('mdiam;udiam;diam25;diam75;diam_min;diam_max;' + 
                          'mdepth;udepth;depth25;depth75;depth_min;depth_max;' + 
                          'mh;uh;h25;h75;hmin;hmax;' +
                          'mhr;uhr;hr25;hr75;hrmin;hrmax;' +
                          'm_mcw;u_mcw;mcw_25;mcw_75;mcw_min;mcw_max;' +
                          'm_ucw;u_ucw;ucw_25;ucw_75;ucw_min;ucw_max;' +
                          'mcse;ucse;cse_25;cse_75;cse_min;cse_max;dD;'
                          'mrupcw;urupcw;rupcw_25;rupcw_75;rupcw_min;rupcw_max;'
                          'mrufrc;urufrc;rufrc_25;rufrc_75;rufrc_min;rufrc_max;'
                          'mlrs;ulrs;lrs_25;lrs_75;lrs_min;lrs_max;'
                          'murs;uurs;urs_25;urs_75;urs_min;urs_max;'
                          'mfsa;ufsa;fsa_25;fsa_75;fsa_min;fsa_max;'
                          'mcrdl;ucrdl;crdl_25;crdl_75;crdl_min;crdl_max;'
                          'mfrdl;ufrdl;frdl_25;frdl_75;frdl_min;frdl_max;vol')

dD = (med_h-med_depth) / med_diam

# calculate volume here
vol = cv.main(path, pathdata, filenamecrater, resolution_DEM=7.4031617)

name_crater_f = 'final_res2D_coldspots.txt'
            
arc = np.column_stack((med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            med_depth, depthf, depth_25, depth_75, depth_min, depth_max,  
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_hr, unc_hr, hr_25, hr_75, hr_min, hr_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max, dD,
            med_rupcw, unc_rupcw, rupcw_25, rupcw_75, rupcw_min, rupcw_max,
            med_rufrc, unc_rufrc, rufrc_25, rufrc_75, rufrc_min, rufrc_max,
            med_lrs, unc_lrs, lrs_25, lrs_75, lrs_min, lrs_max,
            med_urs, unc_urs, urs_25, urs_75, urs_min, urs_max,
            med_fsa, unc_fsa, fsa_25, fsa_75, fsa_min, fsa_max,
            med_crdl, unc_crdl, crdl_25, crdl_75, crdl_min, crdl_max,
            med_frdl, unc_frdl, frdl_25, frdl_75, frdl_min, frdl_max, vol)) # and volume somewhere
        
np.savetxt(pathdata + name_crater_f, arc, delimiter = ";", header=header_txt,fmt='%10.5f', comments='#') 
