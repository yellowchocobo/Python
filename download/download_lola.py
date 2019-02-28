# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 09:04:15 2017

@author: root
"""

#from ftplib import FTP
import os, os.path
import re
import urllib
import urllib2
import numpy as np

# main path
paths = '/var/tmp/mysshfs/bydn/projects/beyondearth/nilscp/lola/automatic/Tray2016/'
'''
path_down = '/work/nilscp/isaleruns/project/paper2/APR/DINNER_DOUTER/'
filename = 'non_valid_concentric.csv'
'''
# I need to load maxlat, minlat, westernlon, easternlon
path_down = '/home/nilscp/Downloads/'

data = np.genfromtxt(path_down + filename,delimiter=",",skip_header=1,dtype='str')

# get name
name = data[:,0]
del data

# get folders name
folders = []
for nn in name:
    nd = nn.replace(" ", "_")
    ne = nd.replace(".", "_")
    folders.append(ne)
    
# create folders

for ix, f in enumerate(folders): 
    path_f = paths + f
    
    if not os.path.exists(path_f):
        os.makedirs(path_f)
    
# load data
data = np.genfromtxt(path_down + filename,delimiter=",",skip_header=1)

maxlat = np.around(data[:,3],decimals=2)
minlat = np.around(data[:,4],decimals=2)
westernlon = np.around(data[:,5],decimals=2)
easternlon = np.around(data[:,6],decimals=2)

easternlon2 = np.around(data[:,2],decimals=2)


'''
# test
maxlat = np.around([30.05,29.34],decimals=2)
minlat = np.around([29.76,29.02],decimals=2)
westernlon = np.around([275.65,280.07],decimals=2)
easternlon = np.around([275.98,280.41],decimals=2)

folders = ['Unnamed_30','Unnamed_31']


'''


# default url
default_url1 = 'http://oderest.rsl.wustl.edu/livegds?query=lolardr&results=vsi&maxlat='
default_url2 = '&minlat='
default_url3 = '&westernlon='
default_url4 = '&easternlon='



# get the lola data url for all craters 
lst_url = []

for ix, val in np.ndenumerate(maxlat):
    
    tmp_url = (default_url1 + str(maxlat[ix]) + default_url2 + str(minlat[ix])
    + default_url3 + str(westernlon[ix]) + default_url4 + str(easternlon[ix]))

    lst_url.append(tmp_url)
    

# download the shape files for each crater
for ix, t_url in enumerate(lst_url):
    
    patht = paths + folders[ix] + '/'
    os.chdir(patht)

    response = urllib2.urlopen(t_url)
    html = response.read()

    data = html.split('</URL>')

    lst_shp = []

    word = '_shp'
    pattern = re.compile(word)

    for line in data:
    
        if (pattern.search(line)):
            tmp_line = line.split('<URL>')[1]
            lst_shp.append(tmp_line)
            
    for shpf in lst_shp:
        fname = shpf.split('/')[-1]
        urllib.urlretrieve(shpf,patht+fname)
        

    
### defining the geographic coordinate systems

def proj(LON,LAT):
    wkt = 'PROJCS["Equirectangular_MOON",GEOGCS["GCS_MOON",DATUM["D_MOON",SPHEROID["MOON_localRadius",1737400.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Equidistant_Cylindrical"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],'
    wkt1 = 'PARAMETER["central_meridian",' + str(LON) + '],'
    wkt2 = 'PARAMETER["standard_parallel_1",' + str(LAT) + '],'
    wkt3 = 'UNIT["Meter",1.0]]'
    wktf = wkt + wkt1 + wkt2 + wkt3
    return wktf

filename_tmp = "reproj_tmp.prj"

for ix, f in enumerate(folders):
    os.chdir(paths + f)
    
    if np.isnan(easternlon2[ix]):
        None
    else:
        wkt = proj(easternlon2[ix],maxlat[ix])	
        file = open(filename_tmp,"w")
        file.write(wkt)
        file.close()
        
'''
***********************************************************************
'''        

### playing with the data. only 1 and 2 (from the final results)
path = '/var/tmp/mysshfs/bydn/projects/beyondearth/nilscp/lola/data/'
filename = 'lunar_conc_10042017.csv'

data = np.genfromtxt(path + filename ,delimiter=",",skip_header=1,dtype='str')

# get name
name = data[:,0]
del data

# get folders name
folders = []
for nn in name:
    nd = nn.replace(" ", "_")
    ne = nd.replace(".", "_")
    folders.append(ne)

#get data
data = np.genfromtxt(path + filename ,delimiter=",",skip_header=1)

ORBIT = data[:,1]
CONC = data[:,2]

del data

'''
***********************************************************************
'''

def loadData3(path,filename):
    
    os.chdir(path)
    
    data = np.genfromtxt(filename,delimiter=",",skip_header=1)

    X = data[:,-2]
    Y = data[:,-1]
    Z = data[:,4]
    ORB = data[:,16]
    ID = data[:,14]
    
    return X, Y, Z, ORB, ID
    
'''
***********************************************************************
'''

import matplotlib.pyplot as plt


def plotConcProf3(drim,alt,slope,header_title,figtitle,paths,(drim_l,alt_l),(drim_r,alt_r),(dinner_l,inn_l),(dinner_r,inn_r)):
    
    '''
    would be nice to plot the values found from DINNER, DOUTER AND SO ON...
    '''
        
    fig = plt.figure(figsize=(9,5))
    plt.scatter(drim,alt,c=np.abs(slope),s=40,vmin=0,vmax=30)
    plt.xlim((0,np.max(drim)))
    plt.colorbar()
    plt.tick_params('both', labelsize=16,length=10, width=1., which='major')
    plt.xlabel('Distance (m)',fontsize=20)
    plt.ylabel('Elevation (m)',fontsize=20)
    dr = str(np.around((dinner_l-dinner_r)/(drim_l-drim_r),decimals=2))
    plt.title(header_title + ' ' + dr,color='k',fontsize=18)
    plt.plot(drim_l,alt_l,"ko",ms=10)
    plt.plot(drim_r,alt_r,"ko",ms=10)
    plt.plot(dinner_l,inn_l,"ko",ms=10)
    plt.plot(dinner_r,inn_r,"ko",ms=10)
    fig.tight_layout()    
    fig.subplots_adjust(top=0.9)    
    fig.savefig(paths+figtitle,dpi=300)
    plt.legend()
    #plt.close()
    
'''
***********************************************************************
'''

import matplotlib.pyplot as plt


def plotConcProf2(drim,alt,slope,header_title,figtitle,paths):
    
    '''
    would be nice to plot the values found from DINNER, DOUTER AND SO ON...
    '''
        
    fig = plt.figure(figsize=(9,5))
    plt.scatter(drim,alt,c=np.abs(slope),s=40,vmin=0,vmax=30)
    plt.xlim((0,np.max(drim)))
    plt.colorbar()
    plt.tick_params('both', labelsize=16,length=10, width=1., which='major')
    plt.xlabel('Distance (m)',fontsize=20)
    plt.ylabel('Elevation (m)',fontsize=20)
    plt.title(header_title,color='k',fontsize=18)
    fig.tight_layout()    
    fig.subplots_adjust(top=0.9)    
    fig.savefig(paths+figtitle,dpi=300)
    plt.legend()
    plt.close()    
'''
***********************************************************************
'''
# We need to calculate DOUTER, DINNER, HINNER, HOUTER

'''
***********************************************************************
'''
def slope_crater(Y,X,window):

    def movingaverage(values,window):
        
        weights = np.repeat(1.0,window)/window
        ma = np.convolve(values,weights,'valid')
        return ma
    
    altma = movingaverage(Y,window)
    dma = movingaverage(X,window)
    
    alt2 = altma[1:]
    alt1 = altma[:-1]
    drim2 = dma[1:]
    drim1 = dma[:-1]
    
    slope = np.arctan((alt2-alt1)/(drim2-drim1)) * (180. / np.pi)
        
    return (drim2+drim1)/2.,(alt2+alt1)/2., slope

'''
***********************************************************************
''' 

def findplateau(X,Y,Z):
    
    X0 = X[np.argmax(Y)]
    Y0 = Y[np.argmax(Y)]
        
    XD = X - X0
    YD = Y - Y0
    D = np.sqrt((XD**2.)+(YD**2.))
        
    alt2 = Z[1:]
    alt1 = Z[:-1]
    drim2 = D[1:]
    drim1 = D[:-1]
        
    alt = (alt1 + alt2) /2.
    drim = (drim1 + drim2) /2.
        
    slope = np.arctan((alt2-alt1)/(drim2-drim1)) * (180. / np.pi)
        
    argmi = np.argmin(alt)            
        
    mindep = alt[argmi]
    minx = drim[argmi]
        
    argma1 = np.argmax(alt[:argmi])            
    argma2 = np.argmax(alt[argmi:])
        
    # drim, alt of the northern rim
    rim1x = drim[argma1]
    rim1y = alt[argma1]
        
    # drim, alt of the southern rim
    rim2x = drim[argmi + argma2]
    rim2y = alt[argmi + argma2]
        
    ixa = np.argmax(slope[argma1+5:argmi-5])
    ixb = np.argmin(slope[argmi+5:argmi+argma2-5])
        
    ixfa = argma1+ixa+5
    ixfb = argmi + 5 + ixb
    
    inn1x = drim[ixfa]
    inn2x = drim[ixfb]
    
    inn1y = alt[ixfa]
    inn2y = alt[ixfb]
    #drim1 = rim2x - rim1x
    #dinner = drim[ixfa] - drim[ixfb]
        
    return drim, alt, slope, (rim1x,rim2x), (rim1y,rim2y), (inn1x,inn2x), (inn1y,inn2y)


'''
***********************************************************************
'''

import glob

path1 = '/var/tmp/mysshfs/bydn/projects/beyondearth/nilscp/lola/automatic/Wood1978/'
path2 = '/var/tmp/mysshfs/bydn/projects/beyondearth/nilscp/lola/automatic/Tray2016/'

paths = '/var/tmp/mysshfs/bydn/projects/beyondearth/nilscp/lola/data/plots_non_conc/'

#abc = (len(range(1,6)),len(folders))

# definition of variables
#DINNER_L = np.ones(abc)
#DOUTER_L = np.ones(abc)
#HINNER_L = np.ones(abc)
#HOUTER_L = np.ones(abc)
#
#DINNER_R = np.ones(abc)
#DOUTER_R = np.ones(abc)
#HINNER_R = np.ones(abc)
#HOUTER_R = np.ones(abc)
#
#DINNER = np.ones(abc)
#DOUTER = np.ones(abc)
#HINNER = np.ones(abc)
#HOUTER = np.ones(abc)

#n = np.ones(abc)

#DRIM = dict()
#ALT = dict()
#SLOPE = dict()


## if only want to plot the best concentric


for ix, f in enumerate(folders):
      
    name = f
    
    if os.path.exists(path1 + f + '/'):
        pathn = path1 + f + '/'
        os.chdir(pathn)
        filenames = glob.glob('*_reproj.csv')
        for filename in filenames:
            if filename.endswith('N_reproj.csv'):
                X, Y, Z, ORB, ID = loadData3(pathn,filename)
                print filename
            elif filename.endswith('S_reproj.csv'):
                X, Y, Z, ORB, ID = loadData3(pathn,filename)
                print filename
    else:
        pathn = path2 + f + '/'
        os.chdir(pathn)
        filenames = glob.glob('*_reproj.csv')
        X, Y, Z, ORB, ID = loadData3(pathn,filenames[0])
        print filenames[0]

                   
    for i in range(1,6):

        idx = np.where((ID == i) & (ORB == ORBIT[ix]))
        
        if (idx[0].size == 0):
            None
            #DINNER_L[i-1,ix] = DOUTER_L[i-1,ix] = HINNER_L[i-1,ix] = HOUTER_L[i-1,ix] = n[i-1,ix] = np.nan
            #DINNER_R[i-1,ix] = DOUTER_R[i-1,ix] = HINNER_R[i-1,ix] = HOUTER_R[i-1,ix] = np.nan
            #DINNER[i-1,ix] = DOUTER[i-1,ix] = HINNER[i-1,ix] = HOUTER[i-1,ix] = np.nan
        elif ((idx[0].size == 1) & (len(idx[0]) < 20)):
            None
            #DINNER_L[i-1,ix] = DOUTER_L[i-1,ix] = HINNER_L[i-1,ix] = HOUTER_L[i-1,ix] = n[i-1,ix] = np.nan
            #DINNER_R[i-1,ix] = DOUTER_R[i-1,ix] = HINNER_R[i-1,ix] = HOUTER_R[i-1,ix] = np.nan
            #DINNER[i-1,ix] = DOUTER[i-1,ix] = HINNER[i-1,ix] = HOUTER[i-1,ix] = np.nan
        else:
            XS = X[idx]
            YS = Y[idx]
            ZS = Z[idx]
            
            X0 = XS[np.argmax(YS)]
            Y0 = YS[np.argmax(YS)]
        
            XD = XS - X0
            YD = YS - Y0
            D = np.sqrt((XD**2.)+(YD**2.))
        
            alt2 = ZS[1:]
            alt1 = ZS[:-1]
            drim2 = D[1:]
            drim1 = D[:-1]
        
            alt = (alt1 + alt2) /2.
            drim = (drim1 + drim2) /2.
        
            slope = np.arctan((alt2-alt1)/(drim2-drim1)) * (180. / np.pi)
            
            #n[i-1,ix] = len(XS)
            #drim, alt, slope, (DOUTER_L[i-1,ix],DOUTER_R[i-1,ix]), (HOUTER_L[i-1,ix],HOUTER_R[i-1,ix]), (DINNER_L[i-1,ix],DINNER_R[i-1,ix]), (HINNER_L[i-1,ix],HINNER_R[i-1,ix]) = findplateau(XS,YS,ZS)
            
            #DINNER[i-1,ix] = DINNER_R[i-1,ix] - DINNER_L[i-1,ix]
            #DOUTER[i-1,ix] = DOUTER_R[i-1,ix] - DOUTER_L[i-1,ix]
            
            #yya = (DOUTER_L[i-1,ix],HOUTER_L[i-1,ix])
            #yyb = (DOUTER_R[i-1,ix],HOUTER_R[i-1,ix])
            #yyc = (DINNER_L[i-1,ix],HINNER_L[i-1,ix])
            #yyd = (DINNER_R[i-1,ix],HINNER_R[i-1,ix])       
            
            plotConcProf2(drim,alt,slope,name + ' #' + str(np.int(ORBIT[ix])) + ' #' + str(i),name+'_' + str(np.int(ORBIT[ix])) + '_' + str(i) + '.png',paths)
            
            
############
            
def flush(aaa):
    for i in range(aaa):
        plt.close()

            
            
            
ix = 86 
f = folders[ix]
name = f
           
if os.path.exists(path1 + f + '/'):
    pathn = path1 + f + '/'
    os.chdir(pathn)
    filenames = glob.glob('*_reproj.csv')
    for filename in filenames:
        if filename.endswith('N_reproj.csv'):
            X, Y, Z, ORB, ID = loadData3(pathn,filename)
            print filename
        elif filename.endswith('S_reproj.csv'):
            X, Y, Z, ORB, ID = loadData3(pathn,filename)
            print filename
else:
    pathn = path2 + f + '/'
    os.chdir(pathn)
    filenames = glob.glob('*_reproj.csv')
    X, Y, Z, ORB, ID = loadData3(pathn,filenames[0])
    print filenames[0]

               
for i in range(1,6):
    
    idx = np.where((ID == i) & (ORB == ORBIT[ix]))
    
    if (idx[0].size == 0):
        None
    elif ((idx[0].size == 1) & (len(idx[0]) < 20)):
        None
    else:
        XS = X[idx]
        YS = Y[idx]
        ZS = Z[idx]
        
        drim, alt, slope, (DOUTER_L,DOUTER_R), (HOUTER_L,HOUTER_R), (DINNER_L,DINNER_R), (HINNER_L,HINNER_R) = findplateau(XS,YS,ZS)
        yya = (DOUTER_L,HOUTER_L)
        yyb = (DOUTER_R,HOUTER_R)
        yyc = (DINNER_L,HINNER_L)
        yyd = (DINNER_R,HINNER_R)
        
        n = len(XS)
        print n
        print DOUTER_L, DOUTER_R, HOUTER_L, HOUTER_R
        plotConcProf3(drim,alt,slope,name + ' #' + str(np.int(ORBIT[ix])) + ' #' + str(i),
                      name+'_' + str(np.int(ORBIT[ix])) + '_' + str(i) + '.png',paths, yya, yyb, yyc, yyd)