import os
import numpy as np

'''
***********************************************************************
'''

path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/UAVG/SANDCOH/'
paths = '/global/work/nilscp/iSALE/benchmarkFI/UAVG/SANDCOH/'
NAME = 'UAVG_SANDCOH'

os.chdir(path)

'''
***********************************************************************
'''

U = 11380.                  #  Impact velocity in m/s
L1 = np.array([10, 20, 40, 60, 80])
L2 = np.arange(100,1001,50)
L = np.concatenate((L1,L2))  
        # Diameter of the projectile (m)
g = 1.62                    # gravity of the targeted body (SI)
rho = 2850. * 0.9           # target's density (kg/m3)
teta = 2850.                # projectile's density (kg/m3)
CPPR = 20

variable = np.array(["PATH","MODEL","TDUMP","GRIDH","GRIDV","GRIDSPC","DT","TEND","DTSAVE"])


GRIDH = np.array([650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650.])

GRIDV = np.array([650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650., 650.,
                  650., 650., 650., 650.])

'''
***********************************************************************
'''

#param = ["YDAM0", "FRICDAM","YLIMDAM", "ALPHA0"]
#arr = ["5.0D+1","0.6D+0","1.0D+9", "1.25D+0"]
MODEL = NAME
#INPUT_MATERIAL(param,arr,MODEL)

'''
***********************************************************************
'''

param = ["PATH","MODEL","TDUMP","GRIDH","GRIDV","GRIDSPC","LAYPOS","DT","TEND","DTSAVE"]
nvar = [1,1,1,2,2,1,1,1,1,1] # 
arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta,CPPR,paths,NAME, GRIDH, GRIDV)
INPUT_ASTEROID(param,arr,nvar)


'''
***********************************************************************
'''

ESTTIME = ['24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00',
           '24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00',
           '24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00',
           '24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00','24-00:00:00']

keyast = 'asteroid*L*'
keymat = 'mat*'
generateBashfilenames(path,ESTTIME, keyast, keymat)



'''
***********************************************************************
'''


# For several cases for Uavg

MODELNAME = ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2',
             'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']

for M in MODELNAME:
    
    path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/UAVG/' + M + '/'
    paths = '/global/work/nilscp/iSALE/benchmarkFI/UAVG/'  + M + '/'
    NAME = M
    
    os.chdir(path)
    
    arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta,CPPR, paths, NAME, GRIDH, GRIDV)
    INPUT_ASTEROID(param,arr,nvar)
    generateBashfilenames(path, ESTTIME, keyast, keymat)
    
'''
***********************************************************************
'''


# For several cases for ULAR

U = 21210.   #  Impact velocity in m/s


MODELNAME = ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2',
             'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']

for M in MODELNAME:
    
    path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/ULAR/' + M + '/'
    paths = '/global/work/nilscp/iSALE/benchmarkFI/ULAR/'  + M + '/'
    NAME = M
    
    os.chdir(path)
    
    arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta,CPPR, paths, NAME, GRIDH, GRIDV)
    INPUT_ASTEROID(param,arr,nvar)
    generateBashfilenames(path, ESTTIME, keyast, keymat)

'''
***********************************************************************
'''

# for USEC
U = 2380.   #  Impact velocity in m/s

MODELNAME = ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2',
             'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']


GRIDH = np.array([500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500.])

GRIDV = np.array([500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500., 500.,
                  500., 500., 500., 500.])

for M in MODELNAME:
    
    path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/USEC/' + M + '/'
    paths = '/global/work/nilscp/iSALE/benchmarkFI/USEC/'  + M + '/'
    NAME = M
    
    os.chdir(path)
    
    arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta,CPPR, paths, NAME, GRIDH, GRIDV)
    INPUT_ASTEROID(param,arr,nvar)
    generateBashfilenames(path, ESTTIME, keyast, keymat)