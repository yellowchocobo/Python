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
L = np.array([100, 1000])   # Diameter of the projectile (m)
g = 1.62                    # gravity of the targeted body (SI)
rho = 2850. * 0.9           # target's density (kg/m3)
teta = 2850.                # projectile's density (kg/m3)
CPPR = 20

variable = np.array(["PATH","MODEL","TDUMP","GRIDH","GRIDV","GRIDSPC","DT","TEND","DTSAVE"])


GRIDH = np.array([650., 650.])

GRIDV = np.array([650., 650.])

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

ESTTIME = ['24-00:00:00','24-00:00:00']

keyast = 'asteroid*L*'
keymat = 'mat*'
generateBashfilenames(path,ESTTIME, keyast, keymat)



'''
***********************************************************************
'''


# For several cases for Uavg

U = 11380.                  #  Impact velocity in m/s


MODELNAME = ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2',
             'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']

CP = [5.0, 10.0, 20.0]
CPT = ['CPPR5', 'CPPR10', 'CPPR20']

GRIDH = np.array(([180., 180.], [325., 325.], [650., 650.]))

GRIDV = np.array(([180., 180.], [325., 325.], [650., 650.]))

for i, CPPR in enumerate(CP):
   
    for M in MODELNAME:
    
        path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/CPPR/UAVG/' + CPT[i] + '/' + M + '/'
        paths = '/global/work/nilscp/iSALE/benchmarkFI/CPPR/UAVG/' + CPT[i] + '/' + M + '/'
        NAME = M
        
        #print (path, paths, NAME, GRIDH[i])
        os.chdir(path)
        
        arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta, CPPR, paths, NAME, GRIDH[i], GRIDV[i])
        INPUT_ASTEROID(param,arr,nvar)
        generateBashfilenames(path, ESTTIME, keyast, keymat)
    

'''
***********************************************************************
'''

# for USEC
U = 2380.                  #  Impact velocity in m/s


MODELNAME = ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2',
             'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']

CP = [5.0, 10.0, 20.0]
CPT = ['CPPR5', 'CPPR10', 'CPPR20']

GRIDH = np.array(([180., 180.], [250., 250.], [500., 500.]))

GRIDV = np.array(([180., 180.], [250., 250.], [500., 500.]))

for i, CPPR in enumerate(CP):
   
    for M in MODELNAME:
    
        path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/projects/benchmarkFI/CPPR/USEC/' + CPT[i] + '/' + M + '/'
        paths = '/global/work/nilscp/iSALE/benchmarkFI/CPPR/USEC/' + CPT[i] + '/' + M + '/'
        NAME = M
        
        #print (path, paths, NAME, GRIDH[i])
        os.chdir(path)
        
        arr = INPUT_ASTEROID_ARRAY(U,L,g,rho,teta, CPPR, paths, NAME, GRIDH[i], GRIDV[i])
        INPUT_ASTEROID(param,arr,nvar)
        generateBashfilenames(path, ESTTIME, keyast, keymat)