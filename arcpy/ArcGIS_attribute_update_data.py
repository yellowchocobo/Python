import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np

#path
path = "X:/Moon/downloading/STEPMED/LAYERS"

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = path + "/layers2.gdb"

#load results
pathdata = 'X:/Moon/downloading/STEPMED/ndata2_05112018/'

alldata = np.loadtxt(pathdata + 'final_res.txt', delimiter = ";", comments="#")

# change nan to -9999
ixnan = np.isnan(alldata)
alldata[ixnan] = -9999

(med_diam, diamf, diam_25, diam_75, diam_min, diam_max, 
            depthf, 
            med_h, unc_h, h_25, h_75, h_min, h_max,
            med_mcw, unc_mcw, mcw_25, mcw_75, mcw_min, mcw_max,
            med_ucw, unc_ucw, ucw_25, ucw_75, ucw_min, ucw_max,
            med_cse, unc_cse, cse_25, cse_75, cse_min, cse_max, dD) = (alldata[:,0],alldata[:,1], alldata[:,2], alldata[:,3], 
																	   alldata[:,4],alldata[:,5], alldata[:,6], alldata[:,7],
																	   alldata[:,8],alldata[:,9], alldata[:,10], alldata[:,11],
																	   alldata[:,12],alldata[:,13], alldata[:,14], alldata[:,15],
																	   alldata[:,16],alldata[:,17], alldata[:,18], alldata[:,19],
																	   alldata[:,20],alldata[:,21], alldata[:,22], alldata[:,23],
																	   alldata[:,24],alldata[:,25], alldata[:,26], alldata[:,27],
																	   alldata[:,28],alldata[:,29], alldata[:,30], alldata[:,31])
infile = ['CENTER', 'fresh_craters_updt_stepmed']

fld_to_be_created = [arcpy.ValidateFieldName("MED_DIAM"), arcpy.ValidateFieldName("UNC_DIAM"), arcpy.ValidateFieldName("DIAM_25"), 
					 arcpy.ValidateFieldName("DIAM_75"), arcpy.ValidateFieldName("DIAM_MIN"), arcpy.ValidateFieldName("DIAM_MAX"),
					 arcpy.ValidateFieldName("DEPTH"),
					 arcpy.ValidateFieldName("MED_H"),arcpy.ValidateFieldName("UNC_H"),arcpy.ValidateFieldName("H_25"),
					 arcpy.ValidateFieldName("H_75"),arcpy.ValidateFieldName("H_MIN"),arcpy.ValidateFieldName("H_MAX"),
					 arcpy.ValidateFieldName("MED_MCW"),arcpy.ValidateFieldName("UNC_MCW"),arcpy.ValidateFieldName("MCW_25"),
					 arcpy.ValidateFieldName("MCW_75"),arcpy.ValidateFieldName("MCW_MIN"),arcpy.ValidateFieldName("MCW_MAX"),
					 arcpy.ValidateFieldName("MED_UCW"),arcpy.ValidateFieldName("UNC_UCW"),arcpy.ValidateFieldName("UCW_25"),
					 arcpy.ValidateFieldName("UCW_75"),arcpy.ValidateFieldName("UCW_MIN"),arcpy.ValidateFieldName("UCW_MAX"),
					 arcpy.ValidateFieldName("MED_CSE"),arcpy.ValidateFieldName("UNC_CSE"),arcpy.ValidateFieldName("CSE_25"),
					 arcpy.ValidateFieldName("CSE_75"),arcpy.ValidateFieldName("CSE_MIN"),arcpy.ValidateFieldName("CSE_MAX"),
					 arcpy.ValidateFieldName("dD")]
					 

# add fields
for shapef in infile:
	for field in fld_to_be_created:
		arcpy.AddField_management(shapef, field, "FLOAT", "", "")

# update all the data
for shapef in infile:
	with arcpy.da.UpdateCursor(shapef, fld_to_be_created) as cursor:
		ix = 0
		for row in cursor:
			for i in range(32):
				row[i] = alldata[ix,i]	
			cursor.updateRow(row)
			ix = ix + 1
			
fld_to_be_created2 = [arcpy.ValidateFieldName("MARE")]
for shapef in infile:
	for field in fld_to_be_created2:
		arcpy.AddField_management(shapef, field, "LONG", "", "")