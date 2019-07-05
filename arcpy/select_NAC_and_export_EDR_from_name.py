import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np


'''
We are writing a routine to select images that meet an already known id name
'''	

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# otherwise it makes the script crash and we delete the pyramid right away
# so there is no point on using time on that
#arcpy.env.pyramid = "NONE"
arcpy.env.addOutputsToMap = "FALSE"

'''
**************************************************************************************************
'''
path = "D:/NAC/"
pathdab = "D:/FRESH_IMPACT_WILLIAMS_2018/database.gdb/"

env.workspace = env.scratchWorkspace = pathdab

NAC_lst = np.genfromtxt(path + "nac_images_AMES_to_download.txt", dtype="|S39")
NAC_lstN = []

for NAC in NAC_lst:
    NAC_lstN.append(NAC.decode("utf-8"))

arcpy.MakeFeatureLayer_management("NAC_ALL", "NAC_lyr")


infile = "NAC_lyr"
edr_source_selected = []
       
fields = ["SHAPE@", "productid", "edr_source"]

for NAC in NAC_lst:
    
    NACnames_tmp = NAC.split(".TIF")[0].split("NAC_DTM_")[1].split("_DEM")[0]
    
    # first NAC image
    NACname1 = NACnames_tmp.split("_")[0] + "LE"
    NACname2 = NACnames_tmp.split("_")[0] + "RE"
    
    # second NAC image
    NACname3 = NACnames_tmp.split("_")[1] + "LE" 
    NACname4 = NACnames_tmp.split("_")[1] + "RE"
    
    for NACname in [NACname1, NACname2, NACname3, NACname4]:
    
        # we use a query to select the NAC image of interest
        query = "productid = '" + NACname + "'"
        arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query)
            
        arcpy.CopyFeatures_management(infile, "tempo")
            
        with arcpy.da.UpdateCursor("tempo", fields) as cursor:
            
            for row in cursor:
                EDRname = str(row[2])
                edr_source_selected.append(EDRname)