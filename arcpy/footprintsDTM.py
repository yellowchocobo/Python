import glob, os
import numpy as np
import shutil
import arcpy
from arcpy import env
from arcpy.sa import *

'''
We are writing a routine to create the footprints of NAC DTMs

I think this should work but maybe it does not work because I am runinng
already another arcpy

'''	

'''
**************************************************************************************************
'''

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# otherwise it makes the script crash and we delete the pyramid right away
# so there is no point on using time on that
#arcpy.env.pyramid = "NONE"
arcpy.env.addOutputsToMap = "FALSE"
    

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "NAC_DTM_ALDROVANDI.TIF"

'''
**************************************************************************************************
'''

def moving_DTM(pathastra):
    
    '''
    move DTMs around the astra folder
    '''
    
    # change to main directory
    os.chdir(pathastra)
    
    # list of all folders in main directory
    folders = glob.glob('NAC_DTM*')
    
    # loop through all folders
    for folder in folders:
        os.chdir(pathastra + folder)
        
        DTM = glob.glob('*.TIF')
        
        #move to the DTM folder if there is a detected DTM
        if not DTM:
            None
            
        else:
            # if the DTM folder exists, simply move the DTM
            if os.path.isdir(pathastra + folder + '/DTM/'):
                shutil.move(pathastra + folder + '/' + DTM, pathastra + folder + '/DTM/' + DTM)
                
            # if not create also the DTM folder    
            else:
                os.makedirs(pathastra + folder + '/DTM/')
                shutil.move(pathastra + folder + '/' + DTM, pathastra + folder + '/DTM/' + DTM)
    

'''
**************************************************************************************************
'''

def list_DTM(pathastra):
    
    '''
    make a list of the absolute path to the DTMs 
    plan to use the output in the footprints_DTM function
    '''
    
    # empty list containing all absolute path to DTMs
    absolute_DTM_path = []
    
    # change to main directory
    os.chdir(pathastra)
    
    # list of all folders in main directory
    folders = glob.glob('*')
    
    # loop through all folders
    for folder in folders:
        
        os.chdir(pathastra + folder + "DTM")
        DTMs = glob.glob("*.TIF")
        
        # in case, there are several hits
        for DTM in DTMs:
            absolute_DTM_path.append(os.path.abspath(DTM))
        
        
    return (absolute_DTM_path)
    
    
'''
**************************************************************************************************
'''

def footprints_DTM(path, pathdab, absolute_DTM_paths):
    
    '''
    path  = "C:/Users/nilscp/Desktop/testarcpy/DTM/"
    pathdab = "C:/Users/nilscp/Desktop/testarcpy/DTM/database.gdb/"
    footprints_DTM(path, pathdab)
    
    pathastra = "/net/astra/astra-01/nilscp/Moon/downloading/NAC_DTM_ALL/"
    absolute_DTM_path = ['Y:/nilscp/Moon/NAC_DTM/NAC_DTM_M102522406_M102529564_DEM.TIF', 
    'Y:/nilscp/Moon/NAC_DTM/NAC_DTM_M104877210_M104884367_DEM.TIF']
    
    '''

    # change to directory of interest
    os.chdir(path)
    
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
        
    # find all the DTMS
    #DTMs = glob.glob('*.tif') #absolute_DTM_path_now
    
    # Spatial reference (plate carree, global spatial reference)
    spatialReference_new = "PROJCS['Equirectangular_Moon',GEOGCS['GCS_Moon',DATUM['D_Moon',SPHEROID['Moon_localRadius',1737400.0,0.0]],PRIMEM['Reference_Meridian',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Plate_Carree'],PARAMETER['false_easting',0.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',0.0],UNIT['Meter',1.0]]"
        
    # loop throgh DTMs
    for ix, pathDTM in enumerate(absolute_DTM_paths):
        
        # name of the DTM
        DTM = pathDTM.split("/")[-1]
        print (DTM)
        
        #path to the DTM
        path_dtm_folder = pathDTM.split(DTM)[0]
        
        # change directory to absolute path
        os.chdir(path_dtm_folder)
        
        # get the footprint
        arcpy.RasterDomain_3d(in_raster=path + DTM, out_feature_class="tmp_polygon_footprints", out_geometry_type="POLYGON")
        
        # projection of the polygon
        arcpy.Project_management("tmp_polygon_footprints", "tmp_polygon_footprints_proj", spatialReference_new)
                
        if arcpy.Exists("polygon_footprints_new"):
            
            # copy tmp all
            arcpy.CopyFeatures_management("polygon_footprints_new", "polygon_footprints_new_tmp")
            
            # add to main polygon and overwrite the main polygon
            arcpy.Merge_management(["polygon_footprints_new_tmp", "tmp_polygon_footprints_proj"] , "polygon_footprints_new")
            
        else:
            arcpy.CopyFeatures_management("tmp_polygon_footprints_proj", "polygon_footprints_new")
        
        # delete old raster and tmp polygon footprints (only one layer will be kept at a time)
        arcpy.Delete_management("tmp_polygon_footprints")
        arcpy.Delete_management("tmp_polygon_footprints_proj")
        arcpy.Delete_management("polygon_footprints_new_tmp")    
        
'''
**************************************************************************************************


def merge_footprints(pathdab, flag_delete):
    

    #pathdab = "C:/Users/nilscp/Desktop/testarcpy/DTM/database.gdb/"
    #flag_delete = False
    
    #merge_footprints(pathdab, flag_delete)

        
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
    
    # list of polygon footprints
    list_polygon_footprints = []
        
    # check if it has been run before
    fcList = arcpy.ListFeatureClasses()
    
    # get only shapefiles corresponding to polygon footprints
    for fc in fcList:
        if fc.endswith("_polygon_footprints"):
            list_polygon_footprints.append(fc)
        else:
            None
            
    
    arcpy.Merge_management(list_polygon_footprints, "polygon_footprints_all")
    
    # if we want to delete all the single footprints, then flag_delete = True
    if flag_delete:
        for fp in list_polygon_footprints:
            arcpy.Delete_management(pathdab + fp)
    else:
        None
        
        
    return list_polygon_footprints



**************************************************************************************************
'''

def add_footprints_attribute(pathdab, infile, list_polygon_footprints, absolute_DTM_paths):
    
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
    
    fieldname1 = arcpy.ValidateFieldName("DTM_name")
    fieldname2 = arcpy.ValidateFieldName("abspath")
    
    # 
    arcpy.AddField_management(infile, fieldname1, "TEXT","","",60)
    arcpy.AddField_management(infile, fieldname2, "TEXT","","",90)
      
    with arcpy.da.UpdateCursor(infile, [fieldname1, fieldname2]) as cursor:    	
        ix = 0
        for row in cursor:
            print (list_polygon_footprints[ix])
            row[0] = list_polygon_footprints[ix]
            row[1] = absolute_DTM_paths[ix]
            cursor.updateRow(row)
            ix = ix + 1
            
    print ("DONE")
    
    
'''
**************************************************************************************************
'''


def intersect_ROI_DTM(pathdab, location_completely_within, outASCII, crater_id):
    
    '''
    simple version of the script above. Only DTMs covering completely 4 times
    the radius of the crater has been selected. (points, buffer, envelope and the
    absolute path to the  DTMs are provided)
    
    
    location_completely_within = "location_overlapDTM_completely_within"
    pathdab = "D:/NAC_DTM/NAC_AMES/arcpy/database.gdb/"
    outASCII = "D:/ANALYSIS/NAC_DTM/ASCII/"
    crater_id = np.genfromtxt("D:/ANALYSIS/NAC_DTM/ASCII/crater_id.txt", dtype=str)
    
    intersect_ROI_DTM(pathdab, location_completely_within, outASCII, crater_id)
    '''
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
       
    # Make a layer from the feature class of the location of the centre of the crater
    arcpy.MakeFeatureLayer_management(location_completely_within, location_completely_within + "_lyr")
    
    # run buffer tool
    out = 'buffer'
    bufferField = "Diameter4txt"
    sideType = "FULL"
    endType = "ROUND"
    dissolveType = "NONE"
    dissolveField = "Distance"
    
    with arcpy.da.UpdateCursor(location_completely_within + "_lyr", ["x_coord", "y_coord", "abspath"]) as cursor:
        
        ix = 0
        
        for row in cursor:
            print (ix)
            if os.path.isfile(outASCII + crater_id[ix] + '.asc'):
                ix = ix + 1
            else:               
                #query selection CENTER
                query = "CRATER_ID = '" + crater_id[ix] + "'"
                arcpy.SelectLayerByAttribute_management(location_completely_within + "_lyr", "NEW_SELECTION", query)
                
                # make a layer of the selection
                arcpy.CopyFeatures_management(location_completely_within + "_lyr", "CENTER_TMP")
                
                # old coordinate systems
                desc = arcpy.Describe("CENTER_TMP")
                spatialReference = desc.spatialReference
                        
                # projection of the right coordinate systems (same as DTM?)
                desc_new = arcpy.Describe(row[2]) 
                spatialReference_new = desc_new.spatialReference
                
                # projection of location
                arcpy.Project_management(in_dataset="CENTER_TMP", out_dataset="CENTER_PROJ", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
                
                # buffer creation 
                arcpy.Buffer_analysis("CENTER_PROJ", out, bufferField, sideType, endType, dissolveType)
                
                # run feature to envelope tool
                arcpy.FeatureEnvelopeToPolygon_management(out,
													      "envelope_proj",
													       "SINGLEPART")
                
                desc_proj = arcpy.Describe("envelope_proj")
                extent = desc_proj.extent
                top = extent.YMax
                bottom = extent.YMin
                left = extent.XMin
                right = extent.XMax
                
                ExtStr = "{} {} {} {}".format(left, bottom, right, top)
                
                # The following inputs are layers or table views: "dtm", "square_test"
                arcpy.Clip_management(in_raster=row[2], rectangle= ExtStr, out_raster= pathdab + "dtm_clip", in_template_dataset="envelope_proj", clipping_geometry="NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
    			
    			# Get input Raster properties
                inRas = arcpy.Raster('dtm_clip')
    						
    			# or I could convert it to ascii
                arcpy.RasterToASCII_conversion(inRas, outASCII + crater_id[ix] + '.asc')
                
                ix = ix + 1
                arcpy.Delete_management("dtm_clip")
                arcpy.Delete_management("envelope_proj")
                arcpy.Delete_management(out)
                arcpy.Delete_management("CENTER_PROJ")
                arcpy.Delete_management("CENTER_TMP")
       
'''
**************************************************************************************************
'''

path  = "D:/NAC_DTM/NAC_DTM_RDR/"
pathdab = "D:/NAC_DTM/NAC_DTM_RDR/arcpy/database.gdb/"
pathastra = ""

os.chdir(path)
filenames = glob.glob("*.TIF")

#NAC_DTM_M1097015769_M1097022916_DEM_polygon_footprints

#filenames2 = filenames[844:]

absolute_DTM_paths = []
for f in filenames:
    absolute_DTM_paths.append(path + f)
    
flag_delete = False

# this need to be tested
#absolute_DTM_paths = list_DTM(pathastra) # get list of the absolute paths to the DTMs
#footprints_DTM(path, pathdab, absolute_DTM_paths)
#list_polygon_footprints = merge_footprints(pathdab, flag_delete) # test it once with flag_delete = False
add_footprints_attribute(pathdab, "polygon_footprints_all", filenames, absolute_DTM_paths)

# for some reasons footprints_DTM stop to work after 70 or so iteration
# I guess it has to do about how much memory
