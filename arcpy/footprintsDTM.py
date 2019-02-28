import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np
import shutil

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
    for pathDTM in absolute_DTM_paths:
                
        # name of the DTM
        DTM = pathDTM.split("/")[-1]
        print (DTM)
        
        #path to the DTM
        path_dtm_folder = pathDTM.split(DTM)[0]
        
        # change directory to absolute path
        os.chdir(path_dtm_folder)
        
        # get the extent of miniarea_square (but this is in the new coordinates! So there is a problem in the next steps
        desc_DTM = arcpy.Describe(pathDTM)
        
        # old spatial reference
        spatialReference_DTM = desc_DTM.spatialReference
        
        # projection
        arcpy.ProjectRaster_management(in_raster=pathDTM, out_raster= path + DTM + "_project.TIF", out_coor_system=spatialReference_new, resampling_type="NEAREST", geographic_transform="", Registration_Point="", in_coor_system=spatialReference_DTM)
                
        # conversion to polygon
        arcpy.RasterDomain_3d(in_raster=path + DTM + "_project.TIF", out_feature_class=DTM.split(".")[0] + "_polygon_footprints", out_geometry_type="POLYGON")

        # delete old raster
        arcpy.Delete_management(path + DTM + "_project.TIF")
        
        
'''
**************************************************************************************************
'''

def merge_footprints(pathdab, flag_delete):
    
    '''
    pathdab = "C:/Users/nilscp/Desktop/testarcpy/DTM/database.gdb/"
    flag_delete = False
    
    merge_footprints(pathdab, flag_delete)
    '''
        
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


'''
**************************************************************************************************
'''

def add_footprints_attribute(pathdab, infile, list_polygon_footprints, absolute_DTM_paths):
    
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
    
    fieldname1 = arcpy.ValidateFieldName("DTM_name")
    fieldname2 = arcpy.ValidateFieldName("abspath")
    
    # 
    arcpy.AddField_management(infile, fieldname1, "TEXT","","",30)
    arcpy.AddField_management(infile, fieldname2, "TEXT","","",60)
      
    with arcpy.da.UpdateCursor(infile, [fieldname1, fieldname2]) as cursor:    	
        ix = 0
        for row in cursor:	
            row[0] = list_polygon_footprints[ix].split("_polygon_footprints")[0]
            row[1] = absolute_DTM_paths[ix]
            cursor.updateRow(row)
            ix = ix + 1
            
    print ("DONE")
    
    
'''
**************************************************************************************************
'''

def intersect_ROI_DTM(pathdab, infile_ROI, infile_DTM, outASCII):
    
    '''
    give absolute path
    need to be careful with the differences in projection coordinates infile_ROI 
    and infile_DTM should be in plate carree
    
    Need to double check that
    
    infile_ROI should be maybe only 4 times the diameter! Need to fix that
    and  be careful with the projection and so on...
    
    '''
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
    
    #
    arcpy.Intersect_analysis(in_features=[infile_ROI, infile_DTM], out_feature_class="intersect", join_attributes="NO_FID", cluster_tolerance="-1 Unknown", output_type="INPUT")
    
               
    # Make a layer from the feature class
    arcpy.MakeFeatureLayer_management("intersect", "intersect_lyr")
    
    with arcpy.da.UpdateCursor("intersect_lyr", ["CRATER_ID", "DTM_name", "abspath"]) as cursor:
        
        ix = 0
        
        for row in cursor:
            
            if ix > 0:
                arcpy.SelectLayerByAttribute_management("intersect_lyr", "CLEAR_SELECTION")
            
            else:
                None
                
            # query
            query = "CRATER_ID = '" + row[0] + "'"
            print (query)
            arcpy.SelectLayerByAttribute_management("intersect_lyr", "NEW_SELECTION", query)
                 
            # make a layer of the selection
            arcpy.CopyFeatures_management("intersect_lyr", "intersect_tmp")
            
            # old coordinate systems
            desc = arcpy.Describe("intersect_tmp")
            spatialReference = desc.spatialReference

                
            # project to the right coordinate systems (same as DTM?)
            desc_new = arcpy.Describe(row[2]) #abspath
            spatialReference_new = desc_new.spatialReference

                
            # projection of footprints to correct coordinate system
            arcpy.Project_management(in_dataset="intersect_tmp", out_dataset="intersect_tmp_proj", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")	
            
            desc_proj = arcpy.Describe("intersect_tmp_proj")
            extent = desc_proj.extent
            top = extent.YMax
            bottom = extent.YMin
            left = extent.XMin
            right = extent.XMax
            
            ExtStr = "{} {} {} {}".format(left, bottom, right, top)
            
            # The following inputs are layers or table views: "dtm", "square_test"
            arcpy.Clip_management(in_raster=row[2], rectangle= ExtStr, out_raster= pathdab + "dtm_clip", in_template_dataset="intersect_tmp_proj", nodata_value="-3.402823e+038", clipping_geometry="NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
			
			# Get input Raster properties
            inRas = arcpy.Raster('dtm_clip')
						
			# or I could convert it to ascii
            arcpy.RasterToASCII_conversion(inRas, outASCII + row[0] + '.asc')
            
    
       
'''
**************************************************************************************************
'''

path  = "C:/Users/nilscp/Desktop/testarcpy/DTM/"
pathdab = "C:/Users/nilscp/Desktop/testarcpy/DTM/database.gdb/"
pathastra = ""
absolute_DTM_paths = ['Y:/nilscp/Moon/NAC_DTM/NAC_DTM_M102522406_M102529564_DEM.TIF', 
                     'Y:/nilscp/Moon/NAC_DTM/NAC_DTM_M104877210_M104884367_DEM.TIF']
flag_delete = False

# this need to be tested
#absolute_DTM_paths = list_DTM(pathastra) # get list of the absolute paths to the DTMs
footprints_DTM(path, pathdab, absolute_DTM_paths)
list_polygon_footprints = merge_footprints(pathdab, flag_delete) # test it once with flag_delete = False
add_footprints_attribute(pathdab, "polygon_footprints_all", list_polygon_footprints, absolute_DTM_paths)
