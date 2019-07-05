import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np


'''
We are writing a routine to calculate the distance to the min, 25p, median, 75p
and max distance to the edge of the continuous ejecta (or other types of ejecta)

'''	

'''
**************************************************************************************************
'''

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")
    

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "NAC_DTM_ALDROVANDI.TIF"

def ejecta_distance(path, pathdab, infile_center_crater, infile_ejecta_polygon):
    
    '''
    path  = "C:/Users/nilscp/Desktop/testarcpy/DTM/"
    pathdab = "C:/Users/nilscp/Desktop/testarcpy/DTM/database.gdb/"
    infile_center_crater = "CENTER_crater001"
    infile_ejecta_polygon = "continuousej_crater001"
    
    
    I need somehow to get the diameter of the crater and infile_center_crater
    and infile_ejecta_polygon should have the same coordinates system (preferentially
    equirectangular with lat, lon of the location of the centre of the crater)
    
    '''

    # change to directory of interest
    os.chdir(path)
    
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
    
    # first we need to densify the number vertices along the polygon
    arcpy.Densify_edit(in_features=infile_ejecta_polygon, densification_method="ANGLE", max_angle="1.0")
    #max angle of 0.5 was making way too many points
    
    # create name for new point shapefile
    tmpstr = infile_ejecta_polygon.split("_")
    infile_ejecta_vertices = tmpstr[0] + "_vertices_" + tmpstr[1]
    
    # Feature vertice to points
    arcpy.FeatureVerticesToPoints_management(in_features=infile_ejecta_polygon, 
                                         out_feature_class=infile_ejecta_vertices,
                                         point_location="ALL")

    # add xy coordinates for the crater centre and ejecta vertices
    arcpy.AddXY_management(in_features=infile_ejecta_vertices)
    arcpy.AddXY_management(in_features=infile_center_crater)
    
    # For the crater center
    with arcpy.da.SearchCursor(infile_center_crater, ["POINT_X", "POINT_Y"]) as cursor:    	
        for row in cursor:
            xcenter = row[0]
            ycenter = row[1]
    
    # For the vertices from the polygon
    n = int(arcpy.GetCount_management(infile_ejecta_vertices)[0])
    xvertices = np.zeros(n)
    yvertices = np.zeros(n)
    
    with arcpy.da.SearchCursor(infile_ejecta_vertices, ["POINT_X", "POINT_Y"]) as cursor:    	
        ix = 0
        for row in cursor:
            xvertices[ix] = row[0]
            yvertices[ix] = row[1]
            ix = ix + 1
    
    # calculating the distance to each vertice
    a = (yvertices-ycenter)**2.0
    b = (xvertices-xcenter)**2.0
    dist = np.sqrt(a + b)
    
    # get the min, 25p, median, 75p and max distance
    ej_min_distance = np.min(dist)
    ej_25p_distance = np.percentile(dist,25.0)
    ej_median_distance = np.percentile(dist,50.0)
    ej_75p_distance = np.percentile(dist,75.0)
    ej_max_distance = np.max(dist)
    
    return (ej_min_distance, ej_25p_distance, ej_median_distance, 
            ej_75p_distance, ej_max_distance)

       
'''
**************************************************************************************************
'''




