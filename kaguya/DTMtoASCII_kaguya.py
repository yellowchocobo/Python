import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np


'''
We are writing a routine to select images that overlap, have similar incidences,
phase angles, ground azimuth but different phases
'''	

# main path
path_to_raster = "D:/Kaguya/SLDEM2013/"

# file to copy (locations of preliminary rims)
#infile = "D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/layers/database.gdb/CRATER_validation"

infile = "D:/ANALYSIS/SIMPLECRATERS_MOON/LAYERS/layers2019.gdb/rayed_craters_UPD_NILS"

# Set overwrite option
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = 0

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")



pathdab = r"D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/layers/database.gdb/"

outASCII = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii/" #"D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/ASCII/" # change

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = pathdab

# cellsize of the kaguya DTMs
cellsize = 7.4031617

'''
**************************************************************************************************
'''
           
def convert_to_simple_cylindrical(x_degree, y_degree, cellsize = 7.4031617):
    
    # latitude
    y_simple_cylindrical = y_degree*(4096*cellsize)
        
    # longitude
    x_simple_cylindrical = x_degree*(4096*cellsize)
    
    return (x_simple_cylindrical, y_simple_cylindrical)

'''
**************************************************************************************************
'''
def convert_sc360_to_sc180(ulx_scy_360, lrx_scy_360, cellsize = 7.403161724669900):
    
    ulx_degree = ulx_scy_360 / (4096*cellsize)
    
    lrx_degree = lrx_scy_360 / (4096*cellsize)
    
    if ulx_degree >= 180.0:
        ulx_degree_180 = -(360.0 - ulx_degree)
    else:
        ulx_degree_180 = ulx_degree
        
    if lrx_degree >= 180.0:
        lrx_degree_180 = -(360.0 - lrx_degree)
    else:
        lrx_degree_180 = lrx_degree
        
    ulx_sc180 = ulx_degree_180 * (4096*cellsize)
    
    lrx_sc180 = lrx_degree_180 * (4096*cellsize)
    
    return (ulx_sc180, lrx_sc180)
    


'''
**************************************************************************************************
'''

def convert_simple_cylindrical_to_degree(x_scy, y_scy, cellsize = 7.4031617):
    
    # latitude
    y_degree = y_scy / (4096*cellsize)
    
    # longitude 
    x_degree = x_scy / (4096*cellsize)
    
    # convert x_degree to -180+180 degrees
    if x_degree > 180.0:
        x_converted = -(360.0 - x_degree)
    else:
        x_converted = x_degree
    
    return (x_converted, y_degree)

'''
**************************************************************************************************
'''

def convert_simple_cylindrical_to_xdegree(x_scy, cellsize = 7.4031617):
        
    # longitude 
    x_degree = x_scy / (4096*cellsize)
        
    return (x_degree)

'''
**************************************************************************************************
'''

def convert_simple_cylindrical_to_ydegree(y_scy, cellsize = 7.4031617):
    
    # latitude
    y_degree = y_scy / (4096*cellsize)
    
    return (y_degree)

'''
**************************************************************************************************
'''

def get_square(x_simple_cylindrical, y_simple_cylindrical, distance):
    
    top = y_simple_cylindrical + distance
    bottom = y_simple_cylindrical - distance
    left = x_simple_cylindrical - distance
    right = x_simple_cylindrical + distance
     
    return (left, bottom, right, top)

'''
**************************************************************************************************
'''

def select_DTMs(bottom, top, left, right):
    
    '''
    equivalent to bottom, top, left, right
    
    
    
    '''    
    dtm_selection = []
        
    min_latitude = np.int(np.floor(bottom))
    max_latitude = np.int(np.ceil(top))
    min_longitude = np.int(np.floor(left))
    max_longitude = np.int(np.ceil(right))
    
    for i in range(min_longitude, max_longitude):
        
        lon2 = 'E' + str(int(i)).zfill(3)
        
        if i + 1 == 360:
            lon3 = 'E' + str(int(0)).zfill(3)
        else:    
            lon3 = 'E' + str(int(i+1)).zfill(3)
        
        '''
        Must modify if we want data all the way up to 90 degrees
        
        Not correct as it is here
        '''
        for j in range(min_latitude, max_latitude):
            
            if j  < -1:
                lat1 = 'S' + str(int(abs(j+1))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            elif j == -1:
                lat1 = 'N' + str(int(abs(j+1))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            else:
                lat1 = 'N' + str(int(abs(j+1))).zfill(2)
                lat2 = 'N' + str(int(abs(j))).zfill(2)

            # do
            dtm_selection.append('DTM_MAP_01_' + 
                                   lat1 + lon2 + lat2 + lon3 + 'SC_fix.tif')            
            
    return (dtm_selection)

'''
**************************************************************************************************
'''

# extract the centers of craters (OK regardless of the projection)

arcpy.FeatureToPoint_management(infile, pathdab + 'CENTERN', 
								"CENTROID")

infile2 = pathdab + "CENTERN"								
								
# crater name and buffer extent
fieldname1 = arcpy.ValidateFieldName("CRATER_ID")
fieldname2 = arcpy.ValidateFieldName("BUFFER_TXT")

# add fields
arcpy.AddField_management(infile, fieldname1, "TEXT","","",30)

# get the number of rows in infile
n = int(arcpy.GetCount_management(infile)[0])

# prepare empty arrays
diam = np.ones(n)
x_coord = np.ones(n)
y_coord = np.ones(n)
crater_id = np.chararray(n, itemsize=30)
buffer_txt = np.chararray(n, itemsize=30)


#crater_id_list = ['flamsteed_s', 'herigonius_k']

#crater_id_list = ['flamsteed_s', 'herigonius_k', 'unnamed_0000' ,'encke_x',
#                  'lassell_d','unnamed_0001','samir','unnamed_0002','unnamed_0003',
#                  'unnamed_0004','unnamed_0005','unnamed_0006','unnamed_0007','unnamed_0008',
#                  'unnamed_0009','unnamed_0010','unnamed_0011','unnamed_0012','unnamed_0013']

#crater_id = np.array(crater_id_list)

with arcpy.da.UpdateCursor(infile, ["Diam_km", "CRATER_ID", "x_coord", "y_coord"]) as cursor:
    ix = 0
    for row in cursor:
        a = 'crater' + str(int(ix)).zfill(4)
        buffer_value = np.round((row[0]/2.) * 8.0, decimals=4) # changed to 8
        b = str(buffer_value) + ' Kilometers'
        row[1] = a
        #row[1] = crater_id[ix]
        cursor.updateRow(row)
        diam[ix] = row[0]
        crater_id[ix] = a
        buffer_txt[ix] = b
        x_coord[ix] = row[2]
        y_coord[ix] = row[3]
        ix = ix + 1

# add two fields
arcpy.AddField_management(infile2, fieldname1, "TEXT","","",30)
arcpy.AddField_management(infile2, fieldname2, "TEXT","","",30)


with arcpy.da.UpdateCursor(infile2, ["CRATER_ID", "BUFFER_TXT"]) as cursor:
	ix = 0
	for row in cursor:
		row[0] = crater_id[ix]
		row[1] = buffer_txt[ix]
		cursor.updateRow(row)
		ix = ix + 1
										  

# run buffer tool
out = 'miniarea_TMP'
bufferField = "BUFFER_TXT"
sideType = "FULL"
endType = "ROUND"
dissolveType = "NONE"
dissolveField = "Distance"

										  
# This will later be done to all features
arcpy.env.addOutputsToMap = 0

# make a lyr
#arcpy.MakeFeatureLayer_management(infile2, infile2 + "_lyr")


'''
**************************************************************************************************
'''

with arcpy.da.UpdateCursor("CENTERN", ["Shape@", "x_coord", "y_coord"]) as cursor:
    ix = 0
    for row in cursor:
        
        print (ix)
        if os.path.isfile(outASCII + crater_id[ix] + '.asc'):
            ix = ix + 1
        else:
            
            distance_8r = ((diam[ix]/2.0) * 8.0) * 1000.0
            
            #query selection CENTER
            query = "CRATER_ID = '" + crater_id[ix] + "'"
            arcpy.SelectLayerByAttribute_management("CENTERN", "NEW_SELECTION", query)
            
            # make a layer of the selection
            arcpy.CopyFeatures_management("CENTERN", "CENTER_TMP")
            
            # old coordinate systems
            desc = arcpy.Describe("CENTER_TMP")
            spatialReference = desc.spatialReference
            
            # project to the right coordinate systems
            # central meridian should be replaced by the longitude
            # standard parallel_1 by the latitude
            cent_med = np.round(row[1],decimals=0)
            std_parall = np.round(row[2],decimals=0)
            
            str_bef = "PROJCS['Equirectangular_Moon',GEOGCS['GCS_Moon',DATUM['D_Moon',SPHEROID['Moon_localRadius',1737400.0,0.0]],PRIMEM['Reference_Meridian',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Equidistant_Cylindrical'],PARAMETER['false_easting',0.0],PARAMETER['false_northing',0.0],"
            str_cent_med = "PARAMETER['central_meridian'," + str(cent_med) + "],"
            str_parall = "PARAMETER['standard_parallel_1'," + str(std_parall) + "],"
            str_after = "UNIT['Meter',1.0]]"
            
            # the whole string is
            spatialReference_new = str_bef + str_cent_med + str_parall + str_after
            
            # projection
            arcpy.Project_management(in_dataset="CENTER_TMP", out_dataset="CENTER_PROJ", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
            
            # buffer creation 
            arcpy.Buffer_analysis("CENTER_PROJ", out, bufferField, sideType, endType, dissolveType)
            
            # run feature to envelope tool
            arcpy.FeatureEnvelopeToPolygon_management("miniarea_TMP",
													  "miniarea_square",
													  "SINGLEPART")           
            
            #reproject in normal (basic for SLDEM2013)
            #arcpy.Project_management(in_dataset="miniarea_square", out_dataset="miniarea_square2", out_coor_system=spatialReference_DTM_first, transform_method="", in_coor_system=spatialReference_new, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
            
            # get the extent of miniarea_square (but this is in the new coordinates! So there is a problem in the next steps
            desc_extent = arcpy.Describe("miniarea_square")
            extent = desc_extent.extent
            top = extent.YMax
            bottom = extent.YMin
            left = extent.XMin
            right = extent.XMax
            
            
            ExtStr = "{} {} {} {}".format(left, bottom, right, top)
            
                        # the extent of the squared area is converted back to degrees
            leftn = convert_simple_cylindrical_to_xdegree(left)
            bottomn = convert_simple_cylindrical_to_ydegree(bottom)
            rightn = convert_simple_cylindrical_to_xdegree(right)
            topn = convert_simple_cylindrical_to_ydegree(top)
            
            # and it is feeded to the select DTMs function
            # it will return the DTMs that cover the squared area
            
            
            leftn = leftn + np.round(row[1],decimals=0)
            rightn = rightn + np.round(row[1],decimals=0)
            
            if leftn < 0:
                leftn = 360 + leftn
                
            if rightn < 0:
                rightn = 360 + rightn
            
            # in the case where you have a left boundary > 180.0 and right boundary < 180.0
            #example left = 359.0 and right = 1.0
            if ((leftn > 180) and (rightn < 180)):
                DTMs_kaguya_list1 = select_DTMs(bottomn, topn, leftn, 360.0)
                DTMs_kaguya_list2 = select_DTMs(bottomn, topn, 0.0, rightn)
                DTMs_kaguya_list = DTMs_kaguya_list1 + DTMs_kaguya_list2
            else:
                DTMs_kaguya_list = select_DTMs(bottomn, topn, leftn, rightn)
            
            # get the projection of the SLDEM2013 (take the spatial coordinates of the first DTM)
            desc_DTM_first = arcpy.Describe(path_to_raster + DTMs_kaguya_list[0])
            spatialReference_DTM_first = desc_DTM_first.spatialReference
            
            # convert all DTMS to the new projection
            for DTM in DTMs_kaguya_list:
                
                desc_tmp = arcpy.Describe(path_to_raster + DTM)
                spatialRef_tmp = desc_tmp.spatialReference
                name_tmp = DTM.split(".")[0]
                
                arcpy.ProjectRaster_management(in_raster= path_to_raster + DTM, 
                                           out_raster= name_tmp, 
                                           out_coor_system = spatialReference_new,
                                           resampling_type ="NEAREST", 
                                           cell_size="7.4031617 7.4031617", 
                                           geographic_transform="", 
                                           Registration_Point="", 
                                           in_coor_system=spatialRef_tmp) 
            
            # create a mosaic of all imaging covering the area
            mosaic = ""
            
            for ik, DTM in enumerate(DTMs_kaguya_list):
                
                if ik < (len(DTMs_kaguya_list) - 1):
                    mosaic +=  DTM.split(".")[0] + ';'
                    
                else:
                    mosaic +=  DTM.split(".")[0]
                    
            #os.chdir(path_to_raster)
            
            arcpy.MosaicToNewRaster_management(input_rasters=mosaic, output_location=pathdab, 
                                       raster_dataset_name_with_extension="mosaic", 
                                       coordinate_system_for_the_raster= spatialReference_new, 
                                       pixel_type="16_BIT_SIGNED", cellsize="7.4031617", 
                                       number_of_bands="1", mosaic_method="LAST", mosaic_colormap_mode="FIRST")
            
            # The following inputs are layers or table views: "dtm", "square_test"
            arcpy.Clip_management(in_raster="mosaic", rectangle= ExtStr, out_raster=pathdab + "dtm_clip", in_template_dataset=pathdab + "miniarea_square", nodata_value="-3.402823e+038", clipping_geometry="NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
            
            # Get input Raster properties
            inRas = arcpy.Raster('dtm_clip')
            
            # or I could convert it to ascii
            arcpy.RasterToASCII_conversion(inRas, outASCII + crater_id[ix] + ".asc")
            
            ix = ix + 1
            arcpy.Delete_management("dtm_clip")
            arcpy.Delete_management("miniarea_square")
            arcpy.Delete_management("miniarea_TMP")
            arcpy.Delete_management("CENTER_PROJ")
            arcpy.Delete_management("CENTER_TMP")
            arcpy.Delete_management("mosaic")
            
            for DTM_tmp in DTMs_kaguya_list:
                arcpy.Delete_management(DTM_tmp.split(".")[0])
           
# save data with crater_id, X, Y, diameter
data_to_export = np.column_stack((x_coord, y_coord, diam))
header_txt = 'X;Y;Diam;'
np.savetxt(outASCII + "data.txt", data_to_export, delimiter = ";", header=header_txt, fmt='%10.5f', comments='#')

header_txt = 'crater_id'
np.savetxt(outASCII + "crater_id.txt", crater_id, delimiter = ";", header=header_txt, fmt='%s', comments='#')

