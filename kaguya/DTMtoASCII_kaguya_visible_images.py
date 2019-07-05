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
path_to_raster = "E:/Kaguya/ORTHOIMAGES/"

# file to copy (locations of preliminary rims)
#infile = "D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/layers/database.gdb/CRATER_validation"

infile = "X:/Moon/ANALYSIS/COLDSPOTS/database/database2.gdb/coldspotsCopy"

# Set overwrite option
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = 0

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")



pathdab = r"D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/layers/database2.gdb/"

outASCII = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_visible16R/" #"D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/SLDEM2013_Kaguya/ASCII/" # change

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = pathdab

# cellsize of the kaguya DTMs
cellsize = 7.4031617

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

def select_visible_im(bottom, top, left, right):
    
    '''
    equivalent to bottom, top, left, right
    
    
    
    '''    
    dtm_selection = []
    
    # ok this works for every degrees, but every three degrees, need to think a bit more
    min_latitude = np.int(np.floor(bottom))
    max_latitude = np.int(np.ceil(top))
    min_longitude = np.int(np.floor(left))
    max_longitude = np.int(np.ceil(right))
    
    all_longitude = np.arange(0,360,3)
    all_latitude = np.arange(-90, 90, 3)
    
    # get the closest min and max longitudes and latitudes
    min_latitude = all_latitude[all_latitude <=  bottom].max()
    max_latitude = all_latitude[all_latitude >=  top].min()

    min_longitude = all_longitude[all_longitude <=  left].max()
    max_longitude = all_longitude[all_longitude >=  right].min() 
    
    for i in range(min_longitude, max_longitude, 3):
        
        lon2 = 'E' + str(int(i)).zfill(3)
        lon3 = 'E' + str(int(i+3)).zfill(3) # not sure about this one
        
        '''
        Must modify if we want data all the way up to 90 degrees
        
        Not correct as it is here
        '''
        for j in range(min_latitude, max_latitude, 3):
            
            if j  < 0:
                lat1 = 'S' + str(int(abs(j+3))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            elif j == 0:
                lat1 = 'N' + str(int(abs(j+3))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            else:
                lat1 = 'N' + str(int(abs(j+3))).zfill(2)
                lat2 = 'N' + str(int(abs(j))).zfill(2)
                
            #print (lat1,lat2)

            # do
            dtm_selection.append('TCO_MAP_02_' + 
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

#with arcpy.da.UpdateCursor(infile, ["Diam_km", "CRATER_ID", "x_coord", "y_coord"]) as cursor:

with arcpy.da.UpdateCursor(infile, ["Diameter", "CRATER_ID", "Lon", "Lat"]) as cursor:
    ix = 0
    for row in cursor:
        a = 'cpcrater' + str(int(ix)).zfill(4)
        buffer_value = np.round((row[0]/2000.) * 16.0, decimals=4) # changed to 8
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

with arcpy.da.UpdateCursor("CENTERN", ["Shape@", "Lon", "Lat"]) as cursor:
    ix = 0
    for row in cursor:
        
        try:
        
            print (ix)
            if os.path.isfile(outASCII + crater_id[ix] + '_visible.asc'):
                ix = ix + 1
            else:
                
                #distance_8r = ((diam[ix]/2.0) * 16.0) #* 1000.0 # this had not been changed
                
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
                
                # 11.06.2019 convert from local ref to lat/lon with the help of arcgis
                xpoints = np.array([left, left, right, right])
                ypoints = np.array([top, bottom, top, bottom])
                datapoints = np.column_stack((xpoints, ypoints))
                
                commentstxt = "x;y"
                np.savetxt(path_to_raster + "tmp_square_coordinates.txt", datapoints, header=commentstxt, delimiter=";", comments="")
                            
                arcpy.MakeXYEventLayer_management(path_to_raster + "tmp_square_coordinates.txt","x","y","square_points_tmp",spatialReference_new)
                arcpy.FeatureToPoint_management("square_points_tmp", "square_points")
                
                # Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
                arcpy.Project_management(in_dataset="square_points", out_dataset="square_points_latlon", out_coor_system="GEOGCS['GCS_Moon_2000',DATUM['D_Moon_2000',SPHEROID['Moon_2000_IAU_IAG',1737400.0,0.0]],PRIMEM['Reference_Meridian',0.0],UNIT['Degree',0.0174532925199433]]", 
                             transform_method="", in_coor_system=spatialReference_new, 
                             preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
                
                # add x, y coordinates
                arcpy.AddXY_management(in_features="square_points_latlon")
                
                # extract data
                xpoints_converted = []
                ypoints_converted = []
                with arcpy.da.UpdateCursor("square_points_latlon", ["POINT_X", "POINT_Y"]) as cursorpoints:
                    for rowp in cursorpoints:
                        xpoints_converted.append(rowp[0])
                        ypoints_converted.append(rowp[1])
                        
                leftn = xpoints_converted[0]
                rightn = xpoints_converted[2]
                topn = ypoints_converted[0]
                bottomn = ypoints_converted[1]
                           
                # and it is feeded to the select DTMs function
                # it will return the DTMs that cover the squared area
                
                if leftn < 0:
                    leftn = 360 + leftn
                    
                if rightn < 0:
                    rightn = 360 + rightn
                
                # in the case where you have a left boundary > 180.0 and right boundary < 180.0
                #example left = 359.0 and right = 1.0
                if ((leftn > 180) and (rightn < 180)):
                    DTMs_kaguya_list1 = select_visible_im(bottomn, topn, leftn, 360.0)
                    DTMs_kaguya_list2 = select_visible_im(bottomn, topn, 0.0, rightn)
                    DTMs_kaguya_list = DTMs_kaguya_list1 + DTMs_kaguya_list2
                else:
                    DTMs_kaguya_list = select_visible_im(bottomn, topn, leftn, rightn)
                
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
                arcpy.RasterToASCII_conversion(inRas, outASCII + crater_id[ix] + "_visible.asc")
                
                ix = ix + 1
                arcpy.Delete_management("dtm_clip")
                arcpy.Delete_management("miniarea_square")
                arcpy.Delete_management("miniarea_TMP")
                arcpy.Delete_management("CENTER_PROJ")
                arcpy.Delete_management("CENTER_TMP")
                arcpy.Delete_management("mosaic")
                
                for DTM_tmp in DTMs_kaguya_list:
                    arcpy.Delete_management(DTM_tmp.split(".")[0])
        except:
            ix = ix + 1
            
# save data with crater_id, X, Y, diameter
data_to_export = np.column_stack((x_coord, y_coord, diam))
header_txt = 'X;Y;Diam;'
np.savetxt(outASCII + "data.txt", data_to_export, delimiter = ";", header=header_txt, fmt='%10.5f', comments='#')

header_txt = 'crater_id'
np.savetxt(outASCII + "crater_id.txt", crater_id, delimiter = ";", header=header_txt, fmt='%s', comments='#')

