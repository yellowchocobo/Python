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
path = "D:/ANALYSIS/SIMPLECRATERS_MOON/VALIDATION/LAYERS/"
pathdab = r"D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2015_LARGE_CRATERS/layers/database.gdb/"

# file to copy (locations of preliminary rims)
infile = r"D:/ANALYSIS/SIMPLECRATERS_MOON/LAYERS/layers2019.gdb/rayed_craters_UPD_NILS"

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = pathdab

'''
**************************************************************************************************
'''

# extract the centers of craters (OK regardless of the projection)

arcpy.FeatureToPoint_management(infile, 'CENTERN', 
								"CENTROID")

infile2 = "CENTERN"								
								
# crater name and buffer extent
fieldname1 = arcpy.ValidateFieldName("CRATER_ID")
fieldname2 = arcpy.ValidateFieldName("BUFFER_TXT")

# add fields
arcpy.AddField_management(infile, fieldname1, "TEXT","","",30)
#arcpy.AddField_management(infile, fieldname2, "TEXT","","",30)

# get the number of rows in infile
n = int(arcpy.GetCount_management(infile)[0])

# prepare empty arrays
diam = np.ones(n)
crater_id = np.chararray(n, itemsize=30)
buffer_txt = np.chararray(n, itemsize=30)

#crater_id_list = ['Flamsteed_S', 'Herigonius_K', 'Unnamed_0000' ,'Encke_X',
#                  'Lassell_D','Unnamed_0001','Samir','Unnamed_0002','Unnamed_0003',
#                  'Unnamed_0004','Unnamed_0005','Unnamed_0006','Unnamed_0007','Unnamed_0008',
#                  'Unnamed_0009','Unnamed_0010','Unnamed_0011','Unnamed_0012','Unnamed_0013']

#crater_id = np.array(crater_id_list)

with arcpy.da.UpdateCursor(infile, ["Diam_km", "CRATER_ID"]) as cursor:
    ix = 0
    for row in cursor:
        a = 'crater' + str(int(ix)).zfill(4)
        buffer_value = np.round((row[0]/2.) * 2.0, decimals=4) # changed to 2
        b = str(buffer_value) + ' Kilometers'
        row[1] = a
        #row[1] = crater_id[ix]
        cursor.updateRow(row)
        diam[ix] = row[0]
        crater_id[ix] = a
        buffer_txt[ix] = b
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
outASCII = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2015_LARGE_CRATERS/ascii_olivine/"

# make a lyr
#arcpy.MakeFeatureLayer_management(infile2, infile2 + "_lyr")

'''
**************************************************************************************************
'''

olivine_raster = os.path.join('X:/Moon/Kaguya_mineralogy/',
                              'Lunar_Kaguya_MIMap_MineralDeconv_OlivinePercent_50N50S_colorize.tif')

with arcpy.da.UpdateCursor(infile2, ["Shape@", "x_coord", "y_coord"]) as cursor:
	ix = 0
	for row in cursor:
		print (ix)
		if os.path.isfile(outASCII + crater_id[ix] + '.asc'):
			ix = ix + 1		
		else:
			#query selection CENTER
			query = "CRATER_ID = '" + crater_id[ix] + "'"
			arcpy.SelectLayerByAttribute_management(infile2, "NEW_SELECTION", query)
			
			# make a layer of the selection
			arcpy.CopyFeatures_management(infile2, "CENTER_TMP")
			
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
			arcpy.FeatureEnvelopeToPolygon_management(pathdab + "miniarea_TMP",
													  pathdab + "miniarea_square",
													  "SINGLEPART")
													  
			#reproject in normal
			arcpy.Project_management(in_dataset="miniarea_square", out_dataset= pathdab + "miniarea_square2", out_coor_system="PROJCS['Equirectangular_Moon',GEOGCS['GCS_Moon',DATUM['D_Moon',SPHEROID['Moon_localRadius',1737400.0,0.0]],PRIMEM['Reference_Meridian',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Plate_Carree'],PARAMETER['false_easting',0.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',0.0],UNIT['Meter',1.0]]", transform_method="", in_coor_system=spatialReference_new, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
											
											
			# get the extent of miniarea_square (but this is in the new coordinates! So there is a problem in the next steps
			desc_extent = arcpy.Describe("miniarea_square2")
			extent = desc_extent.extent
			top = extent.YMax
			bottom = extent.YMin
			left = extent.XMin
			right = extent.XMax

			
			ExtStr = "{} {} {} {}".format(left, bottom, right, top)
											
											
			# The following inputs are layers or table views: "dtm", "square_test"
			arcpy.Clip_management(in_raster=olivine_raster, rectangle= ExtStr, out_raster=pathdab + "dtm_clip", in_template_dataset="miniarea_square2", nodata_value="-3.402823e+038", clipping_geometry="NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
			
			desc_DTM = arcpy.Describe("dtm_clip")
			spatialReference_DTM = desc_DTM.spatialReference
			
			# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
			# The following inputs are layers or table views: "Lunar_LRO_LrocKaguya_DEMmerg13"
			arcpy.ProjectRaster_management(in_raster="dtm_clip", out_raster= pathdab + "dtm_clipn", out_coor_system=spatialReference_new, resampling_type="NEAREST", cell_size="59.2252937999999 59.2252937999983", geographic_transform="", Registration_Point="", in_coor_system=spatialReference_DTM)
			
			#low filter pass?
			
			# Get input Raster properties
			inRas = arcpy.Raster('dtm_clipn')
			
			# Project raster
			
			# or I could convert it to ascii
			arcpy.RasterToASCII_conversion(inRas, outASCII + crater_id[ix] + ".asc")
			
			ix = ix + 1
            
			arcpy.Delete_management(pathdab + "dtm_clip")
			arcpy.Delete_management(pathdab + "dtm_clipn")
			arcpy.Delete_management("miniarea_square")
			arcpy.Delete_management("miniarea_square2")
			arcpy.Delete_management("miniarea_TMP")
			arcpy.Delete_management("CENTER_PROJ")
			arcpy.Delete_management("CENTER_TMP")
			
#maybe delete the asc file that are around from before
#f = glob.glob('X:/Moon/downloading/NAC_DTM_ALL/SOSEGENES/ASCII/crater*.asc')
#for ff in f:
#	os.remove(ff)


# Get the centre of the crater
#infile_centroid = 'CENTER'
X = np.ones(n)
Y = np.ones(n)

with arcpy.da.UpdateCursor(infile2, ["SHAPE@XY"]) as cursor:
	ix = 0
	for row in cursor:
		X[ix] = row[0][0]
		Y[ix] = row[0][1]
		ix = ix + 1

# save data with crater_id, X, Y, diameter
data_to_export = np.column_stack((X, Y, diam))
header_txt = 'X;Y;Diam;'
np.savetxt(outASCII + "data.txt", data_to_export, delimiter = ";", header=header_txt, fmt='%10.5f', comments='#')

header_txt = 'crater_id'
np.savetxt(outASCII + "crater_id.txt", crater_id, delimiter = ";", header=header_txt, fmt='%s', comments='#')

