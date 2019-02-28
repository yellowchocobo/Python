import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *

'''
We are writing a routine to select images that overlap, have incidences between 30 and 70 degrees, and doing so by using the least number of images
'''	

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = "Z:/simple_craters/Censorinus/LAYERS/LAYERS.gdb"

# infile
infile_nacl = "nacl"
infile_nacr = "nacr"

# area boundary
area_bnd = "AREA_CENSO_40radii"
rows = arcpy.da.SearchCursor(area_bnd, ["SHAPE@"])

# get the area of the 7 radii circle
for row in rows:
	feat = row[0]
	area_bnd_value = feat.area

  
def select_NAC(infile, area_bnd, output_name, inc_min, inc_max):
	'''
	# select nac images for a region of interest and an interval of incidences
	'''
	arcpy.SelectLayerByLocation_management(infile, "INTERSECT", area_bnd)

	# search nac images with incidences between 30 and 70 degrees (add to selection)
	arcpy.SelectLayerByAttribute_management(infile, "REMOVE_FROM_SELECTION", '"ctr_incang" < ' + str(inc_min))
	arcpy.SelectLayerByAttribute_management(infile, "REMOVE_FROM_SELECTION", '"ctr_incang" > ' + str(inc_max))
	
	# make a layer of the selection
	arcpy.CopyFeatures_management(infile, output_name)

# select nac images (left and right satellite images) between 30 and 70 degrees of incidence
# I could try something else
select_NAC(infile_nacl, area_bnd, "nacl_selection_20_40", 20., 40.)
select_NAC(infile_nacr, area_bnd, "nacr_selection_20_40", 20., 40.)

# merge selected NACL and NACR in a single layer
arcpy.Merge_management( ["nacl_selection_20_40", "nacr_selection_20_40"], "nac_selection_censo2_20_40")

#sort data by area covered by each nac (descending)
arcpy.Sort_management("nac_selection_censo2_20_40", "nac_selection_censo2_20_40_sorted", [["Shape_Area", "DESCENDING"]])

# need to loop through the data
rows = arcpy.da.SearchCursor("nac_selection_censo2_20_40_sorted", ["SHAPE@", "productid", "edr_source"])

#run it while the area is covered to 90%, and only add file that gives at least 50% larger (no overlap)

NAC_selected = []
edr_source_selected = []
area_merged2 = 0	
infile = "nac_selection_censo2_20_40_sorted"

for row in rows:		
	feat = row[0]
	area_merged2 = feat.area
	NACname = str(row[1])
	EDRname = str(row[2])
	print len(NAC_selected)
	print area_merged2 / area_bnd_value

	
	if area_merged2 / area_bnd_value < 0.9:
		if len(NAC_selected) == 0:
			# this will be run only one time, the first time !!!!
			NAC_selected.append(NACname)
			edr_source_selected.append(EDRname)
			query = "productid = '" + NACname + "'"
			arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query)
			arcpy.CopyFeatures_management(infile, "merged_area_all")
			rows2 = arcpy.da.SearchCursor("merged_area_all", ["SHAPE@"])
			for row2 in rows2:
				# calculate the area of the first NAC
				feat2 = row2[0]
				area_merged2 = feat2.area
		else:
			# copy previous sum of NAC images
			arcpy.Copy_management("Z:/simple_craters/Byrgius_A/LAYERS/LAYERS.gdb/merged_area_all", 
			"Z:/simple_craters/Byrgius_A/LAYERS/LAYERS.gdb/merged_area")
			# the selected NAC is copied into a layer
			query = "productid = '" + NACname + "'"
			arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query)
			arcpy.CopyFeatures_management(infile, "temporary")
            
			# then the area of this NAC is calculated
			rows3 = arcpy.da.SearchCursor("temporary", ["SHAPE@"])
			for row3 in rows3:
				feat3 = row3[0]
				area_merged3 = feat3.area
			# the intersect of the NAC and all selected NAC are done
			arcpy.Intersect_analysis(["merged_area", "temporary"], "temporary_intersect")
			
			# if layer is empty
			if arcpy.management.GetCount("temporary_intersect")[0] == '0':
				area_merged4 = 0.0
			else:
				# the area of the intersect is saved
				rows4 = arcpy.da.SearchCursor("temporary_intersect", ["SHAPE@"])
				area_merged4 = 0.0
				for row4 in rows4:
					feat4 = row4[0]
					area_merged4 += feat4.area
			
			# we only add it up if the area of the intersect covers less than 15% of the entire NAC
			if (area_merged4/area_merged3) < 0.5:
				print (area_merged4/area_merged3)
				# we merge the temporary NAC to the "merged_area" NAC file
				arcpy.Merge_management( ["merged_area", "temporary"], "merged_area_all")
				# we add the name of the NAC and edr source in the lists
				NAC_selected.append(NACname)
				edr_source_selected.append(EDRname)
				# we recalculate the total area covered by the selected NAC
				area_merged2 += (area_merged3 - area_merged4)
				
			else:
				# Nothing is done
				None
			# We delete temporary and temporary_intersect
			arcpy.Delete_management("temporary")
			arcpy.Delete_management("temporary_intersect")
			arcpy.Delete_management("merged_area")
	else:
		None
	

	