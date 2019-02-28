import arcpy
import numpy as np
from arcpy import env
import glob, os
from arcpy.sa import *

'''
We are writing a routine to select images that overlap, have similar incidences,
phase angles, ground azimuth but different phases
'''	

# Set overwrite option
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# define paths and workspace (I need to create the gdb at some points)
env.workspace = env.scratchWorkspace = "Z:/NAC_DTM/COPERNICUS/LAYERS/layers.gdb"

'''
*********************************************************************************
'''

def get_extent(infile):

	# number of rows in the table 
	n = np.int(arcpy.GetCount_management(infile)[0])

	top = np.ones(n)
	bottom  = np.ones(n)
	left = np.ones(n)
	right = np.ones(n)
	# get all the rows in the table
	rows = arcpy.da.SearchCursor(infile, ["OID@", "SHAPE@"])
	
	# index
	ix = 0

	# loop through the rows and get extent of each polygon
	for row in rows:		
		feat = row[1]
		extent = feat.extent
		top[ix] = extent.YMax
		bottom[ix] = extent.YMin
		left[ix] = extent.XMin
		right[ix] = extent.XMax
		ix = ix + 1
	
	return (top, bottom, left, right)
	

'''
*********************************************************************************
'''	

def load_table(filename):
	
	data = np.loadtxt(filename, skiprows = 1,delimiter = ",")
	lat = data[:,3]
	lon = data[:,4]
	res = data[:,5]
	emi = data[:,6]
	inc = data[:,7]
	pha = data[:,8]
	azigrd = data[:,14]
	
	return (lat, lon, res, emi, inc, pha, azigrd)



# I have trouble because the 32-bit Python is installed with ArcGIs for some reasons
# how to read the characters

'''
*********************************************************************************
'''	

def get_characters(infile):

	productid = []
	edr_source = []
	start_time = []
	pilot_url = []
	
	rows = arcpy.da.SearchCursor(infile, ["productid", "edr_source", "start_time", "pilot_url"])
	for row in rows:
		productid.append(row[0])
		edr_source.append(row[1])
		start_time.append(row[2])
		pilot_url.append(row[3])
	
	productid = np.array(productid)
	edr_source = np.array(edr_source)
	start_time = np.array(start_time)
	pilot_url = np.array(pilot_url)
	
	return (productid, edr_source, start_time, pilot_url)

'''
*********************************************************************************
'''	

def tmp_array(table, idx):

	'''return tmp array with one less entry'''
	
	n = len(table)

	if idx == 0:
		var = table[1:]
	elif idx == n - 1:
		var = table[:-1]
	else:
		var = np.concatenate((table[:idx],table[idx+1:]))
	return var 

'''
*********************************************************************************
'''	
def convert_coordinates(coordinates):
	'''
	For some the reasons, the coordinates we get from the extent is kind of fucked up.
	I need to dig a bit more in that. Here is a quick fix. The errors may propagate 
	along.
	'''
	modcoord = 360.0 + (coordinates)
	return modcoord

'''
*********************************************************************************
'''	
	
# inputs
infile = "nac_all"
numbers = "C:/Users/nilscp/Downloads/PHASE/numbers.csv" 
texts =  "C:/Users/nilscp/Downloads/PHASE/string.csv" 

# load extents
top, bottom, left, right = get_extent(infile)
	
# load numbers 	
(lat, lon, res, emi, inc, pha, azigrd) = load_table(numbers)
	
# load strings
(productid, edr_source, start_time, pilot_url) = get_characters(infile)

# change to 0-360 degrees nomenclature (something weird with the coordinates)
left_360 = convert_coordinates(left)
right_360 = convert_coordinates(right)

'''
*********************************************************************************
'''	

# loop through each 

#input files
IMAGE1_F = np.array([])
IMAGE2_F = np.array([])
EMI1_F = np.array([])
EMI2_F = np.array([])
INC1_F = np.array([])
INC2_F = np.array([])
PHA1_F = np.array([])
PHA2_F = np.array([])
AZI1_F = np.array([])
AZI2_F = np.array([])
RES1_F = np.array([])
RES2_F = np.array([])
EDR_PRODUCT1_F = np.array([])
EDR_PRODUCT2_F = np.array([])

for ix_tmp, number in np.ndenumerate(lat):
	
	idx = ix_tmp[0]
	
	# creating temporary array with a length (n - 1)
	tmp_lat = tmp_array(lat, idx)
	tmp_lon = tmp_array(lon, idx)
	tmp_top = tmp_array(top, idx)
	tmp_bottom = tmp_array(bottom, idx)
	tmp_left = tmp_array(left_360, idx)
	tmp_right = tmp_array(right_360, idx)
	tmp_emi = tmp_array(emi, idx)
	tmp_inc = tmp_array(inc, idx)
	tmp_pha = tmp_array(pha, idx)
	tmp_azigrd = tmp_array(azigrd, idx)
	tmp_res = tmp_array(res, idx)
	
	# actual value for the entry
	act_lat = lat[idx]
	act_lon = lon[idx]
	act_emi = emi[idx]
	act_inc = inc[idx]
	act_pha = pha[idx]
	act_azigrd = azigrd[idx]
	act_top = top[idx]
	act_bottom = bottom[idx]
	act_left = left_360[idx]
	act_right = right_360[idx]
	act_res = res[idx]
	
	# threshold values
	degree_inc = 3.0
	degree_azi = 10.0
	
	# calculating boundaries
	min_inc = act_inc - degree_inc
	max_inc = act_inc + degree_inc
	min_azi = act_azigrd - degree_azi
	max_azi = act_azigrd + degree_azi
	
	# search for footprints that overlap, have an emission within +- 2 degrees and azigrd within +-10 degrees or something like that
	#ix = np.where( (act_lat > tmp_bottom) & (act_lat < tmp_top) & (act_lon > tmp_left) & (act_lon < tmp_right))
	
	# is the upper left, upper right, lower left or lower right (in order) of the footprint overlapping or containing other images?	
	ix = np.where(((np.logical_and(act_top > tmp_bottom, act_top < tmp_top)) & (np.logical_and(act_left > tmp_left, act_left < tmp_right))) |
				  ((np.logical_and(act_top > tmp_bottom, act_top < tmp_top)) & (np.logical_and(act_right > tmp_left, act_right < tmp_right))) |
				  ((np.logical_and(act_bottom > tmp_bottom, act_bottom < tmp_top)) & (np.logical_and(act_left > tmp_left, act_left < tmp_right))) |
				  ((np.logical_and(act_bottom > tmp_bottom, act_bottom < tmp_top)) & (np.logical_and(act_right > tmp_left, act_right < tmp_right)))|
				  ((np.logical_and(tmp_top > act_bottom, tmp_top < act_top)) & (np.logical_and(tmp_left > act_left, tmp_left < act_right))) |
				  ((np.logical_and(tmp_top > act_bottom, tmp_top < act_top)) & (np.logical_and(tmp_right > act_left, tmp_right < act_right))) |
				  ((np.logical_and(tmp_bottom > act_bottom, tmp_bottom < act_top)) & (np.logical_and(tmp_left > act_left, tmp_left < act_right))) |
				  ((np.logical_and(tmp_bottom > act_bottom, tmp_bottom < act_top)) & (np.logical_and(tmp_right > act_left, tmp_right < act_right))))
	
	# candidate footprints
	iy = ix[0] + 1
	cand_lat = lat[iy]
	cand_lon = lon[iy]
	cand_emi = emi[iy]
	cand_inc = inc[iy]
	cand_pha = pha[iy]
	cand_azigrd = azigrd[iy]
	cand_top = top[iy]
	cand_bottom = bottom[iy]
	cand_left = left_360[iy]
	cand_right = right_360[iy]
	cand_productid = productid[iy]
	cand_edr_source = edr_source[iy]
	cand_start_time = start_time[iy]
	cand_pilot_url = pilot_url[iy]
	cand_res = res[iy]
	
	iz = np.where(((np.logical_and(cand_inc >= min_inc, cand_inc <= max_inc)) & (np.logical_and(cand_azigrd >= min_azi, cand_azigrd < max_azi))))
	izz = iz[0]
	
	# footprint that respond to criteria
	for z in izz:
		IMAGE1_F = np.append(IMAGE1_F, productid[idx])
		IMAGE2_F = np.append(IMAGE2_F, cand_productid[z])
		RES1_F = np.append(RES1_F, res[idx])
		RES2_F = np.append(RES2_F, cand_res[z])
		EMI1_F = np.append(EMI1_F, act_emi)
		EMI2_F = np.append(EMI2_F, cand_emi[z])
		INC1_F = np.append(INC1_F, act_inc)
		INC2_F = np.append(INC2_F, cand_inc[z])
		PHA1_F = np.append(PHA1_F, act_pha)
		PHA2_F = np.append(PHA2_F, cand_pha[z])
		AZI1_F = np.append(AZI1_F, act_azigrd)
		AZI2_F = np.append(AZI2_F, cand_azigrd[z])
		EDR_PRODUCT1_F = np.append(EDR_PRODUCT1_F, edr_source[idx])
		EDR_PRODUCT2_F = np.append(EDR_PRODUCT2_F, cand_edr_source[z])

# I need to find the union of the two rasters based on the four corners of each footprint
		

# save everything in a .csv file
pathsave = "C:/Users/nilscp/Downloads/"
os.chdir(pathsave)
fname_nbr = "summary_phase_nbr.csv"
fname_txt = "summary_phase_txt.csv"

header_txt = "image1;image2;download1;download2"
header_nbr = "res1;res2;emi1;emi2;inc1;inc2;pha1;pha2;azi1;azi2"

output_txt = np.column_stack((IMAGE1_F, IMAGE2_F, EDR_PRODUCT1_F, EDR_PRODUCT2_F))
output_nbr = np.column_stack((RES1_F, RES2_F, EMI1_F, EMI2_F, INC1_F, INC2_F, PHA1_F, PHA2_F, AZI1_F, AZI2_F))

# save data
np.savetxt(fname_nbr, output_nbr, header=header_nbr, delimiter=";", fmt = ["%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e"])
np.savetxt(fname_txt, output_txt, header=header_txt, delimiter=";", fmt = ["%s","%s","%s","%s"])

# not taking into account places matches for the same productID
IMAGE1_F_mod = np.array([])
IMAGE2_F_mod = np.array([])
EMI1_F_mod = np.array([])
EMI2_F_mod = np.array([])
INC1_F_mod = np.array([])
INC2_F_mod = np.array([])
PHA1_F_mod = np.array([])
PHA2_F_mod = np.array([])
AZI1_F_mod = np.array([])
AZI2_F_mod = np.array([])
RES1_F_mod = np.array([])
RES2_F_mod = np.array([])
EDR_PRODUCT1_F_mod = np.array([])
EDR_PRODUCT2_F_mod = np.array([])

fname_nbr_mod = "summary_phase_nbr_mod.csv"
fname_txt_mod = "summary_phase_txt_mod.csv"

for im, match in np.ndenumerate(IMAGE1_F):
	imm = im[0]
	
	if match[:-2] == IMAGE2_F[imm][:-2]:
		None
	else:
		IMAGE1_F_mod = np.append(IMAGE1_F_mod, IMAGE1_F[imm])
		IMAGE2_F_mod = np.append(IMAGE2_F_mod, IMAGE2_F[imm])
		EMI1_F_mod = np.append(EMI1_F_mod, EMI1_F[imm])
		EMI2_F_mod = np.append(EMI2_F_mod, EMI2_F[imm])
		INC1_F_mod = np.append(INC1_F_mod, INC1_F[imm])
		INC2_F_mod = np.append(INC2_F_mod, INC2_F[imm])
		PHA1_F_mod = np.append(PHA1_F_mod, PHA1_F[imm])
		PHA2_F_mod = np.append(PHA2_F_mod, PHA2_F[imm])
		AZI1_F_mod = np.append(AZI1_F_mod, AZI1_F[imm])
		AZI2_F_mod = np.append(AZI2_F_mod, AZI2_F[imm])
		RES1_F_mod = np.append(RES1_F_mod, RES1_F[imm])
		RES2_F_mod = np.append(RES2_F_mod, RES2_F[imm])
		EDR_PRODUCT1_F_mod = np.append(EDR_PRODUCT1_F_mod, EDR_PRODUCT1_F[imm])
		EDR_PRODUCT2_F_mod = np.append(EDR_PRODUCT2_F_mod, EDR_PRODUCT2_F[imm])
		
output_txt_mod = np.column_stack((IMAGE1_F_mod, IMAGE2_F_mod, EDR_PRODUCT1_F_mod, EDR_PRODUCT2_F_mod))
output_nbr_mod = np.column_stack((RES1_F_mod, RES2_F_mod, EMI1_F_mod, EMI2_F_mod, INC1_F_mod, INC2_F_mod, PHA1_F_mod, PHA2_F_mod, AZI1_F_mod, AZI2_F_mod))

# save data
np.savetxt(fname_nbr_mod, output_nbr_mod, header=header_nbr, delimiter=";", fmt = ["%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e"])
np.savetxt(fname_txt_mod, output_txt_mod, header=header_txt, delimiter=";", fmt = ["%s","%s","%s","%s"])


'''
*********************************************************************************
'''	

def select_match(IMAGE1_F_mod, IMAGE2_F_mod, EDR_PRODUCT1_F_mod, EDR_PRODUCT2_F_mod):

	'''
	# inputs
	infile = "NACL_ALL"
	textpath =  "C:/Users/nilscp/Downloads/PHASE/filter_nac_string.csv"
	lines = []
	with open(textpath) as file: 
		for line in file:
			lines.append(line)
	IMAGE1_F_mod = []
	IMAGE2_F_mod = []
	EDR_PRODUCT1_F_mod = []
	EDR_PRODUCT2_F_mod = []
	for line in lines:
		IMAGE1_F_mod.append(line.split(',')[0])
		IMAGE2_F_mod.append(line.split(',')[1])
		EDR_PRODUCT1_F_mod.append(line.split(',')[2])
		EDR_PRODUCT2_F_mod.append(line.split(',')[3])
	IMAGE1_F_mod = np.array(IMAGE1_F_mod)
	IMAGE2_F_mod = np.array(IMAGE2_F_mod)
	EDR_PRODUCT1_F_mod = np.array(EDR_PRODUCT1_F_mod)
	EDR_PRODUCT2_F_mod = np.array(EDR_PRODUCT2_F_mod)
	UNIQUE, EDR_ALL = select_match(IMAGE1_F_mod, IMAGE2_F_mod,EDR_PRODUCT1_F_mod, EDR_PRODUCT2_F_mod)
	'''

	ALL = np.concatenate((IMAGE1_F_mod, IMAGE2_F_mod))
	EDR_ALL = np.concatenate((EDR_PRODUCT1_F_mod, EDR_PRODUCT2_F_mod))
	UNIQUE, indices = np.unique(ALL,return_index=True)
	
	for iy, id in np.ndenumerate(UNIQUE):
		iyy = iy[0]
		query1 = "productid = '" + UNIQUE[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "ADD_TO_SELECTION", query1)
	arcpy.CopyFeatures_management(infile,"UNIQUE_MATCH")
	
	return UNIQUE, EDR_ALL[indices]
'''
*********************************************************************************
'''	

def get_union2(infile, IMAGE1_F_mod, IMAGE2_F_mod):

	'''
	This is to calculate the final intersection of all filtered NAC images. In other words, the area we will be able to work with!!
	empty layers might lead to some errors
	infile2 = "UNIQUE_MATCH"
	get_union2(infile2, IMAGE1_F_mod, IMAGE2_F_mod)
	'''
	
	n = len(IMAGE1_F_mod)
	
	int_top = np.ones(n)
	int_bottom = np.ones(n)
	int_left = np.ones(n)
	int_right = np.ones(n)
	
	for iy, id in np.ndenumerate(IMAGE1_F_mod):
	
		iyy = iy[0]
	
		query1 = "productid = '" + IMAGE1_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query1)
		arcpy.CopyFeatures_management(infile,"test_query1")


		query2 = "productid = '" + IMAGE2_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query2)
		arcpy.CopyFeatures_management(infile,"test_query2")

		inFeatures = ["test_query1","test_query2"]
		outFeatures = IMAGE1_F_mod[iyy] + "_" + IMAGE2_F_mod[iyy] + "_intersect"
		arcpy.Intersect_analysis(inFeatures, outFeatures)

		arcpy.Delete_management("test_query1")
		arcpy.Delete_management("test_query2")

# after that I can crop the images with maptrim in ISIS3 
(int_top, int_bottom, int_left, int_right) = get_union(infile, IMAGE1_F_mod, IMAGE2_F_mod)
output = np.column_stack((int_top, int_bottom, int_left, int_right))
os.chdir(pathsave)
np.savetxt("maptrim2.txt", output, header="int_top;int_bottom;int_left;int_right", delimiter=";", fmt = ["%1.6e","%1.6e","%1.6e","%1.6e"])

'''
*********************************************************************************
'''	
	

def update_values_with_correct(infile, IMAGE1_F_mod, IMAGE2_F_mod, pathsave):

	'''
	empty layers might lead to some errors
	infile = "NACL_ALL"
	pathsave = "C:/Users/nilscp/Downloads/PHASE/"
	update_values_with_correct(infile, IMAGE1_F_mod, IMAGE2_F_mod, pathsave)
	
	'''
	
	n = len(IMAGE1_F_mod)
	
	res1 = np.ones(n)
	res2 = np.ones(n)
	emi1 = np.ones(n)
	emi2 = np.ones(n)
	inc1 = np.ones(n)
	inc2 = np.ones(n)
	pha1 = np.ones(n)
	pha2 = np.ones(n)
	azi1 = np.ones(n)
	azi2 = np.ones(n)
	
	for iy, id in np.ndenumerate(IMAGE1_F_mod):
	
		iyy = iy[0]
	
		query = "productid = '" + IMAGE1_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query)
		rows = arcpy.da.SearchCursor(infile, ["ctr_pixres", "ctr_emiang", "ctr_incang", "ctr_phase", "sunazgrnd"])
		
		for row in rows:
			res1[iyy] = row[0]
			emi1[iyy] = row[1]
			inc1[iyy] = row[2]
			pha1[iyy] = row[3]
			azi1[iyy] = row[4]

		query = "productid = '" + IMAGE2_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query)
		rows = arcpy.da.SearchCursor(infile, ["ctr_pixres", "ctr_emiang", "ctr_incang", "ctr_phase", "sunazgrnd"])

		for row in rows:
			res2[iyy] = row[0]
			emi2[iyy] = row[1]
			inc2[iyy] = row[2]
			pha2[iyy] = row[3]
			azi2[iyy] = row[4]

	output = np.column_stack((res1, res2, emi1, emi2, inc1, inc2, pha1, pha2, azi1, azi2))
	os.chdir(pathsave)
	np.savetxt(pathsave + "newdata.txt", output, header="res1;res2;emi1;emi2;inc1;inc2;pha1;pha2;azi1;azi2", delimiter=";", fmt = ["%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e","%1.6e"])

'''
*********************************************************************************
'''	

def get_union(infile, IMAGE1_F_mod, IMAGE2_F_mod):

	'''
	empty layers might lead to some errors
	infile2 = "NAC_FILTERED"
	get_union(infile2, IMAGE1_F_mod, IMAGE2_F_mod)
	'''
	
	n = len(IMAGE1_F_mod)
	
	int_top = np.ones(n)
	int_bottom = np.ones(n)
	int_left = np.ones(n)
	int_right = np.ones(n)
	
	for iy, id in np.ndenumerate(IMAGE1_F_mod):
	
		iyy = iy[0]
	
		query1 = "productid = '" + IMAGE1_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query1)
		arcpy.CopyFeatures_management(infile,"test_query1")


		query2 = "productid = '" + IMAGE2_F_mod[iyy] + "'"
		arcpy.SelectLayerByAttribute_management(infile, "NEW_SELECTION", query2)
		arcpy.CopyFeatures_management(infile,"test_query2")

		inFeatures = ["test_query1","test_query2"]
		outFeatures = "intersect"
		arcpy.Intersect_analysis(inFeatures, outFeatures)

		rows = arcpy.da.SearchCursor("intersect", ["Shape@"])

		for row in rows:
			feat = row[0]
			extent = feat.extent
			int_top[iyy] = extent.YMax
			int_bottom[iyy] = extent.YMin
			int_left[iyy] = extent.XMin
			int_right[iyy] = extent.XMax
			#int_top.append(extent.YMax)
			#int_bottom.append(extent.YMin)
			#int_left.append(extent.XMin)
			#int_right.append(extent.XMax)


		arcpy.Delete_management("test_query1")
		arcpy.Delete_management("test_query2")
		arcpy.Delete_management("intersect")
	
	return (int_top, int_bottom, int_left, int_right)

# after that I can crop the images with maptrim in ISIS3 
(int_top, int_bottom, int_left, int_right) = get_union(infile, IMAGE1_F_mod, IMAGE2_F_mod)
output = np.column_stack((int_top, int_bottom, int_left, int_right))
os.chdir(pathsave)
np.savetxt("maptrim2.txt", output, header="int_top;int_bottom;int_left;int_right", delimiter=";", fmt = ["%1.6e","%1.6e","%1.6e","%1.6e"])

'''
*********************************************************************************
'''	
def degree2rad(degree):
	rad = degree * (np.pi/180.0)
	return rad
	
'''
*********************************************************************************
'''	

def rad2degree(rad):
	degree = rad * (180.0/np.pi)
	return degree
'''
*********************************************************************************
'''		
def lommel_Seeliger(inc,emi):

	inc_rad = degree2rad(inc)
	emi_rad = degree2rad(emi)
	
	var = np.cos(inc_rad) / (np.cos(emi_rad) + np.cos(inc_rad))
	
	return var

'''
*********************************************************************************
'''	

'''
I have to download the images I am interested in, possibly automatically and then to process them with ISIS and a in-house bash script
'''

def raster_lommel_Seeliger(raster_path, new_raster_name, inc, emi): 
	raster = arcpy.Raster(raster_path)
	var = lommel_Seeliger(inc, emi)
	outraster = Divide(raster,var)
	outraster.save(new_raster_name)
	print "DONE"
	
'''
*********************************************************************************
'''	
raster_path = "Y:/nilscp/Moon/PHASE/"
os.chdir(raster_path)
raster_files = glob.glob("*tif")
inc_example = [51.72,50.77]
emi_example = [27.07,8.95]
pha_example = [50.77,25.29]

for iu, file in enumerate(raster_files):
	path_raster_files = raster_path + file
	new_raster_name = raster_path + file.split(".t")[0] + "_lommel.img"	
	raster_lommel_Seeliger(path_raster_files, new_raster_name, inc_example[iu], emi_example[iu])
	
'''
*********************************************************************************
'''
def phase_ratio_imagery(raster1, raster2, pha1, pha2, new_raster): 
	rast1 = arcpy.Raster(raster1)
	rast2 = arcpy.Raster(raster2)
	
	if pha1 < pha2:
		outraster = Divide(rast1,rast2)
		outraster.save(new_raster)
	else:
		outraster = Divide(rast2,rast1)
		outraster.save(new_raster)
		
	print "The phase-ratio imagery has been processed"
	
'''
*********************************************************************************
'''
raster1 = "Y:/nilscp/Moon/PHASE/M144748567RE_final_lommel.img"
raster2 = "Y:/nilscp/Moon/PHASE/M144755352RE_final_lommel.img"
pha1 = pha_example[0]
pha2 = pha_example[1]
if pha1 < pha2:
	new_raster = raster1.split('_lommel')[0] + '_' +(raster2.split('PHASE/')[1]).split('_lommel')[0] + '_phase_ratio3.img'
else:
	new_raster = raster2.split('_lommel')[0] + '_' +(raster1.split('PHASE/')[1]).split('_lommel')[0] + '_phase_ratio3.img'

phase_ratio_imagery(raster1, raster2, pha1, pha2, new_raster)