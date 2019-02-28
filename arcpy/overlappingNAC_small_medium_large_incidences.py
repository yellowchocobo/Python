import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np


'''
We are writing a routine to select images that overlap, have similar incidences,
phase angles, ground azimuth but different phases
'''	

'''
**************************************************************************************************
'''

def select_NAC(infile, area_bnd, output_name, inc_min, inc_max):
    
    '''
	# select nac images for a region of interest and an interval of incidences
	'''
    
    arcpy.MakeFeatureLayer_management(infile, "infile_lyr")
    
    arcpy.SelectLayerByLocation_management("infile_lyr", "INTERSECT", area_bnd)
    
    # search nac images with incidences between 30 and 70 degrees (add to selection)
    arcpy.SelectLayerByAttribute_management("infile_lyr", "REMOVE_FROM_SELECTION", '"ctr_incang" < ' + str(inc_min))
    arcpy.SelectLayerByAttribute_management("infile_lyr", "REMOVE_FROM_SELECTION", '"ctr_incang" > ' + str(inc_max))
    
    # make a layer of the selection
    arcpy.CopyFeatures_management("infile_lyr", output_name)
    
    # delete
    arcpy.Delete_management("infile_lyr")


'''
**************************************************************************************************
'''

def calculate_newly_added_area(pathdab, infile, infile_area, fresh_crater_area, spatialReference, spatialReference_new):
    
    '''
    infile = "nac_selection_20_40_proj"
    infile_area = "area_bnd"
    pathdab = "Y:/nilscp/GeologyMoon/FRESH_IMPACT_WILLIAMS_2018/database.gdb/"
    
    rows = arcpy.da.SearchCursor(area_bnd, ["SHAPE@"])
        
    for row in rows:
        feat = row[0]
        fresh_crater_area = feat.area
        
    ((NAC_selected, edr_source_selected, areac) = 
      calculate_newly_added_area(pathdab, infile, infile_area, fresh_crater_area))
    
    '''
    # if no images have been detected for the specified incidence
    if arcpy.management.GetCount(infile)[0] == '0':
        
        # we do nothing
        return ([], [], 0.0)
            
        
    else:
        
        # we define the NAC selected, the EDR source of selected images (to download)
        # and the area_covered by NAC images (area_merged2), which is 0 at the start
        NAC_selected = []
        edr_source_selected = []
        area_covered = 0
               
        fields = ["SHAPE@", "productid", "edr_source", "COVERED_AREA"]
        
        # create a layer to avoid error due to the query
        arcpy.MakeFeatureLayer_management(infile, infile + "_lyr")
        
        with arcpy.da.UpdateCursor(infile + "_lyr", fields) as cursor:
        
            # we loop through all the candidates NAC
            for row in cursor:		
                #feat = row[0]
                #area_covered = feat.area
                NACname = str(row[1])
                EDRname = str(row[2])
                
                # we use a query to select the NAC image of interest
                query = "productid = '" + NACname + "'"
                arcpy.SelectLayerByAttribute_management(infile + "_lyr", "NEW_SELECTION", query)
                #print (query)
                            
                # and a copy of this file is done
                arcpy.CopyFeatures_management(infile + "_lyr", "merged_area_shp")
                              
            	# the intersect of the NAC and all selected NAC are done
                arcpy.Intersect_analysis(["merged_area_shp", infile_area], "temporary_intersect")
            
                # calculate the area that is overlapping
            			
    			# the area of the intersect is saved
                area_intersect = 0.0
                with arcpy.da.SearchCursor("temporary_intersect", ["SHAPE@"]) as rows_intersect:
                         
                    # In case it is divided in several polygons?
                    for row_i in rows_intersect:
                        feat_i = row_i[0]
                        area_intersect += feat_i.area
                    
                # then we save that to fieldname_area
                row[3] = (area_intersect / (fresh_crater_area))
                #print ((area_intersect / (fresh_crater_area)))
                
                #update row
                cursor.updateRow(row)
                
                # delete both temporary shapefiles
                arcpy.Delete_management("temporary_intersect")
                arcpy.Delete_management("merged_area_shp")
            
        # save to new file but sorted this time (depending on the area covered)
        arcpy.SelectLayerByAttribute_management(infile + "_lyr", "CLEAR_SELECTION") #otherwise it keeps from previous selection
        
        # infile + "_sorted" is not in local reference
        arcpy.Sort_management(infile, infile + "_sorted", [["COVERED_AREA", "DESCENDING"]])
        
        # create a layer to avoid error due to the query
        arcpy.MakeFeatureLayer_management(infile + "_sorted", infile + "_sorted" + "_lyr")
           
        # define the shapefile we are going to loop through and the columns that
        # will be read (we use the sorted one now)
        with arcpy.da.SearchCursor(infile + "_sorted" + "_lyr", fields) as rows:    
    
            for row in rows:		
                #feat = row[0]
                #area_covered = feat.area # area of the NAC image
                NACname = str(row[1])
                EDRname = str(row[2])
            	
                # if not completely covered (continue)
                if (area_covered / (fresh_crater_area)) < 0.99:
                
                    # if nothing is selected, this is the first time it is run
                    if len(NAC_selected) == 0:
                
                        # NAC images and EDR source strings are added (because first time)
                        NAC_selected.append(NACname)
                        edr_source_selected.append(EDRname)
                        
                        # clear selection
                        arcpy.SelectLayerByAttribute_management(infile + "_sorted" + "_lyr", "CLEAR_SELECTION")
    
                        # we use a query to select the NAC image of interest
                        query = "productid = '" + NACname + "'"
                        arcpy.SelectLayerByAttribute_management(infile + "_sorted" + "_lyr", "NEW_SELECTION", query)
                
                        # and a copy of this file is done to merged_area_all (because first time)
                        arcpy.CopyFeatures_management(infile + "_sorted" + "_lyr", "temporary_NAC_selection")
                        
                        # but I should still take the intersect? Everything in local
                        arcpy.Intersect_analysis(["temporary_NAC_selection", infile_area], "covered_area")
                                     
                        # This is the first time so all the area that is covered can be added
                        # and it has been calculated based on local projection
                        # this might be wrong
                        # reset area_covered to 0?
                        area_covered = 0
                        
                        # the area of the intersect is saved
                        with arcpy.da.SearchCursor("covered_area", ["SHAPE@"]) as rowsm:
                        
                            # loop through the new addition                    
                            for row_m in rowsm:
                                feat_m = row_m[0]
                                area_covered += feat_m.area
                        
                        # We delete temporary and temporary_intersect
                        arcpy.Delete_management("temporary_NAC_selection")
                                                          
                    else:
                        
                        # need to unselect previous selection
                        arcpy.SelectLayerByAttribute_management(infile + "_sorted" + "_lyr", "CLEAR_SELECTION")
    			
                        # copy area covered to a new shapefile we can work with
                        arcpy.Copy_management(pathdab + "/covered_area", pathdab + "/covered_area_tmp")
                
                        # the selected NAC is copied into a layer (it has been done on the sorted)
                        query = "productid = '" + NACname + "'"
                        arcpy.SelectLayerByAttribute_management(infile + "_sorted" + "_lyr", "NEW_SELECTION", query)
                        arcpy.CopyFeatures_management(infile + "_sorted" + "_lyr", "temporary_NAC_selection")
                        
                        # I should actually merge merged_area_all and the new temporary NAC
                        #arcpy.Merge_management(["covered_area_tmp", "temporary_NAC_selection"], "temporary_merged")
                        
                        # or union
                        arcpy.Union_analysis(["covered_area_tmp", "temporary_NAC_selection"], "temporary_merged")
                
                        # the intersect of the NAC of interest and the already covered
                        # area is checked (all previously selected images)
                        
                        # should overwrite covered_area
                        arcpy.Intersect_analysis(["temporary_merged", infile_area], "covered_area")
                                                
                        # reset area_covered to 0?
                        area_covered = 0
                        
                        # the area of the intersect is saved
                        with arcpy.da.SearchCursor("covered_area", ["SHAPE@"]) as rowsm:
                                               
                            # loop through the new addition                    
                            for row_m in rowsm:
                                feat_m = row_m[0]
                                area_covered += feat_m.area
                            
                        # we add the name of the NAC and edr source in the lists
                        NAC_selected.append(NACname)
                        edr_source_selected.append(EDRname)
                        
                        # We delete temporary and temporary_intersect
                        arcpy.Delete_management("temporary")
                        arcpy.Delete_management("temporary_merged")

                else:
                    # if we have a coverage of more than 99% than we don't have to 
                    # add any other NAC images in the list
                    None
        
        return (NAC_selected, edr_source_selected, (area_covered / fresh_crater_area))              


'''
**************************************************************************************************
'''

def main(pathf, pathdab, infile, restart, ixstop, fname):
    
    '''
    pathf = "Y:/nilscp/GeologyMoon/FRESH_IMPACT_WILLIAMS_2018/"
    infile = pathf + 'XYdownloadSupplement.shp'
    pathdab = "Y:/nilscp/GeologyMoon/FRESH_IMPACT_WILLIAMS_2018/database.gdb/"
    fname = '0_200'
    '''
    
    bufferField = "BUFFER_TXT"
    sideType = "FULL"
    endType = "ROUND"
    dissolveType = "NONE"

    # path to folder
    os.chdir(pathf)
    
    # Set overwrite option
    arcpy.env.overwriteOutput = True
    
    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("Spatial")
    
    # define paths and workspace (I need to create the gdb at some points)
    env.workspace = env.scratchWorkspace = pathdab
        
    # extract the centers of craters (OK regardless of the projection)
    arcpy.FeatureToPoint_management(infile, 'CENTER', "CENTROID")
    
    infile2 = pathdab + 'CENTER'								
    								
    # crater name and buffer extent
    fieldname1 = arcpy.ValidateFieldName("CRATER_ID")
    fieldname2 = arcpy.ValidateFieldName("BUFFER_TXT")
    
    # add fields
    arcpy.AddField_management(infile, fieldname1, "TEXT","","",30)
    
    # get the number of rows in infile
    n = int(arcpy.GetCount_management(infile2)[0])
    
    # prepare empty arrays
    diam = np.ones(n)
    crater_id = np.chararray(n, itemsize=30)
    buffer_txt = np.chararray(n, itemsize=30)
    
    ## some constant variables to define
    if restart:
        areac = np.loadtxt(pathf + "numbers" + fname + ".txt", delimiter=";")
        areaa_20_40 = areac[:,0]#np.concatenate((areac[:,0], np.zeros(1781))) # work around (only need to do that once)
        areaa_40_55 = areac[:,1] #np.concatenate((areac[:,1], np.zeros(1781)))
        areaa_55_80 = areac[:,2] #np.concatenate((areac[:,2], np.zeros(1781)))
        
        NAC_selected_tmp_20_40 = np.genfromtxt(pathf + "NAC_i20_40_" + fname + ".txt", dtype='str')
        NAC_selected_tmp_40_55 = np.genfromtxt(pathf + "NAC_i40_55_" + fname + ".txt", dtype='str')
        NAC_selected_tmp_55_80 = np.genfromtxt(pathf + "NAC_i55_80_" + fname + ".txt", dtype='str')
        edr_selected_tmp_20_40 = np.genfromtxt(pathf + "EDR_i20_40_" + fname + ".txt", dtype='str')
        edr_selected_tmp_40_55 = np.genfromtxt(pathf + "EDR_i40_55_" + fname + ".txt", dtype='str')
        edr_selected_tmp_55_80 = np.genfromtxt(pathf + "EDR_i55_80_" + fname + ".txt", dtype='str')
        
        NAC_selected_all_20_40 = []
        edr_selected_all_20_40 = []
        NAC_selected_all_40_55 = []
        edr_selected_all_40_55 = []
        NAC_selected_all_55_80 = []
        edr_selected_all_55_80 = []
        
        # needs to be list
        for ip, var in np.ndenumerate(NAC_selected_tmp_20_40):
            NAC_selected_all_20_40.append(NAC_selected_tmp_20_40[ip])
            edr_selected_all_20_40.append(edr_selected_tmp_20_40[ip])
            
        for ip, var in np.ndenumerate(NAC_selected_tmp_40_55):            
            NAC_selected_all_40_55.append(NAC_selected_tmp_40_55[ip])
            edr_selected_all_40_55.append(edr_selected_tmp_40_55[ip])
            
        for ip, var in np.ndenumerate(NAC_selected_tmp_55_80):            
            NAC_selected_all_55_80.append(NAC_selected_tmp_55_80[ip])
            edr_selected_all_55_80.append(edr_selected_tmp_55_80[ip])
    else:
        NAC_selected_all_20_40 = []
        edr_selected_all_20_40 = []
        NAC_selected_all_40_55 = []
        edr_selected_all_40_55 = []
        NAC_selected_all_55_80 = []
        edr_selected_all_55_80 = []
        areaa_20_40 = np.zeros(n)
        areaa_40_55 = np.zeros(n)
        areaa_55_80 = np.zeros(n)
    
    
    # we add info about the name of the craters here
    
    with arcpy.da.UpdateCursor(infile, ["Diameter", "CRATER_ID"]) as cursor:    	
        ix = 0
        for row in cursor:
            a = 'crater' + str(int(ix)).zfill(4)
            buffer_value = np.round((row[0]) * 10.0, decimals=4)
            b = str(buffer_value) + ' Meters'	
            row[1] = a
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
            
    #arcpy.AddField_management(infile2, fieldname3, "DOUBLE")
    #arcpy.AddField_management(infile2, fieldname4, "DOUBLE")
    #arcpy.AddField_management(infile2, fieldname5, "DOUBLE")
    
    # Make a layer from the feature class
    arcpy.MakeFeatureLayer_management("CENTER", "CENTER_lyr")
           
    with arcpy.da.UpdateCursor("CENTER_lyr", ["Shape@", "Lon", "Lat"]) as cursor:
        ix = 0
        
        for row in cursor:
            
            #print index
            #print (ix)
            
            # could reset here too (just to be sure)
            if ix > 0:
                arcpy.SelectLayerByAttribute_management("CENTER_lyr", "CLEAR_SELECTION")
            
            else:
                None
            
            # if first time or timesteps over last restarted 
            if ((restart & (ix > ixstop)) | (restart == False)):
                #query selection CENTER         
                query = "CRATER_ID = '" + crater_id[ix] + "'"
                print (query)
                arcpy.SelectLayerByAttribute_management("CENTER_lyr", "NEW_SELECTION", query)
                #print ("YESx2")
                
                # make a layer of the selection
                arcpy.CopyFeatures_management("CENTER_lyr", "CENTER_TMP")
                
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
                arcpy.Buffer_analysis("CENTER_PROJ", "miniarea_TMP", bufferField, sideType, endType, dissolveType)
                
                # run feature to envelope tool
                arcpy.FeatureEnvelopeToPolygon_management(pathdab + "miniarea_TMP",
        												  pathdab + "miniarea_square",
        												  "SINGLEPART")
                
                
                # get area of the 10 diameters squared polygon
                area_bnd = "miniarea_square"
                with arcpy.da.SearchCursor(area_bnd, ["SHAPE@"]) as rows:
                
                    # get the area of the 10 diameters squared polygon
                    for row in rows:
                        feat = row[0]
                        fresh_crater_area = feat.area
                    
                # select nac images (left and right satellite images) between x and y degrees of incidence
                select_NAC("nac_all", area_bnd, "nac_selection_20_40", 20., 40.)
                select_NAC("nac_all", area_bnd, "nac_selection_40_55", 40., 55.)
                select_NAC("nac_all", area_bnd, "nac_selection_55_80", 55., 80.)
                
                # define covered area here!
                fieldname_area = arcpy.ValidateFieldName("COVERED_AREA")
                
                for shpfile in ["nac_selection_20_40", "nac_selection_40_55", "nac_selection_55_80"]:
                    if len(arcpy.ListFields(shpfile,"COVERED_AREA"))>0:
                        None
                    else:
                        arcpy.AddField_management(shpfile, fieldname_area, "DOUBLE") 
                
                #I could project the nac_selection here
                arcpy.Project_management(in_dataset="nac_selection_20_40", out_dataset="nac_selection_20_40_proj", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")	
                arcpy.Project_management(in_dataset="nac_selection_40_55", out_dataset="nac_selection_40_55_proj", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")	
                arcpy.Project_management(in_dataset="nac_selection_55_80", out_dataset="nac_selection_55_80_proj", out_coor_system= spatialReference_new, transform_method="", in_coor_system=spatialReference, preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")	
                
                
                # selected NAC + EDR + Area covered
                #print ("last function")
                (NAC_selected_20_40, edr_source_selected_20_40, areac_20_40) = calculate_newly_added_area(pathdab, "nac_selection_20_40_proj", area_bnd, fresh_crater_area, spatialReference, spatialReference_new)
                #print ("2")
                (NAC_selected_40_55, edr_source_selected_40_55, areac_40_55) = calculate_newly_added_area(pathdab, "nac_selection_40_55_proj", area_bnd, fresh_crater_area, spatialReference, spatialReference_new)
                #print ("3")
                (NAC_selected_55_80, edr_source_selected_55_80, areac_55_80) = calculate_newly_added_area(pathdab, "nac_selection_55_80_proj", area_bnd, fresh_crater_area, spatialReference, spatialReference_new)
                
                
                #print (ix)
                arcpy.Delete_management("miniarea_square")
                arcpy.Delete_management("miniarea_TMP")
                arcpy.Delete_management("CENTER_PROJ")
                arcpy.Delete_management("CENTER_TMP")
                arcpy.Delete_management("nac_selection_20_40")
                arcpy.Delete_management("nac_selection_20_40_proj")
                arcpy.Delete_management("nac_selection_40_55")
                arcpy.Delete_management("nac_selection_40_55_proj")
                arcpy.Delete_management("nac_selection_55_80")
                arcpy.Delete_management("nac_selection_55_80_proj")            
                
                #NAC_selected_tmp = NAC_selected_20_40 + NAC_selected_40_55 + NAC_selected_55_80
                #edr_selected_tmp = edr_source_selected_20_40 + edr_source_selected_40_55 + edr_source_selected_55_80
                
                NAC_selected_all_20_40 = NAC_selected_all_20_40 + NAC_selected_20_40
                edr_selected_all_20_40 = edr_selected_all_20_40 + edr_source_selected_20_40
                areaa_20_40[ix] =  areac_20_40
                
                NAC_selected_all_40_55 = NAC_selected_all_40_55 + NAC_selected_40_55
                edr_selected_all_40_55 = edr_selected_all_40_55 + edr_source_selected_40_55
                areaa_40_55[ix] =  areac_40_55
                
                NAC_selected_all_55_80 = NAC_selected_all_55_80 + NAC_selected_55_80
                edr_selected_all_55_80 = edr_selected_all_55_80 + edr_source_selected_55_80
                areaa_55_80[ix] =  areac_55_80
                
                # save every 5 timesteps
                if (ix % 25 == 0):
                    
                    
                    output_nbr = np.column_stack((np.array(areaa_20_40), np.array(areaa_40_55), np.array(areaa_55_80)))
    
                    np.savetxt("numbers0_" + str(int(ix)) + ".txt", output_nbr, delimiter=";")
                    np.savetxt("NAC_i20_40_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_20_40)), delimiter=";",fmt="%s")
                    np.savetxt("NAC_i40_55_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_40_55)), delimiter=";",fmt="%s")
                    np.savetxt("NAC_i55_80_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_55_80)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i20_40_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_20_40)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i40_55_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_40_55)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i55_80_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_55_80)), delimiter=";",fmt="%s")
                    
                elif (ix == 2281):
                    
                    output_nbr = np.column_stack((np.array(areaa_20_40), np.array(areaa_40_55), np.array(areaa_55_80)))
    
                    np.savetxt("numbers0_" + str(int(ix)) + ".txt", output_nbr, delimiter=";")
                    np.savetxt("NAC_i20_40_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_20_40)), delimiter=";",fmt="%s")
                    np.savetxt("NAC_i40_55_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_40_55)), delimiter=";",fmt="%s")
                    np.savetxt("NAC_i55_80_0_" + str(int(ix)) + ".txt", np.array((NAC_selected_all_55_80)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i20_40_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_20_40)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i40_55_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_40_55)), delimiter=";",fmt="%s")
                    np.savetxt("EDR_i55_80_0_" + str(int(ix)) + ".txt", np.array((edr_selected_all_55_80)), delimiter=";",fmt="%s")
                    
                    
            else:
                None
                
            ix = ix + 1
            
    	
    return (NAC_selected_all_20_40, edr_selected_all_20_40, areaa_20_40, 
            NAC_selected_all_40_55, edr_selected_all_40_55, areaa_40_55, 
            NAC_selected_all_55_80, edr_selected_all_55_80, areaa_55_80) 
'''
**************************************************************************************************
'''

pathf = "C:/Users/nilscp/Desktop/testarcpy/"
infile = pathf + 'XYdownloadSupplement.shp'
pathdab = "C:/Users/nilscp/Desktop/testarcpy/database.gdb/"
restart = True
ixstop = 2275
fname = '0_2275'

(NAC_selected_all_20_40, edr_selected_all_20_40, areac_20_40, 
 NAC_selected_all_40_55, edr_selected_all_40_55, areac_40_55, 
 NAC_selected_all_55_80, edr_selected_all_55_80, areac_55_80) = main(pathf, pathdab, infile, restart, ixstop, fname)


