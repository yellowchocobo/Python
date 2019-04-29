# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:00:14 2019

@author: nilscp
"""

infile = "location_overlapDTM_completely_within"

import arcpy
from arcpy import env
import glob, os
from arcpy.sa import *
import numpy as np

env.workspace = env.scratchWorkspace = "D:/NAC_DTM/NAC_AMES/arcpy/database.gdb/"


arcpy.AddField_management(infile, "abspath", "TEXT","","",100)

data = np.genfromtxt("D:/ANALYSIS/NAC_DTM/" + "abspath.txt", dtype=None)

with arcpy.da.UpdateCursor(infile, ["abspath"]) as cursor:
	ix = 0
	for row in cursor:
		row[0] = data[ix]
		cursor.updateRow(row)
		ix = ix + 1