#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:48:03 2019

@author: nilscp
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

res_SLDEM2013_coldspots = "/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/results/data_cp2.csv"
cpcrater_id = "/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/results/crater_id_cp.csv"
add_data = "/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/results/mare_coldspots.txt" # data about if it's on the MARE or not
add_data2 = "/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/results/geology_coldspots2.csv"
add_data3 = "/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/results/flag_quality.csv"

# read data
data_cp = pd.read_csv(res_SLDEM2013_coldspots, delimiter=",")
data_mare = np.loadtxt(add_data, delimiter = ",", skiprows=1)
crater_id_cp = np.loadtxt(cpcrater_id, comments="#", dtype='str') 
data_geology = pd.read_csv(add_data2, delimiter = ",",dtype='str')
flag_quality = np.loadtxt(add_data3, delimiter=",", skiprows=1)
# need to load the data for Mars (Watters et al.) and ask him how the data can be filtered

# extract data from mare
lon = data_mare[:,0]
lat = data_mare[:,1]
mare = data_mare[:,2]
flag = flag_quality[:,0]
quality = flag_quality[:,1]


# create list of complete cp
crater_id_complete_list = []

for i in range(2282):
    
    crater_id_complete_list.append('cpcrater' + str(i).zfill(4))

# byte problem

crater_id_lst = []

for i, c in enumerate(crater_id_cp):
    crater_id_lst.append(crater_id_cp[i][2:-1])
    
# add flag column in addition
flagG = []

for i in crater_id_lst:
    ix = int(i.split(sep="cpcrater")[-1])
    flagG.append(flag[ix])
    
# add a column with crater_id
data_cp['CRATER_ID'] = crater_id_lst
data_cp['FLAG'] = flagG

# only picking the good data
ndata = data_cp[data_cp.QUALITY == 1]
ndata1 = ndata[ndata.FLAG == 1]

ndata2 = ndata1[ndata1.mdiam >= 20.*7.4031617] # take only craters described by at least 20 points

# extract mare values for selection
mare_selection = []
lon_selection = []
lat_selection = []
epoch = []
unit = []
unitS = []
majorG = []


for i in ndata2.CRATER_ID:
    ix = int(i.split(sep="cpcrater")[-1])
    mare_selection.append(mare[ix])
    lon_selection.append(lon[ix])
    lat_selection.append(lat[ix])
    epoch.append(data_geology.EpochG[ix])
    unit.append(data_geology.UnitNameG[ix])
    unitS.append(data_geology.UnitSymbol[ix])
    majorG.append(data_geology.MajorGroup[ix])

'''
***********************************************************************
'''
    
# add column mare
ndata2['MARE'] = mare_selection[:]
ndata2['LON'] = lon_selection[:]
ndata2['LAT'] = lat_selection[:]
ndata2['epoch'] = epoch[:]
ndata2['unit'] = unit[:]
ndata2['unitS'] = unitS[:]
ndata2['majorG'] = majorG[:]

'''
***********************************************************************
'''
# select the 10 smallest and largest dD for mare/highlands 
idx_smallest = np.where(np.nanargsort(ndata2.dD) < 10)
idx = np.argpartition(ndata2.dD.values, 0.05)
 
# 100 smallest
idxs = np.argsort(ndata2.dD.values)[:200]
i = -5
ndata2.CRATER_ID.values[idxs][i], ndata2.mdiam.values[idxs][i], ndata2.mdepth.values[idxs][i],ndata2.mh.values[idxs][i],ndata2.dD.values[idxs][i], ndata2.MARE.values[idxs][i], ndata2.LON.values[idxs][i], ndata2.LAT.values[idxs][i]

# 20 largest
idxl = np.argsort(ndata2.dD.values)[-20:]
i = 1
ndata2.CRATER_ID.values[idxl][i], ndata2.mdiam.values[idxl][i], ndata2.mdepth.values[idxl][i],ndata2.mh.values[idxl][i],ndata2.dD.values[idxl][i], ndata2.MARE.values[idxl][i], ndata2.LON.values[idxl][i], ndata2.LAT.values[idxl][i]

'''
***********************************************************************
'''
# copernican against the rest
ndata2.groupby('epoch').count()


plt.plot(ndata2.mdiam[ndata2.epoch != 'Nectarian System'], ndata2.dD[ndata2.epoch != 'Nectarian System'], 'bo')
plt.plot(ndata2.mdiam[ndata2.epoch == 'Imbrian System'], ndata2.dD[ndata2.epoch == 'Imbrian System'], 'go')
plt.plot(ndata2.mdiam[ndata2.epoch == 'Nectarian System'], ndata2.dD[ndata2.epoch == 'Nectarian System'], 'ro')
plt.plot(ndata2.mdiam[ndata2.epoch == 'Eratosthenian System'], ndata2.dD[ndata2.epoch == 'Eratosthenian System'], 'bo')
plt.plot(ndata2.mdiam[ndata2.epoch == 'pre-Nectarian System'], ndata2.dD[ndata2.epoch == 'pre-Nectarian System'], 'yo')

'''
***********************************************************************
'''

# material
ndata2.groupby('majorG').count()
plt.plot(ndata2.mdiam[ndata2.majorG != 'Crater Materials'], ndata2.dD[ndata2.majorG != 'Crater Materials'], 'ro')

plt.plot(ndata2.mdiam[ndata2.majorG == 'Crater Materials'], ndata2.dD[ndata2.majorG == 'Crater Materials'], 'bo') #usually deeper
plt.plot(ndata2.mdiam[ndata2.majorG == 'Dark Materials'], ndata2.dD[ndata2.majorG == 'Dark Materials'], 'go') #usually shallower
plt.plot(ndata2.mdiam[ndata2.majorG == 'Other Terra Materials'], ndata2.dD[ndata2.majorG == 'Other Terra Materials'], 'yo')
plt.plot(ndata2.mdiam[ndata2.majorG == 'Mare and Other Materials'], ndata2.dD[ndata2.majorG == 'Mare and Other Materials'], 'ro')
plt.plot(ndata2.mdiam[ndata2.majorG == 'Orientale Group'], ndata2.dD[ndata2.majorG == 'Orientale Group'], 'co')

'''
***********************************************************************
'''

# same as watters
# depth and diameter figure 5.
plt.loglog(ndata2.mdiam, ndata2.mhr - ndata2.mdepth, 'ro')
plt.loglog(ndata2.mdiam[ndata2.dD > 0.05], ndata2.mhr[ndata2.dD > 0.05] - ndata2.mdepth[ndata2.dD > 0.05], 'ro')

#plt.loglog(ndata2.mdiam[ndata2.mdiam > 300.0], ndata2.mhr[ndata2.mdiam > 300.0] - ndata2.mdepth[ndata2.mdiam > 300.0], 'ro')
DD = np.linspace(10,1000)
dd = 0.205*(DD**1.012)
dd2 = 0.10*DD
dd3 = 0.05*DD

plt.loglog(DD,dd,'b')
plt.loglog(DD,dd2,'b')
plt.loglog(DD,dd3,'b')

fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

HHH = ndata2.mhr[ndata2.dD > 0.05] - ndata2.mdepth[ndata2.dD > 0.05]
H25 = ndata2.hr25[ndata2.dD > 0.05] - ndata2.depth25[ndata2.dD > 0.05]
H75 = ndata2.hr75[ndata2.dD > 0.05] - ndata2.depth75[ndata2.dD > 0.05]

markers, caps, bars = plt.errorbar(ndata2.mdiam[ndata2.dD > 0.05], ndata2.mhr[ndata2.dD > 0.05] - ndata2.mdepth[ndata2.dD > 0.05], xerr=np.array([ndata2.mdiam[ndata2.dD > 0.05] - ndata2.diam25[ndata2.dD > 0.05],ndata2.diam75[ndata2.dD > 0.05] - ndata2.mdiam[ndata2.dD > 0.05]]), yerr=np.array([HHH - H25, H75 - HHH]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]


###

fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

HHH = ndata2.mhr - ndata2.mdepth
H25 = ndata2.hr25 - ndata2.depth25
H75 = ndata2.hr75 - ndata2.depth75

markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.mhr - ndata2.mdepth, xerr=np.array([ndata2.mdiam - ndata2.diam25,ndata2.diam75 - ndata2.mdiam]), yerr=np.array([HHH - H25, H75 - HHH]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]


# histogram of dD
plt.hist(ndata2.dD.values, bins=20)

plt.hist(ndata2.mcse.values, bins=20)
'''
***********************************************************************
'''

# depth and diameter additional

fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')

# not sure about that as dD25 and dD75 can smaller or larger
dD25 = (ndata2.hr25 - ndata2.depth25) / ndata2.diam25
dD75 = (ndata2.hr75 - ndata2.depth75) / ndata2.diam75

markers, caps, bars = plt.errorbar(ndata2.mdiam[ndata2.mdiam > 100.0], ndata2.dD[ndata2.mdiam > 100.0],fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]

'''
***********************************************************************
'''

# Cavity shape exponent figure 6.
plt.semilogx(ndata2.mdiam, ndata2.mcse,"ko")

# good -∕+ 
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.mcse, xerr=np.array([ndata2.mdiam - ndata2.diam25,ndata2.diam75 - ndata2.mdiam]), yerr=np.array([ndata2.mcse - ndata2.cse_25,ndata2.cse_75 - ndata2.mcse]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
plt.ylim((0.4,3.6))

'''
***********************************************************************
'''

# rim decay length
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.mcrdl, yerr=np.array([ndata2.mcrdl - ndata2.crdl_25,ndata2.crdl_75 - ndata2.mcrdl]),fmt='o')
markers2, caps2, bars2 = plt.errorbar(ndata2.mdiam, ndata2.mfrdl, yerr=np.array([ndata2.mfrdl - ndata2.frdl_25,ndata2.frdl_75 - ndata2.mfrdl]),fmt='ro')

[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
[bar.set_alpha(0.25) for bar in bars2]
[cap.set_alpha(0.25) for cap in caps2]
crdl_synth = 0.124*(DD**0.838)
plt.loglog(DD, crdl_synth,'g')

'''
***********************************************************************
'''

# radius of curvature of upper cavity
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.mrupcw, yerr=np.array([ndata2.mrupcw - ndata2.rupcw_25,ndata2.rupcw_75 - ndata2.mrupcw]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]

'''
***********************************************************************
'''

# middle cavity versus diameter vs flank slope - Figure 12
fig, ax = plt.subplots()
markers, caps, bars = plt.errorbar(ndata2.mfsa, ndata2.m_mcw, xerr=np.array([ndata2.mfsa - ndata2.fsa_25,ndata2.fsa_75 - ndata2.mfsa]), yerr=np.array([ndata2.m_mcw - ndata2.mcw_25,ndata2.mcw_75 - ndata2.m_mcw]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
plt.ylim((0,40))

'''
***********************************************************************
'''

# middle cavity versus diameter (there are some negative values! Need to check those craters) - Figure 13
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.m_mcw, yerr=np.array([ndata2.m_mcw - ndata2.mcw_25,ndata2.mcw_75 - ndata2.m_mcw]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
plt.ylim((0,40))

'''
***********************************************************************
'''

# upper cavity versus diameter (there are some negative values! Need to check those craters) - Figure 13
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.m_ucw, yerr=np.array([ndata2.m_ucw - ndata2.ucw_25,ndata2.ucw_75 - ndata2.m_ucw]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]
plt.ylim((0,20))

'''
***********************************************************************
'''

# lower rim span (there are some negative values! Need to check those craters) - Figure 14
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.mlrs, yerr=np.array([ndata2.mlrs - ndata2.lrs_25,ndata2.lrs_75 - ndata2.mlrs]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]

'''
***********************************************************************
'''

# Upper rim span (there are some negative values! Need to check those craters) - Figure 14
fig, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
markers, caps, bars = plt.errorbar(ndata2.mdiam, ndata2.murs, yerr=np.array([ndata2.murs - ndata2.urs_25,ndata2.urs_75 - ndata2.murs]),fmt='o')
[bar.set_alpha(0.25) for bar in bars]
[cap.set_alpha(0.25) for cap in caps]

plt.scatter(ndata2.mdiam, ndata2.dD, c=ndata2.LAT)
plt.colorbar()

plt.scatter(ndata2.mdiam, ndata2.dD, c=ndata2.LON)
plt.colorbar()

plt.plot(ndata2.mdiam[ndata2.MARE == 0], ndata2.dD[ndata2.MARE == 0], 'ro')
plt.plot(ndata2.mdiam[ndata2.MARE == 1], ndata2.dD[ndata2.MARE == 1], 'bo')


plt.plot(ndata2.m_mcw[ndata2.MARE == 0], ndata2.dD[ndata2.MARE == 0], "ro")
plt.plot(ndata2.m_mcw[ndata2.MARE == 1], ndata2.dD[ndata2.MARE == 1], "bo")


plt.plot(ndata2.mdiam, ndata2.dD, "ro")
plt.plot(ndata2.mdiam, ndata2.m_mcw, "ro")
plt.plot(ndata2.m_mcw,ndata2.mdiam, "ro")
plt.plot(ndata2.m_ucw,ndata2.mdiam, "bo")
plt.plot(ndata2.mdiam,ndata2.mhr, "bo")

plt.plot(ndata2.m_mcw, ndata2.dD, "ro")

plt.plot(ndata2.mdiam, np.abs(ndata2.mdepth), "ro")
diamsynth = np.linspace(10,10000,1000)
depthsynth = diamsynth*0.2
depthsynth15 = diamsynth*0.15

depthsynth2 = diamsynth*0.1
depthsynth3 = diamsynth*0.05
depthsynth4 = diamsynth*0.01


plt.plot(diamsynth, depthsynth, "b")
plt.plot(diamsynth, depthsynth15, "k")
plt.plot(diamsynth, depthsynth2, "g")
plt.plot(diamsynth, depthsynth3, "k")
plt.plot(diamsynth, depthsynth4, "c")
plt.xlim((0,1000))
plt.ylim((0,200))
