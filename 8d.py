# Plot OD on first cycle, for both cocultures
# Date: 3/10/22

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

import math
import re
import os
import glob as glob
import itertools
import pickle

####### KEY VARIABLES ######
#### Key variables
# 1. lambdas (ls in function calls): list of wavelengths at which data was collected
# 2. wells (ws): list of all wells on a 96 well plate (e.g. 'A1','A2'...'H12')
# 3. t (t_p): time points for measurement
# 4. all_data (ad): 96 wells x number of wavelengths x number of time points array containing all data from
#                   experiment. The wells list is indexed according to the 1st dimension of this array.

####### BASIC FUNCTIONS #########

def get_odt(well,wavelength,sub,ls,ws,ad):
	'''Get the OD(t) at wavelength for well (str)'''
	turb = np.where(ls == wavelength)	# get idx for wavelength of interest
	turb = turb[0][0] 
	idx = ws.index(well)				# find idx corresponding to well
	od = ad[idx][turb][:]			# get data at selected wavelength for selected well
	if sub == 1: 
		od = od - od[0]					# subtract off first timepoint if you say so
	return od

def cm2inch(*tupl):
	'''Get the OD(t) at wavelength for well (str)'''
	inch = 2.54
	if isinstance(tupl[0], tuple):
		return tuple(i/inch for i in tupl[0])
	else:
		return tuple(i/inch for i in tupl)


######## Extract Data #################

with open('ecpp.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_ecpp = all_days
all_times_ecpp = t

with open('cfpf.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_cfpf = all_days
all_times_cfpf = t

######## Plot Properties #################

# OD properties (y axis)
od_lim = [np.log(0.01),np.log(1.28)]   # OD_600 limit
od_ticks_vals = np.array([0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28])
od_ticks = np.log(od_ticks_vals)
od_ticklabels = ['0.01','0.02','0.04','0.08','0.16','0.32','0.64','1.28']


# t axis properties
time_lim = [0,24]
time_ticks = np.array([0,6,12,18,24])


SMALL_SIZE = 7
MEDIUM_SIZE = 8
BIGGER_SIZE = 9

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

######## Plot #################
numrows = 1
numcols = 1
fig,axs = plt.subplots(numrows,numcols,figsize=cm2inch(6,4.5),dpi=400)
fig.subplots_adjust(bottom=0.34,top=0.85,right=0.71,left=0.25)

well_co = 'A1'
well_mono = 'D1'
# ad_ecpp = all_days_ecpp[0]
# t_ecpp = all_times_ecpp[0]
# odt_mono = get_odt(well_mono,420,1,lambdas,wells,ad_ecpp)
# odt_co = get_odt(well_co,420,1,lambdas,wells,ad_ecpp)
ad_cfpf = all_days_cfpf[0]
t_cfpf = all_times_cfpf[0]
odt_mono = get_odt(well_mono,420,1,lambdas,wells,ad_cfpf)
odt_co = get_odt(well_co,420,1,lambdas,wells,ad_cfpf)

# axs.set_title('Ec+Pp coculture',fontname='Arial')
axs.set_ylabel(r'OD$_{600}$',fontname='Arial')
axs.set_xlabel('Time in 1st cycle (h)',fontname='Arial')

# Select which coculture to plot"
# Ec/pp
# axs.plot(t_ecpp,np.log(odt_mono+0.02),'--',color='k',alpha=0.6,label='Ec')
# axs.plot(t_ecpp,np.log(odt_co+0.02),'-',color='k',label='Ec+Pp')

# Cf/Pf
axs.plot(t_cfpf,np.log(odt_mono+0.02),'--',color='k',alpha=0.6,label='Cf')
axs.plot(t_cfpf,np.log(odt_co+0.02),'-',color='k',label='Cf+Pf')

axs.set_ylim(od_lim)
axs.set_yticks(od_ticks)
axs.set_yticklabels(od_ticklabels)
axs.set_xlim(time_lim)
axs.set_xticks(time_ticks)

for tick in axs.get_xticklabels():
    tick.set_fontname("Arial")
for tick in axs.get_yticklabels():
    tick.set_fontname("Arial")        

L = axs.legend(loc='lower right')
plt.setp(L.texts, family='Arial')

fig

fig.savefig('CfPf_MonoCo_OD.png',format='png')



