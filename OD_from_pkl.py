# Obtain OD from pickle file

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

######## Extract Data #################

with open('ecpp.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_ecpp = all_days
all_times_ecpp = t

with open('cfpf.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_cfpf = all_days
all_times_cfpf = t

######## Plot #################
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



