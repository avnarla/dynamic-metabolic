# Obtain pH from pickle file

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
	od = all_data[idx][turb][:]			# get data at selected wavelength for selected well
	if sub == 1: 
		od = od - od[0]					# subtract off first timepoint if you say so
	return od

def get_pH(well_dye,well_nodye,ad,ws,t_p):
	'''Get the pH(t) for well (str)'''
	pH = np.zeros((len(t_p),1)) # initiate pH vector
	for i in range(0,len(t_p)):

		idx = ws.index(well_dye) # find idx corresponding to well
		od_wl_dye = ad[idx][:,i] # get spectrum for time point i

		idx2 = ws.index(well_nodye) # find idx corresponding to well
		od_wl_nodye = ad[idx2][:,i] # get spectrum for time point i

		od_wl = od_wl_dye - od_wl_nodye # subtract off cells

		pH[i] = od_wl[19]-od_wl[25] # subtract background absorption at 650 nm from peak at 590 nm OR

	pH_val = 8.75*pH**0.158 # convert absorption difference to pH using standard curve

	return pH_val

######## Extract Data #################
day = 3 # day 3 is stable cycle

with open('cfpf.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_cfpf = all_days
all_times_cfpf = t

t_cp = np.array([0, 6.08, 8.1, 9.57, 11.1, 24])
gn_cp = np.array([10, 8.209,5.758,4.568,0,0])
ac_cp = np.array([0, 0.929, 1.289, 1.250, 1.524, 0])
pyr_cp = np.array([0, 0.137,0.949, 1.104, 0.570, 0])

with open('ecpp.pkl','rb') as f:
    t,wells,all_days,lambdas = pickle.load(f)
all_days_ecpp = all_days
all_times_ecpp = t

######## Plot #################
fig,axs = plt.subplots(1,1,figsize=(6,4.5),dpi=400)
fig.subplots_adjust(bottom=0.2,top=0.85,right=0.8,left=0.2)


###### Plot metabolite dynamics of Cf + Pf, D3 (comment out if want to plot Ec/Pp #######
axs.set_ylabel('Concentration (mM)',fontname='Arial')
axs.set_xlabel('Time in Cycle (h)')
# # axs[1,0].plot(t_cp,gn_cp,'b^-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='b',label='GlcNAc')
axs.plot(t_cp,ac_cp,'rs-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='r',label='Acetate')
axs.plot(t_cp,pyr_cp,'^-',color='xkcd:purple',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:purple',label='Pyruvate')


# ###### Plot metabolite dynamics of Ec + Pp, D3 (comment out if want to plot Cf/Pf) #######
# axs.set_ylabel('Concentration (mM)',fontname='Arial')
# axs.set_xlabel('')
# axs[1,1].plot(t_ep,g_ep,'b^-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='b',label='GlcNAc')
# axs.plot(t_ep,ac_ep,'rs-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='r',label='Acetate')
# axs.plot(t_ep,form_ep+etoh_ep+succ_ep,'s-',color='xkcd:pink',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:pink',label='Formate+Ethanol+Succinate')
# axs.plot(t_ep,pyr_ep+lac_ep,'^-',color='xkcd:purple',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:purple',label='Pyruvate+Lactate')
