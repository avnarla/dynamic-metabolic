# Plot pH and metabolites for CfPf and EcPp cocultures.
# Date: 2/11/22 -> plot metabolites for cocultures separately

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

def cm2inch(*tupl):
	'''Get the OD(t) at wavelength for well (str)'''
	inch = 2.54
	if isinstance(tupl[0], tuple):
		return tuple(i/inch for i in tupl[0])
	else:
		return tuple(i/inch for i in tupl)


######## Extract Data #################
day = 3 # third day of gr/dil expts from 1/15/22 and 1/26/22

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

t_ep = np.array([0,5.85,8.92,11.27,20.98,24])
g_ep = np.array([10,7.234,3.801,1.819,0,0])
ac_ep = np.array([0,1.559,2.815,2.741,0,0])
succ_ep = np.array([0,0.020,0.087,0.066,0,0])
form_ep = np.array([0,1.946,4.433,3.815,0.253,0.089])
etoh_ep = np.array([0,0.410,1.038,1.300,0.05,0])
lac_ep = np.array([0,0.005,0.518,1.307,0.127,0])
pyr_ep = np.array([0,0.002,0.330,0.640,0,0])


######## Plot Properties #################

# pH properties (y axis)
pH_lim = [4.9,7.1]   # OD_600 limit
pH_ticks = np.array([5,5.5,6,6.5,7])
pH_ticklabels = pH_ticks.astype(str)

# conc prop (CfPf)
conc_lim = [-0.2,2.5]
conc_ticks = np.array([0,1,2])
conc_ticklabels = conc_ticks.astype(str)

# conc prop (EcPp)
conc_lim_2 = [-0.2,6.2]
conc_ticks_2 = np.array([0,2,4,6])
conc_ticklabels_2 = conc_ticks_2.astype(str)

# t axis properties
time_lim = [-0.5,24.5]
time_ticks = np.array([0,6,12,18,24])
time_ticklabels = time_ticks.astype(str)

SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 9

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=7)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

ms = 4 # marker size
######## Plot #################
fig,axs = plt.subplots(1,1,figsize=cm2inch(6,4.5),dpi=400)
fig.subplots_adjust(bottom=0.2,top=0.85,right=0.8,left=0.2)


###### Plot metabolite dynamics of Cf + Pf, D3 (comment out if want to plot Ec/Pp #######
axs.set_ylabel('Concentration (mM)',fontname='Arial')
axs.set_xlabel('Time in Cycle (h)')
# # axs[1,0].plot(t_cp,gn_cp,'b^-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='b',label='GlcNAc')
axs.plot(t_cp,ac_cp,'rs-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='r',label='Acetate')
axs.plot(t_cp,pyr_cp,'^-',color='xkcd:purple',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:purple',label='Pyruvate')
axs.set_xlim(time_lim)
axs.set_ylim(conc_lim)
axs.set_yticks(conc_ticks)
axs.set_xticks(time_ticks)
axs.set_yticklabels(conc_ticklabels)
axs.set_xticklabels(time_ticklabels)

# ###### Plot metabolite dynamics of Ec + Pp, D3 (comment out if want to plot Cf/Pf) #######
# axs.set_ylabel('Concentration (mM)',fontname='Arial')
# axs.set_xlabel('')
# axs[1,1].plot(t_ep,g_ep,'b^-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='b',label='GlcNAc')
# axs.plot(t_ep,ac_ep,'rs-',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='r',label='Acetate')
# axs.plot(t_ep,form_ep+etoh_ep+succ_ep,'s-',color='xkcd:pink',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:pink',label='Formate+Ethanol+Succinate')
# axs.plot(t_ep,pyr_ep+lac_ep,'^-',color='xkcd:purple',markersize=ms,markerfacecolor='w', markeredgewidth=1, markeredgecolor='xkcd:purple',label='Pyruvate+Lactate')
# axs.set_xlim(time_lim)
# axs.set_ylim(conc_lim_2)
# axs.set_yticks(conc_ticks_2)
# axs.set_xticks(time_ticks)
# axs.set_yticklabels(conc_ticklabels_2)
# axs.set_xticklabels([])

for t,a,f,p in zip(t_ep,ac_ep,form_ep+etoh_ep+succ_ep,pyr_ep+lac_ep):
        print(str(t)+","+str(a)+","+str(f)+","+str(p))


for tick in axs.get_xticklabels():
    tick.set_fontname("Arial")
for tick in axs.get_yticklabels():
    tick.set_fontname("Arial")

# L=axs.legend(loc='upper right',framealpha=1,fontsize=7) # for Cf/Pf

# L = axs.legend(bbox_to_anchor=(0.55, 0.5),framealpha=1,fontsize=7,title='Metabolites') # (first bbox number 0.64 for EcPp)
# L = axs.legend(loc='upper right',framealpha=1) # CfPf
L = axs.legend(loc='upper right',framealpha=1) # EcPp
plt.setp(L.texts, family='Arial')
# plt.setp(L.get_title(), fontname='Arial')


### Adjust spacing and save ###
# plt.subplots_adjust(wspace=0.3, hspace=0.12)    

fig

fig.savefig('metabolites_D3_EcPp.pdf',format='pdf')



