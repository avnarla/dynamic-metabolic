# Km dependence
# Author: KA
# Date: 1/23/19

import numpy as np
import matplotlib.pyplot as plt


##### Function definitions #####


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

#### INPUTS ####

# parameters for the two species #
Y_a = 1/6.1                     # yield of 1A01 on GlcNAc, [OD/mM]
E_a = 9.3                       # Ac excretion by 1A01, [mM/OD]
lambda_a0 = 0.7                # steady-state growth rate of 1A01 

Y_b = 1/32                    # Ac uptake by 3B05, [OD/mM]
lambda_b0 = 0.35                # steady-state growth rate of 3B05, [hr^-1]
K_ac_range = np.arange(0.001,0.1,0.005)                    # Km for 3B05 growth on acetate (mM)
num_Km = len(K_ac_range)

X = lambda_b0/(lambda_a0*Y_b*E_a)	# constant containing growth rates and yields

# r0 = 1/X							# B0/A0

c_crash = 3.2                   # [Ac] at which growth stops, [mM]

# csv_fn = 'coculture_fullsim_varyKmAc.csv'    # filename for csv output file

# varied parameter; params to calculate as a function of varied parameter #
# rhoB0_range = np.arange(0.001,0.1,0.002)                   # initial total OD: A0+B0
od_crash = np.zeros_like(K_ac_range)   # array containing stopping ODs for different r0
t_crash = np.zeros_like(K_ac_range)    # array containing stopping times for different r0 
A_crash = np.zeros_like(K_ac_range)    # array containing stopping OD for 1A01
B_crash = np.zeros_like(K_ac_range)    # array containing stopping OD for 3B05
nag_crash = np.zeros_like(K_ac_range)  # array containing stopping [GlcNAc]
# num_rho0 = len(rhoB0_range)               # number of r0's tried

# for 1 parameter (r0) value #
tfin = 25      # The length of the simulation
dt = 0.01       # dt of time
t = np.arange(0,tfin,dt)    # Time (h)

A = np.zeros((len(t),num_Km))   # 1A01 growth dynamics, each column is the dynamics for 1 r0 value 
B = np.zeros_like(A)            # 3B05 growth dynamics, each column is the dynamics for 1 r0 value
C = np.zeros_like(A)            # [Acetate] dynamics, each column is the dynamics for 1 r0 value

################### RUN SIMULATION #######################

for vix in range(num_Km):    # run ODEs for a range of some 'v'aried paramater
    # rho_0 = rho0_range[vix]       # set all parameters for loop iteration, starting with initial ratio

    # A[0,vix] = rho_0/(1+r0)					# initial density of A
    # B[0,vix] = rho_0 - A[0,vix] 					# initial density of B

    # B[0,vix] = rhoB0_range[vix]    # B0
    # A[0,vix] = B[0,vix]*X  # A0 for a particular r0

    K_ac = K_ac_range[vix]
    A[0,vix] = 0.01
    B[0,vix] = 0.01

    for tix in range(len(t)-1):        # run ODEs for a single coculture with some set of params

        if C[tix,vix] < c_crash:       # set growth rates of A and B depending on [Ac] (pH)
            lambda_a = lambda_a0
            lambda_b = lambda_b0
        else:
            lambda_a = 0
            lambda_b = 0

        A[tix+1,vix] = A[tix,vix] + lambda_a*A[tix,vix]*dt
        B[tix+1,vix] = B[tix,vix] + lambda_b*C[tix,vix]/(C[tix,vix]+K_ac)*B[tix,vix]*dt
        C[tix+1,vix] = C[tix,vix] + (lambda_a*E_a*A[tix,vix] - lambda_b/Y_b*C[tix,vix]/(C[tix,vix]+K_ac)*B[tix,vix])*dt
        if C[tix+1,vix] < 0: # The [Ac] cannot be < 0 (which could happen if some numerical error)
            C[tix+1,vix] = 0

    sdx = np.argmin(abs(C[:,vix] - c_crash))    # get loop at which [Ac] got to c_crash (stopping index)
    A_crash[vix] = A[sdx,vix]
    B_crash[vix] = B[sdx,vix]
    od_crash[vix] = A[sdx,vix] + B[sdx,vix]
    nag_crash[vix] = 5-((A_crash[vix]-0.01)*(1/Y_a))
    t_crash[vix] = sdx*dt

############# PLOT ############
od_tot = A + B

### Figure of OD, growth of A and B ###
# Time limits and properties
time_lim = [0,15]     
time_ticks = np.array([0,2,4,6,8,10,12,14])
time_ticklabels = time_ticks.astype(str)

# OD properties
od_lim = [np.log(0.01),np.log(1.5)]   # OD_600 limit
od_ticks_vals = np.array([0.05,0.1,0.2,0.4,0.6,1.0,1.5])
od_ticks = np.log(od_ticks_vals)
od_ticklabels = od_ticks_vals.astype(str)

SMALL_SIZE = 7
MEDIUM_SIZE = 8
BIGGER_SIZE = 8

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# xlabel = 'Time (hr)' # set x label 

ms = 4 # markersize


## Plot dependence of crashing concentraitons on 3B05 Km of acetate ##
# Km properties
km_lim = [np.log(0.0008),np.log(0.12)]   # OD_600 limit
km_ticks_vals = np.array([0.001,0.01,0.1])
km_ticks = np.log(km_ticks_vals)
km_ticklabels = [r'$\mathdefault{10^{-3}}$',r'$\mathdefault{10^{-2}}$',r'$\mathdefault{10^{-1}}$']



fig,axs = plt.subplots(1,1,figsize=cm2inch(6,4.5),dpi=400)
fig.subplots_adjust(bottom=0.25,top=0.9,right=0.8,left=0.2)

panel=0
p1, = axs.plot(np.log(K_ac_range),t_crash,color='k')
axs.set_ylabel(r'$t_{\times}$ (hr)',fontname='Arial')
axs.set_xlabel(r'$K_{B,E}}$ (mM)',fontname='Arial')
axs.set_xlim(km_lim)
axs.set_xticks(km_ticks)
axs.set_xticklabels([])
axs.yaxis.label.set_color(p1.get_color())
axs.tick_params(axis="y",colors=p1.get_color())

for tick in axs.get_xticklabels():
    tick.set_fontname("Arial")
for tick in axs.get_yticklabels():
    tick.set_fontname("Arial")


ax1 = axs.twinx()
ax1.plot(np.log(K_ac_range),nag_crash,color='b')
ax1.set_xlim(km_lim)
ax1.set_xticks(km_ticks)
ax1.set_xticklabels(km_ticklabels)
ax1.set_ylabel(r'$G(t_{\times})$ (mM)',fontname='Arial')
# ax1.legend(loc='center left')
ax1.spines["left"].set_color(p1.get_color())
ax1.spines["right"].set_color('b')
ax1.yaxis.label.set_color('b')
ax1.tick_params(axis="y",colors='b')

for tick in ax1.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax1.get_yticklabels():
    tick.set_fontname("Arial")

fig
fig.savefig('Km_dep.eps',format='eps')

