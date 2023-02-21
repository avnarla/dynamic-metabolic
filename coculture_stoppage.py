# Showing stoppage of coculture
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
K_ac = 0.01                  # Km for 3B05 growth on acetate (mM)
# num_Km = len(K_ac_range)

X = lambda_b0/(lambda_a0*Y_b*E_a)	# constant containing growth rates and yields

# r0 = 1/X							# B0/A0
c_crash_B = 3.2                   # [Ac] at which growth stops, B [mM]
c_crash_A = 3.5                   # [Ac] at which growth stops, A [mM]

# csv_fn = 'coculture_fullsim_varyKmAc.csv'    # filename for csv output file

# varied parameter; params to calculate as a function of varied parameter #
# rhoB0_range = np.arange(0.001,0.1,0.002)                   # initial total OD: A0+B0
# od_crash = np.zeros_like(K_ac_range)   # array containing stopping ODs for different r0
# t_crash = np.zeros_like(K_ac_range)    # array containing stopping times for different r0 
# A_crash = np.zeros_like(K_ac_range)    # array containing stopping OD for 1A01
# B_crash = np.zeros_like(K_ac_range)    # array containing stopping OD for 3B05
# nag_crash = np.zeros_like(K_ac_range)  # array containing stopping [GlcNAc]
# num_rho0 = len(rhoB0_range)               # number of r0's tried

# for 1 parameter (r0) value #
tfin = 25      # The length of the simulation
dt = 0.01       # dt of time
t = np.arange(0,tfin,dt)    # Time (h)

A = np.zeros((len(t)))   # 1A01 growth dynamics, each column is the dynamics for 1 r0 value 
B = np.zeros_like(A)            # 3B05 growth dynamics, each column is the dynamics for 1 r0 value
C = np.zeros_like(A)            # [Acetate] dynamics, each column is the dynamics for 1 r0 value
G = np.zeros_like(A)            # [GlcNAc] dynamics
################### RUN SIMULATION #######################

A[0] = 0.01
B[0] = 0.01
G[0] = 5

for tix in range(len(t)-1):        # run ODEs for a single coculture with some set of params

    if C[tix] < c_crash_B:       # set growth rates of A and B depending on [Ac] (pH)
        lambda_b = lambda_b0
    else:
        lambda_b = 0

    if C[tix] < c_crash_A:
        lambda_a = lambda_a0
    else:
        lambda_a = 0

    A[tix+1] = A[tix] + lambda_a*A[tix]*dt
    B[tix+1] = B[tix] + lambda_b*C[tix]/(C[tix]+K_ac)*B[tix]*dt
    C[tix+1] = C[tix] + (lambda_a*E_a*A[tix] - lambda_b/Y_b*C[tix]/(C[tix]+K_ac)*B[tix])*dt
    if C[tix+1] < 0: # The [Ac] cannot be < 0 (which could happen if some numerical error)
        C[tix+1] = 0
    G[tix+1] = G[tix] - lambda_a/Y_a*A[tix]*dt
    if G[tix+1] < 0:  # The [GlcNAc] cannot be negative
        G[tix+1] = 0

sdx = np.argmin(abs(C[:] - c_crash_A))    # get loop at which [Ac] got to c_crash (stopping index)
A_crash = A[sdx]
B_crash = B[sdx]
od_crash = A[sdx] + B[sdx]
nag_crash = 5-((A_crash-0.01)*(1/Y_a))
t_crash = sdx*dt

############# PLOT ############
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

# ##### Figure  #####
fig,axs = plt.subplots(2,1,figsize=cm2inch(6,9),dpi=400)
fig.subplots_adjust(bottom=0.18,top=0.85,right=0.9,left=0.28)

time_lim = [0,10]
time_ticks = np.array([0,2,4,6,8,10])
time_ticklabels = ['0','2','4','6','8','10']

conc_ticks = np.array([0,1,2,3,4,5])
conc_ticklabels = ['0','1','2','3','4','5']

# OD properties
od_lim = [np.log(0.005),np.log(1)]   # OD_600 limit
od_ticks_vals = np.array([0.01,0.025,0.05,0.1,0.2,0.4,0.8])
od_ticks = np.log(od_ticks_vals)
od_ticklabels = od_ticks_vals.astype(str)

# top plot: A and B growth
clix=1
axs[clix].plot(t,np.log(A),'k-',label='1A01')
axs[clix].plot(t,np.log(B),'k--',label='3B05')
axs[clix].set_xlim(time_lim)
axs[clix].set_xticks(time_ticks)
#axs[clix].set_xticklabels([])
axs[clix].set_ylim(od_lim)
axs[clix].set_yticks(od_ticks)
axs[clix].set_yticklabels(od_ticklabels)
axs[clix].set_ylabel(r'$\mathdefault{OD_{600}}$',fontname='Arial')
axs[1].set_xlabel('Time (hr)')
L1 = axs[0].legend(loc='upper left')

for tick in axs[0].get_xticklabels():
    tick.set_fontname("Arial")
for tick in axs[0].get_yticklabels():
    tick.set_fontname("Arial")

plt.setp(L1.texts, family='Arial')

# concentrations of GlcNAc/Ac

# axs[1].plot(tar,gar_d[:,cycle],'b-',label='GlcNAc')
clix=0
axs[clix].plot(t,G,'b-',label='GlcNAc')
axs[clix].plot(t,C,'r-',label='Acetate')
axs[clix].set_xlim(time_lim)
axs[clix].set_ylim(-0.3,5.1)
axs[clix].set_yticks(conc_ticks)
axs[clix].set_xticks([])#time_ticks)
axs[clix].set_yticklabels(conc_ticklabels)
axs[clix].set_ylabel('Concentration (mM)',fontname='Arial')
L2 = axs[1].legend(loc='lower right')
axs[1].set_xlabel('Time (h)')

for tick in axs[1].get_xticklabels():
    tick.set_fontname("Arial")
for tick in axs[1].get_yticklabels():
    tick.set_fontname("Arial")

plt.setp(L2.texts, family='Arial')


plt.subplots_adjust(wspace=0, hspace=0.1)


fig
fig.savefig('stoppage.pdf',format='pdf')

