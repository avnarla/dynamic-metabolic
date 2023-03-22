# Figures for plotting coculture simulation

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pickle
import itertools

################################################

### Data ####
t = [-10,0, 0.27, 1.30, 2.67, 4.28, 6.30, 7.55, 8.85, 11.20, 24]

glcnac = [0, 0, 4.970, 4.951, 4.881, 4.374, 1.742, 0, 0, 0, 0]
ac = [0, 0, 0.028, 0.042, 0.142, 0.625, 3.558, 5.287, 3.547, 0, 0]
a_cfu = [1.3E+09, 1.3E+09, 3.8E+07, 3.7E+07, 3.6E+07, 7.6E+07, 3.2E+08, 9.9E+08, 1.4E+09, 1.4E+09, 1.4E+09]
b_cfu = [4.1E+08, 4.1E+08, 8.4E+06, 8.0E+06, 9.9E+06, 2.0E+07, 5.5E+07, 1.2E+08, 1.7E+08, 4E8, 4E8]
pH = [7.81, 7.81, 7.80, 7.71, 7.55, 7.46, 7.53, 7.74, 7.78]


##### Set params of figure #####
time_ticks = np.array([0,6,12,18,24])
time_ticklabels = ['0','6','12','18','24']

conc_ticks = np.array([0,1,2,3,4,5])
conc_ticklabels = conc_ticks.astype(str)

# CFU properties
cfu_lim = [np.log(10**6),np.log(2*10**9)]   # OD_600 limit
cfu_ticks_vals = np.array([10**6,10**7,10**8,10**9])
cfu_ticks = np.log(cfu_ticks_vals)
cfu_ticklabels = [r'$10^6$',r'$10^7$',r'$10^8$',r'$10^9$']

pH_lim = [6.8,8.6]
pH_ticks = np.array([7,7.5,8,8.5])
pH_ticklabels = pH_ticks.astype(str)

### PLOT ###
numrows = 2
numcols = 1
fig,axs = plt.subplots(numrows,numcols,figsize=(3,4),dpi=400)
fig.subplots_adjust(bottom=0.15,top=0.95,right=0.85,left=0.25)

# plot top plot: 1A01 and 3B05 CFUs
clix = 0
# axs[clix].title.set_text('Stable cycle')
axs[clix].set_ylabel(r'CFU/mL')
# axs[clix].set_xlabel('Time (hr)')
axs[clix].plot(t[0:2],np.log(a_cfu[0:2]),'o-',color='magenta')
axs[clix].plot(t[0:2],np.log(b_cfu[0:2]),'o--',color='magenta',markerfacecolor='w', markeredgewidth=1, markeredgecolor='magenta')
axs[clix].plot(t[2:],np.log(a_cfu[2:]),'o-',color='magenta',label='1A01')
axs[clix].plot(t[2:],np.log(b_cfu[2:]),'o--',color='magenta',markerfacecolor='w', markeredgewidth=1, markeredgecolor='magenta',label='3B05')
axs[clix].set_xlim(-2,25)
axs[clix].set_ylim(cfu_lim)
axs[clix].set_xticks(time_ticks)
axs[clix].set_xticklabels([])
axs[clix].set_yticks(cfu_ticks)
axs[clix].set_yticklabels(cfu_ticklabels)
axs[clix].legend(loc='lower right',fontsize=7)

clix = 1
axs[clix].set_ylabel(r'Concentration (mM)')
axs[clix].set_xlabel('Time (hr)')
axs[clix].plot(t[0:2],glcnac[0:2],'b^-',markerfacecolor='w', markeredgewidth=1, markeredgecolor='b')
axs[clix].plot(t[0:2],ac[0:2],'rs-',markerfacecolor='w', markeredgewidth=1, markeredgecolor='r')
axs[clix].plot(t[2:],glcnac[2:],'b^-',markerfacecolor='w', markeredgewidth=1, markeredgecolor='b',label='GlcNAc')
axs[clix].plot(t[2:],ac[2:],'rs-',markerfacecolor='w', markeredgewidth=1, markeredgecolor='r',label='Acetate')
axs[clix].set_xlim(-2,25)
axs[clix].set_ylim(-0.5,5.8)
axs[clix].set_yticks(conc_ticks)
axs[clix].set_xticks(time_ticks)
axs[clix].set_yticklabels(conc_ticklabels)
axs[clix].set_xticklabels(time_ticklabels)
axs[clix].legend(loc='upper right',fontsize=7)

ax1 = axs[clix].twinx()
ax1.plot(t[2:],pH,'o-',color='xkcd:orange')
ax1.set_ylim(pH_lim)
ax1.set_yticks(pH_ticks)
ax1.set_yticklabels(pH_ticklabels)
ax1.set_ylabel(r'pH')
ax1.spines["right"].set_color('xkcd:orange')
ax1.yaxis.label.set_color('xkcd:orange')
ax1.tick_params(axis="y",colors='xkcd:orange')

plt.subplots_adjust(wspace=0, hspace=0.1)

fig
fig.savefig('figS7_hepes.eps',format='eps')



