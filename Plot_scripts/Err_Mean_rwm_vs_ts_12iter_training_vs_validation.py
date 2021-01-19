#!/usr/bin/env python2i
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 10:50:28 2019

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
    
rwm_mean_training = [1.8507754,  1.3737906,  1.03473818, 1.01393374, 1.02244668, 1.03256836,  1.03255869, 1.05723132, 1.06402569, 1.06403518, 1.06402687]

dT_mean_training = [-3.12,       -2.24,       -1.6,        -0.8,        -0.8,        -0.8,  -0.8,        -0.64,       -0.56,       -0.56,       -0.56 ]

rwm_mean_val =  [2.15324507, 1.57872213, 1.29707394, 1.0917369,  1.10137437, 1.10910221,  1.10909995, 1.08382754, 1.07858917, 1.07858923, 1.07858901]

dT_mean_val = [-3.98693878, -2.54460641, -1.76046647, -0.94880466, -0.98612245, -1.00618076,  -1.00618076, -0.9296793,  -0.93574344, -0.93574344, -0.93574344]
#
iN=10
xint=range(0,iN+4)
#xint=range(10,iN+9)

rwm_mean_training2 = [1.07108062, 1.08839586, 1.11471405, 1.14673575, 1.14697491, 1.14801913]
##
dT_mean_training2 = [-1.37106719, -1.37391304, -1.2113834,  -1.12853755, -1.13106719, -1.12853755]
#
#plt.figure(figsize=(10,10))
# create figure and axis objects with subplots()
fig,ax = plt.subplots(figsize=(5,5))
# make a plot
#ax.plot(xint, Err[0:len(xint),0], color="red", marker="o")
#ax.scatter(0, 2.4623305, color="red", marker="*")
#ax.plot(xint, rwm_mean_training[0:len(xint)], color="red", marker="o", label='Mean of $RWM$: training set')
#ax.plot(xint, rwm_mean_val[0:len(xint)],'--', color="red", marker="x", label='Mean of $RWM$: validation set')
yy=np.array(np.linspace(0, 1))*2
ax.plot(9*np.ones(yy.shape),yy,c='k',linestyle='dashed')
ax.plot(xint[0:10], rwm_mean_training[0:10], color="red", marker="o", label='$\mu_{RWM}$: training set')
ax.plot(xint[0:10], rwm_mean_val[0:10],'--', color="red", marker="x", label='$\mu_{RWM}$: validation set')
#ax.plot(xint[0:10], rwm_mean_training[0:10], color="red", marker="o")
#ax.plot(xint[0:10], rwm_mean_val[0:10],'--', color="red", marker="x")

ax.set_xlabel("Iteration, $\mathbf{m}$",fontsize=14)
# set y-axis label
ax.set_ylabel("Mean of relative waveform misfit, $\mu_{RWM}$",color="red",fontsize=14)
#plt.ylim([2.0, 2.5])
ax.legend(loc='upper right',bbox_to_anchor=(0.95, 0.95))
# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
#ax2.plot(xint, Err[0:len(xint),1],color="blue",marker="o")
#ax2.scatter(0, 396, color="blue", marker="*")
ax2.plot(xint[0:10], np.abs(dT_mean_training[0:10]),color="blue",marker="o",label='$\mu_{\Delta T}$: training set')
ax2.plot(xint[0:10], np.abs(dT_mean_val[0:10]),'--',color="blue",marker="x",label='$\mu_{\Delta T}$: validation set')
#ax2.plot(xint[0:10], np.abs(dT_mean_training[0:10]),color="blue",marker="o")
#ax2.plot(xint[0:10], np.abs(dT_mean_val[0:10]),'--',color="blue",marker="x")

ax2.plot(xint[9:-1], np.abs(rwm_mean_training2[0:4]),color="red",marker="o")
ax2.plot(xint[9:-1], np.abs(dT_mean_training2[0:4]),color="blue",marker="o")
         
ax2.set_ylabel("Mean of travel-time shift, $\mu_{\Delta T}$(s)",color="blue",fontsize=14)
ax2.legend(loc='upper right',bbox_to_anchor=(0.95, 0.8))
plt.xticks(xint[0:16])
plt.xlim([-0.2, 12.2])
ax.set_ylim([0, 2.2])
ax2.set_ylim([0, 5])
#plt.xticks(xint[0:10])
#plt.xticks(xint)
plt.show()
# save the plot as a file
fig.savefig('Misfit_vs_window.png',
            format='png',
            dpi=300,
            bbox_inches='tight')