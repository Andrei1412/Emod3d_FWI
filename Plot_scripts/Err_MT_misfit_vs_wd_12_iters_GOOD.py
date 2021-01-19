#!/usr/bin/env python2i
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 10:50:28 2019

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
    
##iN=16
iN=10
Err=np.array([[8.10509518,304],[6.24488777,295],[4.48133162,292],[4.71763736,293],[4.08310807,295],[3.8712934,286],[3.99048825,286],[4.02149162,289],[4.00253364,287],[3.90183374,288],[3.27640018, 276],[3.20957784, 270],[3.21473381, 267],[3.18271712, 271],[3.214283, 268],[3.2315786, 272]])
#xint=range(0,iN-1)
Err_0 = Err[0,0]
#iN=5
Err2=np.array([[2.60653187,553],[2.55095285,548],[2.49932458,541],[2.47281439,533],[2.48493813, 532]])
Err2_0 = Err2[0,0];Err2_1 = Err2_0;

#Err=np.array([843.13996323,624.72674779,524.31773419,455.72303331,459.29979017,462.16926876,462.16840173,452.55299466,450.55013137,450.55020534,450.55011971])
#Err2=np.array([483.54925189,491.95850618,506.33704642,520.97264321])
#Err2=np.array([450.55011971,436.28985938,444.77741887,484.27325394,485.714463,486.00849363])
#Err2_0 = Err2[0];Err2_1 = Err2_0;
#xint=range(0,iN+3)
xint=range(0,iN+3)
# create figure and axis objects with subplots()
#fig,ax = plt.subplots(figsize=(7.5,5))
fig,ax = plt.subplots(figsize=(5,5))

#Err2_new = Err2[0:6,0]/(Err2_1)*(Err2_1/Err2_0)*(Err[len(xint)-7,0]/Err_0)
Err2_new = Err2[0:6,0]/Err_0*(Err[9,0]/Err2_0)
#Err2_new[0] =  (Err[len(xint)-7,0]/Err_0)

#line2 = ax.plot(xint[len(xint)-4:len(xint)], Err2_new[0:4], color="k", marker="o")
#line1 = ax.plot(xint[0:len(xint)-3], Err[0:len(xint)-3]/Err_0, color="k", marker="o")

line2 = ax.plot(xint[len(xint)-4:len(xint)], Err2_new[0:4], color="r", marker="o")
line1 = ax.plot(xint[0:len(xint)-3], Err[0:len(xint)-3,0]/Err_0, color="r", marker="o")

#Err10 = np.array([])
#Err10_new = Err10/Err2_0*(Err[len(xint)-5,0]/Err_0)
#line10 = ax.plot(xint[len(xint)-5:len(xint)-3], Err10_new[0:2], color="r", marker="o",linestyle='dashed')

yy=np.array(np.linspace(0, 1))
#ax.plot(9*np.ones(yy.shape),yy,c='k',linestyle='dashed')
ax.plot(9*np.ones(yy.shape),yy,c='k',linestyle='dashed')
#line2 = ax.plot(xint[len(xint)-5:len(xint)],[Err[len(xint)-5,0]/Err_0,Err2[0:4,0]/(Err2_1)*(Err2_1/Err2_0)*(Err[len(xint)-5,0]/Err_0)], color="k", marker="o", label='misfit with revised CMT')
#ax.legend(loc='upper right',bbox_to_anchor=(0.95, 0.95))
# set x-axis label
ax.set_xlabel("Iteration, $\mathbf{m}$",fontsize=14)
# set y-axis label
#ax.set_ylabel("Misfit",color="red",fontsize=14)
ax.set_ylabel("Normalized misfit, $\chi_{p}(\mathbf{m})/\chi_{p}(\mathbf{m_{00}})$",color="red",fontsize=14)
#ax.title("Normalized RWM according to 14 validation events")
#plt.ylim([0.4, 1.01])
# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
## make a plot with different y-axis using second axis object
line3 = ax2.plot(xint[0:len(xint)-3], Err[0:len(xint)-3,1],color="blue",marker="o")
line4 = ax2.plot(xint[len(xint)-4:len(xint)], Err2[0:4,1],color="blue",marker="o")
##ax2.scatter(0, 396, color="blue", marker="*")
#line3 = ax2.plot(xint, Err[0:len(xint),1],color="blue",marker="o",label='windows with revised CMT')
#line4 = ax2.scatter(0, 396, color="blue", marker="*",label='windows without revised CMT')
#ax2.legend(loc='upper right',bbox_to_anchor=(0.95, 0.75))
ax2.set_ylabel("Windows",color="blue",fontsize=14)
##plt.ylim([270, 305])
#plt.ylim([370, 400])
plt.xticks(xint[0:16])
#ax.text(5,0.8, '13 events', fontsize=14)
#ax.text(5,0.75, 'inversion', fontsize=14)
##ax.text(5,0.7, 'cc>0.7', fontsize=14)
#ax.text(9.2,0.8, '27 events ', fontsize=14)
#ax.text(9.2,0.75, 'inversion', fontsize=14)
#ax.text(9.4,0.7, 'cc>0.8', fontsize=14)
plt.xlim([-0.2, 12.2])
ax.set_ylim([0.4, 1.01])
#ax2.xlim([-0.1, 12.1])
ax2.set_ylim([100, 700])
plt.show()
# save the plot as a file
fig.savefig('Misfit_vs_window.jpg',
            format='jpeg',
            dpi=300,
            bbox_inches='tight')