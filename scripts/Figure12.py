## Figure 11 should be run first!

import os
import flopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import Manager

modelname = 'sva_tm'
tmodelname = 'sva_tm_t'
dsi_dir = os.path.join('.','history-matching//DSI_cal1a//model//')
ies_dir = os.path.join('.','history-matching//cal1a_ies')
figure_dir = os.path.join('.','figures')
modelfiles_dir = os.path.join('.',modelname)

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

pw_4_arr = np.load(os.path.join(dsi_dir,'DSI_end_pw4.npy'))
pw_8_arr = np.load(os.path.join(dsi_dir,'DSI_end_pw8.npy'))
pw_9_arr = np.load(os.path.join(dsi_dir,'DSI_end_pw9.npy'))
pw_10_arr = np.load(os.path.join(dsi_dir,'DSI_end_pw10.npy'))

cm = 1/2.54
fig, ax = plt.subplots(2,2,figsize=(18*cm,18*cm),dpi=600)

ax[0][0].hist(np.array(pw_4_arr),bins=25,color='blue',density=True,alpha=0.4)
ax[0][0].tick_params(axis='both', which='major', labelsize=7)
ax[0][0].set_title('pw4',size=8)
ax[0][0].axvline(x=0.99,color='blue',linestyle='dashed',linewidth=0.72)
ax[0][0].text(0.94, ax[0][0].get_ylim()[1]*0.67, 'True value = 0.99', rotation=90,size=7)
ax[0][0].axvline(x=0.77,color='red',linestyle='dashed',linewidth=0.72)
ax[0][0].text(0.79, ax[0][0].get_ylim()[1]*0.67, 'RI-EN MAP = 0.77', rotation=90,size=7)
ax[0][0].axvline(x=0.75,color='black',linestyle='dashed',linewidth=0.72)
ax[0][0].text(0.71, ax[0][0].get_ylim()[1]*0.67, 'IES mean = 0.75', rotation=90,size=7)
ax[0][0].set_xlim([0,1])
ax[0][0].set_xlabel('depletion',fontsize=8)
ax[0][0].set_ylabel('density',fontsize=8)

ax[0][1].hist(np.array(pw_8_arr),bins=25,color='blue',density=True,alpha=0.4)
ax[0][1].tick_params(axis='both', which='major', labelsize=8)
ax[0][1].set_title('pw8',size=8)
ax[0][1].axvline(x=0.57,color='blue',linestyle='dashed',linewidth=0.72)
ax[0][1].text(0.59, ax[0][1].get_ylim()[1]*0.67, 'True value = 0.57', rotation=90,size=7)
ax[0][1].axvline(x=0.81,color='red',linestyle='dashed',linewidth=0.72)
ax[0][1].text(0.83, ax[0][1].get_ylim()[1]*0.67, 'RI-EN MAP = 0.81', rotation=90,size=7)
ax[0][1].axvline(x=0.27,color='black',linestyle='dashed',linewidth=0.72)
ax[0][1].text(0.29, ax[0][1].get_ylim()[1]*0.67, 'IES mean = 0.27', rotation=90,size=7)
ax[0][1].set_xlim([0,1])
ax[0][1].set_xlabel('depletion',fontsize=8)
ax[0][1].set_ylabel('density',fontsize=8)

ax[1][0].hist(np.array(pw_9_arr),bins=25,color='blue',density=True,alpha=0.4)
ax[1][0].tick_params(axis='both', which='major', labelsize=8)
ax[1][0].set_title('pw9',size=8)
ax[1][0].axvline(x=0.80,color='blue',linestyle='dashed',linewidth=0.72)
ax[1][0].text(0.82, ax[1][0].get_ylim()[1]*0.67, 'True value = 0.80', rotation=90,size=7)
ax[1][0].axvline(x=0.30,color='red',linestyle='dashed',linewidth=0.72)
ax[1][0].text(0.32, ax[1][0].get_ylim()[1]*0.67, 'RI-EN MAP = 0.30', rotation=90,size=7)
ax[1][0].axvline(x=0.10,color='black',linestyle='dashed',linewidth=0.72)
ax[1][0].text(0.12, ax[1][0].get_ylim()[1]*0.67, 'IES mean = 0.10', rotation=90,size=7)
ax[1][0].set_xlim([0,1])
ax[1][0].set_xlabel('depletion',fontsize=8)
ax[1][0].set_ylabel('density',fontsize=8)

ax[1][1].hist(np.array(pw_10_arr),bins=25,color='blue',density=True,alpha=0.4)
ax[1][1].tick_params(axis='both', which='major', labelsize=8)
ax[1][1].set_title('pw10',size=8)
ax[1][1].axvline(x=0.04,color='blue',linestyle='dashed',linewidth=0.72)
ax[1][1].text(0.06, ax[1][1].get_ylim()[1]*0.67, 'True value = 0.04', rotation=90,size=7)
ax[1][1].axvline(x=0.39,color='red',linestyle='dashed',linewidth=0.72)
ax[1][1].text(0.41, ax[1][1].get_ylim()[1]*0.67, 'RI-EN MAP = 0.39', rotation=90,size=7)
ax[1][1].axvline(x=0.55,color='black',linestyle='dashed',linewidth=0.72)
ax[1][1].text(0.57, ax[1][1].get_ylim()[1]*0.67, 'IES mean = 0.55', rotation=90,size=7)
ax[1][1].set_xlim([0,1])
ax[1][1].set_xlabel('depletion',fontsize=8)
ax[1][1].set_ylabel('density',fontsize=8)

fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure12.jpg'),dpi='figure',format='jpeg')
plt.close()