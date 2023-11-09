#basic libs
import numpy as np
import os
import shutil
import pandas as pd
from pathlib import Path
import random
import time as tm
from scipy.stats import halfnorm

#plotting libs
from matplotlib import pyplot as plt

#flopy
import flopy

import sys
sys.path.append('..//PEST_utils')
import PEST_utils

ies_run_name ='cal1a_ies'
phi_meas_filename = ies_run_name+'.phi.actual.csv'
figure_dir = os.path.join('.','figures')
iesfiles_dir = os.path.join('.//history-matching',ies_run_name)

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

#plot phi
phi_df = pd.read_csv(os.path.join(iesfiles_dir,phi_meas_filename))
phi_df.dropna(how='all', axis=1, inplace=True)
phi_df.drop(columns=['iteration','total_runs','mean','standard_deviation','min','max'],inplace=True)
#phi_df = np.log10(phi_df)
phi_t = phi_df.T

phi_t_values = phi_t.values
niter = len(phi_t.columns)

mask = ~np.isnan(phi_t_values)
filtered_data = [d[m] for d, m in zip(phi_t_values.T, mask.T)]

phi_df = pd.read_csv(os.path.join(iesfiles_dir,phi_meas_filename))
cm = 1/2.54
fig, ax = plt.subplots(figsize=(18*cm,6.0*cm),dpi=600)
ax.boxplot(filtered_data,positions=np.array(phi_df['total_runs']),widths=100)
x_ticks = np.arange(0,12100,2000)
x_ticks_labels = np.arange(0,14000,2000)
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_ticks_labels,fontsize=7)
ax.set_xlim([0,10000])


ax.axhline(y=1000,linestyle='dashed',color='grey')
ax.text(100,500,'target measurement objective function',fontsize=7)

ax.set_yscale('log')
#ax.grid(which='both',axis='y')
ax.set_yticks([1000,10000,100000,1000000,10000000])
ax.tick_params(axis='both', which='major', labelsize=7)
ax.set_ylim([400,10000000])
ax.set_xlabel('Model runs',fontsize=7)
ax.set_ylabel('log10(phi)',fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure5.jpg'),dpi='figure',format='jpeg')

plt.close()