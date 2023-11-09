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

run_name ='cal1a_new_reg6'
phi_meas_filename = run_name+'.phi.actual.csv'
figure_dir = os.path.join('.','figures')
files_dir = os.path.join('.//history-matching',run_name)

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

with open(os.path.join(files_dir,run_name+'.rec'),'r') as f:
    rec_lines = f.readlines()

iterations = []
reg_factors = []
measured_phis= []
model_runs = []

for line in rec_lines:
    if "OPTIMISATION ITERATION NO.        :" in line:
        line_split = line.split(':')
        iteration = line_split[1]
        iterations.append(int(iteration))
    if "Current regularisation weight factor                     :" in line:
        line_split = line.split(':')
        reg_factor = line_split[1]
        reg_factors.append(float(reg_factor))
    if "Current value of measurement objective function          :" in line:
        line_split = line.split(':')
        measured_phi = line_split[1]
        measured_phis.append(float(measured_phi))
    if 'Model calls so far' in line:
        line_split = line.split(':')
        imodel_runs = line_split[1]
        model_runs.append(int(imodel_runs))
    #if "Sum of squared weighted residuals (ie phi)                =" in line:
    #    line_split = line.split('=')
    #    measured_phi = line_split[1]
    #    measured_phis.append(float(measured_phi))


cm = 1/2.54
fig, ax = plt.subplots(figsize=(18*cm,6.0*cm),dpi=600)
ax.plot(np.array(model_runs),np.array(measured_phis),'ko',markersize=3.5)
ax.set_xlabel('Model runs',fontsize=7)
ax.set_ylabel('log10(phi)',fontsize=7)
x_ticks = np.arange(0,12100,2000)
x_ticks_labels = np.arange(0,14000,2000)
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_ticks_labels,fontsize=7)
ax.set_xlim([-100,12000])
ax.set_yscale('log')
ax.set_ylim([400,10000000])
ax.set_yticks([1000,10000,100000,1000000,10000000])
ax.tick_params(axis='both', which='major', labelsize=7)
ax.axhline(y=1000,linestyle='dashed',color='grey')
ax.text(0,500,'target measurement objective function',fontsize=8)

fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure8.jpg'),dpi='figure',format='jpeg')
plt.close()