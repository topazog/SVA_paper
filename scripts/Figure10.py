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

run_name='linear_post_new6_ies'
files_dir = os.path.join('.//history-matching',run_name)
pest_run_name='cal1a_new_reg6'
pest_files_dir = os.path.join('.//history-matching',pest_run_name)
phi_meas_filename = 'linear_post_new6_ies.phi.actual.csv'
modelname = 'sva_tm'
tmodelname = 'sva_tm_t'
modelfiles_dir = os.path.join('.',modelname)

figure_dir = os.path.join('.','figures')
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

#plot hydrographs from non-linear analysis ------------------------------------------------

#load phi_actual csv file

phi_df = pd.read_csv(os.path.join(files_dir,phi_meas_filename))
phi_df.dropna(how='all', axis=1, inplace=True)
phi_df.drop(columns=['iteration','total_runs','mean','standard_deviation','min','max'],inplace=True)
#phi_df = np.log10(phi_df)
phi_t = phi_df.T

phi_t_values = phi_t.values
niter = len(phi_t.columns)

mask = ~np.isnan(phi_t_values)
# = [d[m] for d, m in zip(phi_t_values.T, mask.T)]

int_length = 600.0
ur_interval = int_length/5

r_interval = int(np.round(ur_interval/10**(np.round(np.log10(ur_interval),0)),0)*10**(np.round(np.log10(ur_interval),0)))

x_ticks = np.arange(1,int_length,r_interval)
x_ticks_labels = np.arange(0,int_length,r_interval)

sim_df = pd.read_csv(os.path.join(files_dir,run_name+'.0.obs.csv'))
sim_T = sim_df.T

obs_time_thresh = 200.0

#load observation conc data and add noise to it
np.random.seed(seed=4) #the same seed as the original process
fpth = os.path.join(modelfiles_dir,tmodelname+'_.obs.csv')
df = pd.read_csv(fpth,index_col=0)
df_obs = df[df.index<=obs_time_thresh]
df_pred = df[df.index>obs_time_thresh]
nobs_times = len(df_obs)

#add noise to true values
sd = 0.01
noise = np.random.normal(0,0.01,(len(df.columns),nobs_times))

for i,column in enumerate(df.columns):
    df_obs[column]+=noise[:][i]

with open(os.path.join(files_dir,"concobs_.dat"),'r') as f:
    lines = f.readlines()
f.close()

obs_data = []

for i,line in enumerate(lines):
    obs_name = line.split()[0]
    obs_value = line.split()[3]
    obs_data.append([obs_name+'_'+str(i+1),float(obs_value),float(100.0),'conc'])

len_obs = len(obs_data)

with open(os.path.join(files_dir,"p_concobs_.dat"),'r') as f:
    lines = f.readlines()
f.close()

for i,line in enumerate(lines):
    obs_name = line.split()[0]
    obs_value = line.split()[3]
    obs_data.append([obs_name+'_'+str(len_obs+i+1),float(obs_value),float(0.0),'pconc'])

columns = ['obsnme', 'obsval', 'weight','obgnme']
    
obs_data_df = pd.DataFrame(obs_data, columns=columns)

#load calibrated concentrations
fpth = os.path.join(pest_files_dir,tmodelname+'_.obs.csv')
cal_con_sim = pd.read_csv(fpth,index_col=0)

fig, ax = plt.subplots(2,2,figsize=(5,5),dpi=300)

i=0
j=0

time = df.index

well_list = ['mw4','mw8','mw9','mw10']
well_list_names = ['pw4','pw8','pw9','pw10']

icount=0

for k,well in enumerate(well_list):
    sim_T2 = sim_T[sim_T.index.str.contains(str.lower(well)+'_')]
    for l,real in enumerate(sim_T2.columns):
        if phi_t_values[l]<3797:
            icount+=1
            ax[i][j].plot(time,sim_T2[real],color='black',alpha=0.1,linewidth=1)
    iobs_data = obs_data_df['obsval'][(obs_data_df['obsnme'].str.contains(str.lower(well)+'_')) | (obs_data_df['obsnme'].str.contains(str.upper(well)+'_'))]
    ax[i][j].plot(time,cal_con_sim[str.upper(well)],color='red',alpha=0.8,linewidth=1)
    ax[i][j].plot(time,iobs_data,'bo',markersize=1)
    ax[i][j].set_ylim([0,1.0])
    ax[i][j].set_xlim([0,int_length])
    ax[i][j].set_xticks(x_ticks)
    ax[i][j].set_xticklabels(x_ticks_labels.astype(int))
    ax[i][j].tick_params(axis='both', which='major', labelsize=7)
    ax[i][j].axvline(x=200.0,linestyle='dashed',color='grey')
    ax[i][j].set_title(well_list_names[k],fontsize=7)
    ax[i][j].set_xlabel('time',fontsize=7)
    ax[i][j].set_ylabel('depletion',fontsize=7)
    j+=1
    if j==2:
        i+=1
        j=0
fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure10.pdf'),dpi='figure',format='pdf')
plt.close()













