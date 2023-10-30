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
modelname = 'sva_tm'
tmodelname = 'sva_tm_t'
modelfiles_dir = os.path.join('.',modelname)
nz = 900
nz_hyper = 100

cwd = os.getcwd()

#load modflow models and get geometry details
sim = flopy.mf6.MFSimulation.load(modelname, "mf6", "mf6", iesfiles_dir,verbosity_level=0)
fmdl = sim.get_model(modelname)
tmdl = sim.get_model(tmodelname)
nrow = fmdl.dis.nrow.get_data()
ncol = fmdl.dis.ncol.get_data()
nlay = fmdl.dis.nlay.get_data()
delc = fmdl.dis.delc.get_data()[0]
delr = fmdl.dis.delr.get_data()[0]
Lx = ncol*delr
Ly = nrow*delc
#load the npf
npf = fmdl.get_package("npf")
kh_filename = eval(npf.k.get_file_entry().split()[2])

# get obs information from mf6 observation package
obs_dict = tmdl.conc_obs.continuous.get_data()
obs = []
for value in obs_dict.values():
    for ivalue in value:
        x_coord=delr/2+delr*ivalue[2][2]
        y_coord=Ly-(delc/2+delc*ivalue[2][1])
        obs.append([ivalue[0],x_coord,y_coord,ivalue[2][0]+1])

#load observation conc data to add noise to it
#time observation threshold
t_thresh = 200.0
fpth = os.path.join(modelfiles_dir,tmodelname+'_.obs.csv')
df = pd.read_csv(fpth,index_col=0)
df_obs = df[df.index<=t_thresh]
df_pred = df[df.index>t_thresh]
nobs_times = len(df_obs)
time = df.index
#add noise to true values
sd = 0.01
noise = np.random.normal(0,0.01,(len(df.columns),nobs_times))

for i,column in enumerate(df.columns):
    df_obs[column]+=noise[:][i]

#load ies observed data
obs_data_filename = 'obs_data.csv'
obs_data_df = pd.read_csv(os.path.join(iesfiles_dir,obs_data_filename))

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

#load phi evol file
phi_df = pd.read_csv(os.path.join(iesfiles_dir,phi_meas_filename))
phi_df.dropna(how='all', axis=1, inplace=True)
phi_df.drop(columns=['iteration','total_runs','mean','standard_deviation','min','max'],inplace=True)
#phi_df = np.log10(phi_df)
phi_t = phi_df.T

phi_t_values = phi_t.values
niter = len(phi_t.columns)

#plot hydrographs
int_length = df.index[len(df.index)-1]
ur_interval = int_length/5

r_interval = int(np.round(ur_interval/10**(np.round(np.log10(ur_interval),0)),0)*10**(np.round(np.log10(ur_interval),0)))

x_ticks = np.arange(1,int_length,r_interval)
x_ticks_labels = np.arange(0,int_length,r_interval)

sim_df = pd.read_csv(os.path.join(iesfiles_dir,ies_run_name+'.{0}.obs.csv'.format(str(niter-1))))
sim_T = sim_df.T

fig, ax = plt.subplots(2,2,figsize=(5,5),dpi=300)

i=0
j=0

time = df.index

well_list = ['mw4','mw8','mw9','mw10']
well_list_names = ['pw4','pw8','pw9','pw10']

for k,well in enumerate(well_list):
    sim_T2 = sim_T[sim_T.index.str.contains(str.lower(well)+'_')]
    for real in sim_T2.columns:
        ax[i][j].plot(time,sim_T2[real],color='black',alpha=0.1,linewidth=1)
    iobs_data = obs_data_df['obsval'][(obs_data_df['obsnme'].str.contains(str.lower(well)+'_')) | (obs_data_df['obsnme'].str.contains(str.upper(well)+'_'))]
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
fig.savefig(os.path.join(figure_dir,'Figure7.pdf'),dpi='figure',format='pdf')
plt.close()