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

#4 random calibrated fields are plotted
sim_df = pd.read_csv(os.path.join(iesfiles_dir,ies_run_name+'.{0}.par.csv'.format(str(niter-1))))

ireals = [0,253,44,109]

fig, ax = plt.subplots(2,2,figsize=(10,10),dpi=300)
i=0
j=0
for ireal in ireals:
    sim_df_filt = sim_df.iloc[ireal, :]
    k_m=sim_df_filt['k_m']
    sill_m=sim_df_filt['sill_m']
    range_m=sim_df_filt['range_m']
    af_m=sim_df_filt['af_m']
    aa_m=sim_df_filt['aa_m']
    z = []
    z_aa = []
    for k in range(nz):
        z.append(sim_df_filt['z'+str(k+1)])
    for k in range(nz_hyper):
        z_aa.append(sim_df_filt['z_aa'+str(k+1)])
    with open(os.path.join(iesfiles_dir,'mean.dat'),'w') as f:
        f.write(str(k_m)+'\n')
    f.close()
    with open(os.path.join(iesfiles_dir,'hyperparameter_specs.dat'),'w') as f:
        f.write('sill {0} 0.5 0.4\n'.format(str(sill_m)))
        f.write('range {0} 0.6 0.9\n'.format(str(range_m)))
        f.write('af {0} 0.6 0.9\n'.format(str(af_m)))
        f.write('aa {0} 0.27 0.5\n'.format(str(aa_m)))
    f.close()
    with open(os.path.join(iesfiles_dir,'z.dat'),'w') as f:
        for k in range(nz):
            f.write(str(z[k])+'\n')
    f.close()
    with open(os.path.join(iesfiles_dir,'z_aa.dat'),'w') as f:
        for k in range(nz_hyper):
            f.write(str(z_aa[k])+'\n')
    f.close()

    os.chdir(iesfiles_dir)
    os.system('ncpar2d_sva < ncpar2d_sva.in > nul')
    os.chdir(cwd)

    k = np.log10(np.loadtxt(os.path.join(iesfiles_dir,'k.dat')))

    sim.run_simulation(silent=True)

    ucnobj_mf6 = tmdl.output.concentration()
    conc_mf6 = ucnobj_mf6.get_alldata()
    wel1 = fmdl.get_package("wel1")
    wel2 = fmdl.get_package("wel2")
    mapview = flopy.plot.PlotMapView(model=tmdl,layer=0,ax=ax[i,j])
    k_map = mapview.plot_array(a=k,cmap='jet_r',alpha=0.35,vmin=-2,vmax=2)
    levels = np.array([0.1,0.5,0.9,0.99])
    mapview.contour_array(a=conc_mf6[len(time)-1][0], levels=levels,colors='black',alpha= 0.9,linewidths=0.5)
    mapview.plot_bc('wel1',wel1,color='red', kper=1)
    mapview.plot_bc('wel1',wel1,color='magenta', kper=0)
    mapview.plot_bc('wel2',wel2,color='black')

    for k in range(len(obs)-3,len(obs)):
        label = 'pw'+str(k+1)
        y_label = obs[k][2]
        x_label = obs[k][1]
        ax[i,j].text(x_label,y_label,label, ha='center',va='top')

    label = 'pw'+str(3+1)
    y_label = obs[3][2]
    x_label = obs[3][1]
    ax[i,j].text(x_label,y_label,label, ha='center',va='top')

    ax[i,j].set_aspect('equal', 'box')

    j+=1
    if j==2:
        j=0
        i=1
fig.savefig(os.path.join(figure_dir,'Figure6.pdf'),dpi='figure',format='pdf')
plt.close()