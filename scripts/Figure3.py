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

figure_dir = os.path.join('.','figures')

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

modelname = 'sva_tm'
tmodelname = 'sva_tm_t'
modelfiles_dir = os.path.join('.',modelname)
ciwell = 1.0

#load modflow models and get geometry details
fsim = flopy.mf6.MFSimulation.load(modelname, "mf6", "mf6", modelfiles_dir,verbosity_level=0)
fmdl = fsim.get_model(modelname)
tsim = flopy.mf6.MFSimulation.load(tmodelname, "mf6", "mf6", modelfiles_dir,verbosity_level=0)
tmdl = tsim.get_model(tmodelname)
nrow = fmdl.dis.nrow.get_data()
ncol = fmdl.dis.ncol.get_data()
nlay = fmdl.dis.nlay.get_data()
delc = fmdl.dis.delc.get_data()[0]
delr = fmdl.dis.delr.get_data()[0]
Lx = ncol*delr
Ly = nrow*delc
#load the npf and wells
npf = fmdl.get_package("npf")
wel1 = fmdl.get_package("wel1")
wel2 = fmdl.get_package("wel2")

#s = flopy.mf6.MFSimulation.load('mfsim.nam', "mf6", "mf6",modelfiles_dir,verbosity_level=0)
#tdis = s.get_package('tdis')
#print(tdis)

kh_filename = eval(npf.k.get_file_entry().split()[2])

#post process results
mf6_out_path = tsim.simulation_data.mfpath.get_sim_path()
ucnobj_mf6 = tmdl.output.concentration()
conc_mf6 = ucnobj_mf6.get_alldata()

#load observation conc data
fpth = os.path.join(modelfiles_dir,tmodelname+'_.obs.csv')
df = pd.read_csv(fpth,index_col=0)

# get obs information from mf6 observation package
obs_dict = tmdl.conc_obs.continuous.get_data()
obs = []
obs_rcl = []
for value in obs_dict.values():
    for ivalue in value:
        x_coord=delr/2+delr*ivalue[2][2]
        y_coord=Ly-(delc/2+delc*ivalue[2][1])
        obs.append([ivalue[0],x_coord,y_coord,ivalue[2][0]+1])
        obs_rcl.append([ivalue[2][0],ivalue[2][1],ivalue[2][2]])

fig, ax = plt.subplots(1,2,dpi=200,figsize=(10,5))
ax[0].set_aspect('equal', 'box')
ax[0].tick_params(labelsize=8)
ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
ax[0].text(0,1.05,'a)')
#ax.set_xticks([])
#ax.set_yticks([])
levels = np.array([0.1,0.5,0.9,0.99])*ciwell
mapview = flopy.plot.PlotMapView(model=tmdl,layer=0,ax=ax[0])
a=np.log10(npf.k.array[0])
a_mean = a.mean()
a_std = a.std()
k_map = mapview.plot_array(a=np.log10(npf.k.array[0]),cmap='jet_r',alpha=0.35,vmin=-2.0,vmax=2.0)
mapview.contour_array(a=conc_mf6[len(df)-1][0], levels=levels,colors='black',alpha= 0.9,linewidths=0.5)
mapview.plot_bc('wel1',wel1,color='red', kper=1)
mapview.plot_bc('wel1',wel1,color='magenta', kper=0)
mapview.plot_bc('wel2',wel2,color='black')

for i in range(len(obs)):
    label = 'pw'+str(i+1)
    y_label = obs[i][2]
    x_label = obs[i][1]
    ax[0].text(x_label,y_label,label, ha='center',va='top',fontsize=7)

#cbar = fig.colorbar(k_map,orientation="vertical",shrink=0.6)
#cbar.set_label("True log10K",fontsize=5)
#cbar.ax.tick_params(labelsize=5)

#fig.savefig(os.path.join(figure_dir,modelname+'_true_conc_map.png'),format='png')
#plt.close(fig)

#fig, ax = plt.subplots(1,1,dpi=200,figsize=(5,5))
ax[1].tick_params(labelsize=8)
ax[1].text(-20,1.1,'b)')
ax[1].text(0,-0.03,'history-matching period',fontsize=7)
ax[1].text(210,-0.03,'predictive period',fontsize=7)
ax[1].axvline(x=200.0,linestyle='dashed',color='grey')
sel_wells = [3,7,8,9]
for i in range(len(obs_rcl)):
    label = 'pw'+str(i+1)
    if i in sel_wells:
        ax[1].plot(df.index,df['mw'+str(i+1)],label=label,linewidth=2)
    else:
        ax[1].plot(df.index,df['mw'+str(i+1)],label=label,linewidth=0.7)
ax[1].set_xlabel('Time units')
ax[1].set_ylabel('Brine depletion')
ax[1].legend(loc='upper left',fontsize=6)
fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure3.pdf'),format='pdf')
plt.close(fig)