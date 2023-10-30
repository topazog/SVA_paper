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
modelname = 'sva_tm'
tmodelname = 'sva_tm_t'
phi_meas_filename = run_name+'.phi.actual.csv'
figure_dir = os.path.join('.','figures')
files_dir = os.path.join('.//history-matching',run_name)

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

#plot calibrated field and conc distribution ---------------------------

cwd = os.getcwd()

#run PEST with best parameters
os.chdir(files_dir)
os.system("pstclean.exe "+run_name+".pst temp.pst")
os.system('parrep.exe '+run_name+".par.6 temp.pst temp2.pst 0")
os.system('pest.exe temp2.pst')
os.chdir(cwd)

#load calibrated concentrations
fpth = os.path.join(files_dir,tmodelname+'_.obs.csv')
cal_con_sim = pd.read_csv(fpth,index_col=0)


#load modflow models and get geometry details
fsim = flopy.mf6.MFSimulation.load(modelname, "mf6", "mf6", files_dir,verbosity_level=0)
fmdl = fsim.get_model(modelname)
tsim = flopy.mf6.MFSimulation.load(tmodelname, "mf6", "mf6", files_dir,verbosity_level=0)
tmdl = tsim.get_model(tmodelname)
nrow = fmdl.dis.nrow.get_data()
ncol = fmdl.dis.ncol.get_data()
delc = fmdl.dis.delc.get_data()[0]
delr = fmdl.dis.delr.get_data()[0]
Lx = ncol*delr
Ly = nrow*delc

#post process results
mf6_out_path = tsim.simulation_data.mfpath.get_sim_path()
ucnobj_mf6 = tmdl.output.concentration()
conc_mf6 = ucnobj_mf6.get_alldata()
#load the npf y wells packages
npf = fmdl.get_package("npf")
wel1 = fmdl.get_package("wel1")
wel2 = fmdl.get_package("wel2")

fig, ax = plt.subplots(1,1,dpi=200,figsize=(5,5))
ax.set_aspect('equal', 'box')
#ax.set_xticks([])
#ax.set_yticks([])
levels = [0.1,0.5,0.9,0.99]
mdl_true = tsim.get_model(modelname)
mapview = flopy.plot.PlotMapView(model=tmdl,layer=0,ax=ax)
k_map = mapview.plot_array(a=np.log10(npf.k.array[0]),cmap='jet_r',alpha=0.4,vmin=-2,vmax=2)
contour_set=mapview.contour_array(a=conc_mf6[299][0], levels=levels,colors='black',alpha= 0.9,linewidths=0.5)
mapview.plot_bc('wel1',wel1,color='red', kper=1)
mapview.plot_bc('wel1',wel1,color='magenta', kper=0)
mapview.plot_bc('wel2',wel2,color='black')

# get obs information from mf6 observation package
obs_dict = tmdl.conc_obs.continuous.get_data()
obs = []
for value in obs_dict.values():
    for ivalue in value:
        x_coord=delr/2+delr*ivalue[2][2]
        y_coord=Ly-(delc/2+delc*ivalue[2][1])
        obs.append([ivalue[0],x_coord,y_coord,ivalue[2][0]+1])

well_lst = [3,7,8,9]

for i in well_lst:
    label = 'pw'+str(i+1)
    y_label = obs[i][2]
    x_label = obs[i][1]
    ax.text(x_label,y_label,label, ha='center',va='top')

#cbar = fig.colorbar(k_map,orientation="vertical",shrink=0.6)
#cbar.set_label("Cal log10K",fontsize=7)
#cbar.ax.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure9.pdf'),format='pdf')
plt.close(fig)