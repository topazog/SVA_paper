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

obs_data = pd.read_csv(os.path.join(ies_dir,'obs_data.csv'))

#plot hydrographs

#get times from model
sim = flopy.mf6.MFSimulation.load(modelname, "mf6", "mf6", modelfiles_dir,verbosity_level=0)
tmdl = sim.get_model(tmodelname)
ucnobj_mf6 = tmdl.output.concentration()
time = ucnobj_mf6.get_times()

int_length = time[len(time)-1]
ur_interval = int_length/5

r_interval = int(np.round(ur_interval/10**(np.round(np.log10(ur_interval),0)),0)*10**(np.round(np.log10(ur_interval),0)))

x_ticks = np.arange(1,int_length,r_interval)
x_ticks_labels = np.arange(0,int_length,r_interval)

obs_data = pd.read_csv(os.path.join(ies_dir,'obs_data.csv'),header=0)

cols = np.arange(1,50001)

min_arr = np.ones(3000)
max_arr = np.zeros(3000)

pw_4_end = []
pw_8_end = []
pw_9_end = []
pw_10_end = []


sim_df = pd.read_csv(os.path.join(dsi_dir,'dsi_results.csv'),names=cols,header=None)
sim_df=sim_df.reset_index(drop=True)
obs_data=obs_data.reset_index(drop=True)

sel_cols = (((sim_df.sub(obs_data['obsval'],axis='rows')).mul(obs_data['weight'],axis='rows'))**2).sum(axis=0)<3700
bool_df = pd.DataFrame(sel_cols)
sim_df_f = sim_df.loc[:, bool_df[0]]

pw_4_end.append(sim_df_f.iloc[1799].to_numpy())
pw_8_end.append(sim_df_f.iloc[2599].to_numpy())
pw_9_end.append(sim_df_f.iloc[2799].to_numpy())
pw_10_end.append(sim_df_f.iloc[2999].to_numpy())

min_df_f = sim_df_f.min(axis=1)
max_df_f = sim_df_f.max(axis=1)

min_arr = np.minimum(min_arr,min_df_f.array)
max_arr = np.maximum(max_arr,max_df_f.array)

pw_4_arr = np.hstack(np.array(pw_4_end))
pw_8_arr = np.hstack(np.array(pw_8_end))
pw_9_arr = np.hstack(np.array(pw_9_end))
pw_10_arr = np.hstack(np.array(pw_10_end))

np.save(os.path.join(dsi_dir,'DSI_end_pw4.npy'),pw_4_arr)
np.save(os.path.join(dsi_dir,'DSI_end_pw8.npy'),pw_8_arr)
np.save(os.path.join(dsi_dir,'DSI_end_pw9.npy'),pw_9_arr)
np.save(os.path.join(dsi_dir,'DSI_end_pw10.npy'),pw_10_arr)

sorted_min = np.zeros(len(min_arr))
sorted_max = np.zeros(len(max_arr))
obs_data_sorted = np.zeros(len(obs_data))

for k in range(10):
    sorted_min[300*k:300*k+100]=min_arr[100*k:100*(k+1)]
    sorted_min[300*k+100:300*(k+1)]=min_arr[1000+200*k:1000+200*(k+1)]
    sorted_max[300*k:300*k+100]=max_arr[100*k:100*(k+1)]
    sorted_max[300*k+100:300*(k+1)]=max_arr[1000+200*k:1000+200*(k+1)]
    obs_data_sorted[300*k:300*k+100]=obs_data['obsval'][100*k:100*(k+1)]
    obs_data_sorted[300*k+100:300*(k+1)]=obs_data['obsval'][1000+200*k:1000+200*(k+1)]


fig, ax = plt.subplots(2,2,figsize=(5,5),dpi=300)

i=0
j=0

k_sel = [3,7,8,9]
well_list_names = ['pw4','pw8','pw9','pw10']

int_length = 600.0
ur_interval = int_length/5

r_interval = int(np.round(ur_interval/10**(np.round(np.log10(ur_interval),0)),0)*10**(np.round(np.log10(ur_interval),0)))

x_ticks = np.arange(1,int_length,r_interval)
x_ticks_labels = np.arange(0,int_length,r_interval)

for l,k in enumerate(k_sel):
    ax[i][j].fill_between(time,sorted_min[k*300:(k+1)*300,],sorted_max[k*300:(k+1)*300],color='blue',alpha=0.4,linewidth=0.1)
    ax[i][j].plot(time,obs_data_sorted[k*300:(k+1)*300],'bo',markersize=1)
    ax[i][j].set_ylim([0,1.0])
    ax[i][j].set_title(well_list_names[l],fontsize=7)
    ax[i][j].set_xticks(x_ticks)
    ax[i][j].set_xticklabels(x_ticks_labels.astype(int))
    ax[i][j].set_xlim([0,int_length])
    ax[i][j].tick_params(axis='both', which='major', labelsize=7)
    ax[i][j].axvline(x=200.0,linestyle='dashed',color='grey')
    ax[i][j].set_xlabel('time',fontsize=7)
    ax[i][j].set_ylabel('depletion',fontsize=7)

    j+=1
    if j==2:
        i+=1
        j=0
fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure11.pdf'),dpi='figure',format='pdf')
plt.close()