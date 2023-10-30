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

figure_dir = os.path.join('.','figures')

if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

lx = ly = 1.0
ncol = nrow = 50
nz = 2500
nz_hyper = 100

#hyperparameters details (k, aa and af)
k_m = 1.0
sd_m = 2.0
a_m = 0.5
aa_m = np.pi * (45)/180.0
af_m = 5.0
sd_aa = np.pi *(30)/180.0
a_aa = 0.5

lz = np.sqrt(nz)  # it is the number of standard variates but also the total length in z space
nz = int(lz) # length of m space 
# number of points in m space along the y/x axis

#coordinates of m space
delr = lx/ncol
delc = ly/nrow

x_coord = np.zeros((nrow,ncol))
y_coord = np.zeros((nrow,ncol))

for i in range(nrow):
    for j in range(ncol):
        x_coord[i,j]=delr/2+delr*j
        y_coord[i,j]=ly-(delc/2+delc*i)

#generate random z vector for m space
np.random.seed(seed=8) #reproducibility of true field

z_m = np.random.standard_normal((nz,nz))

#sigma and range values for m space

#define the anisotropy details
#angle respect to north of max range

#we define the number of z standard variates for the anisotropy angle
lz_aa = np.sqrt(nz_hyper)
nz_hyper = int(lz_aa)
#scale factor of anisotropy angle
sf_aa = lz_aa/max(lx,ly)

#sigma and range values for anisotropy angle space

#generate new standard variate vector zaa  
z_aa = np.random.standard_normal((nz_hyper,nz_hyper))

aa = np.zeros((nrow,ncol))

#now we need to scale the kernel (that means scale the range)
a_aa_z = sf_aa * a_aa

sx_coord = x_coord*sf_aa
sy_coord = y_coord*sf_aa

for i in range(nrow):
    for j in range(ncol):
        num=0.0
        xy_m = np.array([sx_coord[i,j],sy_coord[i,j]])
        for k in range(nz_hyper):
            for l in range(nz_hyper):
                d = np.sqrt((xy_m[0]-(l+1))**2.0+(xy_m[1]-(k+1))**2.0)
                num+= (sd_aa * np.exp(-2.*(d/a_aa_z)**2) *(4/((a_aa_z)**2 * 3.1416))**(0.5))*z_aa[k,l]
        aa[i,j]=num+aa_m

#the scale factor
sf = lz/max(lx,ly)

#rescale the m-coord to the z-space
sx_coord = x_coord*sf
sy_coord = y_coord*sf

m1 = np.zeros((nrow,ncol))
m2 = np.zeros((nrow,ncol))

a_z = sf * a_m

for i in range(nrow):
    for j in range(ncol):
        #rotation matrix 
        R = np.array([[np.cos(aa[i,j]), np.sin(aa[i,j])], [-np.sin(aa[i,j]), np.cos(aa[i,j])]])
        #scaling matrix
        S = np.array([[1.0, 0], [0, af_m]])
        H = ((S@R).T)@(S@R)
        num=0.0
        xy_coord = np.array([sx_coord[i,j],sy_coord[i,j]])
        for k in range(nz):
            for l in range(nz):
                #calculate the distance between m and zi
                d = np.sqrt((xy_coord-np.array([l+1,k+1]))@H@(xy_coord-np.array([l+1,k+1])))
                num+= (sd_m * np.exp(-2.*(d/a_z)**2) *(4/((a_z)**2 * 3.1416))**(0.5))*z_m[k,l]
        m2[i,j]=num+np.log10(k_m)

for i in range(nrow):
    for j in range(ncol):
        #rotation matrix 
        R = np.array([[np.cos(aa_m), np.sin(aa_m)], [-np.sin(aa_m), np.cos(aa_m)]])
        #scaling matrix
        S = np.array([[1.0, 0], [0, af_m]])
        H = ((S@R).T)@(S@R)
        num=0.0
        xy_coord = np.array([sx_coord[i,j],sy_coord[i,j]])
        for k in range(nz):
            for l in range(nz):
                #calculate the distance between m and zi
                d = np.sqrt((xy_coord-np.array([l+1,k+1]))@H@(xy_coord-np.array([l+1,k+1])))
                num+= (sd_m * np.exp(-2.*(d/a_z)**2) *(4/((a_z)**2 * 3.1416))**(0.5))*z_m[k,l]
        m1[i,j]=num+np.log10(k_m)

fig, ax = plt.subplots(1,2,dpi=200,figsize=(10,5))
#ax[0].set_aspect('equal', 'box')
k_map = ax[0].imshow(m1,cmap='jet_r',vmin=-2,vmax=2,alpha=0.55)
k_map = ax[1].imshow(m2,cmap='jet_r',vmin=-2,vmax=2,alpha=0.55)
#cbar = fig.colorbar(k_map,orientation="vertical",shrink=0.6)
#cbar.set_label("True log10K",fontsize=7)
#cbar.ax[0].tick_params(labelsize=7)

ax[0].text(0,-2,'a)')
ax[1].text(0,-2,'b)')

#fig.tight_layout()
fig.savefig(os.path.join(figure_dir,'Figure2.pdf'),format='pdf')
plt.close(fig)