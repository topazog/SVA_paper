#General workflow for evaluating spatial-varying anisotropy
#Model is generated in MF6

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

def generate_true_k_sva_field(true_k_dir,
                              nrow,
                              ncol,
                              lx,
                              ly,
                              k_m,
                              sd_m,
                              a_m,
                              aa_m, 
                              af_m,
                              nz,
                              sd_aa,
                              a_aa,
                              nz_hyper,
                              seed_number=8):

    """function that generates a random hydraulic conductivity field assuming spatial variability
    of the anisotropy angle, using a squared gaussian kernel as a spatial averaging function.

    Args:
        true_k_dir (str): directory where k realization will be saved
        nrow (int): number of rows for the k grid
        ncol (int): number of cols for the k grid
        lx (float): length along x axis
        ly (float): length along y axis
        sd_m (float): standard deviation of the k field
        a_m (float): correlation range of k field for the squared gaussian kernel
        aa_m (float): mean anisotropy angle (respect to north of max range) of k field in radians
        af_m (float): mean anisotropy factor (max range/min range) of k field
        nz (int): number of k-standard variates
        sd_aa (float): standard deviation of the anisotropy angle (aa) field
        a_aa (float): correlation range of the anisotropy angle (aa) field
        nz_hyper (int): number of aa-standard variates
        seed_number (int, optional): seed number for reproducibility purposes. Defaults to 8.

    Returns:
        t: output filename with full path
    """

    if not os.path.exists(true_k_dir):
        os.makedirs(true_k_dir)

    seed = np.random.seed(seed_number) #reproducibility

    lz = np.sqrt(nz)  # it is the number of standard variates but also the total length in z space
    nz = int(lz)
    # length of m space
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
    z_m = np.random.standard_normal((nz,nz))

    #generate z parameters (standard variates) for model parameters
    with open(os.path.join(true_k_dir,'z_true.dat'),'w') as f:
        for i in range(nz):
            for j in range(nz):
                f.write(str(z_m[i][j])+'\n')
    f.close()

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

    #generate z parameters (standard variates) for model parameters
    with open(os.path.join(true_k_dir,'z_aa_true.dat'),'w') as f:
        for i in range(nz_hyper):
            for j in range(nz_hyper):
                f.write(str(z_aa[i][j])+'\n')
    f.close()

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

    m = np.zeros((nrow,ncol))

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
            m[i,j]=num+np.log10(k_m)

    true_k_fpth = os.path.join(true_k_dir,'true_k_sva.dat')

    if not os.path.exists(true_k_dir):
        os.makedirs(true_k_dir)

    np.savetxt(true_k_fpth,np.power(10,m))

    #fig, ax = plt.subplots(1,1,dpi=200,figsize=(5,5))
    #ax.set_aspect('equal', 'box')
    #k_map = ax.imshow(m,cmap='jet_r',vmin=np.log10(k_m)-1.5*sd_m,vmax=np.log10(k_m)+1.5*sd_m)
    #cbar = fig.colorbar(k_map,orientation="vertical",shrink=0.6)
    #cbar.set_label("True log10K",fontsize=7)
    #cbar.ax.tick_params(labelsize=7)

    #fig.savefig(os.path.join(true_k_dir,'true_k_field.png'),format='png')
    #plt.close(fig)

    return true_k_fpth

if __name__ == "__main__":
    true_k_dir = os.path.join('.','true_k_field')
    #model details
    nlay = 1
    nrow = 50
    ncol = 50
    lx = ly = 1.0
    #hyperparameters details (k, aa and af)
    k_m = 1.0
    sd_m = 2.0
    a_m = 0.5
    aa_m = np.pi * (45)/180.0
    af_m = 5.0
    sd_aa = np.pi * (30)/180.0
    a_aa = 0.5
    #total number of standard variates for true k generation
    nz = 2500 #50x50
    nz_hyper = 100 #10x10
    generate_true_k_sva_field(true_k_dir,nrow,ncol,lx,ly,
                                        k_m,sd_m,a_m,aa_m, af_m,nz,sd_aa,
                                        a_aa,nz_hyper,seed_number=8)