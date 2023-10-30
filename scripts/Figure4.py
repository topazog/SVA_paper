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

num_reals = 30000

np.random.seed(seed=4) #reproducibility

range_m = np.power(10,np.random.normal(-0.6,0.3, size=num_reals))

af_m = halfnorm.rvs(loc = 1, scale = 2.5, size=num_reals)
aa_m = np.random.normal(0.0,0.6,num_reals)

fig, ax = plt.subplots(1,3,dpi=200,figsize=(12,5))


ax[2].hist(range_m,bins=100,color='grey',density=True)
ax[2].tick_params(axis='both', which='major', labelsize=8)
ax[2].axvline(x=0.5,color='black',linestyle='dashed')
ax[2].text(0.53, 1.5, 'True range = 0.5', rotation=90,size=8)
ax[2].set_xlabel('a value',fontsize=10)
ax[2].set_xlim([0,2])
ax[2].set_ylabel('density',fontsize=10)

ax[1].hist(af_m,bins=50,color='grey',density=True)
ax[1].tick_params(axis='both', which='major', labelsize=8)
ax[1].axvline(x=5.0,color='black',linestyle='dashed')
ax[1].text(5.3, 0.15, 'True anisotropy ratio = 5.0', rotation=90,size=8)
ax[1].set_xlabel('anisotropy ratio',fontsize=10)
ax[1].set_ylabel('density',fontsize=10)

ax[0].hist(aa_m,bins=50,color='grey',density=True)
ax[0].tick_params(axis='both', which='major', labelsize=8)
ax[0].axvline(x=np.pi*45/180,color='black',linestyle='dashed')
ax[0].text(np.pi*45/180+0.1, 0.3, 'True anisotropy angle = '+str(np.round(np.pi*45/180,2)), rotation=90,size=8)
ax[0].set_xlabel('anisotropy angle (radians)',fontsize=10)
ax[0].set_ylabel('density',fontsize=10)

fig.tight_layout()

fig.savefig(os.path.join(figure_dir,'Figure4.pdf'),format='pdf')
plt.close(fig)