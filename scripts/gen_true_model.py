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
sys.path.append('..//PEST_utils')
import PEST_utils

def sva_tm_setup(modelname,
                 modelfiles_dir,
                 true_k_dir,
                 lx,
                 ly,
                 nrow,
                 ncol,
                 nlay,
                 top,
                 bot,
                 tdis_rc,
                 qit,
                 ciwell,
                 al,
                 at,
                 por,
                 iconc,
                 mixelm = -1):
    
    #TODO : add CONTINUE option in MFSIM

    # working folder
    if not os.path.exists(modelfiles_dir):
        os.makedirs(modelfiles_dir)

    #spatial discretization
    delr = lx / ncol
    delc = ly / nrow
    #temp discretization
    nper = len(tdis_rc)

    #pumping wells locations
    pw_0_loc = [(0,int(nrow/3),int(ncol/3)),
          (0,2*int(nrow/3),int(ncol/3)),
          (0,int(nrow/6),int(ncol/2)),
          (0,int(nrow/3),ncol-int(ncol/3)),
          (0,2*int(nrow/3),ncol-int(ncol/3)),
          (0,int(nrow/2),int(ncol/2)),
          (0,nrow-int(nrow/5),int(ncol/2))]

    pw_1_loc = [(0,int(nrow/3),int(ncol/3)),
          (0,2*int(nrow/3),int(ncol/3)),
          (0,int(nrow/6),int(ncol/2)),
          (0,int(nrow/3),ncol-int(ncol/3)),
          (0,2*int(nrow/3),ncol-int(ncol/3)),
          (0,int(nrow/2),int(ncol/2)),
          (0,nrow-int(nrow/5),int(ncol/2)),
          (0,3,ncol-12),
          (0,int(nrow/2)+1,2),
          (0,nrow-5,ncol-5)]

    nwells = [len(pw_0_loc),len(pw_1_loc)]
    #injection rate and conc from each of 2 wells
    qiwell = qit/2.0
    #pumping rate from each of N wells
    qpwell = [qiwell*2.0/nwells[0],qiwell*2.0/nwells[1]]

    iw_loc = [(0,int(nrow/3),int(ncol/2)),(0,2*int(nrow/3),int(ncol/2))]
    pwells_0 = []
    pwells_1 = []
    iwells = []
    for i in range(len(pw_0_loc)):
        pwells_0.append([pw_0_loc[i],-qpwell[0],0.0])
    for i in range(len(pw_1_loc)):
        pwells_1.append([pw_1_loc[i],-qpwell[1],0.0])
    for i in range(len(iw_loc)):
        iwells.append([iw_loc[i],qiwell,ciwell])

    #some general solver parameters
    nouter, ninner = 1500, 10
    hclose, rclose, relax = 1e-6, 1e-7, 1.0

    mf6_exe = flopy.which('mf6')

    #create simulation
    sim = flopy.mf6.MFSimulation(sim_name=modelname, exe_name=mf6_exe,
                                version='mf6', sim_ws=modelfiles_dir)
    #create tdis package
    tdis = flopy.mf6.ModflowTdis (sim,nper=nper, pname='tdis',time_units='DAYS',
                                perioddata=tdis_rc)
    #create the gwf model
    gwf = flopy.mf6.ModflowGwf (sim, modelname = modelname,
                                model_nam_file='{}.nam'.format(modelname),
                            newtonoptions='',save_flows=False)
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, 
                                delr=delr, delc=delc,
                                top=top, botm=bot)
    # Create the Flopy iterative model solver (ims) Package object
    ims = flopy.mf6.ModflowIms(sim, pname='ims', print_option='SUMMARY', complexity='COMPLEX', outer_hclose=hclose, 
                            outer_maximum=nouter, under_relaxation='DBD',under_relaxation_theta=0.7, under_relaxation_kappa=0.07,
                            under_relaxation_gamma=0, under_relaxation_momentum=0,backtracking_number=0, backtracking_tolerance=1.1,
                            backtracking_reduction_factor=0, backtracking_residual_limit=0, inner_maximum=ninner, inner_hclose=rclose, 
                            rcloserecord=0.00001, linear_acceleration='BICGSTAB',preconditioner_levels=7, preconditioner_drop_tolerance=1e-05,
                            number_orthogonalizations=0, scaling_method='NONE', reordering_method='NONE', relaxation_factor=0.97,
                            filename="{}.ims".format(modelname))
    sim.register_ims_package(ims, [gwf.name])

    #create the wel packages
    wel1 = flopy.mf6.ModflowGwfwel(gwf, pname='wel1', print_input=True, print_flows=True,
                                auxiliary="CONCENTRATION", stress_period_data={0:pwells_0,1:pwells_1},
                                save_flows=False,filename="{}.wel1".format(modelname))
    wel2 = flopy.mf6.ModflowGwfwel(gwf, pname='wel2', print_input=True, print_flows=True,
                                auxiliary="CONCENTRATION", stress_period_data={0:iwells,1:iwells},
                                save_flows=False,filename="{}.wel2".format(modelname))

    #create the npf package
    shutil.copy(os.path.join(true_k_dir,'true_k_sva.dat'),os.path.join(modelfiles_dir,'k.dat'))
    k11_fname = ['k.dat']
    npf = flopy.mf6.ModflowGwfnpf(gwf, pname='npf', save_flows=True, k=k11_fname,k33overk=False,
                                k33=k11_fname, icelltype=0,filename="{}.npf".format(modelname))
    #create the sto package
    sto = flopy.mf6.ModflowGwfsto(gwf, pname='sto', save_flows=True, ss=0.0, sy=0.0,
                                filename="{}.sto".format(modelname))
    #create the initial cond package
    ic = flopy.mf6.ModflowGwfic(gwf, pname='ic', strt=1.0,filename='{}.ic'.format(modelname))
    #create the oc package
    budget_file = modelname + ".bud"
    head_file = modelname + ".hds"
    oc = flopy.mf6.ModflowGwfoc(gwf,budget_filerecord=budget_file,head_filerecord=head_file, 
                                saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')],printrecord=[('BUDGET', 'ALL')])

    #now for transport

    #create the gwf model
    tmodelname = modelname+'_t'
    gwt = flopy.mf6.MFModel (sim, modelname = tmodelname,model_type="gwt6",
                                model_nam_file='{}.nam'.format(tmodelname))

    gwt.name_file.save_flows = True

    # Create the Flopy iterative model solver (ims) Package object
    imsgwt = flopy.mf6.ModflowIms(sim, pname='imst', print_option='SUMMARY', complexity='COMPLEX', outer_hclose=hclose, 
                            outer_maximum=nouter, under_relaxation='DBD',under_relaxation_theta=0.7, under_relaxation_kappa=0.07,
                            under_relaxation_gamma=0, under_relaxation_momentum=0,backtracking_number=0, backtracking_tolerance=1.1,
                            backtracking_reduction_factor=0, backtracking_residual_limit=0, inner_maximum=ninner, inner_hclose=rclose, 
                            rcloserecord=0.00001, linear_acceleration='BICGSTAB',preconditioner_levels=7, preconditioner_drop_tolerance=1e-05,
                            number_orthogonalizations=0, scaling_method='NONE', reordering_method='NONE', relaxation_factor=0.97,
                            filename="{}.ims".format(tmodelname))
    sim.register_ims_package(imsgwt, [gwt.name])
    #dis package for transport
    dis = flopy.mf6.ModflowGwfdis(gwt, nlay=nlay, nrow=nrow, ncol=ncol, 
                                delr=delr, delc=delc,idomain=1,
                                top=top, botm=bot,filename="{}.dis".format(tmodelname))
    #ic for transport
    ict = flopy.mf6.ModflowGwfic(gwt, pname='ict', strt=iconc,filename='{}.ic'.format(tmodelname))
    #advection package
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()
    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename="{}.adv".format(tmodelname))
    #dispersion package
    if al != 0:
        flopy.mf6.ModflowGwtdsp(gwt,xt3d_off=True,alh=al,ath1=at,filename="{}.dsp".format(tmodelname))
    #transport mass storage package
    flopy.mf6.ModflowGwtmst(gwt,porosity=por,first_order_decay=False,decay=None,decay_sorbed=None,
                sorption=None,bulk_density=None,distcoef=None,filename="{}.mst".format(tmodelname))
    #transport source-sink mixing package
    sourcerecarray = [("wel2", "AUX", "CONCENTRATION")]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename="{}.ssm".format(tmodelname))
    #transport output control package
    flopy.mf6.ModflowGwtoc(gwt,budget_filerecord="{}.cbc".format(tmodelname),
                concentration_filerecord="{}.ucn".format(tmodelname),
                concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
                saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
                printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "ALL")])
    #flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(sim,
                exgtype="GWF6-GWT6",exgmnamea=modelname,exgmnameb=tmodelname,
                filename="{}.gwfgwt".format(modelname))
    #add conc observation package
    model_c_obs_pts = []
    for i in range(len(pw_1_loc)):
        model_c_obs_pts.append(['mw'+str(i+1),'concentration',pw_1_loc[i]])

    c_continuous = {tmodelname+"_.obs.csv":model_c_obs_pts}
    c_obs = flopy.mf6.ModflowUtlobs(gwt, continuous=c_continuous, digits=4,pname='conc_obs')

    sim.write_simulation(silent=True)
    sim.run_simulation(silent=True)

    ucnobj_mf6 = gwt.output.concentration()
    conc_mf6 = ucnobj_mf6.get_alldata()

    #for some reason the csv that mf6 generates has an error where the value is near zero
    #the csv file is split and save for further PEST/IES usage, into obs and pred files

    columns = ['time']
    for i in range(len(model_c_obs_pts)):
        columns.append(model_c_obs_pts[i][0])

    time = 0
    cum_t = []
    nstp = 0

    for iper in range(nper):
        delta_t = []
        iperlen = tdis_rc[iper][0]
        instp = tdis_rc[iper][1]
        itsmult = tdis_rc[iper][2]
        nstp +=instp
        if itsmult > 1.0:
            delta_t1 = iperlen*(itsmult-1)/(itsmult**instp-1)
        else:
            delta_t1 = iperlen/instp
        delta_t.append(delta_t1)
        time += delta_t1
        cum_t.append(time)
        for i in range(1,instp):
            idelta_t = delta_t[i-1]*itsmult
            delta_t.append(idelta_t)
            time=time+idelta_t
            cum_t.append(time)

    data = []

    for i in range(len(cum_t)):
        elem = [cum_t[i]]
        for j in range(len(model_c_obs_pts)):
            rcl = model_c_obs_pts[j][2]
            elem.append(conc_mf6[i][rcl[0]][rcl[1]][rcl[2]])
        data.append(elem)
        
    df = pd.DataFrame(data, columns=columns)

    df.to_csv(os.path.join(modelfiles_dir,tmodelname+"_.obs.csv"),index=False,float_format='%.4E')

if __name__ == "__main__":
    modelname = 'sva_tm'
    tmodelname = 'sva_tm_t'
    true_k_dir = os.path.join('.','true_k_field')   
    modelfiles_dir = os.path.join('.',modelname)
    bin_dir = os.path.join('.','bin')
    #model details
    nlay = 1;nrow = 50;ncol = 50
    lx = ly = 1.0;top = 0.0;bot = -1.0
    qit = 0.0005
    ciwell = 1.0
    al = 0.005 #long dispersivity
    at = 0.1*al #trans dispersity
    por = 0.25 #porosity
    iconc = 0.0 #initial concentration
    tdis_rc = [(200.0, 100, 1.0),
            (400.0,200,1.0)]
    sva_tm_setup(modelname,modelfiles_dir,true_k_dir,
                lx,ly,nrow,ncol,nlay,top,bot,tdis_rc,
                qit,ciwell,al,at,por,iconc,mixelm = -1)