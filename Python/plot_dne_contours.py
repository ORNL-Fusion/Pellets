#%%
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import PelletSumfile as PS
import os
import pandas as pd

# plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.serif': 'Times New Roman'})
plt.rcParams.update({'font.family': 'serif'})
# plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.fontset': 'stix'})
plt.rcParams.update({'mathtext.default': 'it'})
# plt.rcParams.update({'mathtext.default': 'cm'})

#%%
base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
gfile = base_dir + "st50010dn.eqdsk"

degrees_plot = ['m20','m10','p00','p10','p20']
v_plot = ['v250','v375','v500','v625','v750','v875','v1000','v1125','v1250','v1375','v1500']
d_plot = ['d5.0','d6.0','d7.0']
drift_path = ['drift_na','drift_DIII-D','drift_HPI2']

npts_deg = len(degrees_plot)
npts_v = len(v_plot)
npts_d = len(d_plot)
npts_dr = len(drift_path)

dne_max = np.zeros((npts_d,npts_v,npts_deg))
dne_depth = np.zeros((npts_d,npts_v,npts_deg))
rho_max = np.zeros((npts_d,npts_v,npts_deg))
rho_depth = np.zeros((npts_d,npts_v,npts_deg))
dne_arr = np.zeros((205,npts_d,npts_deg))

# dne_max = np.zeros((npts_d,npts_v))
# dne_depth = np.zeros((npts_d,npts_v))
# rho_max = np.zeros((npts_d,npts_v))
# rho_depth = np.zeros((npts_d,npts_v))
# dne_arr = np.zeros((70,npts_d))

save_path = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/'
names = ["drift","radius (m)","velocity (m/s)","angle (deg)","rho_depth (-)","dne_depth (1/m**3)","rho_max (-)","dne_max (1/m**3)"]
drifts = ["NA", "DIII-D", "HPI2"]
angles = ["-20","-10","0","10","20"]

#%%
# save as numpy arrays (one per drift) or create dataframe to save as csv later
save = 1
dataframe = 1

rows = []

if save:

    for ii in range(npts_d):
        for jj in range(npts_v):
            count_deg = 0
            for kk in range(npts_deg):
                if dataframe:
                    count_drift = 0
                    for ll in range(npts_dr):
                        row = {}
                        sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/' + degrees_plot[kk] \
                        + '_degrees/' + drift_path[ll] + '/' + d_plot[ii] + '_' + v_plot[jj] + '/sum_pellet_pellet_test_.dat'
                        pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
                        if drift_path[ll] == 'drift_na':
                            dne = pelsum.profiles['delta_ne']
                        else:
                            dne = pelsum.profiles['delta_ne_drift']
                        rho = pelsum.profiles['rho_t']
                        idx_dnem = np.where(dne==max(dne))
                        idx_dned = np.where(dne>=1.0)
                        if np.shape(idx_dned)[1]==0:
                            rho_max[ii,jj,kk] = np.nan
                            dne_max[ii,jj,kk] = np.nan
                            rho_depth[ii,jj,kk] = np.nan
                            dne_depth[ii,jj,kk] = np.nan
                            continue
                        else: 
                            idx_dnem = idx_dnem[0][0]
                            idx_dned = idx_dned[0][0]
                        rho_max[ii,jj,kk] = rho[idx_dnem]
                        dne_max[ii,jj,kk] = dne[idx_dnem]
                        rho_depth[ii,jj,kk] = rho[idx_dned]
                        dne_depth[ii,jj,kk] = dne[idx_dned]
                        if jj==9:
                            dne_arr[:,ii,kk] = dne
                        row["radius (m)"] = pelsum.r_pel
                        row["velocity (m/s)"] = pelsum.v_pel
                        row["angle (deg)"] = angles[count_deg]
                        row["drift"] = drifts[count_drift]
                        row["rho_depth (-)"] = rho[idx_dned]
                        row["dne_depth (1/m**3)"] = dne[idx_dned]
                        row["rho_peak (-)"] = rho[idx_dnem]
                        row["dne_peak (1/m**3)"] = dne[idx_dnem]
                        rows.append(row)
                        # print("file: ", sumfile_loc)
                        # print("drift = ", drifts[count_drift], "angle = ", angles[count_deg], "r_pel = ", pelsum.r_pel, \
                        #       "v_pel = ", pelsum.v_pel)
                        # print("ii = ", ii, "jj = ", jj, "kk = ", kk, "ll = ", ll)
                        count_drift += 1
                else:
                    sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/' + degrees_plot[kk] \
                        + '_degrees/drift_na/' + d_plot[ii] + '_' + v_plot[jj] + '/sum_pellet_pellet_test_.dat'
                    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
                    # dne = pelsum.profiles['delta_ne_drift']
                    dne = pelsum.profiles['delta_ne']
                    rho = pelsum.profiles['rho_t']
                    idx_dnem = np.where(dne==max(dne))
                    idx_dned = np.where(dne>=1.0)
                    if np.shape(idx_dned)[1]==0:
                        rho_max[ii,jj,kk] = np.nan
                        dne_max[ii,jj,kk] = np.nan
                        rho_depth[ii,jj,kk] = np.nan
                        dne_depth[ii,jj,kk] = np.nan
                        continue
                    else: 
                        idx_dnem = idx_dnem[0][0]
                        idx_dned = idx_dned[0][0]
                    rho_max[ii,jj,kk] = rho[idx_dnem]
                    dne_max[ii,jj,kk] = dne[idx_dnem]
                    rho_depth[ii,jj,kk] = rho[idx_dned]
                    dne_depth[ii,jj,kk] = dne[idx_dned]
                    if jj==9:
                        dne_arr[:,ii,kk] = dne
                count_deg += 1
    if dataframe:
        data = pd.DataFrame(rows)
    else:
        np.save(save_path+'dne_max_drift_na', dne_max)
        np.save(save_path+'rho_max_drift_na', rho_max)
        np.save(save_path+'dne_depth_drift_na', dne_depth)
        np.save(save_path+'rho_depth_drift_na', rho_depth)
        np.save(save_path+'dne_arr_drift_na', dne_arr)

#%%
# save as numpy arrays (one per drift) or create dataframe to save as csv later

dne_max = np.zeros((npts_d,npts_v))
dne_depth = np.zeros((npts_d,npts_v))
rho_max = np.zeros((npts_d,npts_v))
rho_depth = np.zeros((npts_d,npts_v))
dne_arr = np.zeros((205,npts_d))

save = 1
dataframe = 1

rows = []

if save:

    for ii in range(npts_d):
        for jj in range(npts_v):
            if dataframe:
                count_drift = 0
                for ll in range(npts_dr):
                    row = {}
                    sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/straight/' + \
                        drift_path[ll] + '/' + d_plot[ii] + '_' + v_plot[jj] + '/sum_pellet_pellet_test_.dat'
                    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
                    if drift_path[ll] == 'drift_na':
                        dne = pelsum.profiles['delta_ne']
                    else:
                        dne = pelsum.profiles['delta_ne_drift']
                    rho = pelsum.profiles['rho_t']
                    idx_dnem = np.where(dne==max(dne))
                    idx_dned = np.where(dne>=1.0)
                    if np.shape(idx_dned)[1]==0:
                        rho_max[ii,jj] = np.nan
                        dne_max[ii,jj] = np.nan
                        rho_depth[ii,jj] = np.nan
                        dne_depth[ii,jj] = np.nan
                        continue
                    else: 
                        idx_dnem = idx_dnem[0][0]
                        idx_dned = idx_dned[0][0]
                    rho_max[ii,jj] = rho[idx_dnem]
                    dne_max[ii,jj] = dne[idx_dnem]
                    rho_depth[ii,jj] = rho[idx_dned]
                    dne_depth[ii,jj] = dne[idx_dned]
                    if jj==9:
                        dne_arr[:,ii] = dne
                    row["radius (m)"] = pelsum.r_pel
                    row["velocity (m/s)"] = pelsum.v_pel
                    row["drift"] = drifts[count_drift]
                    row["rho_depth (-)"] = rho[idx_dned]
                    row["dne_depth (1/m**3)"] = dne[idx_dned]
                    row["rho_peak (-)"] = rho[idx_dnem]
                    row["dne_peak (1/m**3)"] = dne[idx_dnem]
                    rows.append(row)
                    # print("file: ", sumfile_loc)
                    # print("drift = ", drifts[count_drift], "angle = ", angles[count_deg], "r_pel = ", pelsum.r_pel, \
                    #       "v_pel = ", pelsum.v_pel)
                    # print("ii = ", ii, "jj = ", jj, "kk = ", kk, "ll = ", ll)
                    count_drift += 1
            else:
                sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/' + degrees_plot[kk] \
                    + '_degrees/drift_na/' + d_plot[ii] + '_' + v_plot[jj] + '/sum_pellet_pellet_test_.dat'
                pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
                # dne = pelsum.profiles['delta_ne_drift']
                dne = pelsum.profiles['delta_ne']
                rho = pelsum.profiles['rho_t']
                idx_dnem = np.where(dne==max(dne))
                idx_dned = np.where(dne>=1.0)
                if np.shape(idx_dned)[1]==0:
                    rho_max[ii,jj,kk] = np.nan
                    dne_max[ii,jj,kk] = np.nan
                    rho_depth[ii,jj,kk] = np.nan
                    dne_depth[ii,jj,kk] = np.nan
                    continue
                else: 
                    idx_dnem = idx_dnem[0][0]
                    idx_dned = idx_dned[0][0]
                rho_max[ii,jj,kk] = rho[idx_dnem]
                dne_max[ii,jj,kk] = dne[idx_dnem]
                rho_depth[ii,jj,kk] = rho[idx_dned]
                dne_depth[ii,jj,kk] = dne[idx_dned]
                if jj==9:
                    dne_arr[:,ii,kk] = dne
    if dataframe:
        data = pd.DataFrame(rows)
    else:
        np.save(save_path+'dne_max_drift_na', dne_max)
        np.save(save_path+'rho_max_drift_na', rho_max)
        np.save(save_path+'dne_depth_drift_na', dne_depth)
        np.save(save_path+'rho_depth_drift_na', rho_depth)
        np.save(save_path+'dne_arr_drift_na', dne_arr)


#%%

degrees = [-20.0,-10.0,0.0,10.0,20.]
velocity = [250.,375.,500.,625.,750.,875.,1000.,1125.,1250.,1375.,1500.]

rho_max_na = np.load(save_path+'rho_max_drift_na.npy')
rho_depth_na = np.load(save_path+'rho_depth_drift_na.npy')
dne_arr_na = np.load(save_path+'dne_arr_drift_na.npy')
rho_max_HPI2 = np.load(save_path+'rho_max_drift_HPI2.npy')
rho_depth_HPI2 = np.load(save_path+'rho_depth_drift_HPI2.npy')
dne_arr_HPI2 = np.load(save_path+'dne_arr_drift_HPI2.npy')
rho_max_DIIID = np.load(save_path+'rho_max_drift_DIII-D.npy')
rho_depth_DIIID = np.load(save_path+'rho_depth_drift_DIII-D.npy')
dne_arr_DIIID = np.load(save_path+'dne_arr_drift_DIII-D.npy')

rho = np.linspace(0.0,1.015,205)

#%%

X,Y = np.meshgrid(velocity,degrees)
# print(X,Y)
plt.figure()
cont = plt.pcolormesh(X,Y,(rho_max_na[0,:,:]).T,cmap='plasma',shading='nearest',vmin=np.nanmin(rho_max_na[0,:,:]),vmax=np.nanmax(rho_max_na[0,:,:]))
cb = plt.colorbar()
cb.set_label(r"$\delta n_e$ (peak)")
plt.ylim([-30,30])
plt.xlim([250,1500])
plt.xticks(velocity)
plt.ylabel(r"Degrees ($\degree$)")
plt.xlabel(r"Velocity (ms$^{-1}$)")
plt.title(r"Diameter d = 5.0 (mm)")
# plt.clabel(cont,cont.levels)
# plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/figs/d5_peak_drift_na.png")
plt.show()

#%%
# rho_max = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_max_drift_baylor.npy')

plt.figure()
cont = plt.pcolormesh(X,Y,rho_depth_na[0,:,:].T,cmap='plasma',vmin=np.nanmin(rho_depth_na[0,:,:]),vmax=np.nanmax(rho_depth_na[0,:,:]))
cb = plt.colorbar()
cb.set_label(r"$\delta n_e$ (depth)")
plt.ylim([-30,30])
plt.xlim([250,1500])
plt.ylabel(r"Degrees ($\degree$)")
plt.xlabel(r"Velocity (ms$^{-1}$)")
plt.title(r"Diameter d = 5.0 (mm)")
# plt.clabel(cont,cont.levels)
# plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/figs/d5_depth_drift_na.png")
plt.show()

#%%

X,Y = np.meshgrid(rho[120:204],degrees)
# print(X,Y)
plt.figure()
cont = plt.pcolormesh(X,Y,dne_arr_na[120:204,0,:].T,shading='nearest',cmap='plasma')
cb = plt.colorbar()
cb.set_label(r"$\delta n_e$")
# plt.ylim([-25,25])
# plt.xlim([250,1500])
plt.ylabel(r"Degrees ($\degree$)")
plt.xlabel(r"$\rho$")
plt.title(r"5 mm Diameter, Velocity 1500 $ms^{-1}$")
# plt.clabel(cont,cont.levels)
# plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/figs/d5_dne_drift_na.png")
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

degrees = [-20.0,-10.0,0.0,10.0,20.0]
velocity = [250.,375.,500.,625.,750.,875.,1000.,1125.,1250.,1375.,1500.]

depth_diff_na = 1.015 - rho_depth_na[0,:,:]
max_diff_na = 1.015 - rho_max_na[0,:,:]
depth_diff_HPI2 = 1.015 - rho_depth_HPI2[0,:,:]
max_diff_HPI2 = 1.015 - rho_max_HPI2[0,:,:]
depth_diff_DIIID = 1.015 - rho_depth_DIIID[0,:,:]
max_diff_DIIID = 1.015 - rho_max_DIIID[0,:,:]
set = [0,1,2,3,4]
depth_set_na = [depth_diff_na[:,ii] for ii in set]
max_set_na = [max_diff_na[:,ii] for ii in set]
depth_set_HPI2 = [depth_diff_HPI2[:,ii] for ii in set]
max_set_HPI2 = [max_diff_HPI2[:,ii] for ii in set]
set = [0,1,2,3,4]
depth_set_DIIID = [depth_diff_DIIID[:,ii] for ii in set]
max_set_DIIID = [max_diff_DIIID[:,ii] for ii in set]

if dataframe:
    depth_full_na = 1.015 - rho_depth_na
    max_full_na = 1.015 - rho_max_na
    depth_full_HPI2 = 1.015 - rho_depth_HPI2
    max_full_HPI2 = 1.015 - rho_max_HPI2
    depth_full_DIIID = 1.015 - rho_depth_DIIID
    max_full_DIIID = 1.015 - rho_max_DIIID
    depth_stacked = np.stack([depth_full_na,depth_full_DIIID,depth_full_HPI2],axis=-1)
    depth_full = depth_stacked.reshape(-1)
    max_stacked = np.stack([max_full_na,max_full_DIIID,max_full_HPI2],axis=-1)
    max_full = max_stacked.reshape(-1)
    data['diff_depth (-)'] = depth_full
    data['diff_peak (-)'] = max_full
    data.to_csv(save_path + 'ST50010DN_straight_dataset.csv',index=False)
    data.to_hdf(save_path + 'ST50010DN_straight_dataset.h5', key="data", mode="w")

datasets = [
    (max_set_na,depth_set_na),
    (max_set_HPI2,depth_set_HPI2),
    (max_set_DIIID,depth_set_DIIID)
]

datasets = [
    (max_stacked.reshape(3,11,5,3)[0,:,:,0].T,depth_stacked.reshape(3,11,5,3)[0,:,:,0].T),
    (max_stacked.reshape(3,11,5,3)[0,:,:,2].T,depth_stacked.reshape(3,11,5,3)[0,:,:,2].T),
    (max_stacked.reshape(3,11,5,3)[0,:,:,1].T,depth_stacked.reshape(3,11,5,3)[0,:,:,1].T)
]

depth_data = [
    (depth_set_na),
    (depth_set_HPI2),
    (depth_set_DIIID)
]

max_data = [
    (max_set_na),
    (max_set_HPI2),
    (max_set_DIIID)
]

X, Y = np.meshgrid(velocity, degrees)
# Global normalization for consistent color scale
vmin_max = min(np.min(d) for d in max_data)
vmax_max = max(np.max(d) for d in max_data)
norm_max = colors.Normalize(vmin=vmin_max, vmax=vmax_max)
vmin_depth = min(np.min(d) for d in depth_data)
vmax_depth = max(np.max(d) for d in depth_data)
norm_depth = colors.Normalize(vmin=vmin_depth, vmax=vmax_depth)

# Define threshold for “largest regions” (e.g., top 20%)
threshold_fraction = 0.8

# Create 2×3 grid of plots
fig, axs = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)

for j, (var1,var2) in enumerate((datasets)):
    # --- Row 1: (rho_max - rho_peak)
    pcm1 = axs[0, j].pcolormesh(X, Y, var1, cmap='plasma', norm=norm_max, shading='auto')

    # --- Row 2: (rho_max - rho_depth)
    pcm2 = axs[1, j].pcolormesh(X, Y, var2, cmap='plasma', norm=norm_depth, shading='auto')

cbar1 = fig.colorbar(pcm1, ax=axs[0, :], orientation='vertical', fraction=0.025, pad=0.02)
cbar1.set_label(r'$\rho_{\text{max}} - \rho_{\text{peak}}$')
# cbar1.set_ticks([0.075,0.125,0.175,0.225])


cbar2 = fig.colorbar(pcm2, ax=axs[1, :], orientation='vertical', fraction=0.025, pad=0.02)
cbar2.set_label(r'$\rho_{\text{max}} - \rho_{\text{depth}}$')
# cbar2.set_ticks([0.07,0.14,0.21,0.28])

for ii in range(0,2):
    axs[ii,1].set_yticks([])
    axs[ii,2].set_yticks([])
    axs[ii,0].set_ylabel(r'Degrees ($\degree$)')
for jj in range(0,3):
    axs[0,jj].set_xticks([])
    axs[1,jj].set_xlabel(r'Velocity (ms$^{-1}$)')

axs[0,0].set_title(r'No Drift')
axs[0,1].set_title(r'$\Delta_{\text{HPI2}}$')
axs[0,2].set_title(r'$\Delta_{\text{DIII-D}}$')
axs[0,0].set_yticks([-20.0,-10.0,0.0,10.0,20.0])
axs[1,0].set_yticks([-20.0,-10.0,0.0,10.0,20.0])

# plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/figs/diff_comp.png",dpi=300)
# plt.savefig("/Users/gz6/Documents/TokE/papers/PFS_figures/RUN7_15100088_highres_diff_comp.png",dpi=300)
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

degrees = [-20.0,-10.0,0.0,10.0,20.0]
velocity = [250.,375.,500.,625.,750.,875.,1000.,1125.,1250.,1375.,1500.]

depth_diff_na = 1.015 - rho_depth_na[0,:]
max_diff_na = 1.015 - rho_max_na[0,:]
depth_diff_HPI2 = 1.015 - rho_depth_HPI2[0,:]
max_diff_HPI2 = 1.015 - rho_max_HPI2[0,:]
depth_diff_DIIID = 1.015 - rho_depth_DIIID[0,:]
max_diff_DIIID = 1.015 - rho_max_DIIID[0,:]
set = [0]
depth_set_na = [depth_diff_na[ii] for ii in set]
max_set_na = [max_diff_na[ii] for ii in set]
depth_set_HPI2 = [depth_diff_HPI2[ii] for ii in set]
max_set_HPI2 = [max_diff_HPI2[ii] for ii in set]
set = [0]
depth_set_DIIID = [depth_diff_DIIID[ii] for ii in set]
max_set_DIIID = [max_diff_DIIID[ii] for ii in set]

if dataframe:
    depth_full_na = 1.015 - rho_depth_na
    max_full_na = 1.015 - rho_max_na
    depth_full_HPI2 = 1.015 - rho_depth_HPI2
    max_full_HPI2 = 1.015 - rho_max_HPI2
    depth_full_DIIID = 1.015 - rho_depth_DIIID
    max_full_DIIID = 1.015 - rho_max_DIIID
    depth_stacked = np.stack([depth_full_na,depth_full_DIIID,depth_full_HPI2],axis=-1)
    depth_full = depth_stacked.reshape(-1)
    max_stacked = np.stack([max_full_na,max_full_DIIID,max_full_HPI2],axis=-1)
    max_full = max_stacked.reshape(-1)
    data['diff_depth (-)'] = depth_full
    data['diff_peak (-)'] = max_full
    data.to_csv(save_path + 'ST50010DN_straight_dataset.csv',index=False)
    data.to_hdf(save_path + 'ST50010DN_straight_dataset.h5', key="data", mode="w")

datasets = [
    (max_set_na,depth_set_na),
    (max_set_HPI2,depth_set_HPI2),
    (max_set_DIIID,depth_set_DIIID)
]

datasets = [
    (max_stacked.reshape(3,11,3)[0,:,0].T,depth_stacked.reshape(3,11,3)[0,:,0].T),
    (max_stacked.reshape(3,11,3)[0,:,2].T,depth_stacked.reshape(3,11,3)[0,:,2].T),
    (max_stacked.reshape(3,11,3)[0,:,1].T,depth_stacked.reshape(3,11,3)[0,:,1].T)
]

depth_data = [
    (depth_set_na),
    (depth_set_HPI2),
    (depth_set_DIIID)
]

max_data = [
    (max_set_na),
    (max_set_HPI2),
    (max_set_DIIID)
]

X, Y = np.meshgrid(velocity, degrees)
# Global normalization for consistent color scale
vmin_max = min(np.min(d) for d in max_data)
vmax_max = max(np.max(d) for d in max_data)
norm_max = colors.Normalize(vmin=vmin_max, vmax=vmax_max)
vmin_depth = min(np.min(d) for d in depth_data)
vmax_depth = max(np.max(d) for d in depth_data)
norm_depth = colors.Normalize(vmin=vmin_depth, vmax=vmax_depth)

# Define threshold for “largest regions” (e.g., top 20%)
threshold_fraction = 0.8

# Create 2×3 grid of plots
fig, axs = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)

for j, (var1,var2) in enumerate((datasets)):
    # --- Row 1: (rho_max - rho_peak)
    pcm1 = axs[0, j].pcolormesh(X, Y, var1, cmap='plasma', norm=norm_max, shading='auto')

    # --- Row 2: (rho_max - rho_depth)
    pcm2 = axs[1, j].pcolormesh(X, Y, var2, cmap='plasma', norm=norm_depth, shading='auto')

cbar1 = fig.colorbar(pcm1, ax=axs[0, :], orientation='vertical', fraction=0.025, pad=0.02)
cbar1.set_label(r'$\rho_{\text{max}} - \rho_{\text{peak}}$')
# cbar1.set_ticks([0.075,0.125,0.175,0.225])


cbar2 = fig.colorbar(pcm2, ax=axs[1, :], orientation='vertical', fraction=0.025, pad=0.02)
cbar2.set_label(r'$\rho_{\text{max}} - \rho_{\text{depth}}$')
# cbar2.set_ticks([0.07,0.14,0.21,0.28])

for ii in range(0,2):
    axs[ii,1].set_yticks([])
    axs[ii,2].set_yticks([])
    axs[ii,0].set_ylabel(r'Degrees ($\degree$)')
for jj in range(0,3):
    axs[0,jj].set_xticks([])
    axs[1,jj].set_xlabel(r'Velocity (ms$^{-1}$)')

axs[0,0].set_title(r'No Drift')
axs[0,1].set_title(r'$\Delta_{\text{HPI2}}$')
axs[0,2].set_title(r'$\Delta_{\text{DIII-D}}$')
axs[0,0].set_yticks([-20.0,-10.0,0.0,10.0,20.0])
axs[1,0].set_yticks([-20.0,-10.0,0.0,10.0,20.0])

# plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/figs/diff_comp.png",dpi=300)
# plt.savefig("/Users/gz6/Documents/TokE/papers/PFS_figures/RUN7_15100088_highres_diff_comp.png",dpi=300)
plt.show()


# %%

DF_path = "/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/TokE_pellet_data/RUN1_15100081_data/data/"
DF_name = "RUN1_15100081_dataset.csv"

data = pd.read_csv(DF_path+DF_name)




# %%
