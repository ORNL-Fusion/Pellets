#%%
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
import PelletSumfile as PS
import hermite_ped as hp

#%%
pellet_func = 0
eped_func = 0
hermite_func = 1
save = 0
plot = 1

#%%
ncplas = 205
te0 = 25.0
te1 = 0.075
tetop = 3.5
den0 = 11.0
dentop = 7.0
pedfrac = 0.7
pedne = den0/2.0
# pedne = dentop
# den1 = 0.25*pedne
den1 = 4.0
pedte = tetop
pedwid = 0.1
print(den0/1.5)
# pedte = 4.5
# print(dentop/den0)
# print(tetop/te0)
# 0.6363636363636364
# 0.12

xped = 1.0 - pedwid
xmid = 1.0 - 0.5*pedwid
rho_r = np.linspace(0,1,ncplas)
# print(rho_r)

#%%

if pellet_func:
    pedte1 = np.zeros(ncplas)
    pedte2 = np.zeros(ncplas)
    pedne1 = np.zeros(ncplas)
    pedne2 = np.zeros(ncplas)

    for ii in range(0,ncplas):
        px_den = 4.5
        qx_den = 0.075
        px_te = 1.5
        qx_te = 1.0     
        thill = te0 - pedte
        nhill = den0 - pedne
        pedtlambda = (te0-te1)/(1.0+np.tanh(1.0))
        pednlambda = (den0-den1)/(1.0+np.tanh(1.0))
        pedte1[ii] = te1 + pedtlambda*(np.tanh(1.0)-np.tanh(2*(rho_r[ii]-xmid)/pedwid))
        pedne1[ii] = den1 + pednlambda*(np.tanh(1.0)-np.tanh(2*(rho_r[ii]-xmid)/pedwid))
        if rho_r[ii] < (1.0 - pedwid):
            pedte2[ii] = thill*((1-(rho_r[ii]/xped)**px_te)**qx_te)
            pedne2[ii] = nhill*((1-(rho_r[ii]/xped)**px_den)**qx_den)
        else:
            pedte2[ii] = 0.0
            pedne2[ii] = 0.0

    te_r = pedte1+pedte2
    den_r = pedne1+pedne2

if eped_func:
    aT = 1.2
    betaT = 1.4
    ane = 1.1
    betane = 1.1
    H1 = 1.0 - (rho_r/xped)
    H2T = (1.0 - (rho_r/xped)**aT)**betaT
    H2ne = (1.0 - (rho_r/xped)**ane)**betane
    H2T = np.nan_to_num(H2T)
    H2ne = np.nan_to_num(H2ne)
    te_prof = te1 + pedtlambda*(np.tanh(2.0*(1.0-xmid)/pedwid) - np.tanh(2.0*(rho_r - xmid)/pedwid)) + te1*np.heaviside(H1,H2T)
    ne_prof = den1 + pednlambda*(np.tanh(2.0*(1.0-xmid)/pedwid) - np.tanh(2.0*(rho_r - xmid)/pedwid)) + 0.1*np.heaviside(H1,H2ne)

if hermite_func:
    ne_prof = hp.hermite_ped(rho_r,xped,valaxis=den0,valsep=den1,fact=1.0,valpedtop=dentop)
    te_prof = hp.hermite_ped(rho_r,xped,valaxis=te0,valsep=te1,fact=0.5,valpedtop=tetop)

#%%

def load_blocks(filename):
    blocks = {}
    current_label = None
    current_data = []

    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue

            # Detect header lines: anything with letters
            if any(c.isalpha() for c in stripped):
                # Save previous block
                if current_label is not None and current_data:
                    blocks[current_label] = np.array(current_data, float)
                current_label = stripped
                current_data = []
            else:
                # numeric row
                row = [float(x) for x in stripped.split()]
                current_data.append(row)

    # Save last block
    if current_label is not None and current_data:
        blocks[current_label] = np.array(current_data, float)

    return blocks


#%%

root_dir = "/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/"
filename = "ST40_15100088_ASTRA_RUN7_24610p00ms.pfile"

blocks = load_blocks(root_dir + filename)

labels = ['65 psinorm ne(10^20/m^3) dne/dpsiN', '65 psinorm te(keV) dte/dpsiN']

rho = blocks[labels[0]][:,0]
ne = blocks[labels[0]][:,1]*1.0e1
te = blocks[labels[1]][:,1]
dnedpsiN = blocks[labels[0]][:,2]
dtedpsiN = blocks[labels[1]][:,2]

npts = len(rho)
npts = 200


from scipy import interpolate

rho_uniform = np.linspace(min(rho),max(rho),npts)
rho_uniform[[0,npts-1]] = rho[[0,-1]]
f_ne = interpolate.interp1d(rho,ne)
f_te = interpolate.interp1d(rho,te)
ne_uniform = f_ne(rho_uniform)
te_uniform = f_te(rho_uniform)


#%%

header_str = ['ST40 profiles: rho ne te']

outdir = "/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/refined/"
outname = "RUN7_p15100088_ref.dat"

with open(outdir+outname, "w") as f:
        # line 1: header string
        f.write(header_str[0] + "\n")

        # line 2: number of points
        f.write(f"{npts}\n")

        # data lines
        for a, b, c in zip(rho_uniform, ne_uniform, te_uniform):
            f.write(f"{a: .6f}  {b: .6f}  {c: .6f}\n")



#%%

# plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'font.serif': 'Times New Roman'})
plt.rcParams.update({'font.family': 'serif'})
# plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.fontset': 'stix'})
plt.rcParams.update({'mathtext.default': 'it'})
# plt.rcParams.update({'mathtext.default': 'cm'})

#%%

plot = 1
if plot:
    fig, ax = plt.subplots(2,1,figsize=(8.5,11),sharex=True,layout="constrained")
    ax[0].plot(rho_uniform,ne_uniform,'k.-')
    ax[1].plot(rho_uniform,te_uniform,'b.-')
    ax[0].set_ylabel(r"Density ($\times 10^{19}$ m$^{-3}$)")
    ax[1].set_ylabel(r"Temperature (keV)")
    ax[1].set_xlabel(r"$\rho$")
    ax[0].set_xlim([0.0,1.0])
    ax[0].grid(True)
    ax[1].grid(True)
    # plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/RUN7_15100088_highres_profiles.png",dpi=150)
    plt.show()

#%%

base_dir = "/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/"

RUN1_15100081_file = "RUN1_15100081/refined/RUN1_p15100081_ref.dat"
RUN3_15100089_file = "RUN3_15100089/RUN3_p15100089in.dat"
RUN7_15100081_file = "RUN7_15100081/RUN7_p15100081in.dat"
RUN7_15100088_file = "RUN7_15100088/RUN7_p15100088in.dat"
original_file = "/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/prof_ped_7.0_4.0.dat"

F1 = np.loadtxt(base_dir+RUN1_15100081_file,skiprows=2)[:,:3]
F2 = np.loadtxt(base_dir+RUN3_15100089_file,skiprows=2)[:,1:]
F3 = np.loadtxt(base_dir+RUN7_15100081_file,skiprows=2)[:,1:]
F4 = np.loadtxt(base_dir+RUN7_15100088_file,skiprows=2)[:,1:]
FO = np.loadtxt(original_file,skiprows=2)

data = np.column_stack((F1,F2,F3,F4))

fig, ax = plt.subplots(2,1,figsize=(8.5,11),layout="constrained")

ax[0].plot(FO[:,0],FO[:,1],label="Original")
ax[0].plot(data[:,0],data[:,1],'k',label="RUN1_81")
ax[0].plot(data[:,0],data[:,3],'r',label="RUN3_89")
ax[0].plot(data[:,0],data[:,5],'b',label="RUN7_81")
ax[0].plot(data[:,0],data[:,7],'m',label="RUN7_88")
ax[0].plot(data[-7,0],data[-7,1],'k*')
ax[0].plot(data[-7,0],data[-7,3],'r*')
ax[0].plot(data[-7,0],data[-7,5],'b*')
ax[0].plot(data[-7,0],data[-7,7],'m*')
ax[1].plot(FO[:,0],FO[:,2],label="Original")
ax[1].plot(data[:,0],data[:,2],'k',label="RUN1_81")
ax[1].plot(data[:,0],data[:,4],'r',label="RUN3_89")
ax[1].plot(data[:,0],data[:,6],'b',label="RUN7_81")
ax[1].plot(data[:,0],data[:,8],'m',label="RUN7_88")
ax[1].plot(data[-7,0],data[-7,2],'k*')
ax[1].plot(data[-7,0],data[-7,4],'r*')
ax[1].plot(data[-7,0],data[-7,6],'b*')
ax[1].plot(data[-7,0],data[-7,8],'m*')

ax[0].set_xticks([])
ax[0].set_ylabel(r"$n_e (\times10^{19}$m$^{-3}$)")
ax[1].set_xlabel(r"$\rho$")
ax[1].set_ylabel(r"$T_e$ (keV)")

ax[0].legend()
ax[1].legend()
plt.show()

#%%

plot = 1
if plot:
    fig, ax = plt.subplots(2,1,figsize=(8.5,11),sharex=True,layout="constrained")
    ax[0].plot(FO[:,0],FO[:,1],'k.-')
    ax[1].plot(FO[:,0],FO[:,2],'b.-')
    ax[0].set_ylabel(r"Density ($\times 10^{19}$ m$^{-3}$)")
    ax[1].set_ylabel(r"Temperature (keV)")
    ax[1].set_xlabel(r"$\rho$")
    ax[0].set_xlim([0.0,1.0])
    ax[0].grid(True)
    ax[1].grid(True)
    # plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/ST50010DN/ST50010DN_profiles.png",dpi=150)
    plt.show()


# %%
