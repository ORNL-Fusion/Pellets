import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import PelletSumfile as PS
import hermite_ped as hp

pellet_func = 0
eped_func = 0
hermite_func = 1
save = 0
plot = 1

ncplas = 105
te0 = 25.0
te1 = 0.075
px_te = 1.5
qx_te = 1.0
tetop = 3.0
den0 = 11.0
pedfrac = 1.5
denmult = 0.85
dentop = 7.0
pedne = dentop/pedfrac
den1 = pedne*denmult
px_den = 4.5
qx_den = 0.075
pedwid = 0.1
pedte = 4.5

xped = 1.0 - pedwid
xmid = 1.0 - 0.5*pedwid
rho_r = np.linspace(0,1,ncplas)
print(rho_r)

if pellet_func:
    pedte1 = np.zeros(ncplas)
    pedte2 = np.zeros(ncplas)
    pedne1 = np.zeros(ncplas)
    pedne2 = np.zeros(ncplas)

    for ii in range(0,ncplas):
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
    ne_prof = hp.hermite_ped(rho_r,xped,valaxis=den0,valsep=den1,valpedtop=dentop)
    te_prof = hp.hermite_ped(rho_r,xped,valaxis=te0,valsep=te1,valpedtop=tetop)

if plot:
    fig, ax = plt.subplots(2,1,sharex=True)
    ax[0].plot(rho_r,ne_prof,'k.-')
    ax[1].plot(rho_r,te_prof,'b.-')
    ax[0].set_ylabel(r"Density ($\times 10^{19}$ m$^{-3}$)")
    ax[1].set_ylabel(r"Temperature (keV)")
    ax[1].set_xlabel(r"$\rho_r$")
    plt.show()

if save:
    df = pd.DataFrame({"xr": rho_r, "denxr": ne_prof, "texr": te_prof})
    pfile_name = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/' + 'p_'+'ped_'+str(den0)+'_'+str(np.round(den1))+'.dat'
    pfile = open(pfile_name, 'a')
    pfile.write("User set pedestal.\n")
    pfile.write(str(ncplas)+'\n')
    pfile.write(df.to_string(header=False,index=False))
    pfile.close()