import numpy as np
import matplotlib.pyplot as plt
import PelletSumfile as PS
from scipy import integrate

base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
gfile = base_dir + "st50010dn.eqdsk"

v_plot = ['v1000','v1500']
d_plot = ['d2.0','d3.0','d4.0']
lines = ['--g','--b','--r']

count = 0
fig, ax = plt.subplots(1,3,sharey=True)
for ii in d_plot:
    NP = 0
    for jj in v_plot:
        sumfile_loc = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/' + d_plot[count] + '_' + v_plot[NP] + '/sum_pellet_pellet_test_.dat'
        pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
        if NP==0:
            ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)']*1e-19,'k-',label="Initial")    
        ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)']*1e-19,lines[NP],label=jj + r' (ms$^{-1}$)')
        NP += 1
    ax[count].grid(True)
    if count==0:
        ax[count].set_ylabel(r'$n_e$ ($1\times10^{19}$ m$^{-3}$)')
    ax[count].set_xlabel('rho')
    ax[count].set_title('d = ' + d_plot[count] + ' (mm)')
    ax[count].legend()
    count += 1
plt.show()

# ne_init = pelsum.profiles['ne(tpel-)']
# ne = pelsum.profiles['ne(tpel+)']
# dne = pelsum.profiles['delta_ne']
# x = pelsum.profiles['rho_t']
# check_tot = np.trapezoid(ne)
# check_del = np.trapezoid(dne)
# check_init = np.trapezoid(ne_init)
# print(check_init)
# print(check_del)
# print(check_tot)
# print(check_init+check_del)
# print((4.0/3.0)*np.pi*(1.5e-3)**3*3.03697419E+28)