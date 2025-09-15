import numpy as np
import matplotlib.pyplot as plt
import PelletSumfile as PS
from scipy import integrate

# plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 12})
# plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})

base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
gfile = base_dir + "st50010dn.eqdsk"

dv_plot = 0
dd_plot = 0
single_plot = 1

v_plot = ['v200','v300']
# v_plot = ['v1000','v1250','v1500']
d_plot = ['d2.0','d3.0','d4.0']
lines = ['--g','-.b','-+r']

den_plot = ['den0_4e+18','den0_6e+18','den0_8e+18','den0_1e+19','den0_1.2e+19','den0_1.4e+19']

if single_plot:
    sumfile_loc = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/build/sum_pellet_DIII-D_drift            .dat'
    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # pelsum.plot_density_change()
    plt.figure()
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)'])
    # plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)'])
    # plt.ylim([0.3e20,1.0e20])
    plt.show()
    plt.figure()
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['delta_ne'])
    # plt.ylim([-1.0e18,5.0e19])
    plt.show()
    print(integrate.trapezoid(pelsum.profiles['ne(tpel+)'],pelsum.profiles['rho_t']))

if dv_plot:
    count = 0
    fig, ax = plt.subplots(3,1,sharex=True)
    for ii in d_plot:
        NP = 0
        for jj in v_plot:
            sumfile_loc = '/Users/gz6/Documents/pellets/tokamak_energy/results/k0_NGS_ref/injection_angles/curved_centre/refined/' + d_plot[count] + '_' + v_plot[NP] + '/sum_pellet_pellet_test_.dat'
            pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
            if NP==0:
                ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)']*1e-19,'k-',label="Initial")    
            ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)']*1e-19,lines[NP],label=jj + r' (ms$^{-1}$)')
            NP += 1
        ax[count].grid(True)
        if count==2:
            ax[count].set_xlabel('rho')
            
        ax[count].set_ylabel(r'$n_e$ ($1\times10^{19}$ m$^{-3}$)')
        ax[count].set_title('d = ' + d_plot[count] + ' (mm)')
        ax[count].legend()
        count += 1
    plt.suptitle("Curved Injection (centre) \n Entry: R=2.9 (m), Z=1.9 (m) \n Exit: R=7.0 (m), Z=-1.25 (m)")
    plt.tight_layout()
    plt.subplots_adjust(left=0.075, bottom=0.075, right=0.975, top=0.85, wspace=0.2, hspace=0.25)
    # fig_width, fig_height = plt.gcf().get_size_inches()
    # print(fig_width, fig_height)
    fig.set_size_inches(10.0,12.0)
    plt.show()

    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # print(pelsum.pelpath['R'], pelsum.pelpath['Z'])
    pelsum.plot_trajectory(contours=np.linspace(0.0,1.0,15))
    plt.gcf()
    fig.set_size_inches(5.0,12.0,forward=True)
    plt.tight_layout()

if dd_plot:
    rho_max = np.zeros(len(den_plot))
    fig, ax = plt.subplots(2,3)
    count = 0
    for ii in 0,1:
        for jj in 0,1,2:
            sumfile_loc = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/' + den_plot[count] + '/sum_pellet_pellet_test_.dat'
            pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
            idx = np.where(pelsum.profiles['ne(tpel+)']==np.max(pelsum.profiles['ne(tpel+)']))[0][0]
            print(sumfile_loc)
            print(idx)
            print(np.shape(pelsum.profiles['ne(tpel+)']))
            rho_max[count] = pelsum.profiles['rho_t'][idx]
            ax[ii,jj].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)']*1e-19,'k-',label="Initial")
            ax[ii,jj].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)']*1e-19,lines[jj],label=den_plot[count] + r' (ms$^{-1}$)')
            ax[ii,jj].axvline(x=rho_max[count],color='k',linestyle='-',alpha=0.5,label='Max density')  
            count += 1
    plt.show()
    print(rho_max)

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