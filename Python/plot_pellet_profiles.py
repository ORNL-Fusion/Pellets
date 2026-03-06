#%%
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import PelletSumfile as PS
from scipy import integrate
import EFIT.equilParams_class as epc

# plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.serif': 'Times New Roman'})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'mathtext.fontset': 'stix'})
# plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'it'})

#%%

base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
gfile = base_dir + "RUN7_15100088.geqdsk"
# gfile = base_dir + "g200201.00000"

dv_plot = 0
dd_plot = 0
single_plot = 1
traj_plot = 1

# v_plot = ['v200','v300']
# v_plot = ['v1000','v1250','v1500']
v_plot = ['v250','v375','v500','v625','v750','v875','v1000','v1125','v1250','v1375','v1500']
d_plot = ['d5.0','d6.0','d7.0']
# degrees_plot = ['m30','m20','m10','p00','p10','p20','p30']
degrees_plot = ['p00','m10','m20','p10','p20']
labels = [r"$\pm 20\degree$",r"$\pm 10\degree$",r"$0\degree$"][::-1]
# lines = ['--g','-.b','-.r']
# lines_2 = ['-xg','-+b','-+r']

den_plot = ['den0_4e+18','den0_6e+18','den0_8e+18','den0_1e+19','den0_1.2e+19','den0_1.4e+19']

#%%

if traj_plot:
    count = [0,1,2,1,2]
    n = 4
    colors = plt.cm.plasma(np.linspace(0,1,n))
    fig, ax = plt.subplots(figsize=(5.0,8.5),layout="constrained")
    # sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/straight/drift_HPI2/d7.0_v1500/sum_pellet_pellet_test_.dat'
    # pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # Rtraj = pelsum.pelpath['R']
    # Ztraj = pelsum.pelpath['Z']
    # idr = np.where(pelsum.pelpath['rho']==min(pelsum.pelpath['rho']))[0][0]

    # for tr in range(0,idr):
    #     plt.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]],'k-',lw = 2.5)
    #     if tr == idr-1:
    #         plt.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]],'k-',lw = 2.5,label='Off-Axis')

    for kk in range(len(degrees_plot)):
    # for kk in range(0):
        sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/refined/' + degrees_plot[kk] + '_degrees/drift_HPI2/d7.0_v1500/sum_pellet_pellet_test_.dat'
        pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
        Rtraj = pelsum.pelpath['R']
        Ztraj = pelsum.pelpath['Z']
        gf = epc.equilParams(gfile)
        psiRZn = gf.g['psiRZn']
        gR = gf.g['R']
        gZ = gf.g['Z']
        wall = gf.g['wall']
        idr = np.where(pelsum.pelpath['rho']==min(pelsum.pelpath['rho']))[0][0]
        # print(len(Rtraj)-1)
        if kk==0:
            contours = np.linspace(0.0,1.0,11)
            cont = ax.contour(gR, gZ, np.sqrt(psiRZn), contours, cmap='plasma', linewidths = [1.0], alpha=0.5)
            ax.clabel(cont,cont.levels)
            ax.contour(gR, gZ, np.sqrt(psiRZn), [1], colors = 'k', linewidths = [2.0])
        for tr in range(len(Rtraj)-1):
            # print(tr)
        # for tr in range(idr):
            if Rtraj[tr] > 6.7:
                continue
            elif (kk==0 or kk==1 or kk==2) and tr==len(Rtraj)-2:
                print(kk)
                ax.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]], '-', color=colors[count[kk]],lw = 1.8, markersize=1,label=labels[kk])
            else:
                ax.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]], '-', color=colors[count[kk]],lw = 1.8, markersize=1)
        if kk==0:
            ax.plot([Rtraj[idr], Rtraj[idr]], [Ztraj[idr], Ztraj[idr]], 'kx', lw = 1, markersize=5)

    plt.legend()
    # plt.title('Injection Angles')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.xlim([2.,8.])
    plt.ylim([-7.5,7.5])
    # plt.tight_layout()
    # plt.savefig("/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/straight/contour_angles.png")
    # plt.savefig("/Users/gz6/Documents/TokE/papers/PFS_figures/RUN7_15100088_contour_all_trajectories.png",dpi=300)
    plt.show()

#%%

fig, ax = plt.subplots(figsize=(5.0,8.5),layout="constrained")
sumfile_loc = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/build/sum_pellet_pellet_test_            .dat'
pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
Rtraj = pelsum.pelpath['R']
Ztraj = pelsum.pelpath['Z']
gf = epc.equilParams(gfile)
psiRZn = gf.g['psiRZn']
gR = gf.g['R']
gZ = gf.g['Z']
wall = gf.g['wall']
contours = np.linspace(0.0,1.0,11)
cont = ax.contour(gR, gZ, np.sqrt(psiRZn), contours, cmap='plasma', linewidths = [1.0], alpha=0.5)
ax.clabel(cont,cont.levels)
ax.contour(gR, gZ, np.sqrt(psiRZn), [1], colors = 'k', linewidths = [2.0])
for tr in range(len(Rtraj)-1):
    ax.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]], '-',lw = 1.8, markersize=1)

plt.show()

#%%

if single_plot:
    sumfile_loc = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/build/sum_pellet_pellet_test_            .dat'
    # sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/p00_degrees/drift_na/d5.0_v1250/sum_pellet_pellet_test_.dat'
    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # pelsum.plot_density_change()
    plt.figure()
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)'])
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)'])
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['ne_d(tpel+)'])
    # plt.show()
    # plt.figure()
    # plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['Te(tpel-)'])
    # plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['Te(tpel+)'])
    # plt.ylim([0.3e20,1.0e20])
    # plt.show()
    plt.figure()
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['delta_ne']*1.0e-20,'k-')
    plt.plot(pelsum.profiles['rho_t'],pelsum.profiles['delta_ne_drift']*1.0e-20,'r--')
    # ax1.set_xlim([min(pelsum.profiles['rho_t']), max(pelsum.profiles['rho_t'])])
    # ax2.set_xlim([min(pelsum.profiles['r_grid']), max(pelsum.profiles['r_grid'])])
    # plt.ylim([-1.0e18,5.0e19])
    plt.show()
    # print(integrate.trapezoid(pelsum.profiles['delta_ne'],pelsum.profiles['rho_t']))
    # print(integrate.trapezoid(pelsum.profiles['delta_ne_drift'],pelsum.profiles['rho_t']))

    # gf = epc.equilParams(gfile)
    # psiRZn = gf.g['psiRZn']
    # gR = gf.g['R']
    # gZ = gf.g['Z']
    # R = pelsum.pelpath['R']
    # Z = pelsum.pelpath['Z']
    # wall = gf.g['wall']
    # idr = np.where(pelsum.pelpath['rho']==min(pelsum.pelpath['rho']))[0][0]
    # # print(len(Rtraj)-1)
    # plt.figure()
    # contours = np.linspace(0.0,1.0,11)
    # cont = plt.contour(gR, gZ, psiRZn, contours, cmap='plasma', linewidths = [1.0], alpha=0.5)
    # plt.clabel(cont,cont.levels)
    # plt.contour(gR, gZ, psiRZn, [1], colors = 'k', linewidths = [2.0])
    # check = gf.getTorPsi()
    # # plt.plot(R[0:200], np.sqrt(check['psitorN1D']), 'k.-')
    # plt.plot(R[0:200], np.sqrt(check['psitor1D']/check['psitor1D']), 'k.-')
    # plt.show()
    
    # fig = plt.figure()
    # pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # idr = np.where(pelsum.pelpath['rho']==min(pelsum.pelpath['rho']))[0][0]
    # # print(idr)
    # # print(pelsum.pelpath['rho'][idr])
    # # print(pelsum.pelpath['R'][idr], pelsum.pelpath['Z'][idr])
    # pelsum.plot_trajectory(idr,r"30$\degree$",contours=np.linspace(0.0,1.0,11))
    # plt.gcf()
    # # plt.plot(pelsum.profiles['R'][idr],pelsum.profiles['Z'][idr])
    # # print(pelsum.profiles['R'][idr],pelsum.profiles['Z'][idr],pelsum.pelpath['R'][idr],pelsum.pelpath['Z'][idr])
    # fig.set_size_inches(5.0,12.0,forward=True)
    # plt.tight_layout()
    # plt.show()

#%%

dv_plot = 1

lines=['--',':','-.']
colours = ['k','orange','magenta']

if dv_plot:
    count = 0
    fig, ax = plt.subplots(3,1,figsize=(8.5,11.0),sharex=True,layout="constrained")
    for ii in d_plot:
        NP = 0
        for jj in v_plot[0:11:4]:
            sumfile_loc = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/refined/p20_degrees/drift_na/' + d_plot[count] + '_' + v_plot[NP] + '/sum_pellet_pellet_test_.dat'
            pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
            sumfile_loc_2 = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/refined/p20_degrees/drift_DIII-D/' + d_plot[count] + '_' + v_plot[NP] + '/sum_pellet_pellet_test_.dat'
            pelsum_2 = PS.PelletSumfile(sumfile_loc=sumfile_loc_2,gfile=gfile)
            sumfile_loc_3 = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100088/refined/p20_degrees/drift_HPI2/' + d_plot[count] + '_' + v_plot[NP] + '/sum_pellet_pellet_test_.dat'
            pelsum_3 = PS.PelletSumfile(sumfile_loc=sumfile_loc_3,gfile=gfile)
            if NP==0:
                ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel-)']*1e-19,'k-',label="Initial")    
            ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)']*1e-19,color=colours[0],linestyle=lines[NP],label="NA") 
            # ax[count].plot(pelsum.profiles['rho_t'],pelsum.profiles['ne(tpel+)']*1e-19,label=jj + r' (ms$^{-1}$)')
            ax[count].plot(pelsum_2.profiles['rho_t'],pelsum_2.profiles['ne_d(tpel+)']*1e-19,color=colours[1],linestyle=lines[NP],label=jj + r' (ms$^{-1}$) (DIII-D)')
            ax[count].plot(pelsum_3.profiles['rho_t'],pelsum_3.profiles['ne_d(tpel+)']*1e-19,color=colours[2],linestyle=lines[NP],label=jj + r' (ms$^{-1}$) (HPI2)')
            NP += 1
        ax[count].grid(True)
        if count==2:
            ax[count].set_xlabel(r'$\rho$')
            
        ax[count].set_ylabel(r'$n_e$ ($1\times10^{19}$ m$^{-3}$)')
        ax[count].set_title('d = ' + d_plot[count] + ' (mm)')
        ax[count].legend(loc=6)
        count += 1
 
    # plt.suptitle("Curved Injection (centre) \n Entry: R=2.9 (m), Z=1.9 (m) \n Exit: R=7.0 (m), Z=-1.25 (m)")
    # plt.suptitle("Curved Injection (off-centre) \n Entry: R=2.85 (m), Z=1.33 (m) \n Exit: R=5.6 (m), Z=-3.75 (m)")
    # plt.suptitle("Straight Injection \n Entry: R=3.1 (m), Z=3.0 (m) \n Exit: R=5.0 (m), Z=-4.25 (m)")
    plt.suptitle(r"20$\degree$")
    # ax[2].set_xlabel(r"$\rho$")
    # plt.tight_layout()
    # plt.subplots_adjust(left=0.075, bottom=0.075, right=0.975, top=0.85, wspace=0.2, hspace=0.25)
    # fig_width, fig_height = plt.gcf().get_size_inches()
    # print(fig_width, fig_height)
    # fig.set_size_inches(10.0,12.0)
    # plt.savefig("/Users/gz6/Documents/TokE/papers/PFS_figures/RUN7_15100088_20deg_dv_highres.png",dpi=200)
    plt.show()
    print(integrate.trapezoid(pelsum.profiles['ne(tpel+)'],pelsum.profiles['rho_t']))
    print(integrate.trapezoid(pelsum_2.profiles['ne_d(tpel+)'],pelsum_2.profiles['rho_t']))

    pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile=gfile)
    # print(pelsum.pelpath['R'], pelsum.pelpath['Z'])
    pelsum.plot_trajectory(r"", contours=np.linspace(0.0,1.0,11))
    plt.gcf()
    # fig.set_size_inches(5.0,12.0,forward=True)
    # plt.tight_layout()

#%%

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

#%%

fig, ax = plt.subplots(3,1,layout="constrained")
fig.set_size_inches(6.0,7.5,forward=True)
sumfile_na = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/p00_degrees/drift_na/d5.0_v1000/sum_pellet_pellet_test_.dat'
pelsum_na = PS.PelletSumfile(sumfile_loc=sumfile_na,gfile=gfile)
v_plot = ["v1000","v1250","v1500"]
x = np.linspace(0,1.015,205)
x1, x2, y1, y2 = 0.6, 1.015, 3.5, 17.0
mask_x = (x <= x2) & (x >= x1)

rho_max_na = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_max_drift_na.npy')
rho_depth_na = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_depth_drift_na.npy')
rho_max_HPI2 = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_max_drift_HPI2.npy')
rho_depth_HPI2 = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_depth_drift_HPI2.npy')
rho_max_baylor = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_max_drift_baylor.npy')
rho_depth_baylor = np.load('/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/data/rho_depth_drift_baylor.npy')

degrees = [-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0]
velocity = [250.,375.,500.,625.,750.,875.,1000.,1125.,1250.,1375.,1500.]
count = [6,8,10]

n = 8
colors = plt.cm.viridis(np.linspace(0,1,n))

print(rho_max_HPI2[0,10,3])

for ii in range(0,len(v_plot)):
    sumfile_na = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/p00_degrees/drift_na/d5.0_'+v_plot[ii]+'/sum_pellet_pellet_test_.dat'
    sumfile_HPI2 = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/p00_degrees/drift_HPI2/d5.0_'+v_plot[ii]+'/sum_pellet_pellet_test_.dat'
    sumfile_baylor = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/injection_angles/updated/angle_scan/p00_degrees/drift_baylor/d5.0_'+v_plot[ii]+'/sum_pellet_pellet_test_.dat'
    pelsum_na = PS.PelletSumfile(sumfile_loc=sumfile_na,gfile=gfile)
    pelsum_HPI2 = PS.PelletSumfile(sumfile_loc=sumfile_HPI2,gfile=gfile)
    pelsum_baylor = PS.PelletSumfile(sumfile_loc=sumfile_baylor,gfile=gfile)
    ax[ii].plot(pelsum_na.profiles['rho_t'],pelsum_na.profiles['ne(tpel-)']*1.0e-19,'k-',label='Initial')
    l1 = ax[ii].plot(pelsum_na.profiles['rho_t'],pelsum_na.profiles['ne(tpel+)']*1.0e-19,color=colors[5],linestyle='-',linewidth='1.5',label='No Drift')
    l2 = ax[ii].plot(pelsum_HPI2.profiles['rho_t'],pelsum_HPI2.profiles['ne_d(tpel+)']*1.0e-19,linewidth='1.5',color=colors[3],linestyle='-',label='HPI2')
    l3 = ax[ii].plot(pelsum_baylor.profiles['rho_t'],pelsum_baylor.profiles['ne_d(tpel+)']*1.0e-19,color=colors[0],linestyle='-',linewidth='1.5',label='DIII-D')
    ax[ii].vlines(rho_depth_na[0,count[ii],3],3.5,15.0,color=colors[4],linestyle='--',alpha=0.5)
    ax[ii].vlines(rho_depth_HPI2[0,count[ii],3],3.5,15.0,color=colors[3],linestyle='--',alpha=0.5)
    ax[ii].vlines(rho_depth_baylor[0,count[ii],4],3.5,15.0,color=colors[0],linestyle='--',alpha=0.5)
    labels = [(r"$\rho$ = %0.2f" % rho_depth_na[0,count[ii],3]),(r"$\rho$ = %0.2f" % rho_depth_HPI2[0,count[ii],3]),
              (r"$\rho$ = %0.2f" % rho_depth_baylor[0,count[ii],3])]
    ax[ii].annotate(
                labels[0],
                xy=(0.02, 0.55), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                color=colors[4],fontsize='large', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))
    ax[ii].annotate(
                labels[1],
                xy=(0.02, 0.4), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                color=colors[3],fontsize='large', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))
    ax[ii].annotate(
                labels[2],
                xy=(0.02, 0.25), xycoords='axes fraction',
                xytext=(+0.5, -0.5), textcoords='offset fontsize',
                color=colors[0],fontsize='large', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))
    # ax[ii].axvspan(rho_max_na[0,count[ii],3],rho_max_HPI2[0,count[ii],3],color='orange',alpha=0.4)
    # ax[ii].axvspan(rho_max_na[0,count[ii],3],rho_max_baylor[0,count[ii],4],color='m',alpha=0.2)
    # ax[ii].vlines(rho_depth_na[0,count[ii],3],3.5,15.0,color='k',linestyle='--',alpha=0.7)
    # ax[ii].vlines(rho_depth_HPI2[0,count[ii],3],3.5,15.0,color='orange',linestyle='--',alpha=0.7)
    # ax[ii].vlines(rho_depth_baylor[0,count[ii],4],3.5,15.0,color='m',linestyle='--',alpha=0.7)
    ax[ii].set_xlim([0.0,1.015])
    ax[ii].set_ylim([3.5,15.0])
    ax[ii].set_ylabel(r"$n_e$ ($\times 10^{19}$ m$^{-3}$)")
    # axins = ax[ii].inset_axes([0.2, 0.2, 0.35, 0.75],xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    # axins.plot(x[mask_x],pelsum_na.profiles['ne(tpel+)'][mask_x]*1.0e-19,'k--',linewidth='2.0')
    # axins.plot(x[mask_x],pelsum_HPI2.profiles['ne_d(tpel+)'][mask_x]*1.0e-19,linewidth='2.0',color='orange',linestyle='--')
    # axins.plot(x[mask_x],pelsum_baylor.profiles['ne_d(tpel+)'][mask_x]*1.0e-19,'m--',linewidth='2.0')
    if ii!=2:
        ax[ii].set_xticks([])

labels = ['(a)','(b)','(c)']

for ii in range(3):
    ax[ii].annotate(
            labels[ii],
            xy=(0, 1), xycoords='axes fraction',
            xytext=(+0.5, -0.5), textcoords='offset fontsize',
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1.0', edgecolor='none', pad=3.0))
    

ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncols=4, mode="expand", borderaxespad=0.,
                      fontsize=14)

# fig.legend(["Initial","No Drift",'HPI2','DIII-D'],fontsize=10)
ax[2].set_xlabel(r'$\rho$')
# plt.tight_layout()
# plt.savefig("/Users/gz6/Documents/TokE/papers/PFS_figures/ne_vs_rho.png")
plt.show()


# %%
