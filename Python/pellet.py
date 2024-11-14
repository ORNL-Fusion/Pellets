import os
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from scipy import interpolate

import EFIT.equilParams_class as epc
import pellets.PelletNml as pnml
import pellets.PelletSumfile as psf


# plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})


# ----------------------------------------------------------------------


def write_profiles_from_pfile(pfile_loc, gfile_loc, efit_tree = None,
                              newfilename = None, plot_mapping = False):
    """
    Read a p file and make a new file in the format that can be read by the pellet code
    
      tested 2018/06/21 -> This reproduces p file data correctly, but doesn't convert psin to rho
    """
    # If not given, name the new file name using the extension on the p file
    
    if newfilename is None:
        folder_ind = pfile_loc.rfind('/')
        if folder_ind == -1:
            pfile_ext = pfile_loc[1:]
        else:
            pfile_ext = pfile_loc[folder_ind + 2]
        
        newfilename = 'prof_' + pfile_ext
    
    pfile_dict = _read_pfile(pfile_loc)
    
    psin = pfile_dict['psinorm']
    ne = pfile_dict['ne(10^20/m^3)'] * 10  # needs to be in units of 10^19/m^3
    te = pfile_dict['te(KeV)']
    try:
        ti = pfile_dict['ti(KeV)']
    except:
        ti = pfile_dict['ti(keV)']
    npts = len(psin)
    
    gfile = epc.equilParams(gfile_loc, tree = efit_tree)
    psin_g = gfile.PSIdict['psiN1D']
    torfluxnorm_g = gfile.getTorPsi()['psitorN1D']
    
    rho = np.sqrt(np.interp(psin, psin_g, torfluxnorm_g))
    
    if plot_mapping:
        plt.figure()
        plt.plot(psin_g, torfluxnorm_g, label = 'Tor flux')
        plt.plot(psin, rho, label = 'rho')
        plt.xlabel('$\psi_N$')
        plt.legend(loc='best')
        plt.show(block=False)

    with open(newfilename, 'w') as f:
        f.write('rho  ne  te  ti\n')
        f.write(str(npts) + '\n')
        for i in range(npts):
            f.write(str(rho[i]) + '  ' + str(ne[i]) + '  ' + str(te[i]) + '  ' + str(ti[i]) + '\n')


# ----------------------------------------------------------------------


def write_profiles_from_txt(folder, newfilename='pelprofs'):
    """
    Read a series of text files and combine into a profiles file readable by the PELLET code
    """
    if folder[-1] != '/':
        folder += '/'

    params = ['rho', 'ne', 'te', 'ti']
    filenames = ['rho_tor_norm', 'electron_density', 'electron_temperature', 'D_temperature']
    norm = [1, 1e19, 1e3, 1e3]

    profiles = {}

    for pi, param in enumerate(params):
        with open(folder + filenames[pi], 'r') as f:
            lines = f.readlines()

        npts = len(lines)

        profiles[param] = np.zeros(npts)
        for i, l in enumerate(lines):
            profiles[param][i]=float(l) / norm[pi]

    with open(newfilename, 'w') as f:
        f.write('rho  ne  te  ti\n')
        f.write(str(npts) + '\n')
        for i in range(npts):
            f.write(str(profiles['rho'][i]) + '  ' + str(profiles['ne'][i]) + '  ' + str(profiles['te'][i]) +
                    '  ' + str(profiles['ti'][i]) + '\n')


# ----------------------------------------------------------------------


def _read_pickled_profiles(profiles_file, verbose = False):
    import pickle

    if verbose: print('Reading profile fit data from saved pickle file: ' + profiles_file)

    with open(profiles_file, 'rb') as f:
        profiles = pickle.load(f)

    return profiles

# ----------------------------------------------------------------------


def write_profiles_from_pickle(pickled_file_loc = '/home/wilcoxr/pellets/SPARC/SPARC_V1E_transp_3.pkl',
                               gfile_loc = '/home/wilcoxr/pellets/SPARC/g100000.01000',
                               newfilename=None, plot_mapping=False):
    """
    Read a pickled profile file and make a new file in the format that
    can be read by the pellet code
    """
    # If not given, name the new file name using the given filename
    if newfilename is None:
        folder_ind = pickled_file_loc.rfind('/')

        if folder_ind == -1:
            filename_ext = pickled_file_loc[:-4]
        else:
            filename_ext = pickled_file_loc[folder_ind + 1:-4]

        newfilename = 'prof_' + filename_ext

    profilefile_dict = _read_pickled_profiles(pickled_file_loc)

    rho_pol = profilefile_dict['rho_pol']  # sqrt poloidal flux
    psin_approx = profilefile_dict['pol_flux'] / np.max(np.abs(profilefile_dict['pol_flux']))
    rho_tor = profilefile_dict['rho_tor']
    te = profilefile_dict['te']
    ti = profilefile_dict['ti']
    ne = profilefile_dict['ne'] * 10  # needs to be in units of 10^19/m^3

    npts = len(rho_pol)

    gfile = epc.equilParams(gfile_loc)
    psin_g = gfile.PSIdict['psiN1D']
    torfluxnorm_g = gfile.getTorPsi()['psitorN1D']

    rho_calc = np.sqrt(np.interp(rho_pol**2, psin_g, torfluxnorm_g))

    if plot_mapping:
        plt.figure()
        plt.plot(psin_g, torfluxnorm_g, lw=2, label='gfile')
        plt.plot(rho_pol**2, rho_tor**2, label='Given')
        plt.plot(rho_pol**2, rho_calc**2, label='rho calc')
        # plt.plot(rho_pol**2, psin_approx, label = '$\psi_N$ calc')  # This matches
        plt.xlabel('$\psi_N$')
        plt.ylabel('Norm tor flux')
        plt.legend(loc='best')
        plt.grid('on')
        plt.show(block=False)

    else:
        with open(newfilename, 'w') as f:
            f.write('rho  ne  te  ti\n')
            f.write(str(npts) + '\n')
            for i in range(npts):
                f.write(str(rho_tor[i]) + '  ' + str(ne[i]) + '  ' +
                        str(te[i]) + '  ' + str(ti[i]) + '\n')

# ----------------------------------------------------------------------


def plot_SPARC_profiles(value,
                        profiles_file = '/home/wilcoxr/pellets/SPARC/SPARC_V1E_transp_3.pkl'):

    profs = _read_pickled_profiles(profiles_file)
    # psin = np.array(profs['pol_flux']) / np.max(np.abs(profs['pol_flux']))

    rho_pol = profs['rho_pol']
    rho_tor = np.array(profs['rho_tor'])

    plt.figure()
    plt.plot(rho_pol, profs[value], lw=2)
    plt.grid('on')
    plt.ylabel(value)
    plt.xlabel('$\psi_n$')
    plt.show(block=False)

# ----------------------------------------------------------------------


def _read_pfile(pfile_loc):
    """
    Read in the kinetic profiles from a p file to be used as inputs (successfully tested 2018/1/3)
    
    Returns a dictionary with a non-intuitive set of keys (units are included)
    """
    with open(pfile_loc, mode = 'r') as pfile:
        lines = pfile.readlines()
    
    profiles = {}
    nprofs = 0  # counter for total number of profiles so far
    linestart = 0  # counter for which line to start at for each profile
    nlines_tot = len(lines)
    
    while True:
        # Read the header line for each profile first
        lin1 = lines[linestart].split()
        npts_prof = int(lin1[0])
        
        xname = lin1[1]
        yname = lin1[2]
        dyname = ''.join(lin1[3:])[:-1]
        
        # Generate and populate the profile arrays
        x = np.zeros(npts_prof)
        y = np.zeros(npts_prof)
        dy = np.zeros(npts_prof)
        for i in range(npts_prof):
            split_line = lines[linestart + i+1].split()
            x[i] = float(split_line[0])
            y[i] = float(split_line[1])
            dy[i] = float(split_line[2][:-1])

        # profiles[xname + '_' + yname] = x  # psinorm
        profiles[xname] = x
        profiles[yname] = y
        profiles[dyname] = dy
            
        nprofs += 1
        linestart += 1 + npts_prof
            
        if linestart >= nlines_tot:
            break
        
    # Check if all psinorms are the same, consolidate them if they are
    # (they are, don't bother separating)
    
    # condense = True
    # psinorm = None
    # for k in profiles.keys():
    #     if k is None or k=='':
    #         continue
    #
    #     if k[:4] == 'psin':
    #         if psinorm is None:
    #             psinorm = profiles[k]
    #
    #         if max(abs(profiles[k] - psinorm)) > 1e-5:
    #             condense = False
    #             break
                
    # if condense:
    #     profiles = {key: value for key, value in profiles.items()
    #                 if key[:4] != 'psin' or key is None or key==''}
    #     profiles['psinorm'] = psinorm
        
    return profiles

# ----------------------------------------------------------------------


def _get_pellet_path(gfile, device = 'd3d', port = 'tang', dist_outside = 0.003,
                     npts_interp = 200, surf_den = 10, plotit = True):
    """
    Get the pellet start and end points given the port and plasma boundary

    Inputs:
      dist_outside   Distance in R to start and end the pellet path just a bit outside the plasma
      npts_interp    Number of points to interpolate along the path to pick start and end points
      surf_den       Factor to increase the number of points on the outer flux surface
                     (needs to be dense enough to find intersection)
    """
    # Get starting and ending points
    R_traj2, Z_traj2, phi_traj2 = get_pellet_traj(device = device, port = port)
    
    # Convert to X, Y, Z to get curved R and Z for the rest of the trajectory
    
    x_traj2 = R_traj2 * np.cos(phi_traj2)
    y_traj2 = R_traj2 * np.sin(phi_traj2)
    
    xtraj = np.linspace(x_traj2[0], x_traj2[1], npts_interp)
    ytraj = np.linspace(y_traj2[0], y_traj2[1], npts_interp)
    Rtraj = np.sqrt(xtraj ** 2 + ytraj ** 2)
    Ztraj = np.linspace(Z_traj2[0], Z_traj2[1], npts_interp)
    phitraj = np.arctan2(ytraj, xtraj)
    
    gf = epc.equilParams(gfile)
    LCFS = gf.get_allFluxSur()[-1]
    npts_LCFS = len(LCFS['Rs'])
    i_surf_R = np.zeros(npts_LCFS)
    i_surf_Z = np.zeros(npts_LCFS)
    
    # Step the LCFS out by "dist_outside" to get i_surf (initializing surface)
    
    Rmid = np.mean(LCFS['Rs'])
    Zmid = np.mean(LCFS['Zs'])
    for i in range(npts_LCFS):
        theta = np.arctan2(LCFS['Zs'][i] - Zmid, LCFS['Rs'][i] - Rmid)
        i_surf_R[i] = LCFS['Rs'][i] + dist_outside * np.cos(theta)
        i_surf_Z[i] = LCFS['Zs'][i] + dist_outside * np.sin(theta)
    
    # Make isurf more dense if necessary
    
    if npts_LCFS < npts_interp:
        isurf_R = np.zeros((npts_LCFS - 1) * surf_den)
        isurf_Z = np.zeros((npts_LCFS - 1) * surf_den)
        for i in range(npts_LCFS - 1):
            isurf_R[i * surf_den: (i + 1) * surf_den] = np.linspace(i_surf_R[i],
                                                                    i_surf_R[i + 1], surf_den)
            isurf_Z[i * surf_den: (i + 1) * surf_den] = np.linspace(i_surf_Z[i],
                                                                    i_surf_Z[i + 1], surf_den)
    
    else:
        isurf_R = i_surf_R
        isurf_Z = i_surf_Z
    
    # Initialize variables that give distance to i_surf from the current starting and ending points
    # then update them as better points are found
    
    start_point_ind = 0
    
    starting_dist_sq_all = (isurf_R - Rtraj[start_point_ind]) ** 2 + \
                           (isurf_Z - Ztraj[start_point_ind]) ** 2
    nearest_LCFS_ind = np.argmin(starting_dist_sq_all)
    starting_dist_sq = (isurf_R[nearest_LCFS_ind] - Rtraj[start_point_ind]) ** 2 + \
                       (isurf_Z[nearest_LCFS_ind] - Ztraj[start_point_ind]) ** 2
    
    # Loop through points along the path, starting at the beginning
    # if it's a better fit for the start than the current choice, then update it
    # Once the distance to i_surf starts getting bigger, the starting point should
    # be found and we can move on
    
    for i in range(1, npts_interp):
        
        # check if the starting point is closer to isurf, update if so
        
        dist_from_isurf_sq = (isurf_R - Rtraj[i]) ** 2 + (isurf_Z - Ztraj[i]) ** 2
        
        # nearest point to the current pellet location
        nearest_ind_start = np.argmin(dist_from_isurf_sq)
        
        if dist_from_isurf_sq[nearest_ind_start] <= starting_dist_sq:
            
            start_point_ind = i
            starting_dist_sq = dist_from_isurf_sq[nearest_ind_start]
        
        else:
            # If it's not getting closer, then it's getting further away
            # and the starting point has been found
            break
    
    # Check if ending point is inside, if so, then we're done
    
    endpath_theta = np.arctan2(Z_traj2[-1] - Zmid, R_traj2[-1] - Rmid)
    closest_isurf_ind = np.argmin(np.arctan2(isurf_Z - Zmid, isurf_R - Rmid) - endpath_theta)
    r2_isurf = (isurf_R[closest_isurf_ind] - Rmid) ** 2 + (isurf_Z[closest_isurf_ind] - Zmid) ** 2
    r2_endpt = (R_traj2[-1] - Rmid) ** 2 + (Z_traj2[-1] - Zmid) ** 2
    
    if r2_endpt < r2_isurf:
        end_point_ind = -1
    
    else:
        
        # If not, then start at the end and work backwards until we find it
        
        end_point_ind = -1
        
        end_dist_sq_all = (isurf_R - Rtraj[end_point_ind]) ** 2 + \
                          (isurf_Z - Ztraj[end_point_ind]) ** 2
        nearest_LCFS_ind = np.argmin(end_dist_sq_all)
        end_dist_sq = (isurf_R[nearest_LCFS_ind] - Rtraj[end_point_ind]) ** 2 + \
                      (isurf_Z[nearest_LCFS_ind] - Ztraj[end_point_ind]) ** 2
        
        # Cycle through points backwards the same way we cycled forward to find the starting point
        for i in [-k for k in range(2, npts_interp - start_point_ind - 2)]:
            
            # check if the starting point is closer to the LCFS + dist_outside, update if so
            
            dist_from_isurf_sq = (isurf_R - Rtraj[i]) ** 2 + (isurf_Z - Ztraj[i]) ** 2
            
            # nearest point to the current pellet location
            nearest_ind = np.argmin(dist_from_isurf_sq)
            
            if dist_from_isurf_sq[nearest_ind] < end_dist_sq:
                
                end_point_ind = i
                end_dist_sq = dist_from_isurf_sq[nearest_ind]
            
            else:
                # If it's not getting closer, then it's getting further away
                # and the starting point has been found
                break
    
    if plotit:
        wall = gf.g['wall']
        
        plt.figure()
        for w in range(len(wall[:, 0]) - 1):
            plt.plot([wall[w, 0], wall[w + 1, 0]], [wall[w, 1], wall[w + 1, 1]], 'k')
        plt.plot(LCFS['Rs'], LCFS['Zs'], 'r')
        plt.plot(isurf_R, isurf_Z, 'm')
        plt.plot(isurf_R[nearest_ind_start], isurf_Z[nearest_ind_start], '.k', mew = 3)
        plt.plot(Rtraj, Ztraj, '.b')
        plt.plot(Rtraj[start_point_ind], Ztraj[start_point_ind], 'kx', markersize = 12, mew = 2)
        plt.plot(Rtraj[end_point_ind], Ztraj[end_point_ind], 'k+', markersize = 12, mew = 2)
        plt.axis('equal')
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
        plt.show(block = False)
        
    phitraj[phitraj < 0] += 2*np.pi
        
    RZ_out = [[Rtraj[start_point_ind], phitraj[start_point_ind], Ztraj[start_point_ind]],
             [Rtraj[end_point_ind], phitraj[end_point_ind], Ztraj[end_point_ind]]]
    
    # The pellet code takes RZ coordinates, no need to calculate this
    # x_start = Rtraj[start_point_ind] * np.cos(phitraj[start_point_ind])
    # y_start = Rtraj[start_point_ind] * np.sin(phitraj[start_point_ind])
    # x_end = Rtraj[end_point_ind] * np.cos(phitraj[end_point_ind])
    # y_end = Rtraj[end_point_ind] * np.sin(phitraj[end_point_ind])
    #
    # XYZ_out = [[x_start, y_start, Ztraj[start_point_ind]],
    #         [x_end, y_end, Ztraj[end_point_ind]]]
    # return {'RZ': RZ_out, 'XYZ': XYZ_out}
    
    return RZ_out

# ----------------------------------------------------------------------


def plot_pellet_trajectory(gfile, device = 'd3d', ports = ['tang'], efit_tree=None, npts = 100,
                           plot_TS_points = False, surf_color = '0.5', 
                           surflabel = None, contours = None):
    """
    Plot the pellet trajectory with the plasma shape
    
    To get it from MDS+, use the file name to give shot and time:
    pel.plot_pellet_trajectory('g175150.02000',efit_tree='EFIT05')
    """
    # If they only gave one port instead of a list, convert it to a 1-element list
    if len(ports[0]) == 1:
        ports = [ports]
    if contours is None:
        contours = np.arange(0.1,1,0.1)
    
    # Get surfaces and wall from g-file
    
    gf = epc.equilParams(gfile, tree = efit_tree)
    
    wall = gf.g['wall']
    gR = gf.g['R']
    gZ = gf.g['Z']
    psiRZn = gf.g['psiRZn']
    
    plt.figure()
    for w in range(len(wall[:,0])-1):
        plt.plot([wall[w, 0], wall[w+1, 0]], [wall[w, 1], wall[w+1, 1]], 'k')
    
    # surfs = gf.get_allFluxSur()  # This breaks for the SPARC gfile
    # R_LCFS = surfs[-1]['Rs']
    # Z_LCFS = surfs[-1]['Zs']    
    # plt.plot(R_LCFS, Z_LCFS, surf_color, lw=2, label = surflabel)
    if contours is not None:
        cont = plt.contour(gR, gZ, psiRZn, contours, colors = surf_color)
    plt.contour(gR, gZ, psiRZn, [1], colors = surf_color, linewidths = [2])
    
    # Loop through and plot all ports requested

    for p in ports:
        # Get trajectory (start and end points)
        R_traj2, Z_traj2, phi_traj2 = get_pellet_traj(device=device, port=p)
    
        # Convert to X, Y, Z to get curved R and Z correct for the rest of the trajectory

        x_traj2 = R_traj2 * np.cos(phi_traj2)
        y_traj2 = R_traj2 * np.sin(phi_traj2)
    
        xtraj = np.linspace(x_traj2[0], x_traj2[1], npts)
        ytraj = np.linspace(y_traj2[0], y_traj2[1], npts)
        Rtraj = np.sqrt(xtraj ** 2 + ytraj ** 2)
        Ztraj = np.linspace(Z_traj2[0], Z_traj2[1], npts)
    
        plt.plot(Rtraj, Ztraj, lw=2, label = p)
        
    if plot_TS_points:
        Z_points_TS = np.array([826,814,803,790,777,764,755,749,744,737,731,725,719,713,706,
                                700,694,687,681,675,668,661,654,648,641,634,624,610,694,580,
                                564,532,476,397,392,356,318,265,164,29]) / 1000.0
        plt.plot(np.ones(len(Z_points_TS))*1.94, Z_points_TS, '.m')
        
        Z_points_TSdiv = np.array([-1.01, -1.055, -1.088, -1.126, -1.162, -1.204, -1.221, -1.241])
        plt.plot(np.ones(len(Z_points_TSdiv)) * 1.490, Z_points_TSdiv, '.m')
        
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.axis('equal')
    plt.xlim([1, 2.5])
    plt.ylim([-1.4, 1.4])
    plt.legend(loc='upper right')
    
    
    # Also plot top-down
    
    phis = np.linspace(0,2*np.pi,500)
    
    x_iw = np.min(wall[:,0]) * np.cos(phis)  # inner wall
    y_iw = np.min(wall[:,0]) * np.sin(phis)
    x_ow = np.max(wall[:,0]) * np.cos(phis)  # outer wall
    y_ow = np.max(wall[:,0]) * np.sin(phis)
    
    plt.figure()
    plt.plot(x_iw, y_iw, 'k')
    plt.plot(x_ow, y_ow, 'k')
    
    for p in ports:
        R_traj2, Z_traj2, phi_traj2 = get_pellet_traj(device=device, port=p)
    
        x_traj2 = R_traj2 * np.cos(phi_traj2)
        y_traj2 = R_traj2 * np.sin(phi_traj2)
    
        # xtraj = np.linspace(x_traj2[0], x_traj2[1], npts)
        # ytraj = np.linspace(y_traj2[0], y_traj2[1], npts)
        
        plt.plot(x_traj2, y_traj2)
    
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.show(block = False)
    
    
    plt.show(block=False)

# ----------------------------------------------------------------------


def add_surf_to_figure(gfile, efit_tree=None, fignum=None, surfnum=-1, zorder = 1,
                       clr=None, lw=2, label = None):
    
    gf = epc.equilParams(gfile, tree = efit_tree)
    surfs = gf.get_allFluxSur()
    
    plt.figure(fignum)
    plt.hold('on')
    plt.plot(surfs[surfnum]['Rs'], surfs[surfnum]['Zs'], color = clr, lw=lw,
             zorder = zorder, label = label)

# ----------------------------------------------------------------------


def plot_shape_scan(shots = [175145, 175148, 175150, 175152, 175155],
                    times = ['02000']*5, plot_TS = True):
    """
    Just use EFIT01s
    """
    
    gfiles = ['g' + str(s) + '.' + times[i] for i, s in enumerate(shots)]
    plot_pellet_trajectory(gfiles[0], efit_tree = 'EFIT01', plot_TS_points = plot_TS,
                           surflabel = shots[0])
    fignum = plt.gcf().number
    
    for i, s in enumerate(shots[1:]):
        add_surf_to_figure(gfiles[i+1], efit_tree = 'EFIT01', fignum = fignum, label = str(s))
        
    plt.legend(loc='best')

# ----------------------------------------------------------------------


def compare_sum_files(sumfileloclist, gfile_loc, keylist=['ne(tpel-)', 'Te(tpel-)', 'delta_ne'],
                      label_ident='sum_pellet', label_len=6,
                      compare_ablation=False):
    """
    Quickly plot a single output from the Pellet code, either from the radial
    or time-dependent grids

    This will break if you try to plot things with different x-axes

    Inputs:
      key      Quantity to plot from the sumfiles, choose from:
                   delta_ne
      label_ident   Character string to search for to identify labels for legend
      label_len     Length of identifying labels for the legend
    """
    labels = []
    lins = ['-', '--', '-.']
    clrs = ['b', 'g', 'r', 'm']

    nplots = len(keylist)
    xlab = None

    # if plot_vs_psin:  # doesn't get anything outside 1
    #     gfile = epc.equilParams(gfile_loc)
    #     psin_g = gfile.PSIdict['psiN1D']
    #     torfluxnorm_g = gfile.getTorPsi()['psitorN1D']
    #     rho_to_psin = interpolate.interp1d(np.sqrt(torfluxnorm_g), psin_g)

    f, ax = plt.subplots(nplots, sharex=True)

    for i, sumfile in enumerate(sumfileloclist):
        sf = psf.PelletSumfile(sumfile, gfile_loc)

        if label_ident is not None:
            ind = sumfile.find(label_ident) + len(label_ident)
            labels.append(sumfile[ind: ind + label_len])
        else:
            labels.append(None)

        for k, key in enumerate(keylist):
            if key in sf.profiles.keys():
                ax[k].plot(sf.profiles['rho_t'], sf.profiles[key],
                           lins[np.mod(i, len(lins))] + clrs[np.mod(i, len(clrs))],
                           lw=2 + np.mod(i, len(lins)), label=labels[i])
                if xlab is None:
                    xlab = 'rho'
                elif xlab != 'rho':
                    print("WARNING: these requested keys aren't all radial profiles")

            elif key in sf.pelpath.keys():
                ax[k].plot(sf.pelpath['t'][:-1] * 1e3, sf.pelpath[key][:-1],
                           lins[np.mod(i, len(lins))] + clrs[np.mod(i, len(clrs))],
                           lw=2 + np.mod(i, len(lins)), label=labels[i])
                if xlab is None:
                    xlab = 'time (ms)'
                elif xlab != 'time (ms)':
                    print("WARNING: these requested keys aren't all time-dependent")

    ax[-1].set_xlabel(xlab)
    for k in range(nplots):
        ax[k].set_ylabel(keylist[k])
        ax[k].grid('on')

    if label_ident is not None:
        ax[0].legend(loc='best')

    if compare_ablation:
        lintyp = ['-', '--']
        plt.figure()

        for i, sumfile in enumerate(sumfileloclist):
            sf = psf.PelletSumfile(sumfile, gfile_loc)
            ablation = sf.profiles['Te(tpel-)']**(11./6) * (sf.profiles['ne(tpel-)']*1e-19)**0.33
            plt.plot(sf.profiles['rho_t'], ablation, lintyp[np.mod(i,2)], lw=2, label = labels[i])

        if label_ident is not None:
            plt.legend(loc = 'best')

        # plt.ylabel('n$_e^{0.33}$ * T$_e^{(11/6)}$')
        plt.ylabel('ne^0.33 * Te^(11/6)')
        plt.grid('on')
        plt.xlabel('rho')

    plt.show(block = False)

# ----------------------------------------------------------------------


def compare_ablation_profiles(sumfileloclist, gfile_loc, figsize_ratio = [3, 3, 5], xlims = None,
                              ylims=None, yticks=None, no_legend = False):
    """
    Plot ablation profiles from a group of sumfiles

    in paper:
        sumfiles=['178555/sum_pellet178555.2815_small.dat','165415/sum_pellet165415.4000.dat']

    d3d_upgrade:
    pel.compare_ablation_profiles(['sum_pelletd3d_upgrade_neped1.0.dat','../HFS45/sum_pelletd3d_upgrade_neped1.0.dat'],
    'g100000.0010',xlims=[0.7,1.05],ylims=[[0,15],[-0.01,3.5],[-1,35]],yticks=[[0,5,10,15],[0,1,2,3],[0,10,20,30]],no_legend=True)

    """
    label_ident = 'sum_pellet'
    # figsize_ratios = [2, 2, 4],
    # label_len = 6,  # these are both breaking for some reason when this is a variable

    labels = []
    lins = ['-', '--', '-.']
    clrs = ['b', 'g', 'r', 'm']
    xlab = 'rho'

    # if plot_vs_psin:  # doesn't get anything outside 1
    #     gfile = epc.equilParams(gfile_loc)
    #     psin_g = gfile.PSIdict['psiN1D']
    #     torfluxnorm_g = gfile.getTorPsi()['psitorN1D']
    #     rho_to_psin = interpolate.interp1d(np.sqrt(torfluxnorm_g), psin_g)

    f, ax = plt.subplots(3, sharex=True, gridspec_kw={'height_ratios': figsize_ratio})

    for i, sumfile in enumerate(sumfileloclist):
        sf = psf.PelletSumfile(sumfile, gfile_loc)

        if label_ident is not None:
            ind = sumfile.find(label_ident) + len(label_ident)
            labels.append(sumfile[ind: ind + 6])
        else:
            labels.append(None)

        ax[0].plot(sf.profiles['rho_t'], sf.profiles['ne(tpel-)'] / 1e19,
                   lins[np.mod(i, len(lins))] + clrs[np.mod(i, len(clrs))],
                   lw=2 + np.mod(i, len(lins)), label=labels[i])

        ax[1].plot(sf.profiles['rho_t'], sf.profiles['Te(tpel-)'],
                   lins[np.mod(i, len(lins))] + clrs[np.mod(i, len(clrs))],
                   lw=2 + np.mod(i, len(lins)), label=labels[i])

        ax[2].plot(sf.profiles['rho_t'], sf.profiles['delta_ne'] / 1e19,
                   lins[np.mod(i, len(lins))] + clrs[np.mod(i, len(clrs))],
                   lw=2 + np.mod(i, len(lins)), label=labels[i])

    ax[-1].set_xlabel(xlab)
    ax[0].set_ylabel('n$_e$ (10$^{19}$ m$^{-3}$)')
    ax[1].set_ylabel('T$_e$ (keV)')
    ax[2].set_ylabel('$\Delta$n$_e$ (10$^{19}$ m$^{-3}$)')

    if xlims is not None:
        plt.xlim(xlims)

    for k in range(3):
        ax[k].grid('on')
        if ylims is not None:
            ax[k].set_ylim(ylims[k])

        if yticks is None:
            yticks_old = ax[k].get_yticks()
            ax[k].set_yticks(yticks_old[::2])
        else:
            ax[k].set_yticks(yticks[k])

    if not no_legend:
        ax[0].legend(loc='lower left')


    plt.show(block = False)
    
# ----------------------------------------------------------------------


def scan_pellet_mass_velocity(top_folder, device = 'SPARC', runid_prefix = 'PFR',
                              diameters_mm = np.around(np.arange(1.8, 4.1, 0.2), 1),
                              velocities = np.arange(100, 1501, 100),
                              pellet_executable_loc = '/Users/wilcox/Codes/Pellets/src/xpellet'):
    """
    Use the same equilibrium, scan the pellet size and velocity
    'top_folder' should contain a g file, a prof file, and a nml file to copy as templates
    This will generate a new folder for each requested pellet radius and velocity
    and run xpellet there
    
    Inputs:
      pellet_executable_loc = '/home/wilcoxr/pellets/src/xpellet' on Iris
    """
    if top_folder[-1] != '/': top_folder += '/'
    os.chdir(top_folder)

    existing_nml = top_folder + 'nml_pellet.dat'
    
    for f in os.listdir(top_folder):
        if f[:2] == 'g1' or f[:6].lower() =='geqdsk':
            equilib_fileloc = f
            
        elif f[:5] == 'prof_':
            profile_fileloc = f
            
    # Read the example
    nml_orig = pnml.PelletNml(equilib_fileloc, profile_fileloc, device = device,
                              runid = 'just read',
                              existing_nml_fileloc = existing_nml)
            
    # Make a new folder for each pellet size and run xpellet in each one
    
    for d in diameters_mm:
        r_pel_in = diameter_to_spherical_radius(d) * 1.0e-3  # needs to print in m
        
        for v in velocities:
            # Need to copy over files to each folder, otherwise it breaks
            new_folder = 'd' + str(d) + '_v' + str(v)
            if runid_prefix is None:
                runid = equilib_fileloc[1:] + '_' + new_folder
            else:
                runid = runid_prefix + '_' + new_folder
        
            os.mkdir(new_folder)
            call(['cp', equilib_fileloc, new_folder])
            call(['cp', profile_fileloc, new_folder])
            os.chdir(new_folder)
        
            nml = pnml.PelletNml(equilib_fileloc, profile_fileloc, device = device,
                                 runid = runid,
                                 existing_nml_fileloc = existing_nml)
            nml.write_nml('nml_pellet.dat', input_dict = {'v_pel': v, 'r_pel' : r_pel_in})
        
            call([pellet_executable_loc])
            rename_output_files_without_spaces()
            os.remove(equilib_fileloc)
            os.remove(profile_fileloc)
            os.chdir(top_folder)
            
# ----------------------------------------------------------------------

def rename_output_files_without_spaces(folder = None):
    """
    Delete the obnoxious spaces in the PELLET output filenames
    Defaults to cwd
    """
    for f in os.listdir(folder):
        if f[:4] == 'sum_' or f[:3] == '1d_' or f[:4] == 'msg_':
            new_name = f.replace(' ', '')
            os.rename(f, new_name)
            
        
# ----------------------------------------------------------------------


def scan_pellet_mass(top_folder, mass_scale = np.arange(0.1, 1.6, 0.1),
                     pellet_executable_loc = '/Users/wilcox/Codes/Pellets/src/xpellet'):
    """
    Use the same equilibrium but scan the pellet size
    'top_folder' should contain a g file, a prof file, and an example nml file as templates
    This will generate a new folder for each requested pellet radius
    and run xpellet there
    
    Inputs:
      mass_scale Array of scale factors to apply to the original pellet mass
                 (in practice, scales the spherical pellet radius)
      pellet_executable_loc = '/home/wilcoxr/pellets/src/xpellet' on Iris
    """
    if top_folder[-1] != '/': top_folder += '/'
    os.chdir(top_folder)

    existing_nml = top_folder + 'nml_pellet.dat'
    
    for f in os.listdir(top_folder):
        if f[:2] == 'g1':
            equilib_fileloc = f
            
        elif f[:5] == 'prof_':
            profile_fileloc = f
            
    # Read the example
    nml_orig = pnml.PelletNml(equilib_fileloc, profile_fileloc, runid = 'just read',
                              existing_nml_fileloc = existing_nml)
    r_pel_orig = nml_orig.indata['r_pel']
            
    # Make a new folder for each pellet size and run xpellet in each one
    
    for m in mass_scale:
        # Need to copy over files to each folder, otherwise it breaks
        new_folder = 'mass' + str(m)
        runid = equilib_fileloc[1:] + '_' + new_folder
        
        os.mkdir(new_folder)
        call(['cp', equilib_fileloc, new_folder])
        call(['cp', profile_fileloc, new_folder])
        os.chdir(new_folder)

        nml = pnml.PelletNml(equilib_fileloc, profile_fileloc, runid = runid,
                             existing_nml_fileloc = existing_nml)
        # Convert from mass scale to radius, adjust r_pel
        r_pel_in = m**(1./3) * r_pel_orig
        nml.write_nml('nml_pellet.dat', input_dict = {'r_pel' : r_pel_in})
        
        call([pellet_executable_loc])
        os.chdir(top_folder)

# ----------------------------------------------------------------------


def scan_pellet_velocity(top_folder, velocities = np.arange(50,301,10),
                         pellet_executable_loc = '/home/shared/bin/xpellet'):
    """
    Use the same equilibrium but scan the pellet velocity
    'top_folder' should contain a g file, a prof file, and an example nml file as templates
    This will generate a new folder for each requested pellet velocity
    and run xpellet there

    Inputs:
      pellet_executable_loc = '/home/wilcoxr/pellets/src/xpellet' on Iris
       '/home/shared/bin/xpellet' on BB8
    """
    if top_folder[-1] != '/': top_folder += '/'
    os.chdir(top_folder)

    existing_nml = top_folder + 'nml_pellet.dat'

    for f in os.listdir(top_folder):
        if f[:2] == 'g1':
            equilib_fileloc = f
    
        elif f[:5] == 'prof_' or f[:8] == 'pelprofs':
            profile_fileloc = f

    nml = pnml.PelletNml(equilib_fileloc, profile_fileloc, runid = 'just read',
                         existing_nml_fileloc = existing_nml)
    # v_pel_orig = nml_orig.indata['v_pel']

    # Make a new folder for each pellet size and run xpellet in each one

    for v in velocities:
        # Need to copy over files to each folder, otherwise it breaks
        new_folder = 'vel' + str(v)
        runid = equilib_fileloc[1:] + '_' + new_folder
    
        os.mkdir(new_folder)
        call(['cp', equilib_fileloc, new_folder])
        call(['cp', profile_fileloc, new_folder])
        os.chdir(new_folder)
    
        # nml = pnml.PelletNml(equilib_fileloc, profile_fileloc, runid = runid,
        #                      existing_nml_fileloc = existing_nml)
        nml.write_nml('nml_pellet.dat', input_dict = {'v_pel': v})
    
        call([pellet_executable_loc])
        os.chdir(top_folder)

# ----------------------------------------------------------------------

def diameter_to_spherical_radius(diameter, length):
    # Assume cylinder with same length as diameter
    # 4/3 pi r**3 = L pi (d/2)**2
    srad = (3.0/16 * length * diameter**2)**(1.0/3)
    return srad

# ----------------------------------------------------------------------

def plot_pendepth_vs_mv(top_folder = '/Users/wilcox/Codes/Pellets/runs/SPARC/LFS_AXIS',
                        scan_ident = 'd', cmap = 'inferno_r', plot_vs_mass = False):
    """
    Gather up all runs from a mass and velocity scan
    """
    runs = PelletRuns(top_folder, scan_ident = scan_ident)
    for f in os.listdir(top_folder):
        if f[:2] == 'g1' or f[:6].lower() =='geqdsk':
            equilib_filename = f
    port = top_folder[top_folder.rfind('/')+1:]
    
    nruns = len(runs.folders)
    peldepth = np.zeros(nruns)
    frac_rad_pen = np.zeros(nruns)  # as reported directly by the code
    diameters = np.zeros(nruns)
    r_eff = np.zeros(nruns)
    velocities = np.zeros(nruns)
    
    for i, sfloc in enumerate(runs.sumfilelist):
        sumfile = psf.PelletSumfile(sfloc, gfile = None)
        # sumfile = psf.PelletSumfile(sfloc, runs.gfile)
        
        vind = runs.folders[i].find('_v')
        diameters[i] = float(runs.folders[i][1:vind])
        velocities[i] = int(sumfile.v_pel)
        r_eff[i] = sumfile.r_pel
        
        peldepth[i] = sumfile.get_penetration_depth()
        frac_rad_pen[i] = sumfile.fracradpen
        
    ndiams = len(set(diameters))
    nvels = len(set(velocities))
    d_grid = np.reshape(diameters, (ndiams, nvels))
    v_grid = np.reshape(velocities, (ndiams, nvels))
    peldepth_grid = np.reshape(peldepth, (ndiams, nvels))
    frac_rad_pen_grid = np.reshape(frac_rad_pen, (ndiams, nvels))
    
    v_order = np.argsort(v_grid[0,:])
    
    
    
    plt.figure()
    
    # meshobj = plt.pcolormesh(v_grid, d_grid, peldepth_grid, cmap = cmap)
    meshobj = plt.contourf(v_grid[:, v_order], d_grid,
                           peldepth_grid[:, v_order], cmap = cmap, levels=50)
    # meshobj = plt.contourf(v_grid[:, v_order], d_grid,
    #                        frac_rad_pen_grid[:, v_order], cmap = cmap, levels=50)
    
    cbar = plt.colorbar(meshobj)
    cbar.set_label('Penetration depth ($\\rho$)', rotation = 90)
        
    # z_grid, r_grid = np.meshgrid(self.z, self.r)
    # plt.pcolormesh(r_grid, z_grid, phase_deg, cmap = cmap)
    
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Diameter (mm)')
    plt.xticks(list(set(velocities)))
    plt.yticks(list(set(diameters)))
    plt.title(equilib_filename + ', trajectory: ' + port)
    # plt.ylabel('Penetration depth (rho)')
    plt.show(block=False)
    
# ----------------------------------------------------------------------


def plot_pendepth_vs_v_all_trajs(top_folder = '/Users/wilcox/Codes/Pellets/runs/SPARC/',
                                 trajs = ['LFS_DIV', 'LFS_AXIS', 'LFS_NORM', 'PFR', 'HFS', 'PFR_UP', 'HFS_UP'],
                                 scan_ident = 'd3.0_v'):
    """
    Given fixed pellet mass, plot velocity scans for each injection location
    """
    for f in os.listdir(top_folder):
        if f[:2] == 'g1' or f[:6].lower() =='geqdsk':
            equilib_filename = f
            
    if top_folder[-1] != '/': top_folder += '/'
            
    velocities = {}
    frac_rad_pen = {}
    peldepth = {}
            
    for t in trajs:
        runs = PelletRuns(top_folder + t + '/', scan_ident = scan_ident)
        
        nruns = len(runs.folders)
        pel_d = np.zeros(nruns)
        frac_pen = np.zeros(nruns)  # as reported directly by the code
        vel = np.zeros(nruns)
        
        for i, sfloc in enumerate(runs.sumfilelist):
            sumfile = psf.PelletSumfile(sfloc, gfile = None)
            # sumfile = psf.PelletSumfile(sfloc, runs.gfile)
            
            vind = runs.folders[i].find('_v')
            vel[i] = int(sumfile.v_pel)
            
            pel_d[i] = sumfile.get_penetration_depth()
            frac_pen[i] = sumfile.fracradpen
        
        v_order = np.argsort(vel)
        velocities[t] = vel[v_order]
        frac_rad_pen[t] = frac_pen[v_order]
        peldepth[t] = pel_d[v_order]
        
    
    plt.figure()
    
    for t in trajs:
        plt.plot(velocities[t], 1 - peldepth[t], lw=2, label = t)
        # plt.plot(velocities[t], frac_rad_pen[t], lw=2, label = t)
    
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Penetration depth (1 - rho)')
    plt.legend(loc='best')
    plt.grid('on')
    plt.show(block=False)
    
# ----------------------------------------------------------------------


class PelletRuns(object):
    """
    Class to gather and manipulate data from a scan of pellet runs
    """
    def __init__(self, top_folder, scan_ident = 'd'):
        if top_folder[-1] != '/': top_folder += '/'
        
        folders = [f for f in os.listdir(top_folder) if f[:len(scan_ident)] == scan_ident]

        sumfilelist = []
        nmllist = []
        for f in folders:
            sumfile = [s for s in os.listdir(top_folder + f) if s[:10] == 'sum_pellet']
            if len(sumfile) > 1:
                print('WARNING: Multiple sum files in folder:')
                print('  ' + f)
                print(' Slecting first one:')
                print(sumfile[0])
            elif len(sumfile) == 0:
                print('ERROR: No sumfile found in folder:')
                print(top_folder + f)
                return
    
            sumfilelist.append(top_folder + f + '/' + sumfile[0])
            
            nmlfile = [s for s in os.listdir(top_folder + f) if s == 'nml_pellet.dat']
            if len(sumfile) > 1:
                print('WARNING: Multiple nml files in folder:')
                print('  ' + f)
                print(' Slecting first one:')
                print(nmlfile[0])
            elif len(nmlfile) == 0:
                print('ERROR: No nml file found in folder:')
                print(top_folder + f)
                return

            nmllist.append(top_folder + f + '/' + nmlfile[0])

        gfiles_in_top = [g for g in os.listdir(top_folder) if g[:2] == 'g1' or g[:6] == 'geqdsk']

        if len(gfiles_in_top) > 1:
            print('WARNING: Multiple g files in top folder:')
            print('  ' + top_folder)
            print(' Slecting first one:')
            print(gfiles_in_top[0])
        elif len(gfiles_in_top) == 0:
            print('ERROR: No g file found in top folder:')
            print(top_folder)
            return
        
        self.gfile = gfiles_in_top[0]
        self.sumfilelist = sumfilelist
        self.nmllist = nmllist
        self.folders = folders
        self.scan_ident = scan_ident

    # ------------------------------------------------------------------

    def plot_key(self, key = 'delta_ne', label_ident = 'mass', label_len = 3):
        # if top_folder[-1] != '/': top_folder += '/'
        #
        # mass_scales = [m for m in os.listdir(top_folder) if m[:4] == 'mass']
        #
        # sumfilelist = []
        # for m in mass_scales:
        #     sumfile = [s for s in os.listdir(top_folder + m) if s[:10] == 'sum_pellet']
        #     if len(sumfile) > 1:
        #         print('WARNING: Multiple sum files in folder:')
        #         print('  ' + m)
        #         print(' Slecting first one:')
        #         print(sumfile[0])
        #     elif len(sumfile) == 0:
        #         print('ERROR: No sumfile found in folder:')
        #         print(top_folder + m)
        #         return
        #
        #     sumfilelist.append(top_folder + m + '/' + sumfile[0])
        #
        # gfiles_in_top = [g for g in os.listdir(top_folder) if g[:2] == 'g1']
        #
        # if len(gfiles_in_top) > 1:
        #     print('WARNING: Multiple g files in top folder:')
        #     print('  ' + top_folder)
        #     print(' Slecting first one:')
        #     print(gfiles_in_top[0])
        # elif len(gfiles_in_top) == 0:
        #     print('ERROR: No g file found in top folder:')
        #     print(top_folder)
        #     return
    
        compare_sum_files(self.sumfilelist, self.gfile, key = key,
                          label_ident = label_ident, label_len = label_len)

    # ----------------------------------------------------------------------
        
    def plot_penetration_depth(self):
        
        peldepth = np.zeros(len(self.sumfilelist))
        for i, sfloc in enumerate(self.sumfilelist):
            sumfile = psf.PelletSumfile(sfloc, self.gfile)
            peldepth[i] = sumfile.get_penetration_depth()
            
        xvar = [np.float(n[len(self.scan_ident):]) for n in self.folders]
        
        # Reorder everything by value of xvar
        xvar, peldepth = zip(*sorted(zip(xvar, peldepth)))
        
        
        plt.figure()
        plt.plot(np.array(xvar), (1-np.array(peldepth)), lw=2)
        plt.xlabel(self.scan_ident)
        plt.ylabel('Penedtration depth (rho)')
        plt.grid('on')
        plt.show(block=False)

# ----------------------------------------------------------------------


def __plot_pellet_radius_scan(top_folder, key = 'delta_ne', label_ident = '_r', label_len = 5):
    """
    **Deprecated since I started scanning pellet mass instead of radius
    """
    if top_folder[-1] != '/': top_folder += '/'
    
    r_pels = [r for r in os.listdir(top_folder) if r[-2:] == 'mm']
    
    sumfilelist = []
    for r in r_pels:
        sumfile = [s for s in os.listdir(top_folder + r) if s[:10] == 'sum_pellet']
        if len(sumfile) > 1:
            print('WARNING: Multiple sum files in folder:')
            print('  ' + r)
            print(' Slecting first one:')
            print(sumfile[0])
        elif len(sumfile) == 0:
            print('ERROR: No sumfile found in folder:')
            print(top_folder + r)
            return
        
        sumfilelist.append(top_folder + r + '/' + sumfile[0])
    
    gfiles_in_top = [g for g in os.listdir(top_folder) if g[:2] == 'g1']
    
    if len(gfiles_in_top) > 1:
        print('WARNING: Multiple g files in top folder:')
        print('  ' + top_folder)
        print(' Slecting first one:')
        print(gfiles_in_top[0])
    elif len(gfiles_in_top) == 0:
        print('ERROR: No g file found in top folder:')
        print(top_folder)
        return
    
    compare_sum_files(sumfilelist, gfiles_in_top[0], key = key,
                      label_ident = label_ident, label_len = label_len)

# ----------------------------------------------------------------------


def _test_profiles(sumfileloc, pfile_loc, gfile_loc):
    """
    Check input and output profiles and make sure they line up (they do)
    """
    
    # Get input profile from p-file
    
    pfile_dict = read_pfile(pfile_loc)
    
    psin = pfile_dict['psinorm']
    ne = pfile_dict['ne(10^20/m^3)'] * 10  # needs to be in units of 10^19/m^3
    te = pfile_dict['te(KeV)']
    ti = pfile_dict['ti(KeV)']
    npts_pf = len(psin)
    
    # Translate psin to rho using g file
    
    gfile = epc.equilParams(gfile_loc)
    psin_g = gfile.PSIdict['psiN1D']
    torfluxnorm_g = gfile.getTorPsi()['psitorN1D']
    rho_pf = np.sqrt(np.interp(psin, psin_g, torfluxnorm_g))  # rho locations of profiles in the p-file
    
    # Get output from the 'sum_pellet...' file
    
    profiles, rho, dVol, Te_minus, ne_minus, Te_plus, ne_plus, delta_ne, delta_ne_ne = read_pellet_output(sumfileloc)
    
    # Plot and compare

    plt.figure()
    plt.plot(rho_pf, te, 'k', label = 'p file')
    plt.plot(profiles['rho_t'], profiles['Te(tpel-)'], 'b', lw=2, label='-')
    plt.plot(profiles['rho_t'], profiles['Te(tpel+)'], 'r', lw=1, label='+')
    plt.ylabel('Te (keV)')
    plt.xlabel('rho')
    plt.legend(loc='best')
    plt.grid('on')

    plt.figure()
    plt.plot(rho_pf, ne, 'k', label = 'p file')
    plt.plot(profiles['rho_t'], profiles['ne(tpel-)']*1e-19, 'b', lw=2, label='-')
    plt.plot(profiles['rho_t'], profiles['ne(tpel+)']*1e-19, 'r', lw=1, label='+')
    plt.ylabel('ne (10^19/m^3)')
    plt.xlabel('rho')
    plt.legend(loc='best')
    plt.grid('on')
    plt.show(block=False)

    # plt.figure()
    # plt.plot(rho, ne_plus-ne_minus, 'b', lw=2, label='manual')
    # plt.plot(rho, delta_ne, 'r', label='delta_ne')
    # plt.ylabel('delta ne')
    # plt.xlabel('rho')
    # plt.legend(loc='best')
    # plt.grid('on')
    #
    # plt.figure()
    # plt.plot(rho, (ne_plus-ne_minus)/ne_plus, 'b', lw=2, label='manual')
    # plt.plot(rho, delta_ne_ne, 'r', label='delta_ne')
    # plt.ylabel('delta ne / ne')
    # plt.xlabel('rho')
    # plt.legend(loc='best')
    # plt.grid('on')
    # plt.show(block=False)

# ----------------------------------------------------------------------


def _write_temp_equilib_file(wout_file, temp_fileloc):
    """
    Write out an abbreviated wout file that the pellet code can read using
    the READEQ subroutine of setup_ajax.f90
    (not tested yet)
    """
    wout_obj = woc.Wout(wout_file)
    
    phitot = np.abs(wout_obj.data['phi'][-1])
    sqrt_torflux = np.sqrt(wout_obj.data['s'])  # radial grid for iota, rz, and lambda arrays (should all be the same)
    npts = len(sqrt_torflux)
    iota_bar = wout_obj.data['iotaf']
    mnmax = wout_obj.data['mnmax']  # Total number of modes to sum over
    lmns = wout_obj.data['lmns']
    lmnc = wout_obj.data['lmnc']
    rmns = wout_obj.data['rmns']
    rmnc = wout_obj.data['rmnc']
    zmns = wout_obj.data['zmns']
    zmnc = wout_obj.data['zmnc']
    
    with open(temp_fileloc, 'w') as f:
        f.write(str(phitot) + '\n')
        f.write(str(npts) + '\n')  # nr_iota
        f.write(str(sqrt_torflux) + '\n')  # rho_iota array
        f.write(str(iota_bar) + '\n')  # iota array
        f.write(str(npts) + '\n')  # nr_rzl
        f.write(str(sqrt_torflux) + '\n')  # rho_rz array
        f.write(str(sqrt_torflux) + '\n')  # rho_lam array
        
        f.write(str(mnmax) + '\n')
        for i in range(mnmax):
            f.write(str(rmns[:,i]) + '\n')
            f.write(str(rmnc[:,i]) + '\n')
            f.write(str(zmns[:,i]) + '\n')
            f.write(str(zmnc[:,i]) + '\n')
            f.write(str(lmns[:,i]) + '\n')
            f.write(str(lmnc[:,i]) + '\n')
    
# ----------------------------------------------------------------------


def _compare_nml_files(file1, file2, verbose = True):
    """
    For testing purposes, compare the original and rewritten files to
    make sure the reading and writing scripts are working properly
    
    As of 4/4/2017, this all works -RSW
    """
    f1 = pnml.PelletNml(file1, verbose = verbose)
    f2 = pnml.PelletNml(file2, verbose = verbose)
    
    no_bad_vals = True
    for k in pnml.input_list:
        if str(f1.indata[k]) != str(f2.indata[k]):
            print(str(k) + ' = ' + str(f1.indata[k]) + ', ' + str(f2.indata[k]))
            no_bad_vals = False
            
    if no_bad_vals:
        print('All values match')
        

# ----------------------------------------------------------------------

def get_pellet_traj(device = 'd3d', port = 'tang'):

    if device.upper() == 'D3D':
        if port.upper() == 'HFS_UPPER':
            R_traj = [1.016, 1.6]
            Z_traj = [0.721, 0]
            phi_traj = [0, 0]
        elif port.upper() == 'HFS_MID':
            R_traj = [1.016, 1.6]
            Z_traj = [0.280, 0]
            phi_traj = [0, 0]
        elif port.upper() == 'OUTBOARD' or port.upper() == 'R0':
            R_traj = [2.365, 1.85]
            Z_traj = [0.089, 0.089]
            phi_traj = [0, 0]
        elif port.upper() == 'SPI':
            R_traj = [2.284, 1.016]
            Z_traj = [0.684, -0.690]
            phi_traj = [0, 0]
        elif port.upper() == 'R2':
            R_traj = [2.246, 1.016]
            Z_traj = [-1.474, -0.007]
            phi_traj = [0, 0]
        elif port.upper() == 'V1':
            R_traj = [1.479, 1.479]
            Z_traj = [1.203, 0.0]
            phi_traj = [0, 0]
        elif port.upper() == 'V3':
            R_traj = [2.098, 2.098]
            Z_traj = [1.017, 0.0]
            phi_traj = [0, 0]
        elif port.upper() == 'TANG':
            x_traj = np.array([1.6639, 1.18177])
            y_traj = np.array([-1.5583, -0.20838])
            Z_traj = [0.7618, 1.2]
            R_traj = np.sqrt(x_traj**2 + y_traj**2)  # [ 2.27966272,  1.2]
            phi_traj = np.arctan2(y_traj, x_traj)  # [-0.75263725, -0.17453462]
            # From(x, y, z) = 1.66390, -1.55831, 0.761823
            # To(x, y, z) = 1.18177, -0.208378, 1.20000
        else:
            print('Requested port not available on device ' + device + ': ' + port)
            return

    elif device.upper() == 'SPARC':
        if port.upper() == 'LFS_AXIS':  # below midplane, directed towards the axis
            # R_traj = np.array([2.0217, 1.8903])
            # Z_traj = np.array([-0.8675, 0.0])
            R_traj = np.array([2.0217, 1.956])
            Z_traj = np.array([-0.8675, -0.43375])
            phi_traj = np.array([0.0, 0.0])
        elif port.upper() == 'LFS_NORM':  # below midplane, directed more perpendicular to surface
            R_traj = np.array([2.0217, 1.8903])
            Z_traj = np.array([-0.8675, -0.6])
            phi_traj = np.array([0.0, 0.0])
        elif port.upper() == 'LFS_DIV':  # Injected from closer to divertor
            R_traj = np.array([1.8100, 1.80772])
            Z_traj = np.array([-1.06402, -0.51541])
            phi_traj = np.array([0.0, 0.27278])
        elif port.upper() == 'PFR':
            R_traj = np.array([1.48944, 1.83647])
            Z_traj = np.array([-1.08008, -0.66252])
            phi_traj = np.array([0.187326, 0.648714])
        elif port.upper() == 'HFS':  
            R_traj = np.array([1.41362, 1.760116])
            Z_traj = np.array([-1.03649, -0.68347])
            phi_traj = np.array([0.0, 0.120725])
        elif port.upper() == 'PFR_UP':  #  Flipped to upper divertor
            R_traj = np.array([1.49035, 1.83647])
            Z_traj = np.array([1.07764, 0.66252])
            phi_traj = np.array([0.19053, 0.648714])
        elif port.upper() == 'HFS_UP':  # goes to ~0.7
            R_traj = np.array([1.43256, 1.77785])
            Z_traj = np.array([1.01651, 0.66594])
            phi_traj = np.array([0.0, 0.11710])
            # R_traj = np.array([1.432561, 1.760116])  # extends to psin~0.7
            # Z_traj = np.array([-1.01651, -0.68347])
            # phi_traj = np.array([0.0, 0.11235])
            # R_traj = np.array([1.432561, 1.921058])
            # Z_traj = np.array([-1.01651, -0.52572])
            # phi_traj = np.array([0.0, 0.15196])
        else:
            print('Requested port not available on device ' + device + ': ' + port)
            return
    
    return R_traj, Z_traj, phi_traj

# ----------------------------------------------------------------------


def find_pellet_traj_drdzdrphi(gfile_loc, npts = 20, extend = 1,
                               rstart = np.array([1.463374, 1.432561]),
                               zstart = np.array([1.21927, 1.01651]),
                               dr = np.array([3.15e-5, 0.632916]),
                               dz = np.array([-0.448493, -0.666078]),
                               drphi = np.array([0.893786, 0.394662]),
                               contours = np.arange(0.1,1.0,0.1)):
    """
    For SPARC (order in defaults is ['LFS_DIV', 'PFR', 'HFS', 'LFS_AXIS'])
    
            rstart = np.array([1.81, 1.463374, 1.432561, 2.0217]),
            zstart = np.array([-1.06402, -1.21927, -1.01651, -0.8675]),
            dr = np.array([-0.09381, 3.15e-5, 0.632916, -0.1314]),
            dz = np.array([0.74454, 0.448493, 0.666078, 0.8675]),
            drphi = np.array([0.660954, 0.893786, 0.394662, 0.0]),
            
            R_traj = np.array([2.0217, 1.8903])
            Z_traj = np.array([-0.8675, -0.40])
    
    Use outputs to grab starting points that are just outside LCFS and
    endpoints so that trajectory doesn't cross a flux surface tangency plane
    
    To go to psin~0.5, selected points [14] (of npts=20, extend=1) for HFS and LFS_DIV
    Both started from given starting points
    
    For PFR, extend=1.5, npts=30, go from index [6] to [24] (gets just inside psin~0.7)
    """
    gfile = epc.equilParams(gfile_loc)
    
    psiRZn = gfile.g['psiRZn']
    gR = gfile.g['R']
    gZ = gfile.g['Z']
    wall = gfile.g['wall']
    
    nlocs = len(rstart)
    # Convert R, phi to x, y (assuming phi_start=0)
    # dx = dr * np.cos(phi) - rstart * np.sin(phi) * dphi
    # dy = dr * np.sin(phi) + rstart * np.cos(phi) * dphi
    # dx = dr
    # dy = drphi
    
    xtraj = np.zeros([nlocs, npts])
    ytraj = np.zeros([nlocs, npts])
    ztraj = np.zeros([nlocs, npts])
    Rtraj = np.zeros([nlocs, npts])
    phitraj = np.zeros([nlocs, npts])
    
    plt.figure()
    for w in range(len(wall[:, 0]) - 1):
        plt.plot([wall[w, 0], wall[w + 1, 0]], [wall[w, 1], wall[w + 1, 1]], 'k')
    if contours is not None:
        cont = plt.contour(gR, gZ, psiRZn, contours, colors = 'r')
    plt.contour(gR, gZ, psiRZn, [1], colors = 'r', linewidths = [2])
    
    for i in range(nlocs):
        xtraj[i,:] = np.linspace(rstart[i], rstart[i] + extend*dr[i], npts)
        ytraj[i,:] = np.linspace(0, extend*drphi[i], npts)
        ztraj[i,:] = np.linspace(zstart[i], zstart[i] + extend*dz[i], npts)
        Rtraj[i,:] = np.sqrt(xtraj[i,:]**2 + ytraj[i,:]**2)
        phitraj[i,:] = np.arctan2(ytraj[i,:], xtraj[i,:])
        
        plt.plot(Rtraj[i,:], ztraj[i,:], 'x')
    
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.axis('equal')
    # plt.xlim([1, 2.5])
    # plt.ylim([-1.4, 1.4])
    
    phis = np.linspace(0,2*np.pi,500)
    
    x_iw = np.min(wall[:,0]) * np.cos(phis)  # inner wall
    y_iw = np.min(wall[:,0]) * np.sin(phis)
    x_ow = np.max(wall[:,0]) * np.cos(phis)  # outer wall
    y_ow = np.max(wall[:,0]) * np.sin(phis)
    
    plt.figure()
    plt.plot(x_iw, y_iw, 'k')
    plt.plot(x_ow, y_ow, 'k')
    
    for i in range(nlocs):
        xtraj = np.linspace(rstart[i], rstart[i] + dr[i], npts)
        ytraj = np.linspace(0, drphi[i], npts)
        
        plt.plot(xtraj, ytraj, 'x')
    
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.show(block = False)
    
    return Rtraj, ztraj, phitraj
