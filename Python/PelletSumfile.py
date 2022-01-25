"""
Class to read pellet output file "sum_pellet*.dat" and plot outputs

R Wilcox 2018
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy import interpolate

import EFIT.equilParams_class as epc


# plt.rcParams.update({'font.weight': 'bold'})
# plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})


# ----------------------------------------------------------------------


class PelletSumfile(object):
    """
    Class to read pellet output file "sum_pellet*.dat" and plot outputs
    """
    def __init__(self, sumfile_loc, gfile = None):
        """
        Inputs:
          sumfile_loc   File location of the output file
        """
        self.sumfile_loc = sumfile_loc
        self._read_pellet_output()
        if gfile is not None:
            self.gfile = epc.equilParams(gfile)
            self.shot_and_time = gfile[1:]
        else:
            self.shot_and_time = None

    # ----------------------------------------------------------------------

    def _read_pellet_output(self):
        """
        Returns a dictionary with all of the output arrays
        """
        with open(self.sumfile_loc, mode = 'r') as sumfile:
            lines = sumfile.readlines()
            
        for line in lines:
            if line[:8] == ' NCPLAS=':
                nc_plas = int(line[8:-2])
            elif line[:7] == ' NCSOL=':
                nc_sol = int(line[7:-2])
            elif line[:7] == ' R_PEL=':
                self.r_pel = float(line[7:-2])
            elif line[:7] == ' V_PEL=':
                self.v_pel = float(line[7:-2])
            elif line[:34] == 'Fractional radius penetrated [-] =':
                self.fracradpen = float(line[35:-2])
        
        
        prof_inds_all = [i for i, s in enumerate(lines[:-10]) if \
            (s[:16] == '               i' and lines[i+1][:16] == '               -')]
        
        last_prof_line = [i for i, line in enumerate(lines) if \
            line[:19]=='*** Pellet path ***'][0]
        
        prof_inds = [i for i in prof_inds_all if i < last_prof_line]

        nrad = nc_plas + nc_sol
        
        nlines_tot = len(lines)
        profiles = {}
        
        for pind in prof_inds:
            prof_names = lines[pind].split()
            prof_units = lines[pind + 1].split()
            
            for p in prof_names[1:]:
                profiles[p] = np.zeros(nrad)
                
            for i in range(nrad):
                lin = lines[pind + 2 + i].split()
            
                for j in range(1, len(lin)):
                    profiles[prof_names[j]][i] = lin[j]
        
        # return rho, dVol, Te_minus, ne_minus, Te_plus, ne_plus, delta_ne, delta_ne_ne
    
    
        # Read the "Pellet path" output array
    
        path_section_start = 'i              t            rho              R              Z'
        path_ind = [i for i, s in enumerate(lines) if path_section_start in s][0]
        # Index where profiles start
    
        src_section_start = 'i            src'  # number of path grid points is 2 lines above this
        src_section_ind = [i for i, s in enumerate(lines) if src_section_start in s][0]
        npts = int(lines[src_section_ind - 2].split()[0])  # number of points along the pellet path
    
        path_varnames = lines[path_ind].split()
        path_units = lines[path_ind + 1].split()
        pelpath = {}
    
        for p in path_varnames[1:]:
            pelpath[p] = np.zeros(npts)
    
        for i in range(npts):
            lin = lines[path_ind + 2 + i].split()
        
            for j in range(1, len(lin)):
                pelpath[path_varnames[j]][i] = lin[j]
    
        path_varnames2 = lines[src_section_ind].split()
        for p in path_varnames2[1:]:
            pelpath[p] = np.zeros(npts)
    
        for i in range(npts):
            lin = lines[src_section_ind + 2 + i].split()
        
            for j in range(1, len(lin)):
                pelpath[path_varnames2[j]][i] = lin[j]
    
        path_varnames3 = lines[src_section_ind + npts + 3].split()
        for p in path_varnames3[1:]:
            pelpath[p] = np.zeros(npts)
    
        for i in range(npts):
            lin = lines[src_section_ind + npts + 5 + i].split()
        
            for j in range(1, len(lin)):
                pelpath[path_varnames3[j]][i] = lin[j]
    
        self.profiles = profiles
        self.pelpath = pelpath

    # ----------------------------------------------------------------------
    
    def get_penetration_depth(self):
        """
        Returns the interpolated position in rho of the pellet where its mass goes to zero
        """
        rho_t = self.profiles['rho_t']  # radial profile locations
        dne = self.profiles['delta_ne']
        finite_inds = np.argwhere(dne > 0)[:2]
        pfit = np.polyfit(dne[finite_inds].flatten(), rho_t[finite_inds].flatten(), 1)
        dne_func = np.poly1d(pfit)
        # dne_func = interpolate.interp1d(dne[finite_inds].flatten(), rho_t[finite_inds].flatten(),
        #                                 fill_value='extrapolate')
        
        return dne_func(0)  # linearly extrapolate to dne=0
        
    # ----------------------------------------------------------------------

    def plot_single_output(self, key = 'src'):
        """
        Quickly plot a single output from the Pellet code, either from the radial
        or time-dependent grids
        """
    
        plt.figure()
        if key in self.profiles.keys():
            plt.plot(self.profiles['rho_t'], self.profiles[key], lw=2)
            plt.xlabel('rho')
        elif key in self.pelpath.keys():
            plt.plot(self.pelpath['t'][:-1] * 1e3, self.pelpath[key][:-1], lw=2)
            plt.xlabel('time (ms)')
    
        plt.ylabel(key)
        plt.grid('on')
        plt.show(block = False)
        
    # ----------------------------------------------------------------------

    def plot_density_change(self):
        """
        Quickly plot a single output from the Pellet code, either from the radial
        or time-dependent grids
        """
    
        plt.figure()
        plt.plot(self.profiles['rho_t'], self.profiles['ne(tpel-)']/1.0e20, 'b', lw=2)
        plt.plot(self.profiles['rho_t'], self.profiles['ne(tpel+)']/1.0e20, 'r', lw=2)
        plt.xlabel('rho')
    
        plt.ylabel('n$_e$ (10$^{20}$ m$^{-3}$)')
        plt.grid('on')
        plt.show(block = False)

    # ----------------------------------------------------------------------

    def reproduce_Raman_fig(self):
        """
        Reproduce the figure from Roger Raman's paper
        """
    
        plt.figure()
        
        for param in ['t', 'rpel1', 'Te0', 'ne0']:
            plt.plot(self.pelpath['R'][1:], self.pelpath[param][1:] / np.max(self.pelpath[param]),
                     label=param)
        
        plt.title(self.shot_and_time)
        plt.legend(loc='best')
        # plt.xlim([1.65,2.1])
        plt.xlabel('R (m)')
        plt.ylabel('normalized value')
        plt.grid('on')
        plt.show(block = False)

    # ----------------------------------------------------------------------
        
    def plot_trajectory(self, contours = [1.0]):
        """
        Plot the pellet trajectory along with flux surfaces
        """
        # Get flux urfaces and wall from g-file
        
        psiRZn = self.gfile.g['psiRZn']
        gR = self.gfile.g['R']
        gZ = self.gfile.g['Z']
    
        # surfs = self.gfile.get_allFluxSur()
        # R_LCFS = surfs[-1]['Rs']
        # Z_LCFS = surfs[-1]['Zs']
    
        wall = self.gfile.g['wall']
        
        # Get trajectory from output file

        Rtraj = self.pelpath['R']
        Ztraj = self.pelpath['Z']

        plt.figure()
        for w in range(len(wall[:, 0]) - 1):
            plt.plot([wall[w, 0], wall[w + 1, 0]], [wall[w, 1], wall[w + 1, 1]], 'k')
        # plt.plot(R_LCFS, Z_LCFS, 'r')
        
        if contours is not None:
            cont = plt.contour(gR, gZ, psiRZn, contours, colors = 'r')
        plt.contour(gR, gZ, psiRZn, [1], colors = 'r', linewidths = [2])
        
        for tr in range(len(Rtraj)-1):
            plt.plot([Rtraj[tr], Rtraj[tr+1]], [Ztraj[tr], Ztraj[tr+1]], 'b', lw = 2)
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
        plt.axis('equal')
        plt.xlim([1, 2.5])
        plt.ylim([-1.4, 1.4])
        plt.show(block = False)
            
# ----------------------------------------------------------------------
