"""
Class to read and write 'nml_pellet.dat' file to run the Pellet code

R. Wilcox 2018
"""

from os import path
import numpy as np
import matplotlib.pyplot as plt

import EFIT.equilParams_class as epc
import Misc.Fnml as fnml


# plt.rcParams.update({'font.weight': 'bold'})
# plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.facecolor': 'w'})
plt.rcParams.update({'mathtext.default': 'regular'})


# This is the standard order that the inputs are listed in the namelist file
# All inputs are explained in the header comments of pellet_dr.f90
input_list = ['cn_runid', 'cn_eq', 'cn_prof', 'cn_device', 'k_equil', 'k_pel', 'k_readd',
              'k_seg_p', 'ncplas', 'ncsol', 'nhorz', 'nvert', 'amu_pel', 'r_pel',
              'v_pel', 'rseg_p', 'r0', 'a0', 'dsol', 's0', 'e0', 'e1', 'd1',
              'bt0', 'cur', 'q0', 'te0', 'te1', 'px_te', 'qx_te', 'rl_te', 'den0',
              'den1', 'px_den', 'qx_den', 'rl_den', 'rvert', 'zhorz', 'pb', 'px_hb',
              'qx_hb', 'amu_b', 'eb0', 'amu_i', 'dn01', 'rl_dn0', 'pa', 'px_ha',
              'qx_ha', 'k_ped', 'k_prl', 'nprlcld', 'prlinjang', 'pedwid', 'pedte',
              'prlq0', 'prlqa', 'prlqf', 'fpelprl', 'iprlcld']

# ----------------------------------------------------------------------


class PelletNml(object):
    """
    Class to read and write 'nml_pellet.dat' file to run the Pellet code
    """
    def __init__(self, equilib_fileloc = None, profile_fileloc= None, runid = None,
                 existing_nml_fileloc = '/home/wilcoxr/pellets/test/nml_pellet.dat',
                 input_dict = None, device = None, verbose = True):
        """
        Inputs:
          equilib_fileloc       Location of the equilibrium file (g or wout)
          input_dict            Dictionary of inputs to modify before setting up the input file
          example_nml_fileloc   Location of an existing namelist input to copy as a start
                                Copies all inputs from this file and only overwrites
                                the ones it knows about
        """
        
        # Initialize the array if there was an input to do so
        if existing_nml_fileloc is not None:
            if verbose:
                print('Reading namelist input from ' + existing_nml_fileloc)
            self.indata = fnml.read_Fnml(existing_nml_fileloc)
            
        else:
            self.indata = {}
            
            # Get equilibrium inputs, either from wout or g file, or use simple 2D approx.
            if equilib_fileloc is not None:
                # self.indata['k_equil'] = 1
                self.indata['cn_eq'] = equilib_fileloc
                filepath, filename = path.split(equilib_fileloc)
    
                if filename[:4].lower() == 'wout':
                    import VMEC.Python.wout_class as woc
                    
                    if filename[-3:] != '.nc':
                        print('Sorry, only .nc wout files currently supported')
                        return
        
                    if verbose:
                        print("Loading VMEC file: " + filename)
        
                    wout = woc.Wout(equilib_fileloc)
                    
                    self.wout_data = wout.data
        
                    self.indata['bt0'] = wout.data['b0']  # be careful, may be negative
                    self.indata['cur'] = wout.data['ctor']
                    self.indata['r0'] = wout.data['Rmajor_p']
                    self.indata['q0'] = 1. / wout.data['iotaf'][0]
                    self.s = wout.data['s']
                    self.phitot = wout.data['phi']
    
                elif filename[0].lower() == 'g':
                    gfile = epc.equilParams(equilib_fileloc)
        
                    self.indata['bt0'] = gfile.bcentr  # maybe should be abs?
                    self.indata['cur'] = gfile.g['Ip']
                    self.indata['r0'] = gfile.g['R0']
                    self.indata['q0'] = gfile.g['q'][0]
                
            else:
                # These paraemters need to be set for simple 2D approximation:
                # r0, a0, s0, e0, e1, d1, bt0, cur, q0
                self.indata['r0'] = 1.7
                self.indata['a0'] = 0.565
                self.indata['dsol'] = 0.01
                self.indata['q0'] = 1.0
                self.indata['k_equil'] = 0
            
            
            # Get kinetic profile data from p file or use input profile approximation factors
                
            if profile_fileloc is not None:
                self.indata['k_readd'] = 1
                self.indata['cn_prof'] = profile_fileloc
                
            else:
                self.indata['k_readd'] = 0
                self.indata['cn_prof'] = ''
        
        if device is not None:
            self.indata['cn_device'] = device
        self.indata['cn_runid'] = runid

        if input_dict is not None:
            for k in input_dict.keys():
                self.indata[k] = input_dict[k]
        
        
    # ---------------------------------------------------------------------


    def write_nml(self, fileloc, input_dict = [], nchars = 14, verbose = True):
        """
        Generate the namelist input file for running PELLET from the loaded equilibrium
        
        Inputs:
          input_dict  Dictionary of variables to overwrite in the input dictionary
                      while leaving existing object unchanged
                      (if a variable is given in 'input_dict' but is not in 'input_list',
                       then it will not be written)
          nchars      Number of characters allotted for the variable name
        """
        if fileloc[-4:] != '.dat':
            fileloc += '.dat'
            if verbose:
                print("'.dat' extension added to requested filename")
                
        
        # Manually overwrite any input values
        for inval in input_dict:
            self.indata[inval] = input_dict[inval]
            
        # Write out to file
        with open(fileloc, 'w') as f:
            f.write(' $indata\n')

            # for i in self.indata:
            for i in input_list:
                if isinstance(self.indata[i], str):
                    f.write('  ' + i.ljust(nchars) + "= '" + str(self.indata[i]) + "'\n")
                    
                elif str(self.indata[i])[0] == '[':
                    if str(self.indata[i])[:2] == '[[':
                        # Multi-dimmensional arrays are transposed in FORTRAN vs Python

                        for j in range(np.size(self.indata[i],0)):
                            var_str = str(i) + '(1:' + str(len(self.indata[i][j,:])) + ',' + \
                                      str(j+1) + ')'
                            val_str = str(self.indata[i][j,:])[1:-1]  # get rid of brackets
    
                            # add zeros after decimals if they aren't there
                            val_str_list = list(val_str)
                            pythonic_decimals = [j for j, c in enumerate(val_str_list) 
                                                 if val_str[j: j + 2] == '. ']
                            if val_str[-1] == '.':
                                pythonic_decimals.append(len(val_str) - 1)
                            for k, pd in enumerate(pythonic_decimals):
                                val_str_list[pd + k + 1: pd + k + 1] = '0'
                                # val_str_list[pd + 1] = '0'
                                # This breaks at the end points, and might
                                # not work if there's only 1 space
                            val_str_out = "".join(val_str_list)

                            # Finally write it all out
                            f.write('  ' + var_str.ljust(nchars) + '= ' + val_str_out + '\n')
                        
                    else:
                        # Make sure indeces are applied correctly
                        var_str = str(i) + '(1:' + str(len(self.indata[i])) + ')'
                        val_str = str(self.indata[i])[2:-1]  # get rid of brackets
                        
                        # add zeros after decimals if they aren't there
                        val_str_list = list(val_str)
                        pythonic_decimals = [j for j, c in enumerate(val_str_list)
                                             if val_str[j : j + 2] == '. ']
                        if val_str[-1] == '.':
                            pythonic_decimals.append(len(val_str) - 1)
                        for k, pd in enumerate(pythonic_decimals):
                            val_str_list[pd+k+1 : pd+k+1] = '0'
                            # val_str_list[pd + 1] = '0'
                            # This breaks at the end points, and might
                            # not work if there's only 1 space
                        val_str_out = "".join(val_str_list)
                        
                        # Finally write it all out
                        f.write('  ' + var_str.ljust(nchars) + '= ' + val_str_out + '\n')
                    
                else:
                    f.write('  ' + i.ljust(nchars) + '= ' + str(self.indata[i]) + '\n')
                
            f.write('/\n')
            
        if verbose:
            newpath, _ = path.split(fileloc)
            if newpath == '':
                newpath = path.abspath('.')
                if newpath[-1] != '/':
                    newpath += '/'
            else:
                newpath = ''
                
            print('File written to ' + newpath + fileloc)

# ----------------------------------------------------------------------
