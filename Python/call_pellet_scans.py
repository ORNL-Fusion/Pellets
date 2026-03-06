#%%
import numpy as np
import pellet as pel
import calc_angle as ca

#%%

base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
nmlfile = base_dir + "nml_pellet.dat"
# gfile = base_dir + "RUN3_15100089.geqdsk"
xfile = base_dir + "build/xpellet"

#%%

entry = [3.4,3.9]
exit = [5.79,-0.00]
shift = 18.5
[orig_angle,shift_plus,shift_minus] = ca.calc_angle(entry,exit,shift)

print(orig_angle,shift_plus,shift_minus)

#%%

pel_old = 1.5
r_pel_input = [2.5,3.0,3.5]
r_pel_input_scaled = np.divide(r_pel_input,pel_old)
mass_scale = r_pel_input_scaled**3
d_input = np.multiply(r_pel_input,2.0)
v_input = [250,375,500,625,750,875,1000,1125,1250,1375,1500]
# v_input = [200,300]

den0 = [0.40,0.60,0.80,1.0,1.2,1.4]
den0_scaled = np.multiply(den0,1.0e19)
den1 = np.multiply(den0_scaled,0.25)

pel.scan_pellet_mass_velocity(base_dir, device = None, runid_prefix = None,
                              diameters_mm = d_input,
                              velocities = v_input,
                              pellet_executable_loc = base_dir + '/build/xpellet')

# pel.scan_density_profile(base_dir,den0 = den0_scaled,
#                          den1=den1,pellet_executable_loc=xfile)