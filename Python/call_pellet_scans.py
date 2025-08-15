import numpy as np
import pellet as pel


base_dir = "/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/"
nmlfile = base_dir + "nml_pellet.dat"
gfile = base_dir + "st50010dn.eqdsk"
xfile = base_dir + "build/xpellet"

pel_old = 1.5
r_pel_input = [1.0,1.5,2.0]
r_pel_input_scaled = np.divide(r_pel_input,pel_old)
mass_scale = r_pel_input_scaled**3
d_input = np.multiply(r_pel_input,2.0)
v_input = [1000,1500]

print(d_input)

pel.scan_pellet_mass_velocity(base_dir, device = None, runid_prefix = None,
                              diameters_mm = d_input,
                              velocities = v_input,
                              pellet_executable_loc = base_dir + '/build/xpellet')