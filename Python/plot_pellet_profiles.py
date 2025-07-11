import numpy as np
import matplotlib.pyplot as plt
import pellet as pel
import PelletSumfile as PS

sumfile_loc = '/Users/gz6/Documents/pellets/code/Pellets/src/build/sum_pelletDIII-D_geqdsk_test_      .dat'

pelsum = PS.PelletSumfile(sumfile_loc=sumfile_loc,gfile="/Users/gz6/Documents/pellets/code/Pellets/src/g200201.00000")
pelsum.plot_density_change()
