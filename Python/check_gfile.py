import numpy as np
import matplotlib.pyplot as plt
import EFIT.equilParams_class as epc

filepath = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/'
filename = 'g200201.00000'
gfile = filepath + filename
gf = epc.equilParams(gfile)

psiRZn = gf.g['psiRZn']
gR = gf.g['R']
gZ = gf.g['Z']

wall = gf.g['wall']

print(wall)

fig = plt.figure()
for w in range(len(wall[:, 0]) - 1):
    plt.plot([wall[w, 0], wall[w + 1, 0]], [wall[w, 1], wall[w + 1, 1]], 'k')


contours=np.linspace(0.0,1.0,15)
if contours is not None:
    cont = plt.contour(gR, gZ, psiRZn, contours, cmap='plasma', linewidths = [1.0])
plt.contour(gR, gZ, psiRZn, [1], colors = 'k', linewidths = [2.0])
# plt.show()

# print("Rmax: ", max(gR))
# print("Rmin: ", min(gR))
# print(dir(gf.g.values))
# print(gf.g.keys)
# print(gf.g)

# gf.plotProfile()
# plt.show()

# print((gf.g['lcfs'][:,0]))
print(gf.help())