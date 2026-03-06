#%%
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import EFIT.equilParams_class as epc
# import EFITUtils as epc

#%%
# filepath = '/Users/gz6/Documents/pellets/code_dev/code/Pellets/src/'
# filename = 'st50010dn.eqdsk'
filepath = '/Users/gz6/Documents/TokE/pellets/results/k6_NGS_Parks_Q/more_runs/RUN7_15100081/'
filename = 'ST40-15100081-ASTRA-RUN7-24610p00ms.geqdsk'

# filepath = '/Users/gz6/Documents/SMARTs-SciDAC-C1/grid/verification/'
# filename = 'geqdsk_cat_reference'
gfile = filepath + filename
gf = epc.equilParams(gfile)

#%%

psiRZn = gf.g['psiRZn']
gR = gf.g['R']
gZ = gf.g['Z']

wall = gf.g['wall']

# print(wall)

# fig = plt.figure()
# for w in range(len(wall[:, 0]) - 1):
#     plt.plot([wall[w, 0], wall[w + 1, 0]], [wall[w, 1], wall[w + 1, 1]], 'k')


# contours=np.linspace(0.0,1.0,11)
# if contours is not None:
#     cont = plt.contour(gR, gZ, psiRZn, contours, cmap='plasma', linewidths = [1.0])
# plt.contour(gR, gZ, psiRZn, [1], colors = 'k', linewidths = [2.0])
# plt.clabel(cont,cont.levels)
# # plt.savefig(filepath+"DIII-D_equil.png",transparent=True)
# plt.show()

psiRZ = gf.g['psiRZ']
Bp_2D = gf.Bp_2D
# print(np.shape(Bp_2D))
# print(np.where(gR<=5.47732)[0][-1])
# print(np.where(gZ<=0.0)[0][-1])
# # print(gf.g['lcfs'])
# print(np.where(gf.g['lcfs'][:,0]==np.max(gf.g['lcfs'][:,0])))
# print(np.where(gf.g['lcfs'][:,1]==np.min(abs(gf.g['lcfs'][:,1]))))
# print(np.shape(gf.g['lcfs']))
# print(gR)
# rows, cols = np.where(gR<=5.47732)
# print(rows,cols)
# idr = np.where(gR<=5.47732)[0][-1]
# idz = 100
# rows = np.where(Bp_2D[:,idz]==np.amax(Bp_2D[:,idz]))[0][0]
# idr_max = rows
# idz_max = idz
# print(np.sqrt((Bp_2D/np.max(Bp_2D))))
# Bp_check = np.sqrt(Bp_2D)*np.amax(Bp_2D)
# Bp_check = (np.sqrt((Bp_2D - Bp_2D[151,99]))/(Bp_2D[171,100] - Bp_2D[151,99]))
# Bp_check = ((Bp_2D - (Bp_2D[idr,idz]))/((Bp_2D[idr_max,idz_max] - (Bp_2D[idr,idz]))))
Bp_check = np.sqrt((psiRZ)/((gf.siBry)))
print(gf.getTorPsi())
# Bp_check = np.sqrt(psiRZn)
print(Bp_2D)
print((Bp_check))
# Bp_check = np.sqrt((Bp_2D - Bp_2D[151,99])/(np.amax(Bp_2D)))
# print(np.shape(Bp_check))
# Bt = gf.g['Bt_2D']
# psiTor = np.sqrt(abs((psiRZ - psiRZ[idr,idz])/(psiRZ[idr_max,idz_max] - psiRZ[idr,idz])))
# psiTor,psiTor = np.meshgrid(Bp_check,Bp_check)
X, Y = np.meshgrid(gR,gZ)
# psiTor = np.tile(Bp_check,(200,1))
# print(psiTor)
print(Bp_check)


#%%

fig, ax = plt.subplots(layout="constrained")
contours=np.linspace(0.0,1.0,21)
if contours is not None:
    cont = ax.contour(X, Y, psiRZn, contours, cmap='plasma', linewidths = [1.0])
ax.contour(gR, gZ, psiRZn, [1], colors = 'k', linewidths = [2.0])
ax.plot(5.7933,0.0,'kx')
# plt.contour(gf.g['lcfs'][:,0], gf.g['lcfs'][:,1], psiTor, [1], colors = 'k', linewidths = [2.0])
ax.clabel(cont,cont.levels)
plt.show()

# print("Rmax: ", max(gR))
# print("Rmin: ", min(gR))
# print(dir(gf.g.values))
# print(gf.g.keys)
# print(gf.g)

# gf.plotProfile()
# plt.show()

# print((gf.g['lcfs'][:,0]))
print(gf.help())
# print(np.shape(psiRZ))
# %%
