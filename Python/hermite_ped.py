import numpy as np

def hermite_ped(rho_r,xped,valaxis,valsep,valpedtop=None):
    prof = np.zeros(len(rho_r))
    t = (rho_r - xped) / (np.max(rho_r) - xped)
    s = rho_r/xped
    mo = -25.0
    if valpedtop is not None:
        mt = 0.5*(valpedtop-valaxis)/xped
    else:
        mt = 0.0
        valpedtop = 0.0
    idx = np.where(rho_r >= xped)[0][0]
    H00,H10,H01,H11 = hermite_poly(rho_r,xped,s)
    prof[0:idx] = valaxis*H00[0:idx] + valpedtop*H01[0:idx] + (xped*mt)*H11[0:idx]
    H00,H10,H01,H11 = hermite_poly(rho_r,xped,t)
    prof[idx:] = valpedtop*H00[idx:] + ((1.0 - xped)*mt)*H10[idx:] + valsep*H01[idx:] + ((1.0 - xped)*mo)*H11[idx:]

    return prof

def hermite_poly(rho_r,xped,var):
    H00 = 2.0*var**3 - 3.0*var**2 + 1.0
    H10 = var**3 - 2.0*var**2 + var
    H01 = -1.0*H00 + 1.0
    H11 = var**3 - var**2

    return H00,H10,H01,H11