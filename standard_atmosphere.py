def standard_atmosphere(p, p_half):
    import numpy as np
    g=9.8066
    kappa=0.2856219
    alpha=-1./5.255877
    beta=-6341.624
    gamma=.0065
    R=287.04
    cp = 1004.7
    p00=1000.
    nlevels = p.size
    nlevels_half = p_half.size

    Z = np.zeros(nlevels)
    phi = np.zeros(nlevels)
    T = np.zeros(nlevels)
    theta = np.zeros(nlevels)

    T_half = np.zeros(nlevels_half)
    theta_half = np.zeros(nlevels_half)
    dthetadp = np.zeros(nlevels_half)
    S = np.zeros(nlevels_half)

    for k in np.arange(nlevels):
        if (p[k] > 226.32):
            Z[k] = (288.15/gamma) * ( 1. - (1013.25/p[k])**alpha)
            T[k] = 288.15 - gamma*Z[k]
            theta[k] = T[k] * (p00/p[k])**(R/cp)
        else:
            T[k] = 216.65
            Z[k] = (11.e3+beta*np.log(p[k]/226.32))
            theta[k] = T[k] * (p00/p[k])**(R/cp)
                
    for k in np.arange(nlevels_half):
        if (p_half[k] > 226.32):
            T_half[k] = 288.15 - gamma*Z[k]
            theta_half[k] = T_half[k] * (p00/p_half[k])**(R/cp)
        else:
            T_half[k] = 216.65
            theta_half[k] = T_half[k] * (p00/p_half[k])**(R/cp)
    
    for k in np.arange(nlevels_half):
        dthetadp[k] = (theta[k+1] - theta[k])/(1.e2*(p[k+1]-p[k]))
        S[k] =  -R*(T_half[k]/theta_half[k])*dthetadp[k]/(p_half[k]*100.)
    return Z, g*Z, T, S
