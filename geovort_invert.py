def geovort_invert(omega,iterations,threshold,num_lats,num_lons,f0, ds, msfm, qgpv,BC):

    import numpy as np
    

#BC_flag = 0 for Dirichlet and 1 for Neumann
    res=np.zeros([19,num_lats,num_lons])
    phip = BC
    

# define coefficients

    mfs=np.zeros([num_lats,num_lons])
    mfs[:] = f0*(ds/(msfm[:]))**2.

    AI_1 = 1.
    AI_2 = 1.
    AI_4 = 1.
    AI_5 = 1.
    AI_3v = -4.

#####
    error = 0.0
    for k  in np.arange(1,19):
        for num_iter in np.arange(iterations):
            for j in np.arange(1,num_lats-1):
                for i in np.arange(1,num_lons-1):
                    RES = AI_1*phip[k,j-1,i] + \
                    AI_2*phip[k,j,i-1] + \
                    AI_3v*phip[k,j,i] + \
                    AI_4*phip[k,j,i+1] + \
                    AI_5*phip[k,j+1,i] - \
                    f0*qgpv[k,j,i]*(ds/msfm[j,i])**2.
                    phip[k,j,i] = phip[k,j,i] - omega*RES/AI_3v
                    res[k,j,i]=RES
                    error += np.abs(RES)
#                if(error/(num_levs*num_lats*num_lons)<1.e-10 and num_iter>0):
        if((np.amax(res[:,])<threshold) and (num_iter>0)): 
            print('stopping at iteration number: ', num_iter)
            print('error, num_levs*num_lats*num_lons: ', num_iter,num_levs*num_lats*num_lons)
            break
    
    return phip,res,num_iter
