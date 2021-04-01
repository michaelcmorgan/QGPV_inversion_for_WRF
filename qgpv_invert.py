def qgpv_invert(omega,iterations,threshold,p_levels,num_levs,num_lats,num_lons,f0, ds, dp,msfm, SS,qgpv, Tb, Tt, BC, BC_flag):

    import numpy as np
    
    R = 287.04
    g = 9.8066

#BC_flag = 0 for Dirichlet and 1 for Neumann
    res=np.zeros([num_levs,num_lats,num_lons])
    phip = BC
    

# define coefficients

    mfs=np.zeros([num_lats,num_lons])
    mfs[:] = ((f0*ds)/(msfm[:]*dp))**2.

    AI_1 = 1.
    AI_2 = 1.
    AI_4 = 1.
    AI_5 = 1.
    AI_3v = -4.

    AI_3s=np.zeros([num_levs,num_lats,num_lons])
    AI_6=np.zeros([num_levs,num_lats,num_lons])
    AI_7=np.zeros([num_levs,num_lats,num_lons])

    for k in np.arange(1,num_levs-1):
        AI_6[k,:]  = mfs[:]/SS[k-1]
        AI_7[k,:]  = mfs[:]/SS[k] 
        AI_3s[k,:] = -(AI_6[k,:]+AI_7[k,:])
#####
    error = 0.0
    for num_iter in np.arange(iterations):
        for k in np.arange(1,num_levs-1):
            for j in np.arange(1,num_lats-1):
                for i in np.arange(1,num_lons-1):
                    RES = AI_1*phip[k,j-1,i] + \
                    AI_2*phip[k,j,i-1] + \
                    (AI_3v+AI_3s[k,j,i])*phip[k,j,i] + \
                    AI_4*phip[k,j,i+1] + \
                    AI_5*phip[k,j+1,i] + \
                    AI_6[k,j,i]*phip[k-1,j,i] + \
                    AI_7[k,j,i]*phip[k+1,j,i] - f0*qgpv[k,j,i]*(ds/msfm[j,i])**2.
                
                    phip[k,j,i] = phip[k,j,i] - omega*RES/(AI_3v+AI_3s[k,j,i])
                    res[k,j,i]=RES
                    error += np.abs(RES)
                    if(BC_flag==1):
                        if(k==1): 
                            phip[0,j,i] = phip[1,j,i] - (R * Tb[j,i]/(100*p_levels[0]))*dp
                        if(k==num_levs-2): 
                            phip[-1,j,i] = phip[-2,j,i] + (R * Tt[j,i]/(100.*p_levels[-1]))*dp
#        if(error/(num_levs*num_lats*num_lons)<1.e-10 and num_iter>0):
        if((np.amax(res)<threshold) and (num_iter>1)): 
            print('stopping at iteration number: ', num_iter)
            print('error, num_levs*num_lats*num_lons: ', num_iter,num_levs*num_lats*num_lons)
            break
    
    return phip,res,num_iter
