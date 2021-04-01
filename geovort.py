def geostrophic_vorticity(num_levs,num_lats,num_lons,msfm,ds,f0,phi_prime):
    import numpy as np

    g=9.8066
    qgpv_geovor = np.zeros([num_levs,num_lats,num_lons])

#    AC_1 = np.zeros([num_lats,num_lons])
#    AC_2 = np.zeros([num_lats,num_lons])
#    AC_3v = np.zeros([num_lats,num_lons])
#    AC_4 = np.zeros([num_lats,num_lons])
#    AC_5 = np.zeros([num_lats,num_lons])

#    for j in np.arange(1,num_lats-1):
#        for i in np.arange(1,num_lons-1): 
    AC_1 =    1.
    AC_2 =    1.
    AC_4 =    1.
    AC_5 =    1.
    AC_3v =  -4.


    for j in np.arange(1,num_lats-1):
        for i in np.arange(1,num_lons-1):
            qgpv_geovor[:,j,i]=(msfm[j,i]**2.)/(f0*ds*ds)*(AC_1*phi_prime[:,j-1,i]+AC_2*phi_prime[:,j,i-1]+ \
            AC_3v*phi_prime[:,j,i]+AC_4*phi_prime[:,j,i+1]+AC_5*phi_prime[:,j+1,i])

    return qgpv_geovor
