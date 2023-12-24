import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

plt.ion()

# frequency vector
freq=["100", "143", "217", "353", "545", "857" ]
nfreq=np.size(freq)

# read each map and add to T (total emission observed at a given frequency)
T=[]
for i in np.arange(nfreq):
    MyMap = hp.read_map("./MyMap_NoDP_Smooth_Degraded_"+freq[i]+".fits")
    T.append(MyMap)

# we want to compute Sc for the CMB given in Hurier et al eq.12

# CT is the covariance matrix of T, CTinv is the inverse of CT
CT = np.cov(T)
CTinv = np.linalg.inv(CT)

# For cmb fc is a column vector of ones col1
col1 = np.ones((6,1))

# calculate Sc
w = (col1.T @ CTinv)/(col1.T @ CTinv @ col1)
Sc = np.dot(w,T)[0]

# write map
hp.write_map("Sc_CMB.fits",Sc,overwrite=True)

# draw and save CMB map
hp.mollview(Sc, norm="hist", title ="Sc for CMB")
plt.savefig("./Sc_CMB.jpg")

# create a mask using the map of max freq
tresh = 0.01
mapmaxfreq = hp.read_map("./MyMap_NoDP_Smooth_Degraded_857.fits")
std_mapmaxfreq = np.std(mapmaxfreq)
mask = mapmaxfreq < tresh*std_mapmaxfreq

hp.mollview(mask, norm="hist", title ="Mask with map of max freq")
plt.savefig("./mask.jpg")

# write the map of the mask
hp.write_map("mask.fits",mask,overwrite=True)

# create the masked map
masked_map =  hp.ma(Sc)
masked_map.mask = np.logical_not(mask)

hp.mollview(masked_map, norm="hist", title ="Masked map")
plt.savefig("./masked_map.jpg")
plt.close()

# write the masked map
hp.write_map("masked_map.fits", masked_map,overwrite=True)

