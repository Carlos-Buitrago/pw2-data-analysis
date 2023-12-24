import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

# obtain the l and corresponding Dl values using the map for CMB
MyMask = hp.read_map("./mask.fits")
MyMapCMB = hp.read_map("./Sc_CMB.fits")
cl_CMB = hp.anafast(MyMapCMB*MyMask)
l_CMB = np.arange(np.size(cl_CMB))
Dl_CMB = ((l_CMB*(l_CMB+1))/(2*np.pi))*cl_CMB 
# plot the Dl vs l
plt.plot(l_CMB,Dl_CMB)
plt.xlabel("l")
plt.ylabel("Dl")
plt.savefig("./angular_power_spectrum_CMB.jpg")
plt.close()

# obtain the l and corresponding Dl values using the map for the lowest freq
MyMap100 = hp.read_map("MyMap_NoDP_Smooth_Degraded_100.fits")
cl_100 = hp.anafast(MyMap100*MyMask)
l_100 = np.arange(np.size(cl_100))
Dl_100 = ((l_100*(l_100+1))/(2*np.pi))*cl_100 
# plot the Dl vs l
plt.plot(l_100,Dl_100)
plt.xlabel("l")
plt.ylabel("Dl")
plt.savefig("./angular_power_spectrum_Map100.jpg")

