import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

# obtain the l and corresponding Dl values using the map for CMB using the mask of highest freq
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

# obtain the l and corresponding Dl values using the map for CMB using the mask of lowest freq
MyMask_low = hp.read_map("./mask_lowestfreq.fits")
MyMapCMB = hp.read_map("./Sc_CMB.fits")
cl_CMB_low = hp.anafast(MyMapCMB*MyMask_low)
l_CMB_low = np.arange(np.size(cl_CMB_low))
Dl_CMB_low = ((l_CMB_low*(l_CMB_low+1))/(2*np.pi))*cl_CMB_low
# plot the Dl vs l
plt.plot(l_CMB_low,Dl_CMB_low)
plt.xlabel("l")
plt.ylabel("Dl")
plt.savefig("./angular_power_spectrum_CMB_lowestfreqmask.jpg")
plt.close()


# obtain the l and corresponding Dl values using the map for the lowest freq using the mask of highest freq
MyMap100 = hp.read_map("MyMap_NoDP_Smooth_Degraded_100.fits")
cl_100 = hp.anafast(MyMap100*MyMask)
l_100 = np.arange(np.size(cl_100))
Dl_100 = ((l_100*(l_100+1))/(2*np.pi))*cl_100 
# plot the Dl vs l
plt.plot(l_100,Dl_100)
plt.xlabel("l")
plt.ylabel("Dl")
plt.savefig("./angular_power_spectrum_Map100.jpg")
plt.close()

# obtain the l and corresponding Dl values using the map for the lowest freq using the mask of lowest freq
MyMap100 = hp.read_map("MyMap_NoDP_Smooth_Degraded_100.fits")
cl_100_low = hp.anafast(MyMap100*MyMask_low)
l_100_low = np.arange(np.size(cl_100_low))
Dl_100_low = ((l_100_low*(l_100_low+1))/(2*np.pi))*cl_100_low 
# plot the Dl vs l
plt.plot(l_100_low,Dl_100_low)
plt.xlabel("l")
plt.ylabel("Dl")
plt.savefig("./angular_power_spectrum_Map100_lowestfreqmask.jpg")

