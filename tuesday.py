import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import astropy.io.fits as fits

plt.ion()

# frequencies used
freq = ["100","143","217","353","545","857"]
nfreq = np.size(freq)

# path to maps
path_maps = "/home/mdouspis/M03_Material/M03_Maps_Planck/"

# dipole effect map
MyDipoleMap = hp.read_map(path_maps+"dipole_nside2048.fits")

# real FWHMs from Planck
infos = fits.getdata("/home/mdouspis/M03_Material/HFI_RIMO_R2.00.fits",2)
real_fwhm = infos[:]["FWHM"]

# when looking at the real FWHM we notice that the one with the worst resolution corresponds to ~10arcmin
# in order to be able to compare all maps we convolve them all to 10arcmin taking into account that
# G(10) = G(real_fhmw)***G(X)
X = np.sqrt(10**2 - real_fwhm**2)

# X in radians
Xrad = X*((np.pi)/(60*180))

# NSIDE value to which the maps will be degraded (512)
final_nside = 512

for i in np.arange(nfreq):
    # read initial map
    MyMap = hp.read_map(path_maps+"HFI_SkyMap_"+freq[i]+"_2048_R4.00_full.fits")
    # remove dipole effect
    MyMapWithoutDipole = MyMap - MyDipoleMap
    # smoothing of the map
    MyMapSmooth = hp.smoothing(MyMapWithoutDipole, fwhm=Xrad[i])
    # degrading of the map
    MyMapDegraded = hp.ud_grade(MyMapSmooth, nside_out = final_nside)
    # write the final map
    hp.write_map("MyMap_NoDP_Smooth_Degraded_"+freq[i]+".fits",MyMapDegraded)
    # draw the maps
    hp.mollview(MyMap,norm="hist",title="Map (Original) " + freq[i] + "GHz")
    plt.savefig("./original_maps/original_"+freq[i]+".jpg")
    hp.mollview(MyMapWithoutDipole,norm="hist",title="Map without dipole eff " + freq[i] + "GHz")
    plt.savefig("./noDP/noDP_"+freq[i]+".jpg")
    hp.mollview(MyMapSmooth,norm="hist",title="Map no DP and smoothed " + freq[i] + "GHz")
    plt.savefig("./noDP_smooth/noDP_smooth_"+freq[i]+".jpg")
    hp.mollview(MyMapDegraded,norm="hist",title="Map no DP, smoothed and degraded " + freq[i] + "GHz")
    plt.savefig("./noDP_smooth_degraded/noDP_smooth_degraded_"+freq[i]+".jpg")
    



