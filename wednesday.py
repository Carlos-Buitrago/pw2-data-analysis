import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.constants as cst

# frequency vector

nu = np.arange(9,1001)
nuG = nu*1e9

# CMB
Acmb = 70 
h = cst.h
kB = cst.k
Tcmb = 2.7255
x = (h*nuG)/(kB*Tcmb)
g = ((np.exp(x) - 1)**2)/((x**2)*np.exp(x))
Scmb = Acmb/g

# synchrotron
As = 30*1e2 
alpha = 1 
nu0 = 4.08*10e8 
ss = As*(((nu0)/(nuG))**3.)

# free-free
EM = 1 
Te = 7000 
T4 = Te/(1e4)
nu9 = nuG/(1e9) 
gff = np.log(np.exp((5.960 - (np.sqrt(3)/math.pi)*np.log(nu9*T4**(-3/2))))+math.e)
tau = 0.05468*(Te**(-3/2))*(nu9**(-2.))*gff
Sff = (10**6)*Te*(1-np.exp(-tau))

# thermal dust
Ad = 100
Betad = 1.55 
Td = 23 
nu0p = 545*1e9
gamma = h/(kB*Td)
Sd = Ad*((nuG/nu0p)**(Betad+1))*((np.exp(gamma*nu0p) - 1)/(np.exp(gamma*nuG) - 1))

# final plot
plt.loglog(nu,Scmb, label = "CMB")
plt.loglog(nu,ss, label = "SS")
plt.loglog(nu,Sff, label = "FF")
plt.loglog(nu,Sd, label = "TD")
plt.legend()
plt.savefig("wednesday_plot.jpg")
plt.show()

########################################

#FWHM for the selected frequencies (THEORETICAL since real telescope has more than one mirror)
freqs = np.array([30,44,70,100,143,217,353,545,857])
D = 1.5 # mirror diameter
fwhm = []
for i in np.arange(np.size(freqs)):
    val = 1.22*cst.c*(1/(freqs[i]*1e9))*(1/D)*((180/(np.pi))*60)
    fwhm.append(val)
print(fwhm)

