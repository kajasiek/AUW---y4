import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


#Constants/Given
Msun = 2 * 10**33 #in g
MBH = 10 * Msun #in g
G = 6.674 * 10**-8 #cgs
c = 3 * 10**10 #cm/s
rin = 20 * G * MBH/c**2 #cgs
rout = 1000 * G * MBH/c**2 #cgs
mp = 1.673 * 10**-24 #cgs
sigmat = 6.652 * 10**-25 #cgs - Thomson
sigmab = 5.6704 * 10**-5 #cgs - SB
eta = 0.01 #accretion efficiency
h = 6.6261 * 10**-27 #cgs <- normal PC
hred = 1.0546 * 10**-27 #cgs <- reduced PC
k = 1.3087 * 10**-16 #cgs
#Eddington shit <- accretion rate, might be useful later lol
Ledd = (4 * np.pi * G * MBH * mp * c)/sigmat
Maccedd = Ledd/(eta * c**2)
const = ((G * MBH)/(8 * np.pi * sigmab) * Maccedd)**0.25 #testing purposes
#Frequency
v = np.arange(15,22,0.01,dtype=float) #this is in logarithmic, remember to convert lol
#leave dtype - required, otherwise overflows
#third number can be reduced/increased to improve/lower the precision
'''


Here starts endless pain and suffering, change my mind


'''
def Temp(r):
    const = ((G * MBH)/(8 * np.pi * sigmab) * Maccedd)**0.25
    return const * 1/r**(0.75)


def bbr(r, v):
    black_body = (2*h*v**3)/c**2 * 1/(np.exp((h*v)/(k*Temp(r))) - 1)
    return 2*np.pi*r*black_body



#Spectrum time - normal



F = np.zeros(v.shape)
for freq in enumerate(v):
    throwaway = quad(bbr, rin, rout, args=(10**freq[1]))
    F[freq[0]] = throwaway[0]

    #print(F[freq[0]]) <- testing lol




#Spectrum time - bremsstrahlung
#Constants/given
Tcor = 3.481 * 10**9 #in K <- corona
tau = np.array([0.01, 0.1, 1])
Z = 1 #charge
gff = 1.2 #dont question this

#Spectrum
Brems = np.zeros(v.shape)
ne = np.zeros(tau.shape)
for t in enumerate(tau):
    #print(t[0]) <- for testing
    ne[t[0]] = t[1]/(sigmat*rout) #Electron density // Assuming bremsstrahlung at max radius - needs clarification // seems like a random number, but I think its correct?
    for freq in enumerate(v):
        Brems[freq[0]] = 6.8 * 10**-38 * Z**2 * ne[t[0]]**2 * Tcor**-0.5 * gff * np.exp(-(h*10**freq[1])/(k*Tcor)) 
    FinalFlux = F + Brems * 4/3 * np.pi * rout**3 
    plt.scatter(v, np.log10(FinalFlux), s=1)
    




'''

To do list:
1) Introducte Compton scattering

'''



plt.show()


