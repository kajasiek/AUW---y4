import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numba import jit
import warnings
#There are a lot of problems with this program
#But to make it look nice, all will be suppressed
warnings.filterwarnings("ignore")


'''

This program was remade by someone else. I do not take any responsibility for this

'''
'''
Improvement of original
'''

#Constants/Given
Msun = 2e33             #[g]
MBH = 10*Msun           #[g]
G = 6.674e-8            #cgs
c = 3e10                #cgs
rin = 20*G*MBH/c**2     #cgs
rout = 1e3*G*MBH/c**2   #cgs
sigmat = 6.652e-25      #cgs - Thomson
sigmab = 5.6704e-5      #cgs - SB
eta = 0.01              #accretion efficiency
h = 6.6261e-27          #cgs - Normal PC
hred = 1.0546e-27       #cgs - Reduced PC
k = 1.3087e-16          #cgs
Tcor = 3.481e9          #[K] - Corona
tau = np.array([0.01,0.1,1.])   #Optical thickness
Z = 1                   #Atomic number
gff = 1.2               #Gaunt factor
me = 9.1094e-28         #[g] - Electron mass
mp = 1.673e-24          #[g] - Proton mass

#Eddington values (accretion rate)
Ledd = (4*np.pi*G*MBH*mp*c)/sigmat
Maccedd  = Ledd/(eta*c**2)

#Frequency range (Log scale)
minv = 3                #Minimum frequency
maxv = 22               #Maximum frequency
#vres = 0.1              #Frequency  resolution
vnum = 200              #Number of frequencies
#v = np.arange(minv,maxv,vres)
v = np.linspace(minv,maxv,vnum)
#Frequency range for Comptonization:
minvc = 8
maxvc = 27
vc = np.linspace(minvc,maxvc,vnum)
#leave dtype - required, otherwise overflows // You should not.


@jit
def T(r):
    '''
Temperature [K] at given distance
r - distance in [cm]
'''
    return (((G*MBH)/(8*np.pi*sigmab)*Maccedd)**0.25)/r**0.75

@jit
def BlackBodyRadiation(r,v):
    '''
Blackbody radiation at given frequency and distance
r - distance in [cm]
v - frequency in [Hz]
'''
    return (4*np.pi*r*h*v**3)/c**2 * 1/(np.exp((h*v)/(k*T(r)))-1)

@jit
def BlackBodyLuminosity(v):
    '''
Observed black body luminosity
v - log frequency in [10^Hz]
'''
    return quad(BlackBodyRadiation,rin,rout,args=(10**v))

@jit
def Bremsstrahlung(v,ne):
    '''
Emissions from Bremsstrahlung in [erg/s/cm^3/Hz]
v - log frequency in [10^Hz]
ne - electron density [cgs]
'''
    rad = 6.8e-38* Z**2 * ne**2 * Tcor**-.5 * gff * np.exp(-h*10**v/(k*Tcor))
    return rad

'''
for taus in bsLum:
    tempFlux = bbLum+taus
    plt.scatter(v,np.log10(tempFlux),s=1)
    plt.scatter(v,np.log10(taus),s=1)
    
plt.scatter(v,np.log10(bbLum),s=1)
plt.show()
'''

@jit
def photflux(n):
    '''Photon flux from black body radiation'''
    return bbLum[n]/(4*np.pi*rout**2)

@jit
def nphot(v,n):
    '''Number of photons'''
    return photflux(n)/(h*v*c)

@jit
def e_distribution(y): #Electron energy distribution
    #print(4*np.pi*c**3 * (me/(2*np.pi*k*Tcor))**1.5 * np.sqrt(1 - 1/y**2) * 1/y**3 * np.exp(-(y*me*c**2)/(k*Tcor)))
    return (4*np.pi*c**3 * (me/(2*np.pi*k*Tcor))**1.5 * np.sqrt(1 - 1/y**2) * 1/y**3 * np.exp(-(y*me*c**2)/(k*Tcor)))

@jit
def x(v0,v1,gamma):
    '''Definition of x in Fc(x)'''
    return v1/(4*gamma**2 *v0)

@jit
def Fc(v0,v1,gamma):
    '''Fc(x) in the gamma integral
v0 - minimal energy
v1 - maximum energy
gamma - gamma'''
    xvar = x(v0,v1,gamma)
    if (xvar>0) and (xvar<1):
        return 1+xvar-2*xvar**2+2*xvar*np.log(xvar)
    else:
        return 0

@jit
def Gamma_int(gamma,v0,v1):
    '''Gamma integral'''
    return Fc(v0,v1,gamma) * e_distribution(gamma)/gamma**2

@jit
def v_int(v0,v1,n):
    '''Frequency/energy integral
v0 - minimal frequency
v1 - maximal frequency
n - loop number'''
    return nphot(v0,n)/v0 * quad(Gamma_int,1,3,args=(v0,v1))[0]

@jit
def emission(minv,maxv,freq,i,ne):
    temp = quad(v_int,10**minv,10**maxv,args=(10**freq,i))[0]
    temp *= 3*c*sigmat*(10**freq)/(16*np.pi)*ne*h*4/3*np.pi*rout**3
    return temp



#Blackbodyradiation
bbLum = []
for freq in v:
    bbLum.append(BlackBodyLuminosity(freq)[0])
bbLum = np.array(bbLum)

#Bremsstrahlung
bsLum = []
for t in tau:
    ne = t/(sigmat*rout)
    temp = []
    for freq in v:
        temp.append(Bremsstrahlung(freq,ne))
    bsLum.append(temp)
bsLum=np.array(bsLum)
bsLum *=(4/3 *np.pi*rout**3)

#Compton Scattering
csLum = []
for t in tau:
    ne = t/(sigmat*rout)
    temp2 = []
    for i,freq in enumerate(vc):
        temp2.append(emission(minv,maxv,freq,i,ne))
    csLum.append(temp2)
csLum=np.array(csLum)
#plt.figure()
plt.xlim(left=15,right=25)
plt.ylim(-40,23)
plt.scatter(v,np.log10(bbLum),s=1)

#plt.figure()
for lum in bsLum:
    plt.scatter(v,np.log10(lum),s=1)

#plt.figure()
for lum in csLum:
    plt.scatter(vc,np.log10(lum),s=1)
plt.show()


'''Structure of Luminosities:
bbLum - independent of tau --> 1-dim (frequencies)
bsLum - 2-dim (tau,frequencies)
csLum - 2-dim (tau,frequencies)
'''
