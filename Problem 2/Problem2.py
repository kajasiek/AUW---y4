import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


'''


PROGRAM IS DONE

Only one thing left to add, but waiting for confirmation.

Everything should be working just fine, no problems whatsoever.

Some values can be changed, will mark them by commenting I guess.

'''

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
v = np.arange(3,22,0.1,dtype=float) #this is in logarithmic, remember to convert lol
#leave dtype - required, otherwise overflows
#third number can be reduced/increased to improve/lower the precision
#IMPORTANT - v.shape and v1.shape (below) need to be the same - cant be bothered to fix

'''
Here starts endless pain and suffering, change my mind
'''


def Temp(r): #Defining temperature at a given distance
    const = ((G * MBH)/(8 * np.pi * sigmab) * Maccedd)**0.25
    return const * 1/r**(0.75)


def bbr(r, v): #defining radiation at given distance and given frequency
    black_body = (2*h*v**3)/c**2 * 1/(np.exp((h*v)/(k*Temp(r))) - 1)
    return 2*np.pi*r*black_body



#Spectrum time - normal



Lum = np.zeros(v.shape)
for freq in enumerate(v): #Calculating luminosity - normal spectrum from the disk
    throwaway = quad(bbr, rin, rout, args=(10**freq[1]))
    Lum[freq[0]] = throwaway[0]






#Spectrum time - bremsstrahlung

#Constants/given
Tcor = 3.481 * 10**9 #in K <- corona
tau = np.array([0.01, 0.1, 1]) #optical depth
Z = 1 #charge
gff = 1.2 #dont question this

#Spectrum
Brems = np.zeros(v.shape)
ne = np.zeros(tau.shape)
plt.figure("Final")
plt.xlim(left = 15, right = 25) #Values can be changed, but those are just fine. THIS NEEDS TO BE HERE - otherwise matplotlib will automatically asign the values and it just breaks the graph
for t in enumerate(tau):
    #print(t[0]) <- for testing
    ne[t[0]] = t[1]/(sigmat*rout) #Electron density // Bremsstrahlung assumed at Rout - confirmed to be fine // You may change the distance, but up to you - expect different result in such case
    for freq in enumerate(v):
        Brems[freq[0]] = 6.8 * 10**-38 * Z**2 * ne[t[0]]**2 * Tcor**-0.5 * gff * np.exp(-(h*10**freq[1])/(k*Tcor)) 
    FinalFlux = Lum + Brems * 4/3 * np.pi * rout**3 #Assuming rout again
    plt.scatter(v, np.log10(FinalFlux), s=1) #Plotting luminosity from disk and bremsstrahlung
    




'''
Oh God, help
'''


#Compton scattering
#Assuming corona is electron based, cause why not - confirmed
#Constants/given
me = 9.1094*10**-28 #g

#Defining quad functions here
def photflux(number): #returns flux at Rout
    return Lum[number]/(4*np.pi*rout**2) 

def nphot(en, number): #returns number of photons of given energy (v) at given distance (I assumed rout - can be changed in photflux)
    return photflux(number)/(h*en) * 1/c 

def electdist(y): #Electron energy distribution
    #return k * Tcor <- This for testing
    return 4*np.pi*c**3 * (me/(2*np.pi*k*Tcor))**1.5 * np.sqrt(1 - 1/y**2) * 1/y**3 * np.exp(-(y*me*c**2)/(k*Tcor)) #This function is fine. Derived this manually, so no mistake (hopefully - tbf emissivity looks just fine, so)

def x(en, en1, y): #definition of x in Fc(x) - lecture
    return en1/(4*y**2 * en) 

def Fc(en, en1, y): #Fc(x) in gamma integral
    if(x(en, en1, y) > 1):
        return 0 #DO NOT REMOVE - IMPORTANT
    else:
        return 1 + x(en, en1, y) - 2*x(en, en1, y)**2 + 2*x(en, en1, y)*np.log(x(en, en1, y))


#Quads themselves
def quadoverv(en, en1, number): #frequency/energy integral
    throwawaygamma = quad(quadovergamma, 1, 3, args=(en,en1))
    return nphot(en, number)/en * throwawaygamma[0]

def quadovergamma(y, en, en1): #gamma integral
    return Fc(en, en1, y) * electdist(y)/y**2


v1 = np.arange(8, 27, 0.1, dtype=float) #epsilon1 from lecture for those interested // REALLY IMPORTANT - v.shape and v1.shape must be the same <- otherwise will fail quite badly // 
#no clue how to fix right now, maybe someone else can? 

#Values in v1 and v can be changed to whatever, but remember to leave the v.shape and v1.shape the same - unless you figure out how to fix it - I cant be bothered

#Compton scattering - spectrum
emissivity = np.zeros(v1.shape) #This is what we wanna end up with
for freq in enumerate(v1):
    throwawayv = quad(quadoverv, 10**3, 10**22, args=(10**freq[1], freq[0]))
    emissivity[freq[0]] = throwawayv[0] * (3 * c * sigmat * 10**freq[1])/(16*np.pi) 



plt.scatter(v1, np.log10(emissivity), s=1)

'''
To do list

1. Somehow convert emissivity to luminosity. Needs confirmation tho, so gonna wait

Otherwise everything else is just fine
'''


plt.show()








