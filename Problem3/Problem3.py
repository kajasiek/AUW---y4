import matplotlib.pyplot as plt
import numpy as np
from numba import jit
from scipy.integrate import quad

#Solution to problem 3

'''

WARNING 1

Coordinates for the letter must be set manually for now - maybe I will fix it at some point

'''

'''

WARNING 2

Cross section/impact parameter/whateveryoucall it has a minimal value of 3*sqrt(3) - otherwise it will be under the photon sphere and the photon will be just invisible

'''

#constants
'''
Msun = 2 * 10**30 #kg
M = 1 * Msun #kg - number can be changed to whatever
c = 3e5 #km/s                                                     Is this even needed? Think not
G = 6.67e-20 #km^3/kg/s^2
r_s = 2 * G * M/c**2 #in km
'''
bmin = 3*np.sqrt(3) #minimal impact parameter
distance_obs = 200 #distance to observer
distance_scr = 20
#Creating the image
def screen_base():
    rad = 2
    theta = np.linspace(0, 2*np.pi, 1000)
    plt.figure("Base screen")

    #Plotting a circle, wow
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    plt.plot(x,y) #Assuming its circular, cause why not?

    #Coordinates - YOU NEED TO MANUALLY FIND THOSE COORDINATES FOR THE LETTER <- maybe will fix later, who knows?
    rad_j = [0.75, 0.75, np.sqrt(1.125), np.sqrt(1.125), np.sqrt(1.125), np.sqrt(1.125), 1.5]
    theta_j = [0, np.pi/2, np.pi/4, 7*np.pi/4, 5*np.pi/4, 3*np.pi/4, 1.5*np.pi]

    #Plotting the letter
    x_j = rad_j * np.cos(theta_j)
    y_j = rad_j * np.sin(theta_j)

    plt.scatter(x_j, y_j, c="gray", s=4)


def black_hole():
    photon_sphere = 1.5
    theta = np.linspace(0, 2*np.pi, 1000)
    fig, ax = plt.subplots()
    ax.set_title("Black hole and photon sphere")
    ax.set_xlim(-50, 50)
    ax.set_ylim(-50, 50)
    ax.set_aspect("equal")
    #Plotting the photon sphere
    x_phot_sph = photon_sphere * np.cos(theta)
    y_phot_sph = photon_sphere * np.sin(theta)
    plt.plot(x_phot_sph, y_phot_sph, c="black")

    #Plotting the black hole
    x_bh =  np.cos(theta) 
    y_bh =  np.sin(theta)
    ax.fill(x_bh, y_bh, "black")


def angle_deviation(): #angle range - the final angle needs to be in this range, otherwise the photon will miss the screen
    dev = [np.pi - np.arctan(0.1), np.pi + np.arctan(0.1)]
    return dev


screen_base()
black_hole()


def rmin_to_b(r):
    return r * np.sqrt(r/(r-1))


def dphi(dr,r,b):
    return dr/(r**2 * np.sqrt(1/b**2 - (1 - 1/r) * 1/r**2))

rmin = 5 #in r_s
N = 10000 #number of steps
dr = (distance_obs - rmin)/N  
r = np.zeros(N)
phi = np.zeros(N)
phi[0] = 0
r[0] = rmin 
b = rmin_to_b(rmin)
for i in range(N-1): 
    r[i+1] = r[i] + dr
    phi[i+1] = phi[i] + dphi(dr, r[i+1], b)
    x = r[i] * np.cos(phi[i])
    y = r[i] * np.sin(phi[i]) 
    plt.scatter(x, y, c='r',s=2)







plt.show()