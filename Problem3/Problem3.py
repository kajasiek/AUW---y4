import matplotlib.pyplot as plt
import numpy as np
from numba import jit

#Solution to problem 3

'''

WARNING 1

Coordinates for the letter must be set manually for now - maybe I will fix it at some point

'''

'''

WARNING 2

Cross section/impact parameter/whateveryoucall it has a minimal value of 1.5 - otherwise it will be under the photon sphere and the photon will be just invisible

'''

#constants
Msun = 2 * 10**30 #kg
M = 1 * Msun #kg - number can be changed to whatever
c = 3e5 #km/s
G = 6.67e-20 #km^3/kg/s^2
r_s = 2 * G * M/c**2 #in km


#Creating the image
def screen_base():
    rad = 2 * r_s
    theta = np.linspace(0, 2*np.pi, 1000)
    plt.figure("Base screen")

    #Plotting a circle, wow
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    plt.plot(x,y)

    #Coordinates - YOU NEED TO MANUALLY FIND THOSE COORDINATES FOR THE LETTER <- maybe will fix later, who knows?
    rad_j = [0.75*r_s, 0.75*r_s, np.sqrt(1.125)*r_s, np.sqrt(1.125)*r_s, np.sqrt(1.125)*r_s, np.sqrt(1.125)*r_s, 1.5*r_s]
    theta_j = [0, np.pi/2, np.pi/4, 7*np.pi/4, 5*np.pi/4, 3*np.pi/4, 1.5*np.pi]

    #Plotting the letter
    x_j = rad_j * np.cos(theta_j)
    y_j = rad_j * np.sin(theta_j)

    plt.scatter(x_j, y_j, c="gray", s=4)


def black_hole():
    photon_sphere = 1.5 * r_s 
    theta = np.linspace(0, 2*np.pi, 1000)
    fig, ax = plt.subplots()
    ax.set_title("Black hole and photon sphere")
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_aspect("equal")
    #Plotting the photon sphere
    x_phot_sph = photon_sphere * np.cos(theta)
    y_phot_sph = photon_sphere * np.sin(theta)
    plt.plot(x_phot_sph, y_phot_sph, c="black")

    #Plotting the black hole
    x_bh = r_s * np.cos(theta) 
    y_bh = r_s * np.sin(theta)
    ax.fill(x_bh, y_bh, "black")


def deviation():
    dev = [np.pi - np.arctan(0.1), np.pi + np.arctan(0.1)]
    return dev

def angular_change(b):
    ang = 2 * r_s * 1/b
    return ang

screen_base()
black_hole()
print(deviation())
print(1.85, 20*r_s*np.tan(angular_change(1.85)), angular_change(1.85))
print(1.95, 20*r_s*np.tan(angular_change(1.95)))
plt.show()