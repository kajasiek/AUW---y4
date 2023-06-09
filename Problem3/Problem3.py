import matplotlib.pyplot as plt
import numpy as np
from numba import jit 
import warnings 
warnings.filterwarnings("ignore")

#Numba can be deleted if you dont want to install it - program works surprisingly fast anyway

#Solution to problem 3

'''

WARNING 1

Coordinates for the letter must be set manually for now \\
Unfixable afaik

'''

'''

WARNING 2

rmin must be higher than 1.5 - otherwise the photon will be \\
under the photon sphere and collide with the black hole

'''

'''

WARNING 3

Schwarzschild radius is 1 - which is why it doesnt show up anywhere

All units are in Schwarzschild radiuses

'''

'''

WARNING 4

The more trajectories you wanna plot at the same time \\
(defined by rmin = np.arange in the trajectory part below), \\
the longer the program will run.
Just be aware of this.

For one trajectory, the program runs surprisingly fast.

'''

'''

Parameters:

One loop around the black hole:
rmin = 1.54 (around this number)

No loop around the black hole:
rmin = 5.5 - 7 (around those numbers)

REMEMBER:
rmin (that you can change) is different from b (impact parameter)

Observer (in the infinity) will see the photons at distance b from the black hole

Also the problem is spherically symmetrical - if a photon coming at b will hit the screen at x, then \\
a photon coming at -b will hit the screen at -x (etc etc)

Center of the coordinate system is in the center of the black hole
'''

#constants
distance_obs = 200 #distance to observer - can be changed, but no need
distance_scr = 20 #distance to screen


#Creating the image
@jit
def screen_base():
    rad = 2
    theta = np.linspace(0, 2*np.pi, 1000)
    plt.figure("Base screen")
    

    #Plotting a circle, wow
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    plt.plot(x,y) #Assuming its circular, cause why not?

    #Coordinates - YOU NEED TO MANUALLY FIND THOSE COORDINATES FOR THE LETTER <- unfixable afaik
    #Coordinates: x: (-1 <-> 1); y: (-1 <-> 1 + kolko)
    y_move = 0.5

    #Step 1 - top line
    x_top = np.arange(-1,1,0.1)
    y_top = np.zeros(x_top.shape)
    for i,j in enumerate(x_top):
        y_top[i] = 1
    
    #Step 2 - side line
    y_side = np.arange(-1,1.01,0.1)
    x_side = np.zeros(y_side.shape)
    for i,j in enumerate(y_side):
        x_side[i] = 1

    #Step 3 - bottom half
    x_center = 0
    y_center = -1
    rad_bottom = 1
    theta_bottom = np.linspace(np.pi, 2*np.pi, 20)
    x_bottom = x_center + rad_bottom * np.cos(theta_bottom)
    y_bottom = y_center + rad_bottom * np.sin(theta_bottom)

    #Moving
    y_top += y_move 
    y_side += y_move 
    y_bottom += y_move

    #Plotting the letter
    plt.scatter(x_top, y_top, c="gray", s=4)
    plt.scatter(x_side, y_side, c="gray", s=4)
    plt.scatter(x_bottom, y_bottom, c="gray", s=4)

    r_top = np.sqrt(x_top**2 + y_top**2)
    r_side = np.sqrt(x_side**2 + y_side**2)
    r_bot = np.sqrt(x_bottom**2 + y_bottom**2)

    return r_top, r_side, r_bot

@jit
def black_hole(): #Function to plot the black hole and a photon sphere around it
    photon_sphere = 1.5
    theta = np.linspace(0, 2*np.pi, 1000)
    fig, ax = plt.subplots()
    ax.set_title("Black hole and photon sphere")
    ax.set_xlim(-25, 25)
    ax.set_ylim(-25, 25)
    ax.set_aspect("equal")
    ax.grid()
    #Plotting the photon sphere
    x_phot_sph = photon_sphere * np.cos(theta)
    y_phot_sph = photon_sphere * np.sin(theta)
    plt.plot(x_phot_sph, y_phot_sph, c="black")

    #Plotting the black hole
    x_bh =  np.cos(theta) 
    y_bh =  np.sin(theta)
    ax.fill(x_bh, y_bh, "black")



r_top, r_side, r_bot = screen_base() #Plotting the screen
black_hole() #Plotting the black hole



'''

Here begins the part, where the program calculates the trajectory of photons around a black hole

Do NOT modify the functions unless you wanna rewrite this entire part

The only parts that can be changed are definitions of rmin, N and e below

'''
@jit
def rmin_to_b(r): #Function to convert rmin to b (impact parameter)
    return r * np.sqrt(r/(r-1))

@jit
def dphi(dr,r,b): #Integrand
    return dr/(r**2 * np.sqrt(1/b**2 - (1 - 1/r) * 1/r**2))

@jit
def deflection(phi): #Rotation
    return phi - np.pi/2

rmin = np.arange(5.5, 7.01, 0.01) #Parameter fitting must be done manually to determine \\
                                #if the photon hits the screen - if I remember, I will add \\
                                #the correct set of parameters to the beginning of the program

print(rmin.shape)
N = 1000000 #number of steps - this can reduced/increased to reduce/increase precision \\
            #This must be a high enough number to actually see the photon curve around a black hole 
            #1e6 is enough to see one loop - higher numbers may result in more loops \\
            #but the program will also run longer

e = 5e-4 #epsilon - increase/decrease to decrease/increase precision
         #Do NOT go lower than the number already put here \\
         #otherwise the program might not find any positions

#Plotting the screen to see if a photon hits it or not
y_screen = [-distance_scr, -distance_scr]
x_screen = [-2, 2]
plt.plot(x_screen, y_screen, c='g')

'''

Here starts the main loop

Integrating to distance_obs (distance to observer). It can be set to whatever you want \\
but 200 works just fine

'''

f = open("wyniki.txt", "w")
f.write(np.array2string(r_top) + "\n")
f.write(np.array2string(r_side) + "\n")
f.write(np.array2string(r_bot) + "\n")
for j, rstart in enumerate(rmin):
    #Constants/initial position
    r = np.zeros(N) 
    phi = np.zeros(N)
    phi[0] = 0
    r[0] = rstart
    dr = (distance_obs - rstart)/N  
    b = rmin_to_b(rstart)

    '''
    Refer to warning 2 at the beginning of the program if you don't understand this statement
    '''
    if(r[0] <= 1.5):
        print("Minimal distance from the black hole is too small:", r[0])
        if(rmin.shape[0] < 2):
            print("Breaking the loop")
            break
        else:
            print("Skipping this iteration")
            continue
    '''

    PART 1

    Integrating trajectory of a photon from rmin to distance_obs

    '''
    for i in range(N-1): 
        r[i+1] = r[i] + dr
        phi[i+1] = phi[i] + dphi(dr, r[i+1], b)

    phi_defl = deflection(phi[N-1]) #Rotating the photons
    phi -= phi_defl 
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.scatter(x, y, c='r', s=3)

    '''

    PART 2

    Continuation of the previous trajectory - integrating from rmin to distance_obs

    '''
    for i in range(N-1):
        r[i+1] = r[i] + dr 
        phi[i+1] = phi[i] - dphi(dr, r[i+1], b) 

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.scatter(x, y, c='b', s=3)

    #Simple loop to determine where the photon hit the screen
    print("Loop number", j)
    for i in range(N-1):
        if(np.abs(y[i]) - 20 < e and np.abs(y[i]) - 20 > 0 and np.abs(x[i]) < 2):
            #print(y[i], x[i])
            f.write(f'{j} {y[i]} {x[i]}\n')
        elif(np.abs(y[i]) - 20 < e and np.abs(y[i]) - 20 > 0 and np.abs(x[i]) > 2):
            print("Photon failed to hit the screen")

    '''

    PART 3

    Plotting a black dot at initial point (at rmin)

    '''

    xstart = r[0] * np.cos(phi[0])
    ystart = r[0] * np.sin(phi[0]) 
    plt.scatter(xstart, ystart, c='black', s=20)



f.close()







plt.show()
