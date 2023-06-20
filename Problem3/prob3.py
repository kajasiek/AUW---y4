import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import pandas as pd
from matplotlib import cm
import matplotlib as mpl
import sys

'''
Trajectory calculations are based on the solution from
David D. Nolte -
Introduction to Modern Dynamics - Chaos, Networks, Space and Time
Oxford University Press (2015)
'''

'''

This code is based on someone else's work
I have not tested it and I do not know how it works

'''

'''Stop writing in my code'''

#constants
RS = 1          #Schwarzschild radius
R_ekran = 2*RS     #
d_ekran = 20*RS    #
d_obs = 500*RS
nFlag = 0       #near flag (particle is within d_ekran)
fFlag = 0       #far flag (particle has left d_ekran radius)
thetas = np.linspace(0,2*np.pi,1000) #angles
#precision
prec = 0.0001

def screen():
    '''Draws the screen and letter'''
    x = R_ekran*np.cos(thetas)
    y = R_ekran*np.sin(thetas)
    plt.figure(0,figsize=(8,8))
    plt.plot(x,y)

#Drawing letter "P"
    tlinex = np.arange(-.75,.5,prec)
    tlinex2 = np.arange(-.5,.25,prec)
    slinex = np.arange(.5,.75,prec)
    slinex2 = np.arange(.25,.5,prec)
    blinex = np.arange(-.5,.25,prec)
    blinex2 = np.arange(-.5,.5,prec)
    vlinex1 = np.array([-.75]*int(2/prec))
    vlinex2 = np.array([-.5]*int(0.75/prec))
    vlinex3 = np.array([.5]*int(.25/prec))
    vlinex4 = np.array([.75]*int(.75/prec))
    tliney = np.array([1]*int(1.25/prec))
    tliney2 = np.array([.75]*int(.75/prec))
    sliney = np.arange(1,.75,-prec)
    sliney2 = np.arange(.75,.5,-prec)
    sliney3 = np.arange(-.25,0,prec)
    sliney4 = np.arange(0,.25,prec)
    bliney = np.array([0]*int(.75/prec))
    bliney2 = np.array([-.25]*int(1/prec))
    vliney1 = np.arange(-1,1,prec)
    vliney2 = np.arange(-1,-.25,prec)
    vliney2b = np.arange(0,.75,prec)
    vliney3 = np.arange(.25,.5,prec)
    vliney4 = np.arange(0,.75,prec)
    bblinex = np.arange(-.75,-.5,prec)
    bbliney = np.array([-1]*int(.25/prec))
    xlist = np.concatenate([tlinex,tlinex2,slinex,slinex2,slinex,slinex2,
                          blinex,blinex2,vlinex1,vlinex2,
                          vlinex2,vlinex3,vlinex4,bblinex])
    ylist = np.concatenate([tliney,tliney2,sliney,sliney2,sliney3,sliney4,
                                bliney,bliney2,vliney1,vliney2,vliney2b,
                                vliney3,vliney4,bbliney])

    #In polar coordinates
    rlist = np.sqrt(xlist**2+ylist**2)
    tlist = np.arctan2(ylist,xlist)

    plt.scatter(xlist,ylist,c='k',s=4)
    plt.xlabel('x [RS]')
    plt.ylabel('y [RS]')
    plt.title('Original letter and screen edge')

    #Debug
    #plt.show()
    return rlist,tlist

def refraction(x,y):
    '''The refraction index
Taken from D.D. Nolte - p. 403'''
    rr0 = np.sqrt(x**2+y**2)
#n(r) = (1-((2*G*M)/(c^2*r)))**-1
#Schwarzschild radius:  RS = 2*G*M/c^2
    fac = np.abs((1-9*(RS/rr0)**2/8))   #Eikonal Correction Factor
    
    n = 1/(1-(RS/rr0))
    nx = -fac*n**2*RS*x/rr0**3
    ny = -fac*n**2*RS*y/rr0**3
    return [n,nx,ny]

def flow_deriv(rr,dt):
    '''Differential equation to solve - Dependent on n'''
    x,y,z,w=rr
    [n,nx,ny] = refraction(x,y)

    yp = np.zeros(shape=(4,))
    yp = z/n,w/n,nx,ny

    return yp

def draw_trajectory():
    '''test function for manually checking if ray hits
and determining range of ystart'''
    for i in range(0,10):
        xstart=-d_obs
        #0loop: 50
        #ymin - 6.825
        #ymax - 8.362
        #1loop: 50
        #ymin - 3.6538205
        #ymax - 3.6539293

        #0loop: 0
        #ymin - 4.0004
        #ymax - 5.5395
        #1loop: 0
        #ymin - 1.50003095
        #ymax - 1.50004367
        ystart=1.50003095+0.000001*i
        print(ystart)
        [n,nx,ny]=refraction(xstart,ystart)
        y0 = [xstart,ystart,n,0]
        tspan=np.linspace(1,700,3500)
        y=integrate.odeint(flow_deriv,y0,tspan)
        plt.figure(1)
        lines=plt.plot(y[:,0],y[:,1])
        plt.setp(lines,lw=1)
        plt.title('Trajectories')

        #Schwarzschild radius of blackhole
    circle= plt.Circle((0,0),radius=RS,color='k')
    ax=plt.gca()
    ax.add_patch(circle)
    plt.axis('scaled')
    axes = plt.gca()
    axes.set_xlim([19,21])
    axes.set_ylim([-2.2,2.2])
    #axes.set_xlim([-21,21])
    #axes.set_ylim([-10,10])
    axes.set_xlabel('x [RS]')
    axes.set_ylabel('y [RS]')

    #Circular photon orbit:
    xstart = 0
    ystart = 1.5*RS
    [n,nx,ny] = refraction(xstart,ystart)
    y0 = [xstart,ystart,n,0]
    tspan=np.linspace(1,12,900)
    y = integrate.odeint(flow_deriv,y0,tspan)
    lines2 = plt.plot(y[1:,0],y[1:,1],lw=.5)
    plt.setp(lines2,lw=2,color='r',ls='--')
    plt.plot([d_ekran,d_ekran],
             [R_ekran,-R_ekran],color='r',lw=1)
    plt.show()
    

def main():
    #Test function below:
    #draw_trajectory()
    #sys.exit()
    radletter,thetaletter = screen()
    #setting up the trajectory:
#loop for several starting points
#Loops: For different observator distances: 500,50,0
    #loop hit lists:
    loop0 = np.linspace(6.958,8.5294,200)
    #loop0 = np.linspace(6.825,8.362,200)
    #loop0 = np.linspace(4.0004,5.5395,200)
    print(loop0[1]-loop0[0])
#Loops: For different observator distances: 500,50,0
    loop1 = np.linspace(3.720735,3.7208456,200)
    #loop1 = np.linspace(3.6538205,3.6539293,200)
    #loop1 = np.linspace(1.50003095,1.50004367,200)
    print(loop1[1]-loop1[0])
    #Second loop trajectories are dismissed due to rounding error
    #loop2 = np.linspace(3.7203734711,3.7203734951,100)
    #defining accuracy of "hit"
    eps = [1e-3,1e-8,5e-10]
    hitlist = []
    startlist = []
    for i,loop in enumerate([loop0,loop1]):
        epsloop = eps[i]
        for ystart in loop:
            xstart=-d_obs
            #no loops:
            #8.5294 - maximum ystart hit
            #6.958 - minimum ystart hit
            #1 loop:
            #3.7208456 - maximum ystart hit
            #3.720735 - minimum ystart hit
            #2 loops: (true for tspan length of 2000!)
            #3.7203734951 - maximum ystart hit
            #3.7203734711 - minimum ystart hit
            #ystart=loop0[i]
            #Debugging:
            #print(ystart)

            [n,nx,ny] = refraction(xstart,ystart)
            
            y0 = [xstart,ystart,n,0]
            #Time amount to consider
            tspan = np.linspace(1,700,3500)
            #integration of differential functions
            y = integrate.odeint(flow_deriv,y0,tspan)
            temp = y[y[:,0]>(d_ekran-epsloop)][0]
            temp = temp[:2]
            #Debug
            #print(f"Hit at:{temp}")
            startlist.append(ystart)
            hitlist.append(temp)
            plt.figure(1)
            lines=plt.plot(y[:,0],y[:,1])
            plt.setp(lines,lw=1)
    plt.title('Trajectories')
            
    #Schwarzschild radius of blackhole
    circle= plt.Circle((0,0),radius=RS,color='k')
    ax=plt.gca()
    ax.add_patch(circle)
    plt.axis('scaled')
    axes = plt.gca()
    axes.set_xlim([-21,21])
    axes.set_ylim([-10,10])
    axes.set_xlabel('x [RS]')
    axes.set_ylabel('y [RS]')

    #Circular photon orbit:
    xstart = 0
    ystart = 1.5*RS
    [n,nx,ny] = refraction(xstart,ystart)
    y0 = [xstart,ystart,n,0]
    tspan=np.linspace(1,12,900)
    y = integrate.odeint(flow_deriv,y0,tspan)
    lines2 = plt.plot(y[1:,0],y[1:,1],lw=.5)
    plt.setp(lines2,lw=2,color='r',ls='--')
    plt.plot([d_ekran,d_ekran],
             [R_ekran,-R_ekran],color='r',lw=1)
    #plt.show()
    #trajectory set
    #Now to calculate the disturbed picture
    hitlist=np.array(hitlist)
    hitlist = np.sqrt((hitlist[:,0]-20)**2+hitlist[:,1]**2)
    startlist = np.array(startlist)
    #print(hitlist)
    #Looking if ray hits point in 0.01 radius (precision of ekran)
    hitletter = []
    hitletterekran = []
    #hitlist - entire list of ray position when hitting ekran
    #radletter - r-coordinates of letter points
    #thetaletter - theta-coordinates of letter points
    #thetas - list of angles
    #startlist - list of ray starting displacement from x axis
    for j,r in enumerate(hitlist):
        for i,rletter in enumerate(radletter):
            if abs(r-rletter)<0.01:
                hitletterekran.append([rletter,thetaletter[i]])
                hitletter.append([startlist[j],thetaletter[i]])
    ekranhit = []
    for j,r in enumerate(hitlist):
        if abs(r-R_ekran)<0.01:
            for theta in thetas:
                ekranhit.append([startlist[j],theta])
    ekranhit = np.array(ekranhit)
    hitletter = np.array(hitletter)
    hitletterekran = np.array(hitletterekran)
    #print(hitletter[:,0])
    #Disturbed letter with edge hitting photons in blue
    plt.figure(3,figsize=(8,8))
    plt.scatter(hitletter[:,0]*np.cos(hitletter[:,1]),
                hitletter[:,0]*np.sin(hitletter[:,1]),s=4,c='k')
    plt.scatter(ekranhit[:,0]*np.cos(ekranhit[:,1]),
                ekranhit[:,0]*np.sin(ekranhit[:,1]),s=0.1,c='b')
    plt.xlim([-10,10])
    plt.ylim([-10,10])
    plt.title(f'Letter and screen edge as seen by observer (at {d_obs=})')
    plt.xlabel('x [RS]')
    plt.ylabel('y [RS]')
    #Letter as "seen" by the projected rays
    plt.figure(4,figsize=(8,8))
    plt.scatter(hitletterekran[:,0]*np.cos(hitletterekran[:,1]),
                hitletterekran[:,0]*np.sin(hitletterekran[:,1]),
                s=4,c='k')
    plt.xlim([-4,4])
    plt.ylim([-4,4])
    plt.title('Letter as seen by rays')
    plt.xlabel('x [RS]')
    plt.ylabel('y [RS]')
    plt.show()
    
    


if __name__ == "__main__":
    main()
    
        
