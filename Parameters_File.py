
import numpy as np
from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt


hbar = 1.0
m = 1.0
s= 1.0
#-------------------------------------------------------------
#Note:
# 1 -  Usual 3 units time double_slit_eckart
# 2 -  12 units of time for the same case above
# 3 -  Deep tunneling case.

# 6 different functions are to be called to set initial values for each run file:
#method_sel,wf_init_set,initial_conditions,pot_select,time_grid,grid_select




#------------------------------------------------------------
def method_sel(method):
    if(method=='FINCO'):
        global Nxr,Nxi,Nyr,Nyi,widthxr,widthxi,widthyr,widthyi
        Nxr = 10
        Nxi = 10
        Nyr = 10
        Nyi = 10
        widthxr = 7
        widthxi = 2
        widthyr = 4
        widthyi = 4

    elif(method=='Split'):
        global dt,tfinal
        dt = 0.01
        tfinal = 12.0

    elif(method=='HHKK'):
        global widthx,widthpx,widthy,widthpy,Nx,Npx,Ny,Npy,deltax,deltapx,deltay,deltapy
        widthx = 8.0
        widthy = 8.0

        deltax = widthx/Nx
        deltay = widthy/Ny
        deltapx = widthpx/Npx
        deltapy = widthpy/Npy

def wf_init_set(method):
    global gammaxi,gammayi,gammaxf,gammayf
    if(method=='FINCO'):
        
        gammaxi = 0.5
        gammayi = 0.5
        gammaxf = 0.5
        gammayf = 0.5

    elif(method=='Split'):
        gammaxi = 0.5
        gammayi = 0.5
        gammaxf = 0.5
        gammayf = 0.5

    elif(method=='HHKK'):
        gammaxi = 1.0
        gammayi = 1.0
        gammaxf = 1.0
        gammayf = 1.0
        
        
def initial_conditions(setindex):
    global x0,y0,px0,py0
    if(setindex==1):
        x0 = -5.0
        y0 = 0.0
        px0 = 3.0
        py0 = 0.0

    elif(setindex==2):
        x0 = -5.0
        y0 = 0.0
        px0 = 3.0
        py0 = 0.0

    elif(setindex==3):            #This is the requisite conditions for deep tunneling conditions
        x0 = -5.0
        y0 = 0.0
        px0 = 0.5
        py0 = 0.0

def grid_select(setindex):
    global div,yi,ly,xi,lx,dx,dy
    div = 512
    if(setindex==1):
        yi = -20.0           #Usual case.
        ly = 40.0
        dy = ly/div

        xi = -20.0
        lx = 40.0
        dx = lx/div

    elif(setindex==2):
        xi = -100.0#              This is the requisite conditions for the reflected and transmitted wavepackets to separate.
        lx = 200.0#
        dx = lx/div
        
        yi = -100.0
        ly = 200.0
        dy = ly/div

    elif(setindex==3):

        xi = -50.0            #This is the requisite conditions for deep tunneling conditions
        lx = 100.0
        yi = -50.0
        ly = 100.0
        
        dx = lx/div
        dy = ly/div

def time_grid(setindex):
    global t,dt
    if(setindex==1):
        dt = 0.5
        t = 3.0 + 0.1*dt

    elif(setindex==2):
        dt = 3.0#              This is the requisite conditions for the reflected and transmitted wavepackets to separate, for the above(usual) case.
        t = 12.0
 
    elif(setindex==3):      #This is the requisite conditions for deep tunneling conditions
        dt = 2.5
        t = 10.0

def pot_select(potkey):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4
    if(potkey=='double_slit'):
        global V0,w,E
        V0 = 16.0
        w = 4.0
        E = 1.0
        def potential(x,y):
            return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*E/(np.cosh(x)**2)

        def dxpotential(x,y):
            return -2*E*(V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*np.sinh(x)/np.cosh(x)**3

        def dypotential(x,y):
            return E*(-1.0*m*w**2*y + m**2*w**4*y**3/(4*V0))/np.cosh(x)**2

        def ddpotential1(x,y):
            return E*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)*(16*V0 - 8.0*m*w**2*y**2 + m**2*w**4*y**4/V0)/(8*np.cosh(x)**2)

        def ddpotential2(x,y):
            return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

        def ddpotential3(x,y):
            return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

        def ddpotential4(x,y):
            return -E*m*w**2*(1 - 3*m*w**2*y**2/(4*V0))/np.cosh(x)**2

    elif(potkey=='harmonic'):
        global k
        k= 1.0
        
        def potential(x,y):
            return 0.5*k*x**2 + 0.5*k*y**2 

        def dxpotential(x,y):
            return k*x

        def dypotential(x,y):
            return k*y

        def ddpotential1(x,y):
            return k

        def ddpotential2(x,y):
            return 0.0

        def ddpotential3(x,y):
            return 0.0

        def ddpotential4(x,y):
            return k

    elif(potkey=='coupled_harmonic'):
        global k
        k=1.0
        def potential(x,y):
            return 0.5*k*x**2 + 0.5*k*y**2 + x*y

        def dxpotential(x,y):
            return k*x + y

        def dypotential(x,y):
            return k*y + x

        def ddpotential1(x,y):
            return k

        def ddpotential2(x,y):
            return 1.0

        def ddpotential3(x,y):
            return 1.0

        def ddpotential4(x,y):
            return k


        
#print(Nxr,Nyr,Nyr,Nyi)
##widthx = 8.0
##widthy = 4.0

#print("width",widthxr,widthxi,widthyr,widthyi) # This is the optimal width. Don't change this anymore
#print("N",Nxr,Nxi,Nyr,Nyi)



##dt = 0.5
##
##dt = 0.5
##t = 3.0 + 0.1*dt

##xi = -50.0            #This is the requisite conditions for deep tunneling conditions
##lx = 100.0
##div = 256
##yi = -50.0
##ly = 100.0
##dx = lx/div
##dy = ly/div
##x0 = -5.0
##y0 = 0.0
##px0 = 0.5
##py0 = 0.0
##dt = 2.5
##t = 10.0

##dt = 3.0#              This is the requisite conditions for the reflected and transmitted wavepackets to separate.
##t = 12.0


##def potential(x,y):
##    return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*E/(np.cosh(x)**2)
##
##def dxpotential(x,y):
##    return -2*E*(V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*np.sinh(x)/np.cosh(x)**3
##
##def dypotential(x,y):
##    return E*(-1.0*m*w**2*y + m**2*w**4*y**3/(4*V0))/np.cosh(x)**2
##
##def ddpotential1(x,y):
##    return E*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)*(16*V0 - 8.0*m*w**2*y**2 + m**2*w**4*y**4/V0)/(8*np.cosh(x)**2)
##
##def ddpotential2(x,y):
##    return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3
##
##def ddpotential3(x,y):
##    return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3
##
##def ddpotential4(x,y):
##    return -E*m*w**2*(1 - 3*m*w**2*y**2/(4*V0))/np.cosh(x)**2


