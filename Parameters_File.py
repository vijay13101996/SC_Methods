
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
        Nxr = 2#(20,80,20,80) needs to be run again; MemoryError
        Nxi = 2
        Nyr = 200
        Nyi = 200
        widthxr = 8 
        widthxi = 8
        widthyr = 6
        widthyi = 6
        print('N',Nxr,Nxi,Nyr,Nyi)
        print('width',widthxr,widthxi,widthyr,widthyi)

    elif(method=='Split'):
        global dt,tfinal
        dt = 0.01
        tfinal = 5.0

    elif(method=='HHKK'):
        global widthx,widthpx,widthy,widthpy,Nx,Npx,Ny,Npy,deltax,deltapx,deltay,deltapy
        Nx = 10
        Npx = 10
        Ny = 10
        Npy = 10
        widthx = 8.0
        widthpx = 5.0
        widthy = 8.0
        widthpy = 5.0

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
        global gammax,gammay
        gammax=0.5
        gammay=0.5
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
        print('initial conditions','x,y,px,py',x0,y0,px0,py0)

    elif(setindex==2):  #This is the requisite conditions for the reflected and transmitted wavepackets to separate.
        x0 = -5.0
        y0 = 0.0
        px0 = 2.0
        py0 = 0.0
        print('initial conditions','x,y,px,py',x0,y0,px0,py0)

    elif(setindex==3):            #This is the requisite conditions for deep tunneling conditions
        x0 = -5.0
        y0 = 0.0
        px0 = 4.0
        py0 = 0.0
        print('initial conditions','x,y,px,py',x0,y0,px0,py0)

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
        print('grid','xi,yi,lx,ly','div',xi,yi,lx,ly,div)

    elif(setindex==2):
        xi = -100.0#              This is the requisite conditions for the reflected and transmitted wavepackets to separate.
        lx = 200.0#
        dx = lx/div
        
        yi = -100.0
        ly = 200.0
        dy = ly/div
        print('grid','xi,yi,lx,ly','div',xi,yi,lx,ly,div)

    elif(setindex==3):

        xi = -50.0            #This is the requisite conditions for deep tunneling conditions
        lx = 100.0
        yi = -50.0
        ly = 100.0
        #div = 128
        dx = lx/div
        dy = ly/div
        print('grid','xi,yi,lx,ly','div',xi,yi,lx,ly,div)

def time_grid(setindex):
    global t,dt,T
    if(setindex==1):
        dt = 1.0
        t = 3.0 + 0.1*dt
        T=arange(0,t+0.1*dt,dt)
        print('dt,t',dt,t)

    elif(setindex==2):
        dt = 3.0              # This is the requisite conditions for the reflected and transmitted wavepackets to separate, for the above(usual) case.
        t = 12.0
        T = [0.0,9.0,12.0]
        T=arange(0,t+0.1*dt,dt)
        print('dt,t',dt,t)
        
 
    elif(setindex==3):      #This is the requisite conditions for deep tunneling conditions
        dt = 2.5
        t = 10.0
        T=arange(0,t+0.1*dt,dt)
        print('dt,t',dt,t)

def pot_select(potkey):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4
    if(potkey=='double_slit'):
        print('Potential:','double_slit')
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
        print('Potential:','harmonic')
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
        print('Potential:','coupled_harmonic')
        
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

    elif(potkey=='eckart'):
        print('Potential:','eckart')
        global V0,w,E
        V0 = 4.0
        w = 4.0
        E = 1.0
        def potential(x,y):
            return (V0)*E/(np.cosh(x)**2)

        def dxpotential(x,y):
            return -2*E*V0*np.sinh(x)/np.cosh(x)**3

        def dypotential(x,y):
            return 0

        def ddpotential1(x,y):
            return 2*E*V0*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2

        def ddpotential2(x,y):
            return 0

        def ddpotential3(x,y):
            return 0

        def ddpotential4(x,y):
            return 0
            
    elif(potkey=='double_slit_tunneling'):
        print('Potential:','double_slit_tunneling')
        global V0,w,E,Vc
        V0 = 16.0
        w = 4.0
        E = 1.0
        Vc = 0.0
        def potential(x,y):
            return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0) + Vc)*E/(np.cosh(x)**2)

        def dxpotential(x,y):
            return -2*E*(V0 + Vc - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*np.sinh(x)/np.cosh(x)**3

        def dypotential(x,y):
            return E*(-1.0*m*w**2*y + m**2*w**4*y**3/(4*V0))/np.cosh(x)**2

        def ddpotential1(x,y):
            return E*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)*(16*V0 + 16*Vc - 8.0*m*w**2*y**2 + m**2*w**4*y**4/V0)/(8*np.cosh(x)**2)

        def ddpotential2(x,y):
            return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

        def ddpotential3(x,y):
            return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

        def ddpotential4(x,y):
            return -E*m*w**2*(1 - 3*m*w**2*y**2/(4*V0))/np.cosh(x)**2

    elif(potkey=='reaction_coordinate'):
        global w,E,Vc,D,alpha
        
        alpha = 1.0
        D = 1.0
        w = 4.0
        E = 1.0
        Vc = 8.0

        def potential(x,y):
            return E*(D*(1 - exp(-alpha*y))**2 + Vc)/np.cosh(x)**2
        
        def dxpotential(x,y):
            return -2*E*(D*(1 - exp(-alpha*y))**2 + Vc)*sinh(x)/np.cosh(x)**3
        
        def dypotential(x,y):
            return 2*D*E*alpha*(1 - exp(-alpha*y))*exp(-alpha*y)/np.cosh(x)**2
        
        def ddpotential1(x,y):
            return 2*E*(D*(1 - exp(-alpha*y))**2 + Vc)*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2
        
        def ddpotential2(x,y):
            return -4*D*E*alpha*(1 - exp(-alpha*y))*exp(-alpha*y)*np.sinh(x)/np.cosh(x)**3
        
        def ddpotential3(x,y):
            return -4*D*E*alpha*(1 - exp(-alpha*y))*exp(-alpha*y)*np.sinh(x)/np.cosh(x)**3
        
        def ddpotential4(x,y):
            return 2*D*E*alpha**2*(-1 + 2*exp(-alpha*y))*exp(-alpha*y)/np.cosh(x)**2
    
        # def potential(x,y):
                    # return ((Vc+w*y**2)*E/np.cosh(x)**2)
                    
        # def dxpotential(x,y):
                    # return (-2*E*(Vc + w*y**2)*np.sinh(x)/np.cosh(x)**3)
                    
        # def dypotential(x,y):
                    # return (2*E*w*y/np.cosh(x)**2)

        # def ddpotential1(x,y):
                    # return (2*E*(Vc + w*y**2)*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2)
                    
        # def ddpotential2(x,y):
                    # return (-4*E*w*y*np.sinh(x)/np.cosh(x)**3)
                    
        # def ddpotential3(x,y):
                    # return (-4*E*w*y*np.sinh(x)/np.cosh(x)**3)

        # def ddpotential4(x,y):
                    # return (2*E*w/np.cosh(x)**2)
                
    elif(potkey=='diffraction'):
        global V0,Vc
        V0 = 16.0
        Vc=0.0
        
        
        def potential(x,y):
                    return (V0/(np.cosh(x)**2*np.cosh(y)**2))
                    
        def dxpotential(x,y):
                    return (-2*V0*np.sinh(x)/(np.cosh(x)**3*np.cosh(y)**2))
                    
        def dypotential(x,y):
                    return (-2*V0*np.sinh(y)/(np.cosh(x)**2*np.cosh(y)**3))

        def ddpotential1(x,y):
                    return (2*V0*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/(np.cosh(x)**2*np.cosh(y)**2))
                    
        def ddpotential2(x,y):
                    return (4*V0*np.sinh(x)*np.sinh(y)/(np.cosh(x)**3*np.cosh(y)**3))
                    
        def ddpotential3(x,y):
                    return (4*V0*np.sinh(x)*np.sinh(y)/(np.cosh(x)**3*np.cosh(y)**3))
                    
        def ddpotential4(x,y):
                    return (2*V0*(3*np.sinh(y)**2/np.cosh(y)**2 - 1)/(np.cosh(x)**2*np.cosh(y)**2))
                    
                    
    elif(potkey=='gaspard'):
    
        global V0,Vc
        V0 = 64.0
        Vc=0.0
        edge = 4.0
        x1 = 0.0
        y1 = -0.5*edge
        x2 = 0.5*edge
        y2 = (0.5*3**0.5 -0.5)*edge
        x3 = -0.5*edge
        y3 = (0.5*3**0.5 -0.5)*edge
        
        
        def potential(x,y):
                    return V0/(cosh(x - x3)**2*cosh(y - y3)**2) + V0/(cosh(x - x2)**2*cosh(y - y2)**2) + V0/(cosh(x - x1)**2*cosh(y - y1)**2)
                    
        def dxpotential(x,y):
                    return -2*V0*sinh(x - x1)/(cosh(x - x1)**3*cosh(y - y1)**2) - 2*V0*sinh(x - x2)/(cosh(x - x2)**3*cosh(y - y2)**2) - 2*V0*sinh(x - x3)/(cosh(x - x3)**3*cosh(y - y3)**2)
                    
        def dypotential(x,y):
                    return -2*V0*sinh(y - y1)/(cosh(x - x1)**2*cosh(y - y1)**3) - 2*V0*sinh(y - y2)/(cosh(x - x2)**2*cosh(y - y2)**3) - 2*V0*sinh(y - y3)/(cosh(x - x3)**2*cosh(y - y3)**3)

        def ddpotential1(x,y):
                    return 2*V0*(3*sinh(x - x1)**2/(cosh(x - x1)**4*cosh(y - y1)**2) + 3*sinh(x - x2)**2/(cosh(x - x2)**4*cosh(y - y2)**2) + 3*sinh(x - x3)**2/(cosh(x - x3)**4*cosh(y - y3)**2) - 1/(cosh(x - x3)**2*cosh(y - y3)**2) - 1/(cosh(x - x2)**2*cosh(y - y2)**2) - 1/(cosh(x - x1)**2*cosh(y - y1)**2))
                    
        def ddpotential2(x,y):
                    return 4*V0*(sinh(x - x1)*sinh(y - y1)/(cosh(x - x1)**3*cosh(y - y1)**3) + sinh(x - x2)*sinh(y - y2)/(cosh(x - x2)**3*cosh(y - y2)**3) + sinh(x - x3)*sinh(y - y3)/(cosh(x - x3)**3*cosh(y - y3)**3))

                    
        def ddpotential3(x,y):
                    return 4*V0*(sinh(x - x1)*sinh(y - y1)/(cosh(x - x1)**3*cosh(y - y1)**3) + sinh(x - x2)*sinh(y - y2)/(cosh(x - x2)**3*cosh(y - y2)**3) + sinh(x - x3)*sinh(y - y3)/(cosh(x - x3)**3*cosh(y - y3)**3))

                    
        def ddpotential4(x,y):
                    return 2*V0*(3*sinh(y - y1)**2/(cosh(x - x1)**2*cosh(y - y1)**4) + 3*sinh(y - y2)**2/(cosh(x - x2)**2*cosh(y - y2)**4) + 3*sinh(y - y3)**2/(cosh(x - x3)**2*cosh(y - y3)**4) - 1/(cosh(x - x3)**2*cosh(y - y3)**2) - 1/(cosh(x - x2)**2*cosh(y - y2)**2) - 1/(cosh(x - x1)**2*cosh(y - y1)**2))
                







                



        
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


