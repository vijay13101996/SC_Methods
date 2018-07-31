import numpy as np
from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt

hbar = 1.0
k= 1.0
xi = -20.0
lx = 40.0
div = 512
dx = lx/div
m = 1.0
D = 0.5

N = 10
Nxr = N
Nxi = N
Nyr = N
Nyi = N
s= 1.0

widthx = 4.0
widthy = 4.0

potkey = 'double_slit'
remark1 = 'eckart'


print("width",widthx,widthy) # This is the optimal width. Don't change this anymore
print("N",Nxr,Nxi,Nyr,Nyi)
yi = -20.0
ly = 40.0
dy = ly/div

x0 = -5.0
y0 = 0.0
px0 = 3.0
py0 = 0.0

gammax = 0.5
gammay = 0.5

gammaxi = 0.5
gammayi = 0.5
gammaxf = 0.5
gammayf = 0.5

V0 = 16.0
w = 4.0
E = 1.0
dt = 0.5

dt = 0.5
t = 3.0 + 0.1*dt

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

##xi = -100.0#              This is the requisite conditions for the reflected and transmitted wavepackets to separate.
##lx = 200.0#
##div = 256
##dx = lx/div
##m = 1.0
##
##yi = -100.0#
##ly = 200.0#
##dy = ly/div
##dt = 4.0
##t = 12.0


def potential(x,y):
    return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*E/(np.cosh(x)**2)#0.5*k*x**2 + 0.5*k*y**2 #V0*E/(np.cosh(x)**2*np.cosh(y)**2)#*(exp(-x**2/s**2)
#(V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0)))*##

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

