import numpy as np
#import trial_code
import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
import pickle, pprint
from matplotlib import pyplot as plt

V0 = 16.0
w = 4.0
E = 16.0
m = 1.0
D = 1.0
alpha = 1

q0 = -5.0
p0 = 2.0

dt=0.01
tre = np.arange(0.0,8+0.1*dt,dt)
tim = np.arange(-4.,4.+0.1*dt,dt)
treal,timag = np.meshgrid(tre,tim)


potkey = 'morse'
potkey = 'quasi_stable_barrier'
potkey = 'pendulum'
potkey = 'double_well'

#potkey = 'inverted_higgs'


if(potkey=='eckart'):
    global potential, dpotential

    def potential(q):
                return E/(np.cosh(q)**2)

    def dpotential(q):
                return -2*E*np.sinh(q)/np.cosh(q)**3
                
if(potkey=='double_well'):
    global potential, dpotential

    def potential(q):
                return (V0 - 0.5*m*w**2*q**2 + m**2*w**4*q**4/(16*V0))

    def dpotential(q):
                return -1.0*m*q*w**2 + m**2*q**3*w**4/(4*V0)
                
if(potkey=='pendulum'):
    global potential,dpotential
    
    def potential(q):
        return (1-np.cos(q))
        
    def dpotential(q):
        return np.sin(q)
   
if(potkey=='gaussian'):
    global potential,dpotential
    
    def potential(q):
        return E*np.exp(-q**2)
        
    def dpotential(q):
        return -2*E*q*np.exp(-q**2)

if(potkey=='inverted_harmonic'):
    global potential,dpotential
    
    def potential(q):
        return E-q**2
        
    def dpotential(q):
        return -2*q

if(potkey=='quartic'):
    global potential, dpotential
    
    def potential(q):
        return E*q**4
        
    def dpotential(q):
        return 4*E*q**3
   
if(potkey=='inverted_higgs'):
    global potential,dpotential
    
    def potential(q):
        return E*(0.5*q**2-0.5*q**4)
        
    def dpotential(q):
        return E*(q-2*q**3)
        
if(potkey=='lorentzian'):
    global potential,dpotential

    def potential(q):
        return V0/(q**6+1)#V0/(q**2 + 1)
        
    def dpotential(q):
        return -6*V0*q**5/(q**6 + 1)**2#-4*V0*q**3/(q**4 + 1)**2     #-2*V0*q/(q**2 + 1)**2

if(potkey=='morse'):
    global potential,dpotential
    
    def potential(q):
        return D*(1-np.exp(-alpha*q))**2
        
    def dpotential(q):
        return 2*D*alpha*(1 - np.exp(-alpha*q))*np.exp(-alpha*q)
        
if(potkey=='quasi_stable_barrier'):
    global potential,dpotential
    q1 = -1.0
    q2 = 1.0
    
    def potential(q):
        return E/np.cosh(q - q2)**2 + E/np.cosh(q - q1)**2
    
    def dpotential(q):
        return -2*E*np.sinh(q - q1)/np.cosh(q - q1)**3 - 2*E*np.sinh(q - q2)/np.cosh(q - q2)**3
   
print('potential at origin',potential(0))  
    
qinit = q0+0j
pinit = p0+0j