from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt
#from matplotlib.figure import figure
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import trial_code
import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
import pickle, pprint
from matplotlib import pyplot as plt
from Parameter_file_Singularity_time_structure import *


print(xinit,yinit)
print(pointarr)

sing_time=0.0
cutoff = 1e3

timear = []



def f(t,r,offsetangle,offsetdistance,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, previous,tPrint, tPrintSpace
    #print('aux')
    #print('r',r)
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    dpdt = -dxpotential(x,y)
    if(t<offsetdistance):
        dtim = exp(1j*offsetangle)
    else:
        dtim=1/dpdt
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtim
    #timear.append(tcomcurr)
    #print('t',t,'tcomcurr',tcomcurr,'t^2px',px)
    #print('t',t,'px',px)#,x,tcomcurr,y)
    
    return drdt.view(dtype=float)


def solout(t,r):
    global sing_time,cutoff
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    timear.append(tcomcurr)
    #expr = tcomcurr-t1
    
    if(abs(px)>cutoff):   # Change at the other location too, when required
        print('Singularity time',tcomcurr,potential(x,y))
        sing_time=tcomcurr
        singtime_arr.append(sing_time)
        return -1 

singtime_arr = []    

    
def finalcondition(f,xinit,yinit,pxinit,pyinit,offsetangle,offsetdistance):
    sol = ode(f)
    sol.set_integrator('dop853',nsteps=1e16)
    sol.set_solout(solout)
    
    
    y0 = array([xinit,pxinit,yinit,pyinit,0.0])
    #print('y0',y0)
    sol.set_initial_value(y0.view(dtype=float),t=0.0)
    sol.set_f_params(offsetangle,offsetdistance,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    sol.integrate(1e12)
    result = sol.y.view(dtype=complex)
    #print('res',result)
    return result
    

offsetdistance = arange(0,3.5+0.1,0.5)
offsetangle = arange(-pi/2,pi/2 +0.01*pi,pi/40)

for dist in offsetdistance:
    for theta in offsetangle:
        finalcondition(f,xinit,yinit,pxinit,pyinit,theta,dist)
#print(singtime_arr)
op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
pickle.dump(singtime_arr,op)
op.close()


    
