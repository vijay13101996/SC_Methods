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
#from Parameter_file_Singularity_time_structure import *
import Complex_plotter
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
import io

plot_full_time = 0



plot = 1
tarr = []
pxarr = []
pyarr = []
xarr = []
yarr=[]
trajon = 0
dt=0.02
tre = np.arange(0.0,8+0.1*dt,dt)
tim = np.arange(-2.,2.+0.1*dt,dt)


#print(tre,tim)
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

height = 1
c = mcolors.ColorConverter().to_rgb
rgb = make_colormap(\
    [c('red'), c('yellow'), 0.15, c('yellow'), c('green'), c('cyan'), 0.5, c('cyan'),
     c('blue'), c('magenta'), 0.85, c('magenta'), c('red')])
     
def f(t,r,dti):
    global current, tPrint, tPrintSpace
    #previous = turningpoints[1]
    #print("Hi")
    q,p =  r
    dxpotential(q) = np.sin(q)
    drdt = array([p,-dxpotential(q)])*dti
    return drdt

    
def faux(t,r,dti):
    global current, tPrint, tPrintSpace,trajon
    r = r.view(dtype=complex)
    #print('r',r)
    #print(len(r))
    q,p=r
    #print('t',t,'px',px)
    #print('changed')
    
    # if(trajon==1):
        # print('here')
        # tarr.append(t)
        # pxarr.append(px)
        # pyarr.append(py)
        
    #print('pole',py,py*(x-1j*pi/2)**0.5)#,(x-1j*pi/2))
    #print('Energy',px**2/2 + py**2/2 + potential(x,y),'Potential',potential(x,y))
    #print('x',x,'px',px)
    dxpotential(q) = np.sin(q)
    
    drdt = array([p,-dxpotential(q)])*dti
    #print('residue',x+1j*pi/2,(x+1j*pi/2)**3)              
    return drdt.view(dtype=float)

def contour_integration(xinit,pxinit,trajdata,T):
    global trajon
    sol = ode(faux)
    sol.set_integrator('dop853',nsteps=1e10)
    sol.set_solout(solout)
            
    y0 = array([xinit,pxinit])
    
    print('y0',y0)
    if(len(T) ==1):
        sol.set_initial_value(y0.view(dtype=float),t=0.0)
        trajdata = sol.y.view(dtype=complex)
    else:
        
        tempsol = y0
        for tc in range(len(T)-1):
            print('Time',T[tc])
            if(tc==len(T)-2):
                trajon=1
                print('started')
                
            #print(tc)
            
            #print(T[tc+1]-T[tc])
            timcom = T[tc+1]-T[tc]
            abstim = abs(timcom)
            dti = timcom/abstim
            sol.set_initial_value(tempsol.view(dtype=float),t=0.0)
            sol.set_f_params(dti)
            #print('y0',tempsol.view(dtype=float))
            sol.integrate(abstim)
            tempsol = sol.y.view(dtype=complex)
        trajdata=sol.y.view(dtype=complex)
        
        #print('Trajdata ',trajdata[:5],'xfinal',trajdata[0])
        return trajdata
        
def solout(t,r):
    global trajon
    r = r.view(dtype=complex)
    x,px=  r
    #if(trajon==1):
        #print('here')
    tarr.append(t)
    pxarr.append(px)
    xarr.append(x)
    
    
def finalcondition(xinit,pxinit,data,tre,tim):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1e3)
            
    y0 = array([xinit,pxinit])
    sol.set_initial_value(y0,t=0.0)

    dtim = tim[1]-tim[0]
    dtre = tre[1]-tre[0]

    print(dtim,dtre)
    

    ti = int(len(tim)/2) + 1
    dti = 1j
    sol.set_f_params(dti)
    trun = dtim
    while(ti<len(tim)):
        #print(Time[0][ti])
        sol.integrate(trun)
        data[0][ti] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        ti+=1
        trun+=dtim
        

    sol.set_initial_value(y0,t=0.0)
    dti = -1j
    sol.set_f_params(dti)
    ti = int(len(tim)/2) - 1
    trun = dtim
    while(ti>=0):
        sol.integrate(trun)
        data[0][ti] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        ti-=1
        trun+=dtim
        

    data[0][int(len(tim)/2)] = y0
    #print(Time[0][int(len(tim)/2)],data[0][int(len(tim)/2)][0] )
    #print(Time[0][int(len(tim)/2)])
    #print(data[:][:][0])
    dti = 1.0
    sol.set_f_params(dti)
    for ti in range(len(tim)):
        sol.set_initial_value(data[0][ti],t=0.0)
        trun=dtre
        for tr in range(1,len(tre)):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti],data[tr][ti][0])
            trun+=dtre

            

xinit = pi/2 + 0j
pxinit = 3*pi + 0j
data = zeros((len(tre),len(tim),2),dtype = complex)

finalcondition(xinit,pxinit,data,tre,tim)

op = open('Singularities_pendulum_point_1_branch_1.pkl','wb')
pickle.dump(data,op)
op.close()

Ttarget = [0.0,3*pi]
trajdata = zeros(2,dtype=complex)
trajdatatunnel = contour_integration(xinit,pxinit,trajdata,Ttarget)

plt.plot(xarr,pxarr)
plt.show()
xfinbar,pxfinbar = data[:,:,:2].transpose((2,0,1))
pxfinbar = matrix.transpose(pxfinbar)
Complex_plotter.plotcomplex(pxfinbar,1,10,tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1])
plt.show()
