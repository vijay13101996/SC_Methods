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
cutoff = 5e2
#print(tre,tim)

#tprintspace=0.0001
timear = []
t1 = 1.3453399115321387+0.2956538655747396j
#t1 = -2.0 + 0.5j
#t1 = (2.9373765435748465-0.019267013663539785j)
#t1 =(1.7600502483227527+0.2454058275120668j)

n=0.5
def faux(t,r,pxinit,tinit,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, previous,tPrint, tPrintSpace
    #print('aux')
    #print('r',r)
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    dtim=1.0
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtim
    timear.append(tcomcurr)
    #print('t',t,'tcomcurr',tcomcurr,'t^2px',px)
    #print('t',t,'px',px)#,x,tcomcurr,y)
    
    return drdt.view(dtype=float)

def f(t,r,pxinit, tinit,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, previous,tPrint, tPrintSpace
    #print('usual')
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    dtim=0.0j
    dtimcheck=0.0j
    # if(t<headstarttime1):
         # dtim = 1.0
    # elif(t>headstarttime1 and t<headstarttime1+headstarttime2):
        # dtim=exp(1j*headstartangle)
    # # elif(t>headstarttime1+headstarttime2):
        # # dtim = -0.5*(1)/(px*dxpotential(x,y)+py*dypotential(x,y))
    # else:    
    expr = tcomcurr-t1
    #print('expr',expr)
    p = px
    dpdt = -1*(dxpotential(x,y))#+py*dypotential(x,y))
    #dtim = 1.0/(2*expr*p + expr**2*dpdt)
    # dtimcheck = px**2/(pxinit*dpdt)
    # dtimcheck1 = px**2*expr**n/(pxinit*(tinit-t1)**n*(dpdt +n*px/expr))
    
    exprt=[]
    expri=[]
    if(len(singtime_arr)==0):
        dtimcheck = px**2/(pxinit*dpdt)
    else:
    
        for singtimes in singtime_arr:
            exprt.append(tcomcurr-singtimes)
            expri.append(tinit-singtimes)
        
        term1 = px**2/pxinit
        #print('t',exprt)
        #print('i',expri)
        for (ext,exi) in zip(exprt,expri):
            #print('ext',ext)
            #print('exi',exi)
            term1*=ext**n/exi**n

        term2 = dpdt
        for ext in exprt:
            term2+=n*px/ext
            
        dtimcheck=term1/term2
    #print(dtimcheck)
    # dtimcheck3 = ((3.63001347383e-08-4.71238854539j)-xinit)/px
    # dtimcheck4 = ((3.63001347383e-08-4.71238854539j)*expr**2-xinit*t1**2)/(expr*(px*expr + 2*x*expr))
    #print('dtimcheck',dtimcheck1,dtimcheck2)
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtimcheck
    #timear.append(tcomcurr)
    #print('t',t,'tcomcurr',tcomcurr,'t^2px',px)
    #print(t,1/px,x,tcomcurr)
    
    return drdt.view(dtype=float)

def solout(t,r):
    global sing_time,cutoff
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    timear.append(tcomcurr)
    expr = tcomcurr-t1
    
    if(abs(px)>cutoff):   # Change at the other location too, when required
        #print('Singularity time',tcomcurr)
        sing_time=tcomcurr
        singtime_arr.append(sing_time)
        return -1 

singtime_arr = []    

    
def finalcondition(f,tfinal,xinit,yinit,pxinit,pyinit,tcomin):
    sol = ode(f)
    sol.set_integrator('dop853',nsteps=1e16)
    #sol.set_solout(solout)
    
    
    y0 = array([xinit,pxinit,yinit,pyinit,tcomin])
    #print('y0',y0)
    sol.set_initial_value(y0.view(dtype=float),t=0.0)
    sol.set_f_params(pxinit,tcomin,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    sol.integrate(tfinal)
    result = sol.y.view(dtype=complex)
    #print('res',result)
    return result
    

result = finalcondition(f,1.0,xinit,yinit,pxinit,pyinit,0.0)
singtime_arr.append(result[4])

print(singtime_arr)
for i in range(3):
    offset = 2.0
    result = finalcondition(faux,offset,xinit,yinit,pxinit,pyinit,0.0)
    result = finalcondition(f,1.0,result[0],result[2],result[1],result[3],offset)
    # if(real(result[4])>0.0):
    singtime_arr.append(result[4])
    print(result[4],result[0],result[2],'pot',potential(0.0,result[2]))
    print(singtime_arr)
print(singtime_arr)


#print(singtime_arr)
op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
pickle.dump(singtime_arr,op)
op.close()


    
