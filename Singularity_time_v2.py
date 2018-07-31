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
t1 =(1.7600502483227527+0.2454058275120668j)

previous =1.0
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
    #print('t',t,'tcomcurr',tcomcurr)#,'t^2px',px)
    #print('t',t,'px',px)#,tcomcurr,y)
    
    return drdt.view(dtype=float)

def f(t,r,pxinit, tinit,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, previous,tPrint, tPrintSpace
    #print('usual')
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    dtim = 0.0j
    expr = tcomcurr-t1
    dpdt = -1*(dxpotential(x,y))#
    exprt=[]
    if(1):#real(tcomcurr)>=0.0):# or abs(real(tcomcurr)-0.0)<1e-12):
        if(len(singtime_arr)==0):
            dtim=1/dpdt#1/(n*expr**(n-1)*px +dpdt*expr*n)
            #print('expr',expr**(n-1))
        else:
            
            for singtimes in singtime_arr:
                #print('tcomcurr',tcomcurr)
                exprt.append(tcomcurr-singtimes)
                
                #expri.append(tinit-singtimes)
            
            term1 = 1.0 + 0j
            
            
            count=0
            #print('previous',previous,'exprt',exprt)
            for ext in exprt:
                #print('ext',ext)
                #print('exi',exi)
                #term1*=ext**n
                if(real(ext)<0.0):# and imag(previous[count])*imag(ext)<1e-5):
                    #print('here,phasejump')
                    term1*=ext**n
                else:
                    term1*=ext**n
                count+=1
                
            #previous=exprt

            term2 = dpdt
            for ext in exprt:
                term2+=n*px/ext
                
            if(term1==-previous):
                print('hereh')
            previous = term1
            
            print('t',tcomcurr,'dVx',dxpotential(x,y),'x cosh',x,cosh(x))#-singtime_arr[len(singtime_arr)-1])
            #print(dxpotential(x,y))
                
            dtim = 1/(term1*term2)#(abs(real(term1)+1j*imag(term1))*term2)
            #print('dtim',1/term1*term2,'px',px)
            #print('term1',term1,previous)
            if(real(tcomcurr)>4.0):
                print('digress')
                dtim = exp(1j*(pi+pi/6))
                #dtim=-dtim
                # #print('negative time','tcomcurr',tcomcurr)
            
            
                
            
    else:
        
        dtim = conj(1/dpdt)
        #print('negative time','tcomcurr',tcomcurr,'dtim',dtim)
    #print('dtim',dtim)
        
    #print('pot',potential(x,y),'dpot',dxpotential(x,y))
    #print('x',x,'px',px,'y',y,'py',py)
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtim
    #timear.append(tcomcurr)
    #print('t',t,'tcomcurr',tcomcurr)#,'t^2px',px)
    #print(t,1/px,x,tcomcurr)
    
    return drdt.view(dtype=float)

def solout(t,r):
    global sing_time,cutoff
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    timear.append(tcomcurr)
    #expr = tcomcurr-t1
    
    if(abs(px)>cutoff):   # Change at the other location too, when required
        #print('Singularity time',tcomcurr)
        #sing_time=tcomcurr
        #singtime_arr.append(sing_time)
        #print('here')
        return -1 

singtime_arr = []    


    
def finalcondition(f,tfinal,xinit,yinit,pxinit,pyinit,tcomin):
    global previous
    sol = ode(f)
    sol.set_integrator('dop853',nsteps=1e6)
    sol.set_solout(solout)
    
    
    y0 = array([xinit,pxinit,yinit,pyinit,tcomin])
    #print('y0',y0)
    sol.set_initial_value(y0.view(dtype=float),t=0.0)
    sol.set_f_params(pxinit,tcomin,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    previous = 1.0#0.0*np.array(singtime_arr)
    #print('previous,start',previous)
    sol.integrate(tfinal)
    result = sol.y.view(dtype=complex)
    #print('res',result)
    return result
    

result = finalcondition(f,1e10,xinit,yinit,pxinit,pyinit,0.0)
singtime_arr.append(result[4])
#print(singtime_arr)

offset=0.0

for i in range(3):
    
        offset = 1.65
        print('new')
        result = finalcondition(faux,offset,xinit,yinit,pxinit,pyinit,0.0)
        result = finalcondition(f,1e4,result[0],result[2],result[1],result[3],offset)
        #if(real(result[4])>0.0):
        singtime_arr.append(result[4])
            #offset+=0.5
        #print(result[4],result[0],result[2],'pot',potential(result[0],result[2]))
        #print(result)
        #print(singtime_arr,'new')
       
    
    
print(singtime_arr)



#print(singtime_arr)
op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
pickle.dump(singtime_arr,op)
op.close()


    
