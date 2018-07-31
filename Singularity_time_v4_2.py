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
import random
import warnings
from matplotlib import pyplot as plt

from Parameter_file_Singularity_time_structure import *

warnings.filterwarnings('error')
print(xinit,yinit)
print(pointarr)

sing_time=0.0
cutoff = 1e3
timear = []

previous =1.0
n=0.5
digress = 0 

phasejumparr=[]
dtim=0.0

xtarget=0.0j
def faux(t,r,digress,theta,pxinit,tinit,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, previous
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    dtim = exp(1j*theta)
    
    # if(t>1e4):
        # print('exceeding')
    # if(abs(t-1e5)<1e-3):
        # print('time past')
    
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtim
    #timear.append(tcomcurr)
    #print(potential(x,y),dxpotential(x,y))
    return drdt.view(dtype=float)


def fq(t,r,digress,pxinit,tinit,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global previous,phasejumparr,dtim
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    #dtim = 0.0j
    dpdt = -1*(dxpotential(x,y))
    
    
       
    
     
    # if(len(singtime_arr)==0):
        # dtim=1/dpdt
        # print('tcomcurr',tcomcurr)
    # else:
        
        # for singtimes in singtime_arr:
            # exprt.append(tcomcurr-singtimes)
            
      
        # term1 = 1.0+0j
        # count=0
        # print('new loop')
        # for ext in exprt:
            # term1*=(ext)**n
            
        # for singtimes in singtime_arr:
            # currentimag = float("%.12f"%imag(tcomcurr-singtimes))
            # current  = real(tcomcurr-singtimes) + 1j*currentimag
            # if(real(current)<0.0 and imag(previous[count])*imag(current)<0.0 and abs(abs(current)-abs(previous[count]))>1e-8):
                # print('current,previous',current,previous[count],'phasejump')
                # phasejumparr[count]=1
            # else:
                # phasejumparr[count]=0
            # previous[count]=current
            # count+=1
          
            
        # print('term1 without phasejump', term1)
        # term1 = 1.0+0j
        # count=0
        # for ext in exprt:
            # term1*=(ext**n)*exp(1j*pi*phasejumparr[count])
            # print('ext',ext)
            # print('phasejump added',phasejumparr[count],exp(1j*pi*phasejumparr[count]))
            # count+=1
            
        # term2 = dpdt
        # count=0
        # for ext in exprt:
            # term2+=n*px/ext
            # count+=1
        
        # dtim = 1/(term1*term2)
        
        # if(len(singtime_arr)==3):
        # print('term1',term1)#,potential(x,y),dxpotential(x,y))
                
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),1.0])*dtim    
    return drdt.view(dtype=float)

def soloutq(t,r):
    global sing_time,cutoff,digress,previous,phasejumparr,dtim,xtarget
    r=r.view(dtype=complex)
    x,px,y,py,tcomcurr =  r
    #timear.append(tcomcurr)
    
    count=0
    
    for singtimes in singtime_arr:
        currentimag = float("%.12f"%imag(tcomcurr-singtimes))
        current  = real(tcomcurr-singtimes) + 1j*currentimag
        #print('current,previous',current,previous[count])
        if(real(current)<0.0 and imag(previous[count])*imag(current)<0.0):
            phasejumparr[count]+=1
            #print('current,previous',current,previous[count],'phasejump',phasejumparr[count])
            
        # else:
            # phasejumparr[count]=0
        previous[count]=current
        
        count+=1
        
    exprt=[]
    dpdt = -1*(dxpotential(x,y))
        
    if(len(singtime_arr)==0):
        dtim=1/dpdt
        #print('tcomcurr',tcomcurr)
    else:
        
        for singtimes in singtime_arr:
            exprt.append(tcomcurr-singtimes)
            
      
        
        #print('new loop')
        # for ext in exprt:
            # term1*=(ext)**n
            
        # for singtimes in singtime_arr:
            # currentimag = float("%.12f"%imag(tcomcurr-singtimes))
            # current  = real(tcomcurr-singtimes) + 1j*currentimag
            # if(real(current)<0.0 and imag(previous[count])*imag(current)<0.0 and abs(abs(current)-abs(previous[count]))>1e-8):
                # print('current,previous',current,previous[count],'phasejump')
                # phasejumparr[count]=1
            # else:
                # phasejumparr[count]=0
            # previous[count]=current
            # count+=1
          
            
        #print('term1 without phasejump', term1)
        term1 = (x-xtarget)**2 
        count=0
        for ext in exprt:
            term1/=(ext**n)*exp(1j*pi*phasejumparr[count])
            #print('ext',ext)
            #print('phasejump added',phasejumparr[count],exp(1j*pi*phasejumparr[count]))
            count+=1
            
        term2 = -px
        count=0
        for ext in exprt:
            term2+=n*(x-xtarget)/ext
            count+=1
        
        dtim = term1/term2
        
        #if(len(singtime_arr)==3):
        #print('term1',term1)#,potential(x,y),dxpotential(x,y))
          
    if(abs(1/(x-xtarget))>cutoff):   
        return -1 
        
    if(real(tcomcurr)<-0.5):
        digress = 1
        #print('digress',1,'tcomcurr',tcomcurr)
        #print(potential(x,y),dxpotential(x,y))
        return -1
        
    if(real(tcomcurr)>4.0):
        digress = 2
        #print('digress',2)
        #print('digress',2,'tcomcurr',tcomcurr)
        return -1
        
    if(imag(tcomcurr)>2.0):
        digress = 3
        #print('digress',3)
        #print('digress',3,'tcomcurr',tcomcurr)
        return -1
        
    if(imag(tcomcurr)<-2.0):
        digress = 4
        #print('digress',4)
        #print('digress',4,'tcomcurr',tcomcurr)
        return -1

def solout1(t,r):
    global digress

    r = r.view(dtype=complex)
    x,px,y,py,tcomcurr=r
    #timear.append(tcomcurr)
    #print('inside solout1',x,y,tcomcurr,potential(x,y),dxpotential(x,y))
   
    if(real(tcomcurr)>0.5 and abs(imag(tcomcurr))<=2.0):
        #print('tcomcurr inside solout1',tcomcurr)
        digress =0
        return -1
    
    elif(real(tcomcurr)>0.5 and imag(tcomcurr)>2.0):
        #print('tcomcurr inside solout1',tcomcurr,'digress',3)
        digress=3
        return -1
     
    elif(real(tcomcurr)>0.5 and imag(tcomcurr)<-2.0):
        #print('tcomcurr inside solout1',tcomcurr,'digress',4)
        digress=4
        return -1
        
def solout2(t,r):
    global digress
    r = r.view(dtype=complex)
    x,px,y,py,tcomcurr=r
    #timear.append(tcomcurr)
    if(real(tcomcurr)<3.5  and abs(imag(tcomcurr))<2.0):
        #print('tcomcurr inside solout2',tcomcurr)
        digress=0
        return -1
    elif(real(tcomcurr)<3.5 and imag(tcomcurr)>2.0):
        #print('tcomcurr inside solout2',tcomcurr,'digress',3)
        digress=3
        return -1
     
    elif(real(tcomcurr)<3.5 and imag(tcomcurr)<-2.0):
        #print('tcomcurr inside solout2',tcomcurr,'digress',4)
        digress=4
        return -1
        
def solout3(t,r):
    global digress
    r = r.view(dtype=complex)
    x,px,y,py,tcomcurr=r
    #timear.append(tcomcurr)
    if(imag(tcomcurr)<1.5 and real(tcomcurr)<4.0 and real(tcomcurr)>0.0):
        #print('tcomcurr inside solout3',tcomcurr)
        digress=0
        return -1
    elif(imag(tcomcurr)<1.5 and real(tcomcurr)<4.0):
        #print('tcomcurr inside solout3',tcomcurr,'digress',1)
        digress=1
        return -1
    elif(imag(tcomcurr)<1.5 and real(tcomcurr)>0.0):
        #print('tcomcurr inside solout3',tcomcurr,'digress',2)
        digress=2
        return -1
        
def solout4(t,r):
    global digress
    r = r.view(dtype=complex)
    x,px,y,py,tcomcurr=r
    #timear.append(tcomcurr)
    if(imag(tcomcurr)>-1.5 and real(tcomcurr)<4.0 and real(tcomcurr)>0.0):
        #print('tcomcurr inside solout4',tcomcurr)
        digress=0
        return -1
    elif(imag(tcomcurr)>-1.5 and real(tcomcurr)<4.0):
        #print('tcomcurr inside solout4',tcomcurr,'digress',1)
        digress=1
        return -1
    elif(imag(tcomcurr)>-1.5 and real(tcomcurr)>0.0):
        #print('tcomcurr inside solout1',tcomcurr,'digress',2)
        digress=2
        return -1
        
singtime_arr = []    

def finalcondition(fq,tfinal,xinit,yinit,pxinit,pyinit,tcomin):
    global digress,previous,phasejumparr
    sol = ode(fq)
    sol.set_integrator('dop853',nsteps=1e6)
    sol.set_solout(soloutq)
    
    
    y0 = array([xinit,pxinit,yinit,pyinit,tcomin])
    sol.set_initial_value(y0.view(dtype=float),t=0.0)
    sol.set_f_params(digress,pxinit,tcomin,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    #previous = 1.0
    previous= tcomin*ones(len(singtime_arr))-np.array(singtime_arr)
    phasejumparr = zeros(len(singtime_arr),dtype=int)
    
    #print(previous)
    sol.integrate(tfinal)
    result = sol.y.view(dtype=complex)
    return result
    
def finalconditionaux(faux,digress,xinit,yinit,pxinit,pyinit,tcomin):
    
    sol = ode(faux)
    sol.set_integrator('dop853',nsteps=1e6)
    theta=0.0
    faceangle = pi/4 +pi/16
    if(digress==1):
        sol.set_solout(solout1)
        theta = random.uniform(0.0-faceangle, 0.0+faceangle)
    elif(digress==2):
        sol.set_solout(solout2)
        theta = random.uniform(pi-faceangle, pi+faceangle)
    elif(digress==3):
        sol.set_solout(solout3)
        theta = random.uniform(3*pi/2-faceangle, 3*pi/2+faceangle)
    elif(digress==4):
        sol.set_solout(solout4)
        theta = random.uniform(pi/2 - faceangle, pi/2 + faceangle)
        
    
    y0 = array([xinit,pxinit,yinit,pyinit,tcomin])
    sol.set_initial_value(y0.view(dtype=float),t=0.0)
    
    sol.set_f_params(digress,theta,pxinit,tcomin, potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    previous = 1.0
    sol.integrate(1e5)
    result = sol.y.view(dtype=complex)
    return result

singpoint_arr = []
result = finalcondition(fq,1e3,xinit,yinit,pxinit,pyinit,0.0)
singtime_arr.append(result[4])
print('singularity point in space',result[0])

# offset=0.0
print(singtime_arr)
xtarget_list = [(0.0+1.5709025174582523j),(0.0-1.5667183672681677j),(0.0-4.7067885713578885j),(0.0+4.7067885713578885j) ]

#xtarget = result[0]
for xtrgt in [xtarget_list[0]]:
    #global xtarget
    xtarget = xtrgt
    print('xtarget',xtarget)
    for i in range(1):
        try:
            print('i',i)
            digress = 0
            result = finalcondition(fq,1e4,xinit,yinit,pxinit,pyinit,0.0)
            if(digress == 0):
                print('digress is zero,appending')
                if(real(result[4])>0.0 and abs(potential(result[0],result[2]))>1e3 and abs(result[0]-xtarget)<1e-3):
                    singtime_arr.append(result[4])
                    print('diff',abs(result[0]-xtarget), 'for point',result[4])
                    print('potential',potential(result[0],result[2]),result[0])
                # op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
                # pickle.dump(timear,op)
                # op.close()
            else:
                
                while(digress!=0):
                        print('digress is not zero,processing,new call',digress)
                        
                        # if(i==5):
                            # op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
                            # pickle.dump(timear,op)
                            # op.close()
                        result =  finalconditionaux(faux,digress,result[0],result[2],result[1],result[3],result[4])
                        
                result = finalcondition(fq,1e10,result[0],result[2],result[1],result[3],result[4])
                print('time',result,abs(result[0]-xtarget))
                    
                
                if(real(result[4])>0.0 and abs(potential(result[0],result[2]))>1e3 and abs(result[0]-xtarget)<1e-3):# and real(result[4])<=4.0 and abs(imag(result[4]))<=2.0):
                    singtime_arr.append(result[4])
                    print('diff',abs(result[0]-xtarget), 'for point',result[4])
                    #print('potential',potential(result[0],result[2]),result[4])
        except RuntimeWarning:
            continue
            
        

    #print(singtime_arr,'new')
       
    
    
print(singtime_arr)


op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
pickle.dump(singtime_arr,op)
op.close()


    
