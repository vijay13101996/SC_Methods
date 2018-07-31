from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import trial_code
from matplotlib.figure import Figure
import pickle, pprint
from matplotlib import pyplot as plt
from Parameter_file_Singularity_time_structure import *
import Complex_plotter


#print(xinit,yinit)
tPrint=0.0
tPrintSpace=0.001

#print(tre,tim)
xi=-20
lx=40
yi=-20
ly=40
div=200
dx=lx/div
dy=ly/div
X = arange(xi,xi+lx,dx)
Y = arange(yi,yi+ly,dy)

x,y = meshgrid(X,Y)

def plotsubmanifold(r1,p1,r2,p2,function,same):
    if(same==0):
        tempr1,tempr2 = meshgrid(r1,r2)
##        plt.imshow(abs(function(tempr1*p1,tempr2*p2)))
##        plt.show()
        #print((tempr1*p1,tempr2*p2))
        return function(tempr1*p1,tempr2*p2)
    elif(same==1):
        #print('here')
        tempr1,tempr2 = meshgrid(r1,r1)
        return function(tempr1*p1 + tempr2*p2,0.0)
    else:
        tempr1,tempr2 = meshgrid(r2,r2)
        return function(0.0,tempr1*p1 + tempr2*p2)
        
        

#print('check1')
##xre = arange(xi,xi+lx,dx)
##xim = arange(xi,xi+lx,dx)
##
##
##
#print (potential(-3.16j,0.0))
##pot = potential(x + y*1j,0.0)
###print(imag(pot))
##plt.imshow(abs(pot))
##plt.show()
#fig = plt.figure(1)
##fig1 = plt.figure(1)   #xr,xi
##fig2 = plt.figure(2)   #yr,yi
##fig3 = plt.figure(3)   #xr,yr
##fig4 = plt.figure(4)   #xi,yi
##fig5 = plt.figure(5)   #xr,yi
##fig6 = plt.figure(6)   #xi,yr

f, axarr = plt.subplots(2,3)

vmin = -10
vmax = 1000



#print('check2')



def f(t,r,dti, momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global currentpx,currentpy, tPrint, tPrintSpace
    #print('r',r)
    #print('tprint',tPrint)
    x,px,y,py, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4 =  r
    ddV1 = ddpotential1(x,y)
    ddV2 = ddpotential2(x,y)
    ddV3 = ddpotential3(x,y)
    ddV4 = ddpotential4(x,y)
    currentpx = px
    currentpy = py
    #print('changed')
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),\
                 -(ddV1*mqp1 + ddV2*mqp3),-(ddV1*mqp2 + ddV2*mqp4),\
                 -(ddV3*mqp1 + ddV4*mqp3),-(ddV3*mqp2 + ddV4*mqp4),\
                 -(ddV1*mqq1 + ddV2*mqq3),-(ddV1*mqq2 + ddV2*mqq4),\
                 -(ddV3*mqq1 + ddV4*mqq3),-(ddV3*mqq2 + ddV4*mqq4),\
                   mpp1,mpp2,mpp3,mpp4, mpq1,mpq2,mpq3,mpq4])*dti

    if (t>tPrint+tPrintSpace):
        #print('here',t,currentpx)
        delp = [momentumx[0]-currentpx,momentumy[0]-currentpy]
        moddelp = (delp[0]**2 + delp[1]**2)**0.5
        if(moddelp>0.25):
            #print('here')
            collision[0]+=1
        tPrint = t
    momentumx[0] = currentpx
    momentumy[0] = currentpy
        
        
    
    
    

    return drdt

def contour_integration(xinit,yinit,pxinit,pyinit,T):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=100000)
            
    y0 = array([xinit,pxinit,yinit,pyinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    #print('y0',y0[0])
    

    if(len(T)==1):
        sol.set_initial_value(y0,t=0.0)
        trajdata = sol.y
    else:
        
        tempsol = y0
        #print(y0)
        for tc in range(len(T)-1):
            
            #print(T[tc+1]-T[tc])
            
            timcom = T[tc+1]-T[tc]
            abstim = abs(timcom)
            dti = timcom/abstim
            sol.set_initial_value(tempsol,t=0.0)
            #tPrint = 0.0
            #print('solup',sol.y)
            momentumx = [3.0]
            momentumy = [0.0]
            sol.set_f_params(dti,momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
            #print('F called')
            sol.integrate(abstim)
            tempsol = sol.y
            #print('soldown',sol.y)
        trajdata=sol.y
        #print('Trajdata',trajdata[:4])

    return trajdata


def finalcondition(xinit,yinit,pxinit,pyinit,timefinal):
    global tPrint,tPrintSpace
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=100000)
    dti =1.0
    #print(xinit,yinit)
            
    y0 = array([xinit,pxinit,yinit,pyinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    momentumx = [pxinit]
    momentumy = [pyinit]
    collision=[0]
    tPrint = -tPrintSpace
    #print('tPrint,here',tPrint)
    
    sol.set_f_params(dti,momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    sol.set_initial_value(y0,t=0.0)
    sol.integrate(timefinal)
    #print(sol.y[:2])
    trajdata[:-1]=sol.y
    trajdata[20] = collision[0]
    #print('trajdata2',trajdata)
    return trajdata
    


DET=lambda M:M[...,0,0]*M[...,1,1]-M[...,0,1]*M[...,1,0]

gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
trajdata = zeros(21,dtype=complex)

collisionarray = zeros(len(pointarr),dtype=int)
collisionarrayindex = []
print('len',len(pointarr))
count=0
for i in range(len(pointarr)):
    xinit = pointarr[i][0]
    yinit = pointarr[i][1]
    pxinit = -1j*(2*gammaxi*x0 + 1j*px0 - 2*gammaxi*xinit)
    pyinit = -1j*(2*gammayi*y0 + 1j*py0 - 2*gammayi*yinit)
    trajdata = finalcondition(xinit,yinit,pxinit,pyinit,3.0)
    collisionarray[i] = trajdata[20]
    if(collisionarray[i]!=0):
        
        collisionarrayindex.append(pointarr[i][2])
        print('i',i,pointarr[i][2])
        #print(collisionarrayindex)
        count+=1

print('count,len',count,len(collisionarrayindex))


op = open('collision_array_20.pkl','wb')
pickle.dump(collisionarray,op)
op.close()

op = open('collision_array_index_20.pkl','wb')
pickle.dump(collisionarrayindex,op)
op.close()






    
    
    



    
