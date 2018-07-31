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


print(xinit,yinit)
tPrint=0
tPrintSpace=0.0005

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

axarr[0,0].imshow(abs(plotsubmanifold(X,1.0,Y,1j,potential,1)),origin='lower',extent=[xi,xi+lx,xi,xi+lx]).set_clim(vmin,vmax)
axarr[0,1].imshow(abs(plotsubmanifold(X,1.0,Y,1j,potential,-1)),origin='lower',extent=[yi,yi+ly,yi,yi+ly]).set_clim(vmin,vmax)
axarr[0,2].imshow(abs(plotsubmanifold(X,1.0,Y,1.0,potential,0)),origin='lower',extent=[xi,xi+lx,yi,yi+ly]).set_clim(vmin,30)
axarr[1,0].imshow(abs(plotsubmanifold(X,1j,Y,1j,potential,0)),origin='lower',extent=[xi,xi+lx,yi,yi+ly]).set_clim(vmin,vmax)
axarr[1,1].imshow(abs(plotsubmanifold(X,1.0,Y,1j,potential,0)),origin='lower',extent=[xi,xi+lx,yi,yi+ly]).set_clim(vmin,vmax)
axarr[1,2].imshow(abs(plotsubmanifold(X,1j,Y,1.0,potential,0)),origin='lower',extent=[xi,xi+lx,yi,yi+ly]).set_clim(vmin,vmax)

#print('check2')



def f(t,r,dti, momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global currentpx,currentpy, tPrint, tPrintSpace
    #print('r',r)
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
##
##    axarr[0,0].scatter(real(x),imag(x))
##    axarr[0,1].scatter(real(y),imag(y))
##    axarr[0,2].scatter(real(x),real(y))
##    axarr[1,0].scatter(imag(x),imag(y))
##    axarr[1,1].scatter(real(x),imag(y))
##    axarr[1,2].scatter(imag(x),real(y))

    if (t>tPrint+tPrintSpace):
        delp = [momentumx[0]-currentpx,momentumy[0]-currentpy]
        moddelp = (delp[0]**2 + delp[1]**2)**0.5
        if(moddelp>0.25):
            #print('collision at ',x,y,'time',t,'potential',potential(x,y))
            axarr[0,0].scatter(real(x),imag(x),color='r')
            axarr[0,1].scatter(real(y),imag(y),color='r')
            axarr[0,2].scatter(real(x),real(y),color='r')
            axarr[1,0].scatter(imag(x),imag(y),color='r')
            axarr[1,1].scatter(real(x),imag(y),color='r')
            axarr[1,2].scatter(imag(x),real(y),color='r')
            collision[0]+=1
        #print('y at time',t,y)
        else:
            axarr[0,0].scatter(real(x),imag(x))
            axarr[0,1].scatter(real(y),imag(y))
            axarr[0,2].scatter(real(x),real(y))
            axarr[1,0].scatter(imag(x),imag(y))
            axarr[1,1].scatter(real(x),imag(y))
            axarr[1,2].scatter(imag(x),real(y))
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
            collision=[0]
            sol.set_f_params(dti,momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
            #print('F called')
            sol.integrate(abstim)
            tempsol = sol.y
            #print('soldown',sol.y)
        trajdata=sol.y
        print('Trajdata',trajdata[:4])

    return trajdata


def finalcondition(xinit,yinit,pxinit,pyinit):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=100000)
    dti =1.0
            
    y0 = array([xinit,pxinit,yinit,pyinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    momentumx = [3.0]
    momentumy = [0.0]
    collision=[0]
    sol.set_f_params(dti,momentumx,momentumy,collision,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    sol.set_initial_value(y0,t=0.0)
    sol.integrate(3.0)
    print(sol.y[:4])
    trajdata[:-1]=sol.y
    trajdata[20] = collision[0]
    print('trajdata2',trajdata)
    return trajdata
    




##xinit = -5.0 
##yinit = 0.0
##pxinit = 3.0
##pyinit = 0.0

DET=lambda M:M[...,0,0]*M[...,1,1]-M[...,0,1]*M[...,1,0]

gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
trajdata = zeros(21,dtype=complex)

trajdata = finalcondition(xinit,yinit,pxinit,pyinit)

#print('trajdata1',trajdata)


##T1 = [0.0,1.7,1.7-0.1j,3.0-0.1j,3.0]#2.8-0.05j,2.8+0.08j,2.7+0.08j,2.7-0.05j,
##
##T2 = [0.0,1.7,1.7-0.05j,2.2-0.05j,2.2+0.08j,2.1+0.08j,2.1-0.05j,2.2-0.05j,2.2+0.08j,2.1+0.08j,2.1-0.05j,3.0-0.05j,3.0]
##T3 =[0.0,1.0,1.0-0.826j,1.989-0.826j]
##T4 =[0.0,0.5,0.5-0.1509j,2.577-0.1509j]
##T5 = [0.0,2.636,2.636+0.201j]
##T6 = [0.0,2.971,2.971+0.056j]
##T7 = [0.0,2.634,2.634-0.247j]
##T8 = [0.0,2.946,2.946-0.017j]
##T9 = [0.0,2.754,2.754+0.02j]
#T10 = [0.0,2.527,2.527-0.071j]
#trajdata=contour_integration(xinit,yinit,pxinit,pyinit,T10)
##
print('potential at singular point:',potential(0.44445405 -1.50603364j  , -5.38276916 +1.05666108j))



mpp = zeros((2,2),dtype=complex)
mqp = zeros((2,2),dtype=complex)
mpq = zeros((2,2),dtype=complex)
mqq = zeros((2,2),dtype=complex)

mpp=trajdata[4:8].reshape((2,2))
mpq=trajdata[8:12].reshape((2,2))
mqp=trajdata[12:16].reshape((2,2))
mqq=trajdata[16:20].reshape((2,2))

#print(mpp,mpq,mqp,mqq,trajdata[4:8])



xsi = 2*np.matmul(gammaf,mqq) + 2*np.matmul(np.matmul(gammaf,mqp),S20) -1j*mpq - 1j*np.matmul(mpp,S20)
j0 = 0.25/(np.linalg.det(gammaf))*DET(xsi)


print('j0',j0)

#print('Check3')


plt.show()

    
    
    



    
