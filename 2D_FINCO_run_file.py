from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
#from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import FINCO_2D
#import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
from datetime import datetime
import time
import pickle, pprint
#import finalconditions as fc
from Parameter_file_double_slit_eckart_10_10_10_10 import *

start_time = time.time()

xbar =  zeros(Nxr*Nxi*Nyr*Nyi) + 0j
ybar = zeros_like(xbar) + 0j
pxbar = zeros_like(xbar) + 0j
pybar = zeros_like(ybar) + 0j

FINCO_2D.ketmanifold(x0,y0, px0,py0, xbar, pxbar, ybar, pybar, gammaxi, gammayi,Nxr,Nxi,Nyr,Nyi,widthx,widthy)


op = open('2D_FINCO_ketmanifold_{}_{}_{}_{}.pkl'.format(Nxr,Nxi,Nyr,Nyi),'wb')
pickle.dump(xbar, op)
pickle.dump(pxbar, op)
pickle.dump(ybar,op)
pickle.dump(pybar,op)

op.close()


print("--- %s seconds ---" % (time.time() - start_time))

xfinbar = zeros_like(xbar) + 0j
yfinbar = zeros_like(ybar) + 0j
pxfinbar = zeros_like(pxbar) + 0j
pyfinbar = zeros_like(pybar) + 0j
Sfinbar = zeros_like(xbar) + 0j

pxfin = zeros_like(pxbar)
xfin = zeros_like(xbar)
pyfin = zeros_like(pybar)
yfin = zeros_like(ybar)
Turningpoints = zeros_like(xbar)

mpp = zeros((Nxr*Nxi*Nyr*Nyi,2,2)) + 0j
mpq = zeros((Nxr*Nxi*Nyr*Nyi,2,2)) + 0j
mqp = zeros((Nxr*Nxi*Nyr*Nyi,2,2)) + 0j
mqq = zeros((Nxr*Nxi*Nyr*Nyi,2,2)) + 0j

gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]

magnitude = zeros_like(xbar)
phaseangle = zeros_like(xbar)

sigma = zeros_like(xbar) + 0j
J0  = zeros_like(xbar) + 0j

successcheck = ones_like(xbar,dtype=int)

maxmomentum =  zeros((Nxr*Nxi*Nyr*Nyi,3),dtype=complex) 

S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])

T = arange(0,t,dt)
data = zeros((Nxr*Nxi*Nyr*Nyi,len(T),22),dtype=complex) 

X = arange(xi,xi+lx,dx)
Y = arange(yi,yi+ly,dy)

x,y = meshgrid(X,Y)

wf = zeros_like(x) + 0j

FINCO_2D.finalconditions(xbar,pxbar,ybar, pybar, x0, px0, y0, py0, gammaxi,gammayi,\
                    gammaxf,gammayf,data,Nxr,Nxi,Nyr,Nyi,potential,dxpotential,dypotential,\
                    ddpotential1, ddpotential2, ddpotential3, ddpotential4, timear,gammaf,gammai,mpp,mpq,mqp,mqq,S20,successcheck,maxmomentum)

#print('xfinbar',xfinbar)

print("--- %s seconds ---" % (time.time() - start_time))

# Code modified to 'pickle'

for tc in range(len(T)):

        xfinbar = data[:,tc,0]
        pxfinbar = data[:,tc,1]
        yfinbar = data[:,tc,2]
        pyfinbar = data[:,tc,3]
        Sfinbar = data[:,tc,4]
        
        mpp[:,0,0] = data[:,tc,5]
        mpp[:,0,1]= data[:,tc,6]
        mpp[:,1,0]= data[:,tc,7]
        mpp[:,1,1]= data[:,tc,8]
        mpq[:,0,0]= data[:,tc,9]
        mpq[:,0,1]= data[:,tc,10]
        mpq[:,1,0]= data[:,tc,11]
        mpq[:,1,1]= data[:,tc,12]
        mqp[:,0,0]= data[:,tc,13]
        mqp[:,0,1]= data[:,tc,14]
        mqp[:,1,0]= data[:,tc,15]
        mqp[:,1,1]= data[:,tc,16]
        mqq[:,0,0]= data[:,tc,17]
        mqq[:,0,1]= data[:,tc,18]
        mqq[:,1,0]= data[:,tc,19]
        mqq[:,1,1]= data[:,tc,20]

        Turningpoints = data[:,tc,21]
                
        FINCO_2D.bramanifold(xfinbar, pxfinbar, yfinbar, pyfinbar, xfin, pxfin, yfin, pyfin, gammaxf, gammayf, Nxr,Nxi,Nyr,Nyi)

        wf = wf*(0.0 + 0j)

##FINCO_2D.FINCO(xfin, yfin, pxfin, pyfin, xfinbar, pxfinbar, yfinbar, pyfinbar,\
##          gammaf,gammaxf,gammayf, Sfinbar, turningpoints, gammaxi, gammayi,\
##          mpp,mpq,mqp,mqq, Nxr,Nxi,Nyr,Nyi, widthx,widthy, x,y, wf,S20, magnitude, phaseangle,sigma,J0,successcheck)

        op = open('2D_FINCO_{}_trajdata_{}_{}_{}_{}_at_time_{}.pkl'.format(potkey,Nxr,Nxi,Nyr,Nyi,tc), 'wb')
        
        #pickle.dump(wf, op)

        pickle.dump(mpp, op)
        pickle.dump(mpq, op)
        pickle.dump(mqp, op)
        pickle.dump(mqq, op)

        pickle.dump(xfinbar,op)
        pickle.dump(yfinbar,op)
        pickle.dump(pxfinbar,op)
        pickle.dump(pyfinbar,op)

        pickle.dump(Sfinbar,op)
        pickle.dump(Turningpoints,op)

        pickle.dump(xfin, op)
        pickle.dump(yfin, op)
        pickle.dump(pxfin, op)
        pickle.dump(pyfin, op)

        pickle.dump(magnitude, op)
        pickle.dump(phaseangle, op)

        pickle.dump(sigma,op)
        pickle.dump(J0,op)
        pickle.dump(successcheck,op)
        pickle.dump(maxmomentum,op)

        op.close()

        print("--- %s seconds ---" % (time.time() - start_time))


#================================================================================


#print(maxmomentum)

print("--- %s seconds ---" % (time.time() - start_time))


#A Plot of the real trajectory propagation with time


print("--- %s seconds ---" % (time.time() - start_time))

# Comparing split operator wavefunction and FINCO wavefunction





