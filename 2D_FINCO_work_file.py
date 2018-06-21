# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 18:46:34 2017

@author: Dell
"""

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
import Complex_plotter
#import finalconditions as fc
from Parameter_file_double_slit_eckart_10_10_10_10 import *


plotcompare=1
plotsection=0
print(Nxr,Nxi,Nyr,Nyi)
#div=256
#----------------------------------------------------------

def xsection(wftemp,ywf,xcut):
   for i in range(len(ywf)):
        ywf[i] = wftemp[i][xcut]#*(pi/2*gammaxi)**0.25
     

def ysection(wftemp,xwf,ycut):
    sum = 0.0
    for i in range(len(xwf)):

        xwf[i] = wftemp[int(ycut)][i]#*(pi/(2*gammayi))**0.25#*exp(1j*0.5*t)
        sum += abs(xwf[i])**2*dx#*(pi)**0.5
    #print(sum)

def printfun(x,wf,outFile):
    for ax,w in zip(x,wf):
        print (real(ax),real(w),imag(w),file=outFile)
    print ("",file=outFile)
    
#---------------------------------------------------------

start_time = time.time()

xbar =  zeros(Nxr*Nxi*Nyr*Nyi) + 0j
ybar = zeros_like(xbar) + 0j
pxbar = zeros_like(xbar) + 0j
pybar = zeros_like(ybar) + 0j

#Interference_further_modified_trimmed.ketmanifold(x0,y0, px0,py0, xbar, pxbar, ybar, pybar, gammaxi, gammayi,Nxr,Nxi,Nyr,Nyi,widthx,widthy)


op = open('2D_FINCO_ketmanifold_{}_{}_{}_{}.pkl'.format(Nxr,Nxi,Nyr,Nyi),'rb')
xbar = pickle.load(op)
pxbar = pickle.load(op)
ybar = pickle.load(op)
pybar = pickle.load(op)

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
data = zeros((Nxr*Nxi*Nyr*Nyi,len(T),28)) + 0j
 

X = arange(xi,xi+lx,dx)
Y = arange(yi,yi+ly,dy)

x,y = meshgrid(X,Y)

wf = zeros_like(x) + 0j

print("--- %s seconds ---" % (time.time() - start_time))

# Code modified to 'pickle'

for tc in range(len(T)):

        op = open('2D_FINCO_{}_trajdata_{}_{}_{}_{}_at_time_{}.pkl'.format(potkey,Nxr,Nxi,Nyr,Nyi,tc), 'rb')
        
        #wf = pickle.load(op)

        mpp = pickle.load(op)
        mpq = pickle.load(op)
        mqp = pickle.load(op)
        mqq = pickle.load(op)

        xfinbar = pickle.load(op)
        yfinbar = pickle.load(op)
        pxfinbar = pickle.load(op)
        pyfinbar = pickle.load(op)

        Sfinbar = pickle.load(op)
        Turningpoints = pickle.load(op)

        xfin = pickle.load(op)
        yfin = pickle.load(op)
        pxfin = pickle.load(op)
        pyfin = pickle.load(op)

        magnitude = pickle.load(op)
        phaseangle = pickle.load(op)

        sigma = pickle.load(op)
        J0 = pickle.load(op)
        successcheck =  pickle.load(op)
        maxmomentum  = pickle.load(op)
        

        op.close()

        print("--- %s seconds ---" % (time.time() - start_time))


#================================================================================



print(maxmomentum[0])
print(Turningpoints)


count=0
for i in range(len(successcheck)):
        if (successcheck[i]!=-1):
                count+=1

print('count',count)

print("--- %s seconds ---" % (time.time() - start_time))


#A Plot of the real trajectory propagation with time

xwf = zeros_like(X) + 0j
ywf = zeros_like(Y) + 0j

print(len(wf))
                
##plt.figure(23)
##plt.imshow(abs(wf)) 
##plt.show()

ycut = 0.65
xcut = 0.5


for i in range(Nxr*Nxi*Nyr*Nyi):
      if ( maxmomentum[i][2]<100):
      #print(sigma[i])
      
         dpxdpy= (widthx*widthy)/(Nxi*Nyi)        # 22/05/18 : Watch out for gamma here, the differential changes....
         dxdy = (widthx*widthy)/(Nxr*Nyr)

         j0 = 0.25/(np.linalg.det(gammaf))*np.linalg.det(2*np.dot(gammaf,mqq[i]) + 2*np.dot(np.dot(gammaf,mqp[i]),S20)\
                          -1j*mpq[i] - 1j*np.dot(mpp[i],S20))
           
         J0[i] = j0
         
         J = -(0.25/np.linalg.det(gammaf))*abs(j0)**2 + 0j
         Integrand = -(-0.25)**2*(4/(pi**2*np.linalg.det(gammaf)))**0.75*J/(j0**0.5)
         
         xgf =  (2*gammaxf/pi)**0.25*exp(-gammaxf*(X[ycut*div]-xfin[i])**2 + 1j*pxfin[i]*(X[ycut*div]-xfin[i]))
         ygf =  (2*gammayf/pi)**0.25*exp(-gammayf*(Y-yfin[i])**2 + 1j*pyfin[i]*(Y-yfin[i]))

         ExponentS = 1j*Sfinbar[i]
         Exponentx = -gammaxf*xfinbar[i]**2 + gammaxf*xfin[i]**2 +\
                      1j*pxfinbar[i]*xfin[i] - 1j*pxfin[i]*xfinbar[i]
         Exponenty = -gammayf*yfinbar[i]**2 + gammayf*yfin[i]**2 +\
                      1j*pyfinbar[i]*yfin[i] - 1j*pyfin[i]*yfinbar[i]
         sigma[i] = Exponentx + Exponenty + ExponentS
         magnitude[i] = abs(Integrand*exp(sigma[i]))
         phaseangle[i] = np.angle(Integrand*exp(sigma[i]))

   ##      if(real(xbar[i])>-4.0 and ybar[i]==0.0 and magnitude[i]>0.0):
   ##          print(magnitude[i])

         # and magnitude[i]>20.0 and magnitude[i]<50.0 and successcheck[i]!=-1 #real(sigma[i])<0.0 and
         # and magnitude[i]>0.0001 and magnitude[i]<50.0 and real(xfinbar[i])>0.0 
          

         if(real(sigma[i])<0.0 and magnitude[i]<50.0 ): # Here the condition for whether magnitude[i] is negative(i.e filtered) is to be taken care.
         
            #if(i not in collisionarrayindex):
            #print('i',i)
               ywf+= Integrand*exp(sigma[i])*xgf*ygf*dxdy*dpxdpy*exp(1j*pi*Turningpoints[i])
               divergearr.append(i)

print('count',count)


pointarr = []
for i in divergearr:
   #print(i,xbar[i],ybar[i])
   pointarr.append([xbar[i],ybar[i],i])

print('divergearr',len(divergearr))

#print(divergearr)
#print(pointarr)

op = open('pointarr_N_20.pkl','wb')
pickle.dump(pointarr,op)
op.close()

print(X[int(ycut*div)],Y[int(xcut*div)]) # : To check for the origin.

if(plotsection==1):

   ysection(wf,xwf,int(xcut*div))
   xsection(wf,ywf,int(ycut*div))

   #print("The Normalizing factor of the wavefunction is:")
   #normal(X,Y,wf)


outFilex = open('2DFINCO_{}_{}_{}_{}_{}_X.txt'.format(potkey,Nxr,Nxi,Nyr,Nyi),'w')
printfun(X,xwf,outFilex)
outFilex.close()

outFiley = open('2DFINCO_{}_{}_{}_{}_{}_Y.txt'.format(potkey,Nxr,Nxi,Nyr,Nyi),'w')
printfun(Y,ywf,outFiley)
outFiley.close()

plotvalue = magnitude*exp(1j*phaseangle)

#print(plotvalue)
submani1 = FINCO_2D.xsubmanifold(plotvalue,Nxr,Nxi,Nyr,Nyi,successcheck)
submani1[real(submani1)>50] = -0.1
##submani1[abs(submani1)<0] = -0.1

##plt.imshow(abs(submani1),origin='lower')
##plt.show()

#Complex_plotter.plotcomplex(submani1,0.1,0.1,x0-widthx/2,x0+widthx/2,-widthx/2,widthx/2)
##fig = plt.figure(3)
##fig.suptitle('xr-xi submanifold', fontsize=14, fontweight='bold')
##ax = fig.add_subplot(111)
##ax.set_xlabel('xreal')
##ax.set_ylabel('ximag')
##cax = ax.imshow(real(submani1), origin = 'lower', extent=[x0-widthx/2,x0+widthx/2,-widthx/2,+widthx/2],interpolation='none')
##cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
##cbar.ax.set_yticklabels(['< -1', '0', '> 1'])
##plt.ioff()
##plt.show()


##plt.suptitle('QM Wavefunction ', fontsize=14, fontweight='bold')
##plt.imshow(abs(wf),origin = 'lower',zorder=0,extent=(xi,xi+lx,yi,yi+ly))#,vmax=0.01)
##plt.plot(X[int(div*ycut)]*ones_like(X),Y, '--', linewidth=1, color='firebrick')
##plt.plot(X,Y[int(div*xcut)]*ones_like(Y), '--', linewidth=1, color='firebrick')
##plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

# Comparing split operator wavefunction and FINCO wavefunction
if (plotcompare ==1):
   for col,fname,m in zip(['blue','red'],['ysection_{}.txt'.format('custom_double_slit_eckart'), '2DFINCO_{}_{}_{}_{}_{}_Y.txt'.format(potkey,Nxr,Nxi,Nyr,Nyi)],['x','.']):
           data=np.loadtxt(fname)
           X=data[:,0]
           Y=data[:,1]
           Z= (data[:,2]**2 + data[:,1]**2)**0.5
           fig = plt.figure(7)
           fig.suptitle('Comparison between Quantum Mechanical and FINCO wavefunctions', fontsize=14, fontweight='bold')
           ax = fig.add_subplot(111)
           #fig.subplots_adjust(top=0.85)
           ax.set_title('Psi(y) vs y')
           ax.set_xlabel('y')
           ax.set_ylabel('psi(y)')
           #plt.plot(X,Y,m,color=col,linestyle='-')
           #plt.draw()
           plt.plot(X,Z,m,color=col,linestyle='-')

   plt.show()
          




