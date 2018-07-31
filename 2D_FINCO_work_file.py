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
#from Parameter_file_FINCO_double_slit_eckart import *
potkey='diffraction'#'reaction_coordinate'###'double_slit_tunneling'     #Code:(2,1,3) for diffraction, (3,3,3) for gaspard
import Parameters_File
Parameters_File.method_sel('FINCO')
Parameters_File.wf_init_set('FINCO')
Parameters_File.initial_conditions(2)
Parameters_File.grid_select(1)
Parameters_File.time_grid(3)
Parameters_File.pot_select(potkey)
from Parameters_File import *
from pylab import *

plotcompare=0
plotsection=0
print(Nxr,Nxi,Nyr,Nyi)
#div=256
print(t)
print('Vc',Vc)
print('T',T)

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
    


def click(event):
   """If the left mouse button is pressed: draw a little square. """
   tb = get_current_fig_manager().toolbar
   if event.button==1 and event.inaxes and tb.mode == '':
       x,y = event.xdata,event.ydata
       plot([x],[y],'rs')
       draw()
    
#---------------------------------------------------------

start_time = time.time()

xbar =  zeros(Nxr*Nxi*Nyr*Nyi) + 0j
ybar = zeros_like(xbar) + 0j
pxbar = zeros_like(xbar) + 0j
pybar = zeros_like(ybar) + 0j

#Interference_further_modified_trimmed.ketmanifold(x0,y0, px0,py0, xbar, pxbar, ybar, pybar, gammaxi, gammayi,Nxr,Nxi,Nyr,Nyi,widthx,widthy)


op = open('2D_FINCO_ketmanifold_{}_{}_{}_{}_width_{}_{}_{}_{}.pkl'.format(Nxr,Nxi,Nyr,Nyi,widthxr,widthxi,widthyr,widthyi),'rb')
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

#T = arange(0,t,dt)
data = zeros((Nxr*Nxi*Nyr*Nyi,len(T),28)) + 0j
 

X = arange(xi,xi+lx,dx)
Y = arange(yi,yi+ly,dy)

x,y = meshgrid(X,Y)

wf = zeros_like(x) + 0j

print("--- %s seconds ---" % (time.time() - start_time))

# Code modified to 'pickle'
print('T',T)
for tc in range(len(T)):
        print('t',T[tc])
        op = open('/home/vijay/Codes/Pickle_files/2D_FINCO_{}_trajdata_{}_{}_{}_{}_tfinal_{}_width_{}_{}_{}_{}_at_time_{}_all_Vc_{}_quadrant.pkl'.format(potkey,Nxr,Nxi,Nyr,Nyi,t,widthxr,widthxi,\
                        widthyr,widthyi,tc,Vc), 'rb')
        
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

        xfin = pickle.load(op)    # If things don't work, check by bra transforming here.
        yfin = pickle.load(op)
        pxfin = pickle.load(op)
        pyfin = pickle.load(op)
        

        #FINCO_2D.bramanifold(xfinbar, pxfinbar, yfinbar, pyfinbar, xfin, pxfin, yfin, pyfin, gammaxf, gammayf, Nxr,Nxi,Nyr,Nyi)

        magnitude = pickle.load(op)
        phaseangle = pickle.load(op)

        sigma = pickle.load(op)
        #print(xfin[:100])
        J0 = pickle.load(op)
        #successcheck =  pickle.load(op)
        maxmomentum  = pickle.load(op)

        #magnitude[magnitude>0.0]=-10.0
        
        x = np.array([0.0,0.0213,0.0260,0.0307,0.0587,0.1053,0.1333,0.1660,0.1846,0.1893,0.2080,0.1800,0.2033,0.1660])
        y = np.array([0.0,0.1446,0.2893,0.4572,0.4619,0.5319,0.6018,0.6112,0.6765,0.7278,0.7698,0.9191,1.1198,1.1804])*-1
        z = np.polyfit(x, y, 7)
        
        plt.figure(1)
        p = np.poly1d(z)
        plt.plot(x,y)
        plt.scatter(x,y)
        
        xc = np.array([0.19-1.03j,0.15-1.05j,0.1-1.08j,0.11-1.12j,0.14-1.15j,0.17-1.18j])
        plt.plot(real(xc),imag(xc),color='g')
        plt.scatter(real(xc),imag(xc),color='g')
        print('Hi')
        
        plotvalue = yfinbar#exp(1j*phaseangle)*magnitude
        submani1 = FINCO_2D.ysubmanifold(plotvalue,Nxr,Nxi,Nyr,Nyi,successcheck)
        submani2 = FINCO_2D.ysubmanifold(sigma,Nxr,Nxi,Nyr,Nyi,successcheck)
        submani3 = FINCO_2D.ysubmanifold(xfin,Nxr,Nxi,Nyr,Nyi,successcheck)
        #submani1[real(submani2)>0.0]=-10.0
        Complex_plotter.plotcomplex(submani1,1,1,y0-widthyr/2,y0+widthyr/2,-widthyi/2,widthyi/2)
        #submani1,1,1,y0,y0+widthyr/2,-widthyi/2,0)
        #submani1,1,1,x0-widthxr/2,x0+widthxr/2,-widthxi/2,widthxi/2)
        #submani1,1,1,y0,y0+widthyr/2,-widthyi/2,0)
        #submani1,1,1,y0,y0+widthyr/2,-widthyi/2,0)
        #submani1,1,1,y0-widthyr/2,y0+widthyr/2,-widthyi/2,widthyi/2)
        
        
        plt.show()

        

        op.close()

        print("--- %s seconds ---" % (time.time() - start_time))


#================================================================================



#print(maxmomentum[0])
#print(Turningpoints)


count=0
##for i in range(len(successcheck)):
##        if (successcheck[i]!=-1):
##                count+=1



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

##3039 (-7-0.8j) (-0.8+1.6j)
##3070 (-7-0.8j) (0.8-2j)
##3071 (-7-0.8j) (0.8-1.6j)
##3655 (-4.6-0.8j) 0j
##4060 (-7-0.4j) (0.4-2j)
##4555 (-5-0.4j) 0j
##4655 (-4.6-0.4j) 0j
##4746 (-4.2-0.4j) (-0.4+0.4j)
##4764 (-4.2-0.4j) (0.4-0.4j)
##5527 (-5+0j) (-1.2+0.8j)
##5555 (-5+0j) 0j
##5583 (-5+0j) (1.2-0.8j)

# 4060 (-7-0.4j) (0.4-2j)
# 5527 (-5+0j) (-1.2+0.8j)
# 5583 (-5+0j) (1.2-0.8j)
# 5652 (-4.6+0j) -1.2j # Stokes divergence
# 5658 (-4.6+0j) 1.2j  # Stokes divergence
# 6148 (-6.6+0.4j) (-0.4+1.2j)
# 6162 (-6.6+0.4j) (0.4-1.2j)
# 6538 (-5+0.4j) (-0.8+1.2j)
# 6572 (-5+0.4j) (0.8-1.2j)
# 7145 (-6.6+0.8j) (-0.4+0j)
# 7165 (-6.6+0.8j) (0.4+0j)
# 7448 (-5.4+0.8j) (-0.4+1.2j)
# 7462 (-5.4+0.8j) (0.4-1.2j)

# 647 (-4.6-2j) (-0.4+0.8j)
# 663 (-4.6-2j) (0.4-0.8j)
# 755 (-4.2-2j) 0j
# 828 (-3.8-2j) (-1.2+1.2j)
# 882 (-3.8-2j) (1.2-1.2j)
# 1847 (-3.8-1.6j) (-0.4+0.8j)
# 1863 (-3.8-1.6j) (0.4-0.8j)
# 2629 (-4.6-1.2j) (-1.2+1.6j)
# 2681 (-4.6-1.2j) (1.2-1.6j)
# 4726 (-4.2-0.4j) (-1.2+0.4j)
# 4784 (-4.2-0.4j) (1.2-0.4j)
# 5717 (-4.2+0j) (-1.6+0.8j)
# 5793 (-4.2+0j) (1.6-0.8j)
# 6045 (-7+0.4j) (-0.4+0j)
# 6065 (-7+0.4j) (0.4+0j)
# 6637 (-4.6+0.4j) (-0.8+0.8j)
# 6673 (-4.6+0.4j) (0.8-0.8j)
# 7255 (-6.2+0.8j) 0j
# 7537 (-5+0.8j) (-0.8+0.8j)
# 7573 (-5+0.8j) (0.8-0.8j)
# 8351 (-5.8+1.2j) -1.6j
# 8359 (-5.8+1.2j) 1.6j


# 1847,1863, - Required, but being eliminated. Also, this region needs additional sampling.


removearr = [3039,3071,3070,3655,4060,4555,4655]

removearr = [5717,5793,1847,1863]

for i in range(Nxr*Nxi*Nyr*Nyi):
    if (maxmomentum[i][2]<100):# and maxmomentum[i][2]<400 ):
      #if i in removearr:
      #print(sigma[i])
      
        dpxdpy= (widthxi*widthyi)/(Nxi*Nyi)        # 22/05/18 : Watch out for gamma here, the differential changes....
        dxdy = (widthxr*widthyr)/(Nxr*Nyr)
        #print(sigma[i],magnitude[i],mpp[i])

        # j0 = 0.25/(np.linalg.det(gammaf))*np.linalg.det(2*np.matmul(gammaf,mqq[i]) + 2*np.dot(np.matmul(gammaf,mqp[i]),S20)\
                      # -1j*mpq[i] - 1j*np.matmul(mpp[i],S20))



        # J0[i] = j0

        # J = -(0.25/np.linalg.det(gammaf))*abs(j0)**2 + 0j
        # Integrand = -(-0.25)**2*(4/(pi**2*np.linalg.det(gammaf)))**0.75*J/(j0**0.5)
                                               # #X[ycut*div]                #X[ycut*div]
        xgf =  (2*gammaxf/pi)**0.25*exp(-gammaxf*(X[ycut*div]-xfin[i])**2 + 1j*pxfin[i]*(X[ycut*div]-xfin[i]))
        ygf =  (2*gammayf/pi)**0.25*exp(-gammayf*(Y-yfin[i])**2 + 1j*pyfin[i]*(Y-yfin[i]))
                                                   # #Y                          #Y

         # ExponentS = 1j*Sfinbar[i]
         # #print(j0)
         # Exponentx = -gammaxf*xfinbar[i]**2 + gammaxf*xfin[i]**2 +\
                      # 1j*pxfinbar[i]*xfin[i] - 1j*pxfin[i]*xfinbar[i]
         # Exponenty = -gammayf*yfinbar[i]**2 + gammayf*yfin[i]**2 +\
                      # 1j*pyfinbar[i]*yfin[i] - 1j*pyfin[i]*yfinbar[i]
         # sigma[i] = Exponentx + Exponenty + ExponentS
         # #print(xfin[i])
         #magnitude[i] = abs(Integrand*exp(sigma[i]))
         #phaseangle[i] = np.angle(Integrand*exp(sigma[i]))
         #print(magnitude[i])

   ##      if(real(xbar[i])>-4.0 and ybar[i]==0.0 and magnitude[i]>0.0):
   ##          print(magnitude[i])

         # and magnitude[i]>20.0 and magnitude[i]<50.0 and successcheck[i]!=-1 #real(sigma[i])<0.0 and
         # and magnitude[i]>0.0001 and magnitude[i]<50.0 and real(xfinbar[i])>0.0

         #dxsi = 2*np.matmul(gammaf,mqq[i]) + 2*np.matmul(np.matmul(gammaf,mqp[i]),S20) -1j*mpq[i] - 1j*np.matmul(mpp[i],S20)

         # pfi = [pxfinbar[i],pyfinbar[i]]
         # pin = [pxbar[i],pybar[i]]
         # #print(dxsi,pfi,pin,1j*np.matmul(pfi,dxsi))

         # Sin = 1j*np.matmul(pfi,dxsi)
         # Sfi = 1j*array(pin)

         # #print(Sfi,Sin)

         # dsigma = Sfi-Sin

         # dsigmax = dsigma[0]
         # dsigmay = dsigma[1]
         # dsigmamag = dsigmax**2 + dsigmay**2
         #print(Turningpoints[i],xbar[i],ybar[i],i)and dsigmamag < 800.0 and magnitude[i]<50.0 and 
          

        if(real(sigma[i])<0.0 and magnitude[i]<10.0): # Here the condition for whether magnitude[i] is negative(i.e filtered) is to be taken care.
         
            ywf+= magnitude[i]*phaseangle[i]*xgf*ygf*dxdy*dpxdpy*exp(1j*pi*Turningpoints[i])#Integrand*exp(sigma[i])
            #print(magnitude[i])
            count+=1
            #if(magnitude[i]>1e-4):
             #   print(i,xbar[i],ybar[i])
               

print('count',count)

#plt.imshow(abs(wf))
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
#submani1 = FINCO_2D.xsubmanifold(plotvalue,Nxr,Nxi,Nyr,Nyi,successcheck)
#submani1[real(submani1)>50] = -0.1
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

##plt.plot(ywf)
##plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

# Comparing split operator wavefunction and FINCO wavefunction
if (plotcompare ==1):
   for col,fname,m in zip(['blue','red'],['ysection_{}.txt'.format('double_slit_tunneling'),\
   '2DFINCO_{}_{}_{}_{}_{}_Y.txt'.format(potkey,Nxr,Nxi,Nyr,Nyi)],['x','.']):
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
          




