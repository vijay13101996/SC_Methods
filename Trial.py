from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt
#from matplotlib.figure import figure
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import Interference
import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
from datetime import datetime
import time
import pickle, pprint
#from Parameter_file_split_custom_double_slit_eckart import *
#from Parameter_file_split import *
potkey = 'reaction_coordinate'
import Parameters_File
Parameters_File.method_sel('Split')
Parameters_File.wf_init_set('Split')
Parameters_File.initial_conditions(3)
Parameters_File.grid_select(3)
Parameters_File.pot_select(potkey)
from Parameters_File import *
print(dt,tfinal)
##import finalconditions as fc

#------ Miller's parameters-----

##s = 50
##w = 600
##V0 = 8000

#(V0 - 2.5*10**(-40)*y**2 + 2.0*10**(-84)*y**4)*(exp(-x**2/s**2))
#print("{0} {1}".format(V0,s))
#-------------------------------

plotsection=1
plotpotential =1

def halfnormal(x,y,wf):
    sum = 0.0
    count = zeros_like(wf)
    for i in range(len(x)):
        for k in range(len(y)):
            if(x[k]>=0.0):
                count[i][k]+=1
                sum+= (abs(wf[i][k])**2)*dx*dy
    
##    plt.imshow(real(count),origin ='lower')
##    plt.show()
    print(sum,count)


def xsection(wf,ywf,xcut):
   for i in range(len(ywf)):
        ywf[i] = wf[i][xcut]#*(pi/gammax)**0.25

def ysection(wf,xwf,ycut):
    for i in range(len(xwf)):
        xwf[i] = wf[ycut][i]#*(pi/gammay)**0.25

def printfun(x,wf,outFile):
    for ax,w in zip(x,wf):
        print(real(ax),real(w),imag(w),file=outFile)
    print("",file=outFile)



xgrid = arange(xi,xi+lx,dx)
ygrid = arange(yi,yi+ly,dy)

x,y = meshgrid(xgrid,ygrid)

pot = potential(x,y) #zeros_like(x)

##for i in range(len(xgrid)):
##    for k in range(len(ygrid)):
##        if(abs(x[i][k])<1.0 and  abs(y[i][k])<3.0):
##            pot[i][k] = V0*E/np.cosh(x[i][k])**2
##        else:
##            pot[i][k] = potential(x[i][k],y[i][k])


if (plotpotential==1):
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.imshow(pot,origin='lower',extent=(xi,xi+lx,yi,yi+ly),vmax=20)
    plt.show()

##fig = plt.figure(2)
##ax = Axes3D(fig)
##ax.set_xlabel('x')
##ax.set_ylabel('y')
###ax.set_zlim(0, 0.01)
###pot[pot>.01]=.01
##ax.plot_surface(x,y,abs(wf))
##plt.show()

wf = (2*gammax/pi)**(0.25)*(2*gammay/pi)**(0.25)*exp(-gammax*(x-x0)**2 + 1j*px0*(x-x0))*exp(-gammay*(y-y0)**2 + 1j*py0*(y-y0))

potop = exp(-0.5j*pot*dt/hbar)

kinx = arange(-len(xgrid)/2,len(xgrid)/2)*2*pi/lxm
             
kiny = arange(-len(ygrid)/2, len(ygrid)/2)*2*pi/ly


kingridx,kingridy = meshgrid(kinx,kiny)

kin = exp((-1j*(kingridx**2)/(2*m)*dt - 1j*(kingridy**2/(2*m)*dt))*hbar)

t = 0 

print(tfinal)
while (t<tfinal):
    wf = wf*potop
    wf = np.fft.fft2(wf)
    wf = np.fft.fftshift(wf)
    wf = wf*kin
    wf = np.fft.ifftshift(wf)
    wf = np.fft.ifft2(wf)
    wf = wf*potop
    t+=dt
    
realwf = real(wf)
rhowf = (real(wf))**2 + (imag(wf))**2
plt.figure(3)

ycut = 0.55
xcut=0.5
plt.suptitle('QM Wavefunction at time t=%s '%tfinal, fontsize=14, fontweight='bold')
plt.imshow(abs(wf),origin = 'lower',zorder=0,extent=[xi,xi+lx,yi,yi+ly],vmax=0.1)
plt.plot(xgrid[int(div*ycut)]*ones_like(xgrid),ygrid, '--', linewidth=1, color='firebrick')
plt.plot(xgrid,ygrid[int(div*xcut)]*ones_like(ygrid), '--', linewidth=1, color='firebrick')
plt.show()

print("Half-normalising factor:")
halfnormal(xgrid,ygrid,wf)

ywf = zeros_like(ygrid) + 0.0j
xwf = zeros_like(xgrid) + 0.0j

##fig = plt.figure(2)
##ax = Axes3D(fig)
##ax.set_xlabel('x')
##ax.set_ylabel('y')
###ax.set_zlim(0, 0.01)
###pot[pot>.01]=.01
##ax.plot_surface(x,y,abs(wf))
##plt.show()

print(xgrid[int(ycut*div)],ygrid[int(xcut*div)]) # : To check for the origin.
xsection(wf,ywf,int(div*ycut))
ysection(wf,xwf,int(div*xcut))

outfile = open('ysection_%s.txt'%potkey,'w')
printfun(ygrid,abs(ywf),outfile)
outfile.close()

outfile1 = open('xsection_%s.txt'%potkey,'w')
printfun(xgrid,xwf,outfile1)
outfile1.close()

if(plotsection == 1):

    plt.figure(10)
    plt.suptitle('ysection_%s'%potkey, fontsize=14, fontweight='bold')
    plt.plotfile('ysection_%s.txt'%potkey, delimiter=' ', cols=(0, 1), 
                 names=('y', 'psi(y)'), marker='_')
    plt.show()

    fig = plt.figure(11)
    fig.suptitle('xsection_%s'%potkey, fontsize=14, fontweight='bold')
    plt.plotfile('xsection_%s.txt'%potkey, delimiter=' ', cols=(0, 1), 
                 names=('x', 'psi(x)'), marker='_')

    plt.show()

##plt.figure(7)
##for col,fname,m in zip(['blue','red'],['ysection_double_slit_gaussian.txt', 'ysection_double_slit_eckart.txt'],['x','.']):
##          data=np.loadtxt(fname)
##          X=data[:,0]
##          Y=data[:,1]
##          Z=data[:,2]#**2+data[:,1]**2
##          plt.plot(X,Y,m,color=col,linestyle='-')
##          plt.plot(X,Z,m,color=col,linestyle='-')
##          
##plt.show()
