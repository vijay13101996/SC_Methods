from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
from datetime import datetime
import time
import pickle, pprint
import Herman_Kluk_2D
#import finalconditions as fc
from mpl_toolkits.mplot3d import Axes3D
import scipy.ndimage
potkey='double_slit_tunneling'
import Parameters_File
Parameters_File.method_sel('HHKK')
Parameters_File.wf_init_set('HHKK')
Parameters_File.initial_conditions(3)
Parameters_File.grid_select(3)
Parameters_File.time_grid(3)
Parameters_File.pot_select(potkey)
from Parameters_File import *

start_time = time.time()
#div = 256

print(Nx)

##Nx=20
##Npx=10
##Ny=20
##Npy=10


plotsection=1
plotsurface=0
plotcompare=1
plotsubmanifold=1
#-----------------------------------------------------------------------------

def normal(x,y,wf):
    sum = 0.0
    for i in range(len(x)):
        for k in range(len(y)):
            sum+= (abs(wf[i][k]))**2*dx*dy
    print(sum)

def xsection(wftemp,ywf,xcut):
   for i in range(len(ywf)):
        ywf[i] = wftemp[i][int(xcut)]
     

def ysection(wftemp,xwf,ycut):
    sum = 0.0
    for i in range(len(xwf)):
        xwf[i] = wftemp[int(ycut)][i]

def printfun(x,wf,outFile):
    for ax,w in zip(x,wf):
        print (real(ax),real(w),imag(w),file=outFile)
    print ("",file=outFile)

#-----------------------------------------------------------------------------


xinit =  zeros((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1)) 
yinit = zeros_like(xinit) 
pxinit = zeros_like(xinit)
pyinit = zeros_like(yinit)

pxfin = zeros_like(pxinit)
xfin = zeros_like(xinit)
pyfin = zeros_like(pyinit)
yfin = zeros_like(yinit)
Turningpoints = zeros_like(xinit)
Sfin = zeros_like(xinit)

mpp = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),2,2)) 
mpq = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),2,2)) 
mqp = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),2,2)) 
mqq = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),2,2))

gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
#gammafin = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),2,2)) + 0j
#Squantum = zeros_like(xinit) + 0j

T = arange(0,t,dt)
data = zeros(((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1),len(T),27)) 

X = arange(xi,xi+lx,dx)
Y = arange(yi,xi+ly,dy)

x,y = meshgrid(X,Y)

wf = zeros_like(x) + 0j

gammax=1
gammay=1


for tc in range(3):

    #wf = wf*(0.0)
    op = open('HHKK_{}_{}_{}_{}_{}_width_{}_{}_{}_{}_tfinal_{}_at_time_{}.pkl'.format(potkey,Nx,Npx,Ny,Npy,widthx,widthpx,widthy,widthpy,t,tc), 'rb')
    #x = pickle.load(op)
    #y = pickle.load(op)
    #wf = pickle.load(op)
    

    xinit = pickle.load(op)
    pxinit = pickle.load(op)
    yinit = pickle.load(op)
    pyinit = pickle.load(op)

##    print('xinit',xinit)
##    print('yinit',yinit)
##    print('pxinit',pxinit)
##    print('pyinit',pyinit)

    xfin = pickle.load(op)
    pxfin = pickle.load(op)
    yfin = pickle.load(op)
    pyfin = pickle.load(op)

    Sfin = pickle.load(op)
    #gammafin = pickle.load(op)

    mpp = pickle.load(op)
    mpq = pickle.load(op)
    mqp = pickle.load(op)
    mqq = pickle.load(op)

    Turningpoints= pickle.load(op)
    #print(np.count_nonzero(Turningpoints))
    #Squantum=pickle.load(op)

    xfin = pickle.load(op)
    yfin = pickle.load(op)
    pxfin = pickle.load(op)
    pyfin = pickle.load(op)

##    plt.figure(23)
##    plt.imshow(abs(wf),origin='lower',extent=[xi,xi+lx,yi,yi+ly])
##    plt.show()

    op.close()
    print("The Normalizing factor of the wavefunction is:")
    normal(X,Y,wf)

##op = open('2D_HHKK_harmonic_N={}_width={}_data.pkl'.format(N,width), 'rb')
##data = pickle.load(op)
##op.close()

#-------------------------------------------------------------------------------
##xanalytic = xinit*np.cos(k*T[1]) + pxinit*np.sin(k*T[1])
##yanalytic = yinit*np.cos(k*T[1]) + pyinit*np.sin(k*T[1])
##
##for i in range(10):
##    #print(xfin[i],xanalytic[i])
##    #print(yfin[i],yanalytic[i])
##    print('mpp',mpp[i])
##    print('mqq',mqq[i])
##    print('mqp',mqp[i])
##    print('mpq',mpq[i])
##
##print(xanalytic,yanalytic)

#--------------------------------------------------------------------------------

X = arange(xi,xi+lx,dx)
Y = arange(yi,yi+ly,dy)

x,y = meshgrid(X,Y)

xwf = zeros_like(X) + 0j
ywf = zeros_like(Y) +0j
print(len(wf))

##tempwf = abs(wf)
###tempwf[tempwf>0.01] = 0.01
##plt.figure(23)
##plt.imshow(tempwf,origin='lower',extent=[xi,xi+lx,yi,yi+ly])
##plt.show()

xoverlap = Herman_Kluk_2D.coherentoverlap(x0,xinit,px0,pxinit,gammax)
yoverlap = Herman_Kluk_2D.coherentoverlap(y0,yinit,py0,pyinit,gammay)

overlap = xoverlap*yoverlap

Rpqt = (np.linalg.det(0.5*(mpp + mqq - 1j*np.matmul(gammaf,mqp) + 1j*np.matmul(np.linalg.inv(gammaf),mpq))))**0.5

xphasespace = Herman_Kluk_2D.xphasespacegrid(Rpqt*exp(1j*Sfin)*overlap,Nx,Npx,Ny,Npy)
yphasespace = Herman_Kluk_2D.yphasespacegrid(Rpqt*exp(1j*Sfin)*overlap,Nx,Npx,Ny,Npy)

#print(xphasespace)
#print(yphasespace)

if(plotsubmanifold==1):
    plt.figure(1)
    plt.imshow(xphasespace,origin='lower', extent=[x0-widthx/2,x0+widthx/2,px0-widthpx/2,px0+widthpx/2],interpolation='None')
    plt.show()

    plt.figure(2)
    plt.imshow(yphasespace,origin='lower', extent=[y0-widthy/2,y0+widthy/2,py0-widthpy/2,py0+widthpy/2],interpolation='None')
    plt.show()

if (plotsurface==1):
    fig = plt.figure(2)
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_zlim(0, 0.01)
    #pot[pot>.01]=.01
    ax.plot_surface(x,y,abs(wf))
    plt.show()


if(plotsection==1):
    ycut = 0.55
    xcut=0.5
    print(X[int(ycut*div)],Y[int(xcut*div)]) # : To check for the origin.

    ysection(wf,xwf,int(xcut*div))
    xsection(wf,ywf,int(ycut*div))
    
    Herman_Kluk_2D.HHKK_section(Nx,Npx,Ny,Npy, X,Y,ywf,deltax,deltay,deltapx,deltapy,mpp,mpq,mqp,mqq,Sfin,\
    x0,y0,px0,py0,xinit,yinit,pxinit,pyinit,xfin,yfin,pxfin,pyfin,gammax,gammay,Turningpoints,int(ycut*div))

    print("The Normalizing factor of the wavefunction is:")
    normal(X,Y,wf)


    outFilex = open('2D_HHKK_{}_{}_{}_{}_X_{}_{}_{}_{}.txt'.format(potkey,Nx,Npx,Ny,Npy,widthx,widthpx,widthy,widthpy),'w')
    printfun(X,xwf,outFilex)
    outFilex.close()

    outFiley = open('2D_HHKK_{}_{}_{}_{}_Y_{}_{}_{}_{}.txt'.format(potkey,Nx,Npx,Ny,Npy,widthx,widthpx,widthy,widthpy),'w')
    printfun(Y,ywf,outFiley)
    outFiley.close()

    # fig = plt.figure(11)
    # #fig.suptitle('2D_HHKK_{}_Y_{}_{}.txt'.format(potkey,N,width), fontsize=14, fontweight='bold')
    # plt.plotfile('2D_HHKK_{}_{}_{}_{}_Y_symmetric_{}_{}.txt'.format(potkey,Nx,Npx,Ny,Npy,width), \
                                            # delimiter=' ', cols=(0, 1), names=('y', 'psi(y)'), marker='_')

    # plt.show()

if (plotcompare==1):
    for col,fname,m in zip(['blue','red'],['ysection_{}.txt'.format('double_slit_tunneling'), \
                                           '2D_HHKK_{}_{}_{}_{}_Y_{}_{}_{}_{}.txt'.format(potkey,Nx,Npx,Ny,Npy,widthx,widthpx,widthy,widthpy)],['x','.']):
            data=np.loadtxt(fname)
            X=data[:,0]
            Y=data[:,1]
            Z=(data[:,2]**2 + data[:,1]**2)**0.5
            fig = plt.figure(7)
            fig.suptitle('Comparison between Quantum Mechanical and HHKK wavefunctions', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            #fig.subplots_adjust(top=0.85)
            ax.set_title('Psi(y) vs y')
            ax.set_xlabel('y')
            ax.set_ylabel('psi(y)')
            #plt.plot(X,Y,m,color=col,linestyle='-')
            #plt.draw()
            plt.plot(X,Z,m,color=col,linestyle='-')

    plt.show()
          
