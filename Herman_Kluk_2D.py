
from scipy import *
from scipy.fftpack import fft,ifft, fftshift, ifftshift
from scipy.integrate import odeint, ode
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool,cpu_count
from functools import partial
import os
#import newfunc as cfunc
#from Parameter_file_HHKK_harmonic import *
##X = arange(xi,xi+lx,dx)
##Y = arange(yi,xi+ly,dy)
tPrint=0
tPrintSpace=0.01
m = 1.0

# Initial phase space grid is initialised correctly 
# Checkpoints: Trajectories : x,y data are going on correctly.

gammaxf = 1.0
gammayf = 1.0
gammaxi = 1.0
gammayi = 1.0
gammai = [[gammaxi,0],[0,gammayi]]
gammaf = [[gammaxf,0],[0,gammayf]]
S20 = ([2*gammaxi*1j,0],[0,2*gammayi*1j])

def initmanifold(x0,y0, px0,py0, xinit, pxinit, yinit, pyinit,Nx,Npx,Ny,Npy,deltax,deltay,deltapx,deltapy):
##    if (N==1):
##        xinit[0] = x0
##        pxinit[0] = px0
##        yinit[0] = y0
##        pyinit[0] = py0
##    else:
        for i in range(Nx+1):
            for j in range(Npx+1):
                for k in range(Ny+1):
                    for l in range(Npy+1):
                        #print('i',((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l)
                        xinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l] = x0  - deltax*Nx/2 + i*deltax                     
                        pxinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l] = px0 - deltapx*Npx/2 + j*deltapx
                        yinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l] = y0 - deltay*Ny/2 + k*deltay
                        pyinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l] = py0 - deltapy*Npy/2 + l*deltapy
                        #print(xinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l],yinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l],\
                         #       pxinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l],pyinit[((i*(Npx+1) + j)*(Ny+1) + k)*(Npy+1) + l])


def f(t,r,turningpoints,gammaf, potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4,gammai):
    global current, tPrint, tPrintSpace
    previous = turningpoints[1]
    #print("Hi")
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4 =  r
    ddV1 = ddpotential1(x,y)
    ddV2 = ddpotential2(x,y)
    ddV3 = ddpotential3(x,y)
    ddV4 = ddpotential4(x,y)
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),px**2/2 + py**2/2 - potential(x,y),\
                 -(ddV1*mqp1 + ddV2*mqp3),-(ddV1*mqp2 + ddV2*mqp4),\
                 -(ddV3*mqp1 + ddV4*mqp3),-(ddV3*mqp2 + ddV4*mqp4),\
                 -(ddV1*mqq1 + ddV2*mqq3),-(ddV1*mqq2 + ddV2*mqq4),\
                 -(ddV3*mqq1 + ddV4*mqq3),-(ddV3*mqq2 + ddV4*mqq4),\
                   mpp1,mpp2,mpp3,mpp4, mpq1,mpq2,mpq3,mpq4])
    mpp = [[mpp1,mpp2],[mpp3,mpp4]]
    mpq = [[mpq1,mpq2],[mpq3,mpq4]]
    mqp = [[mqp1,mqp2],[mqp3,mqp4]]
    mqq = [[mqq1,mqq2],[mqq3,mqq4]]
    
    
    current = (np.linalg.det(0.5*(np.add(mpp,mqq) - 1j*np.matmul(gammai,mqp) + 1j*np.matmul(np.linalg.inv(gammai),mpq))))
    # Check again here..
    #if (t>tPrint+tPrintSpace):
         #print (imag(previous)*imag(current))
    if(real(current) < 0 and imag(previous)*imag(current) < 0):
            #print(turningpoints[0],'hi')
            turningpoints[0] += 1
            #print(turningpoints[0],'hi')
        #tPrint = t
    turningpoints[1] = current
                 
    return drdt

def traj_propagator(init_array,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4,timear):

    global tPrint,current
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1000)
    
    tPrint = -tPrintSpace      #11_05_18 Check this line.
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4=init_array
    mpp = array([[mpp1,mpp2],[mpp3,mpp4]])
    mpq = array([[mpq1,mpq2],[mpq3,mpq4]])
    mqp = array([[mqp1,mqp2],[mqp3,mqp4]])
    mqq = array([[mqq1,mqq2],[mqq3,mqq4]])
    sol.set_initial_value(init_array,t=0.0)
    tempar = [0,(np.linalg.det(0.5*(mpp+mqq - 1j*np.matmul(gammai,mqp) + 1j*np.matmul(np.linalg.inv(gammai),mpq))))]
    current = (np.linalg.det(0.5*(mpp+mqq - 1j*np.matmul(gammai,mqp) + 1j*np.matmul(np.linalg.inv(gammai),mpq))))
    sol.set_f_params(tempar,gammaf,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4,gammai)
    trajdata = zeros((len(timear),len(init_array)+1))
    
    for tc in range(len(timear)):
        
        if tc>0:
          sol.integrate(timear[tc])

        trajdata[tc,:-1] = sol.y.real
        trajdata[tc,21] = tempar[0]
        #print(tempar)
        #print(trajdata)
        
    return trajdata


def finalconditions(xinit,pxinit,yinit, pyinit, x0, px0, y0, py0, gammaxi,gammayi,\
                    gammaxf,gammayf,data,Nx,Npx,Ny,Npy,potential,dxpotential,dypotential,\
                    ddpotential1, ddpotential2, ddpotential3, ddpotential4, timear,gammaf,gammai,mpp,mpq,mqp,mqq):
    global tPrint
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1000)         #11/05/18 Check timesteps, it might be excessive.
    #print("1")
    
    onearr = ones(len(xinit))
    zeroarr = zeros(len(xinit))
    Sarr = zeroarr

    init_conditions_array=array([xinit,pxinit,yinit,pyinit,Sarr,onearr,zeroarr,zeroarr,onearr,\
        zeroarr,zeroarr,zeroarr,zeroarr, zeroarr,zeroarr,zeroarr,zeroarr, onearr,zeroarr,zeroarr,onearr])
        
    print('CPU',len(os.sched_getaffinity(0)))
    p = Pool(cpu_count()-1 or 1 )
    data=p.map(partial(traj_propagator,potential=potential,dxpotential=dxpotential,\
        dypotential=dypotential,ddpotential1=ddpotential1,ddpotential2=ddpotential2,ddpotential3=ddpotential3,ddpotential4=ddpotential4,timear=timear)\
        ,zip(*init_conditions_array))
    data = array(data)
    #print(data[:,3,21])
    return data
    
    
    
    
    # for i in range((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1)):
                    # tPrint = -tPrintSpace      #11_05_18 Check this line.
            
                    # y0 = array([xinit[i],pxinit[i],yinit[i],pyinit[i],\
                               # 0.0,\
                            # 1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
                     
                    # sol.set_initial_value(y0,t=0.0)
                    # tempar = [0,(np.linalg.det(0.5*(mpp[i] + mqq[i] - 1j*np.matmul(gammai,mqp[i]) + 1j*np.matmul(np.linalg.inv(gammai),mpq[i]))))]

                    
                    
                    # sol.set_f_params(tempar,gammaf,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4,gammai)
                  
                    
                    # for tc in range(len(timear)):
                        
                        # if tc>0:
                          # sol.integrate(timear[tc])
# ##                        if all(y0[:4].imag==0):
# ##                            #print(sol.y[:5])
                        # #print(len(sol.y),len(data[i,tc,:-1]))
                        # data[i,tc,:-1] = sol.y.real
                        
                        
                        # data[i,tc,21] = tempar[0]
                        

def finmanifold(xfin, pxfin, yfin, pyfin,N,data,T):
    for i in range(N+1):
        for j in range(N+1):
            for k in range(N+1):
                for l in range(N+1):
                    for tc in range(len(T)):
                        xfin = data[:,tc,0]
                        pxfin = data[:,tc,1]
                        yfin = data[:,tc,2]
                        pyfin = data[:,tc,3]
                        Sfin = data[:,tc,4]
                        #print(xfin,pxfin,yfin,pyfin)


def coherentgaussian(gamma,q,qcoord,pcoord):
    return ((gamma/pi)**0.25)*exp(-gamma*((q-qcoord)**2)/2 + 1j*pcoord*(q-qcoord))

def coherentoverlap(q0,qcoord,p0,pcoord,gamma):
    return exp(-gamma*((q0-qcoord)**2)/4 - ((p0-pcoord)**2)/(4*gamma) + 0.5j*(pcoord+p0)*(qcoord- q0))

##def overlaptemp(w1,w2,dx):
##    overl = 0.0
##    for i in range(len(X)):
##        overl+= conj(w1[i])*w2[i]*dx
##
##    print(overl)
##    
##
####w2  =  coherentgaussian(1.0,X,x0,px0)
####w1 = coherentgaussian(1.0,X,-7.0,1.0)
####
####gammax =1.0
####overlaptemp(w1,w2,dx)
####print(coherentoverlap(x0,-7.0,px0,1.0,gammax))
##
##def overlap(w1,w2):
##    overl = 0.0 + 0j
##    for i in range(len(X)):
##        for j in range(len(Y)):
##            overl+= conj(w1[i][j])*w2[i][j]*dx*dy
##            #print((coherentgaussian(gammax,x[i][j],xcoord,pxcoord).conj()))
##            
##
##    print('overlap',overl)
##    return overl
##
##def normal(w,dx,dy):
##    sum = 0.0
##    for i in range(len(X)):
##        for j in range(len(Y)):
##            sum+= w[i][j]*conj(w[i][j])*dx*dy
            

    print('sum',sum)

def HHKK_section(Nx,Npx,Ny,Npy, x,y,ywf,deltax,deltay,deltapx,deltapy,mpp,mpq,mqp,mqq,S,\
    x0,y0,px0,py0,xinit,yinit,pxinit,pyinit,xfin,yfin,pxfin,pyfin,gammax,gammay,turningpoints,xsection):
    print('gammax',gammax)
    gamma = array([[gammax,0],[0,gammay]])
    for i in range((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1)):
        Rpqt = (np.linalg.det(0.5*(mpp[i] + mqq[i] - 1j*np.matmul(gamma,mqp[i]) + 1j*np.matmul(np.linalg.inv(gamma),mpq[i]))))**0.5
        #print(Rpqt)
        xgf = coherentgaussian(gammax,x[xsection],xfin[i],pxfin[i])#((gammax/pi)**0.25)*exp(-gammax*((x-xfin[i])**2)/2 + 1j*pxfin[i]*(x-xfin[i]))
        ygf = coherentgaussian(gammay,y,yfin[i],pyfin[i])#((gammay/pi)**0.25)*exp(-gammay*((y-yfin[i])**2)/2 + 1j*pyfin[i]*(y-yfin[i]))
        #normal(x,y,xgf*ygf,dx,dy)
##        xgfinit = coherentgaussian(gammax,x,xinit[i],pxinit[i])
##        ygfinit = coherentgaussian(gammay,y,yinit[i],pyinit[i])
##        xgfi = coherentgaussian(1.0,x,x0,px0)
##        ygfi = coherentgaussian(1.0,y,y0,py0)
##        normal(xgfi*ygfi,dx,dy)
##        normal(xgfinit*ygfinit,dx,dy)
        xoverlap = coherentoverlap(x0,xinit[i],px0,pxinit[i],gammax)
        #exp(-gammax*((x0-xinit[i])**2)/4 - ((px0-pxinit[i])**2)/(4*gammax) + 0.5j*(pxinit[i] + px0)*(xinit[i] - x0))
        yoverlap = coherentoverlap(y0,yinit[i],py0,pyinit[i],gammay)
        #exp(-gammay*((y0-yinit[i])**2)/4 - ((py0-pyinit[i])**2)/(4*gammay) + 0.5j*(pyinit[i] + py0)*(yinit[i] - y0))
        ywf += ((0.5/pi)**2)*deltax*deltapx*deltay*deltapy*Rpqt*exp(1j*S[i])*xgf*ygf*yoverlap*xoverlap*exp(1j*pi*turningpoints[i])



def HHKK(Nx,Npx,Ny,Npy, x,y,wf,deltax,deltay,deltapx,deltapy,mpp,mpq,mqp,mqq,S,x0,y0,px0,py0,xinit,yinit,pxinit,pyinit,xfin,yfin,pxfin,pyfin,gammax,gammay,turningpoints):
    gamma = array([[gammax,0],[0,gammay]])
    for i in range((Nx+1)*(Npx+1)*(Ny+1)*(Npy+1)):
        Rpqt = (np.linalg.det(0.5*(mpp[i] + mqq[i] - 1j*np.matmul(gamma,mqp[i]) + 1j*np.matmul(np.linalg.inv(gamma),mpq[i]))))**0.5
        #print(Rpqt)
        xgf = coherentgaussian(gammax,x,xfin[i],pxfin[i])#((gammax/pi)**0.25)*exp(-gammax*((x-xfin[i])**2)/2 + 1j*pxfin[i]*(x-xfin[i]))
        ygf = coherentgaussian(gammay,y,yfin[i],pyfin[i])#((gammay/pi)**0.25)*exp(-gammay*((y-yfin[i])**2)/2 + 1j*pyfin[i]*(y-yfin[i]))
        #normal(x,y,xgf*ygf,dx,dy)
##        xgfinit = coherentgaussian(gammax,x,xinit[i],pxinit[i])
##        ygfinit = coherentgaussian(gammay,y,yinit[i],pyinit[i])
##        xgfi = coherentgaussian(1.0,x,x0,px0)
##        ygfi = coherentgaussian(1.0,y,y0,py0)
##        normal(xgfi*ygfi,dx,dy)
##        normal(xgfinit*ygfinit,dx,dy)
        xoverlap = coherentoverlap(x0,xinit[i],px0,pxinit[i],gammax)
        #exp(-gammax*((x0-xinit[i])**2)/4 - ((px0-pxinit[i])**2)/(4*gammax) + 0.5j*(pxinit[i] + px0)*(xinit[i] - x0))
        yoverlap = coherentoverlap(y0,yinit[i],py0,pyinit[i],gammay)
        #exp(-gammay*((y0-yinit[i])**2)/4 - ((py0-pyinit[i])**2)/(4*gammay) + 0.5j*(pyinit[i] + py0)*(yinit[i] - y0))
        wf += ((0.5/pi)**2)*deltax*deltapx*deltay*deltapy*Rpqt*exp(1j*S[i])*xgf*ygf*yoverlap*xoverlap*exp(1j*pi*turningpoints[i])
        #print(turningpoints[i])

        
##        overlap(xgfinit*ygfinit,xgfi*ygfi)
##        print(xoverlap*yoverlap)
        

                                   # 

    #print(wf)
#*exp(1j*S[i])*((gammax*gammay/pi**2)**0.25)\
#*(np.linalg.det(0.5*(mpp[i] + mqq[i] - 1j*np.matmul(gamma,mqp[i]) + 1j*np.matmul(np.linalg.inv(gamma),mpq[i]))))**0.5\

def xphasespacegrid(arr,Nx,Npx,Ny,Npy):
    subspace = zeros((Nx+1,Npx+1))
    for i in range(Nx+1):
        for j in range(Npx+1):
            subspace[i][j] = arr[((i*(Npx+1) + j)*(Ny+1) + Ny/2)*(Npy+1) + Npy/2]
            #print(arr2[((i + N*j)*N + N/2)*N + N/2])
    return subspace

def yphasespacegrid(arr,Nx,Npx,Ny,Npy):
    subspace = zeros((Ny+1,Npy+1))
    for k in range(Ny+1):
        for l in range(Npy+1):
            subspace[k][l] = arr[((Nx*(Npx+1)/2 + Npx/2)*(Ny+1) + k)*(Npy+1) + l] 
    return subspace




