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

#print(tre,tim)

def f(t,r,dti, potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, tPrint, tPrintSpace
    #previous = turningpoints[1]
    #print("Hi")
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4 =  r
    ddV1 = ddpotential1(x,y)
    ddV2 = ddpotential2(x,y)
    ddV3 = ddpotential3(x,y)
    ddV4 = ddpotential4(x,y)
    #print('changed')
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),px**2/2 + py**2/2 - potential(x,y),\
                 -(ddV1*mqp1 + ddV2*mqp3),-(ddV1*mqp2 + ddV2*mqp4),\
                 -(ddV3*mqp1 + ddV4*mqp3),-(ddV3*mqp2 + ddV4*mqp4),\
                 -(ddV1*mqq1 + ddV2*mqq3),-(ddV1*mqq2 + ddV2*mqq4),\
                 -(ddV3*mqq1 + ddV4*mqq3),-(ddV3*mqq2 + ddV4*mqq4),\
                   mpp1,mpp2,mpp3,mpp4, mpq1,mpq2,mpq3,mpq4])*dti

##    mpp = array([[mpp1,mpp2],[mpp3,mpp4]])
##    mpq = array([[mpq1,mpq2],[mpq3,mpq4]])
##    mqp = array([[mqp1,mqp2],[mqp3,mqp4]])
##    mqq = array([[mqq1,mqq2],[mqq3,mqq4]])
##
##    #print(mpq)
##   
##    #print(mpq)# -1j*np.matmul(mpp,S20))
##    current = (0.25/(np.linalg.det(gammaf)))*np.linalg.det(2*np.matmul(gammaf,mqq) + 2*np.matmul(np.matmul(gammaf,mqp),S20) -1j*mpq -1j*np.matmul(mpp,S20))
##    #exp(-2*1j*Squantum)*4*((gammaf[0][0] + gammaxx)*(gammaf[1][1] + gammayy) - (gammaf[0][1] + gammaxy)*(gammaf[1][0] + gammayx))
##   
##    if (t>tPrint+tPrintSpace):
##         if(real(current) < 0 and imag(previous)*imag(current) < 0):
##            #print(turningpoints[0])
##            turningpoints[0] += 1
##            #print(turningpoints[0])
##         tPrint = t
##    turningpoints[1] = current
                 
    return drdt

    
def faux(t,r,dti, potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, tPrint, tPrintSpace
    #previous = turningpoints[1]
    #print("Hi")
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4 =  r
    ddV1 = ddpotential1(x,y)
    ddV2 = ddpotential2(x,y)
    ddV3 = ddpotential3(x,y)
    ddV4 = ddpotential4(x,y)
    
    #print('t',t,'px',px)
    #print('changed')
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),px**2/2 + py**2/2 - potential(x,y),\
                 -(ddV1*mqp1 + ddV2*mqp3),-(ddV1*mqp2 + ddV2*mqp4),\
                 -(ddV3*mqp1 + ddV4*mqp3),-(ddV3*mqp2 + ddV4*mqp4),\
                 -(ddV1*mqq1 + ddV2*mqq3),-(ddV1*mqq2 + ddV2*mqq4),\
                 -(ddV3*mqq1 + ddV4*mqq3),-(ddV3*mqq2 + ddV4*mqq4),\
                   mpp1,mpp2,mpp3,mpp4, mpq1,mpq2,mpq3,mpq4])*dti
                   
    return drdt

def loop_around(tcoord):
    tlb = (real(tcoord)-0.05) + 1j*(imag(tcoord)-0.05)
    trb = (real(tcoord)+0.05) + 1j*(imag(tcoord)-0.05)
    tlt = (real(tcoord)-0.05) + 1j*(imag(tcoord)+0.05)
    trt = (real(tcoord)+0.05) + 1j*(imag(tcoord)+0.05)
    
    Timecontour = [0.0,tlb,trb,trt,tlt,tlb,0.0]#trb,trt,tlt,tlb,trb,trt,tlt,tlb,
    loopdata = contour_integration(xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,Timecontour,x0,px0,y0,py0)
    return loopdata
    
                   
def finalcondition(loop,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,xin,pxin,yin,pyin):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=100000)
    Sinit = -1j*(0.25*log(2*gammaxi/pi)) -1j*gammaxi*(xinit**2 - xin**2)+ (xinit*pxinit - xin*pxin) +\
                               -1j*(0.25*log(2*gammayi/pi)) -1j*gammayi*(yinit**2 - yin**2)+ (yinit*pyinit - yin*pyin) 
    if(loop==0):        
        y0 = array([xinit,pxinit,yinit,pyinit,Sinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    else:
        y0 = loopdata
    sol.set_initial_value(y0,t=0.0)

    dtim = tim[1]-tim[0]
    dtre = tre[1]-tre[0]

    print(dtim,dtre)
    

    ti = int(len(tim)/2) + 1
    dti = 1j
    sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    trun = dtim
    while(ti<len(tim)):
        #print(Time[0][ti])
        sol.integrate(trun)
        data[0][ti] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        ti+=1
        trun+=dtim
        

    sol.set_initial_value(y0,t=0.0)
    dti = -1j
    sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    ti = int(len(tim)/2) - 1
    trun = dtim
    while(ti>=0):
        sol.integrate(trun)
        data[0][ti] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        ti-=1
        trun+=dtim
        

    data[0][int(len(tim)/2)] = y0
    #print(Time[0][int(len(tim)/2)],data[0][int(len(tim)/2)][0] )
    #print(Time[0][int(len(tim)/2)])
    #print(data[:][:][0])
    dti = 1.0
    sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    for ti in range(len(tim)):
        sol.set_initial_value(data[0][ti],t=0.0)
        trun=dtre
        for tr in range(1,len(tre)):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti],data[tr][ti][0])
            trun+=dtre

def contour_integration(xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,T,xin,pxin,yin,pyin):
    sol = ode(faux)
    sol.set_integrator('zvode',nsteps=100000)
    trajdata=0.0j

    Sinit = -1j*(0.25*log(2*gammaxi/pi)) -1j*gammaxi*(xinit**2 - xin**2)+ (xinit*pxinit - xin*pxin) +\
                               -1j*(0.25*log(2*gammayi/pi)) -1j*gammayi*(yinit**2 - yin**2)+ (yinit*pyinit - yin*pyin) 
            
    y0 = array([xinit,pxinit,yinit,pyinit,Sinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    
    if(len(T) ==1):
        sol.set_initial_value(y0,t=0.0)
        trajdata = sol.y
    else:
        
        tempsol = y0
        for tc in range(len(T)-1):
            
            #print(T[tc+1]-T[tc])
            timcom = T[tc+1]-T[tc]
            abstim = abs(timcom)
            dti = timcom/abstim
            sol.set_initial_value(tempsol,t=0.0)
            sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
            sol.integrate(abstim)
            tempsol = sol.y
        trajdata=sol.y
        print('Trajdata',trajdata[:-16])
        return trajdata


#T = [ [[] for i in range(len(tim))] for j in range(len(tre))]

Time = zeros((len(tre),len(tim))) +0j


for tr in range(len(tre)):
    for ti in range(len(tim)):
        Time[tr][ti] = tre[tr] + 1j*tim[ti]
        
#print(Time)
#print(T)

mpp = zeros((len(tre),len(tim),2,2),dtype=complex)
mqp = zeros((len(tre),len(tim),2,2),dtype=complex)
mpq = zeros((len(tre),len(tim),2,2),dtype=complex)
mqq = zeros((len(tre),len(tim),2,2),dtype=complex)


gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
data = zeros((len(tre),len(tim),21),dtype = complex)


for i in range(len(pointarr)):
    xinit = pointarr[i][0]
    yinit = pointarr[i][1]
    pxinit = -1j*(2*gammaxi*x0 + 1j*px0 - 2*gammaxi*xinit)
    pyinit = -1j*(2*gammayi*y0 + 1j*py0 - 2*gammayi*yinit)
    tcoord = 1.553 +0.075j
    loopdata= loop_around(tcoord)
    finalcondition(1,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,x0,px0,y0,py0)
    op = open('Singularities_double_slit_pointarr_{}_branch_2.pkl'.format(i),'wb')
    pickle.dump(data,op)
    op.close()


op = open('Singularities_double_slit_corrected_point_{}_branch_2.pkl'.format(point),'wb')
pickle.dump(data,op)
op.close()

#j0 = 0.25/(np.linalg.det(gammaf))*np.linalg.det(2*np.dot(gammaf,mqq) + 2*np.dot(np.dot(gammaf,mqp),S20) -1j*mpq - 1j*np.dot(mpp,S20))

    

    
