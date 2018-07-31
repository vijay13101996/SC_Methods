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
import Complex_plotter
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
import io
import matplotlib.cm as cm

plot_full_time = 1
#--------------------------------------------

print(xinit,yinit)
print(pointarr)
plot = 1
tarr = []
pxarr = []
pyarr = []
xarr = []
yarr=[]
colorarr = []
trajon = 0

#print(tre,tim)
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

height = 1
c = mcolors.ColorConverter().to_rgb
rgb = make_colormap(\
    [c('red'), c('yellow'), 0.15, c('yellow'), c('green'), c('cyan'), 0.5, c('cyan'),
     c('blue'), c('magenta'), 0.85, c('magenta'), c('red')])
     
def f(t,r,dti, potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4):
    global current, tPrint, tPrintSpace
    #previous = turningpoints[1]
    #print("Hi")
    #print('r',len(r))
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
    global current, tPrint, tPrintSpace,trajon
    #previous = turningpoints[1]
    #print("Hi")
    r = r.view(dtype=complex)
    #print('r',len(r))
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4 =  r
    ddV1 = ddpotential1(x,y)
    ddV2 = ddpotential2(x,y)
    ddV3 = ddpotential3(x,y)
    ddV4 = ddpotential4(x,y)
    
    #print('t',t,'px',px)
    #print('changed')
    
    # if(trajon==1):
        # print('here')
        # tarr.append(t)
        # pxarr.append(px)
        # pyarr.append(py)
        
    #print('pole',py,py*(x-1j*pi/2)**0.5)#,(x-1j*pi/2))
    #print('Energy',px**2/2 + py**2/2 + potential(x,y),'Potential',potential(x,y))
    #print('x',x,'px',px)
    
    drdt = array([px,-dxpotential(x,y), py, -dypotential(x,y),px**2/2 + py**2/2 - potential(x,y),\
                 -(ddV1*mqp1 + ddV2*mqp3),-(ddV1*mqp2 + ddV2*mqp4),\
                 -(ddV3*mqp1 + ddV4*mqp3),-(ddV3*mqp2 + ddV4*mqp4),\
                 -(ddV1*mqq1 + ddV2*mqq3),-(ddV1*mqq2 + ddV2*mqq4),\
                 -(ddV3*mqq1 + ddV4*mqq3),-(ddV3*mqq2 + ddV4*mqq4),\
                   mpp1,mpp2,mpp3,mpp4, mpq1,mpq2,mpq3,mpq4])*dti
    #print('residue',x+1j*pi/2,(x+1j*pi/2)**3)              
    return drdt.view(dtype=float)

def solout(t,r):
    global trajon,trajcolor
    r = r.view(dtype=complex)
    x,px,y,py,S, mpp1,mpp2,mpp3,mpp4,mpq1,mpq2,mpq3,mpq4,mqp1,mqp2,mqp3,mqp4,mqq1,mqq2,mqq3,mqq4=  r
    #if(trajon==1):
        #print('here')
    tarr.append(t)
    pxarr.append(px)
    pyarr.append(py)
    xarr.append(x)
    yarr.append(y)
    colorarr.append(trajcolor)
    
def contour_integration(xinit,yinit,pxinit,pyinit,gammaf,gammai,trajdata,S20,T,xin,pxin,yin,pyin):
    global trajon,trajcolor
    sol = ode(faux)
    sol.set_integrator('dop853',nsteps=1e10)
    sol.set_solout(solout)

    Sinit = -1j*(0.25*log(2*gammaxi/pi)) -1j*gammaxi*(xinit**2 - xin**2)+ (xinit*pxinit - xin*pxin) +\
                               -1j*(0.25*log(2*gammayi/pi)) -1j*gammayi*(yinit**2 - yin**2)+ (yinit*pyinit - yin*pyin) 
            
    y0 = array([xinit,pxinit,yinit,pyinit,Sinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    
    print('y0',y0)
    if(len(T) ==1):
        sol.set_initial_value(y0.view(dtype=float),t=0.0)
        trajdata = sol.y.view(dtype=complex)
    else:
        
        tempsol = y0
        colors = ['b','c','m','y','w','g']
        
        for tc in range(len(T)-1):
            coloriterindex = tc%len(colors)
            trajcolor = colors[coloriterindex]
            print('trajcolor',trajcolor)
            print('Time',T[tc])
            if(tc==len(T)-2):
                trajon=1
                print('started')
                
            #print(tc)
            
            #print(T[tc+1]-T[tc])
            timcom = T[tc+1]-T[tc]
            abstim = abs(timcom)
            #print('abstim',abstim)
            dti = timcom/abstim
            sol.set_initial_value(tempsol.view(dtype=float),t=0.0)
            sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
            sol.integrate(abstim)
            tempsol = sol.y.view(dtype=complex)
        trajdata=sol.y.view(dtype=complex)
        
        print('Trajdata ',trajdata[:5],'xfinal',trajdata[0])
        return trajdata

def loop_around(tcoord):
    tlb = (real(tcoord)-0.05) + 1j*(imag(tcoord)-0.05)
    trb = (real(tcoord)+0.05) + 1j*(imag(tcoord)-0.05)
    tlt = (real(tcoord)-0.05) + 1j*(imag(tcoord)+0.05)
    trt = (real(tcoord)+0.05) + 1j*(imag(tcoord)+0.05)
    
    Timecontour = [0.0,real(tlb),tlb,tlt,trt,trb,tlb,real(tlb),0.0]#trb,trt,tlt,tlb,trb,trt,tlt,tlb,
    op = open('/home/vijay/Codes/Pickle_files/Loop_contour_for_point_{}.pkl'.format(point),'wb')
    pickle.dump(Timecontour,op)
    op.close()
    trajdata = zeros(21,dtype=complex)
    loopdata = contour_integration(xinit,yinit,pxinit,pyinit,gammaf,gammai,trajdata,S20,Timecontour,x0,px0,y0,py0)
    return loopdata
                   
def finalcondition_longitudinal(xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,xin,pxin,yin,pyin):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1e3)
    Sinit = -1j*(0.25*log(2*gammaxi/pi)) -1j*gammaxi*(xinit**2 - xin**2)+ (xinit*pxinit - xin*pxin) +\
                               -1j*(0.25*log(2*gammayi/pi)) -1j*gammayi*(yinit**2 - yin**2)+ (yinit*pyinit - yin*pyin) 
            
    y0 = array([xinit,pxinit,yinit,pyinit,Sinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
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
            
def finalcondition_transversal(xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,xin,pxin,yin,pyin):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1e3)
    Sinit = -1j*(0.25*log(2*gammaxi/pi)) -1j*gammaxi*(xinit**2 - xin**2)+ (xinit*pxinit - xin*pxin) +\
                               -1j*(0.25*log(2*gammayi/pi)) -1j*gammayi*(yinit**2 - yin**2)+ (yinit*pyinit - yin*pyin) 
            
    y0 = array([xinit,pxinit,yinit,pyinit,Sinit,1.0,0.0,0.0,1.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 1.0,0.0,0.0,1.0])
    sol.set_initial_value(y0,t=0.0)

    dtim = tim[1]-tim[0]
    dtre = tre[1]-tre[0]

    print(dtim,dtre)
    

    tr = 0
    dti = 1.0
    sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    trun = dtre
    while(tr<len(tre)):
        #print('time',Time[tr][int(len(tim)/2)])
        sol.integrate(trun)
        data[tr][int(len(tim)/2)] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        tr+=1
        trun+=dtre
    
    #data[0][int(len(tim)/2)] = y0
    #print(Time[0][int(len(tim)/2)],data[0][int(len(tim)/2)][0] )
    #print(Time[0][int(len(tim)/2)])
    #print(data[:][:][0])
    
    trun=dtim
    for tr in range(len(tre)):
        sol.set_initial_value(data[tr][int(len(tim)/2)],t=0.0)
        dti = 1j
        sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
        
        trun = dtim
        for ti in range(int(len(tim)/2),len(tim)):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti])#,data[tr][ti][0])
            trun+=dtim
        
        sol.set_initial_value(data[tr][int(len(tim)/2)],t=0.0)    
        #print(Time[tr][int(len(tim)/2)],'time')
        dti = -1j
        sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
        trun = dtim
            
        for ti in arange(int(len(tim)/2),-1,-1):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti])#,data[tr][ti][0])
            trun+=dtim
            
def finalcondition_loop_longitudinal(loop,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,xin,pxin,yin,pyin):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1e3)
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

def finalcondition_loop_transversal(loop,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,xin,pxin,yin,pyin):
    sol = ode(f)
    sol.set_integrator('zvode',nsteps=1e3)
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
    

    tr = 0
    dti = 1.0
    sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
    trun = dtre
    while(tr<len(tre)):
        #print('time',Time[tr][int(len(tim)/2)])
        sol.integrate(trun)
        data[tr][int(len(tim)/2)] = sol.y
        #print(Time[0][ti],data[0][ti][0])
        tr+=1
        trun+=dtre
    
    #data[0][int(len(tim)/2)] = y0
    #print(Time[0][int(len(tim)/2)],data[0][int(len(tim)/2)][0] )
    #print(Time[0][int(len(tim)/2)])
    #print(data[:][:][0])
    
    trun=dtim
    for tr in range(len(tre)):
        sol.set_initial_value(data[tr][int(len(tim)/2)],t=0.0)
        dti = 1j
        sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
        
        trun = dtim
        for ti in range(int(len(tim)/2),len(tim)):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti])#,data[tr][ti][0])
            trun+=dtim
        
        sol.set_initial_value(data[tr][int(len(tim)/2)],t=0.0)    
        #print(Time[tr][int(len(tim)/2)],'time')
        dti = -1j
        sol.set_f_params(dti,potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4)
        trun = dtim
            
        for ti in arange(int(len(tim)/2),-1,-1):
            sol.integrate(trun)
            data[tr][ti] = sol.y
            #print(Time[tr][ti])#,data[tr][ti][0])
            trun+=dtim


        
#----------------------------------------------------------------------
Time = zeros((len(tre),len(tim))) +0j


for tr in range(len(tre)):
    for ti in range(len(tim)):
        Time[tr][ti] = tre[tr] + 1j*tim[ti]


mpp = zeros((len(tre),len(tim),2,2),dtype=complex)
mqp = zeros((len(tre),len(tim),2,2),dtype=complex)
mpq = zeros((len(tre),len(tim),2,2),dtype=complex)
mqq = zeros((len(tre),len(tim),2,2),dtype=complex)


gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
data = zeros((len(tre),len(tim),21),dtype = complex)

#----------------------------------------------------------------------

if(plot_full_time ==1):

    for i in range(len(pointarr)):
        xinit = pointarr[i][0]
        yinit = pointarr[i][1]
        pxinit = -1j*(2*gammaxi*x0 + 1j*px0 - 2*gammaxi*xinit)
        pyinit = -1j*(2*gammayi*y0 + 1j*py0 - 2*gammayi*yinit)
        branch = 'transversal'
        
        if(branch=='transversal'):
            loop = 0
            loopdata = zeros(21,dtype=complex)

            if(loop==0):
                finalcondition_loop_transversal(0,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,x0,px0,y0,py0)
                op = open('/home/vijay/Codes/Pickle_files/Singularities_double_slit_pointarr_{}_branch_2.pkl'.format(i),'wb')
                pickle.dump(data,op)
                op.close()
            else:
                tcoord = 2.72+1.28j
                loopdata= loop_around(tcoord)
                finalcondition_loop_transversal(1,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,x0,px0,y0,py0)
                op = open('/home/vijay/Codes/Pickle_files/Singularities_double_slit_pointarr_{}_branch_4.pkl'.format(i),'wb')
                pickle.dump(data,op)
                op.close()
            
        if(branch=='longitudinal'):
            loop = 0
            loopdata = zeros(21,dtype=complex)
            
            if(loop==0):
                finalcondition_loop_longitudinal(0,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,x0,px0,y0,py0)
                op = open('/home/vijay/Codes/Pickle_files/Singularities_double_slit_pointarr_{}_branch_1.pkl'.format(i),'wb')
                pickle.dump(data,op)
                op.close()
            else:
                tcoord = 2.6934+0.6056j
                loopdata= loop_around(tcoord)
                finalcondition_loop_longitudinal(1,loopdata,xinit,yinit,pxinit,pyinit,gammaf,gammai,S20,data,tre,tim,x0,px0,y0,py0)
                op = open('/home/vijay/Codes/Pickle_files/Singularities_double_slit_pointarr_{}_branch_3.pkl'.format(i),'wb')
                pickle.dump(data,op)
                op.close()
#-----------------------------------------------------------------------------------------

trajdata = zeros(21,dtype=complex)


Tusual = [0.0,10.0]
Tinterf = [0.0,3.0,3.0+0.05j,4+0.05j,4.0]
Ttunnel = [0.0,2.26,2.26+0.983j,2.45+0.983j,2.45+0.902j,2.26+0.902j,2.26,4.0]
singpoint = 2.6925+1.3642j
Tnonclassical = [0.0,3.30,3.30+0.04j,3.08+0.04j,3.08,10.0]
Tnonclassical = [0.0,3.787,3.787+0.746j,4.122+0.746j,4.122,10.0]
Tnonclassical = [0.0,2.83959,2.83959-1j]


Tnonclassical = [0.0,2.36581,2.36581+1.85322j,2.45798+1.85322j,2.8+1.85322j]
Tnonclassical = [0.0,2.3375,2.3375+1.9610j,2.8+1.9610j]
Tnonclassical = [0.0,2.8612,2.8612+0.0980j,2.9374+0.0980j,4.0+0.0980j]#,2.9374,10.0]

Tnonclassical = [0.0,2.62786,2.62786+1.1710j,2.8207+1.1710j]#,2.8207+0.9j]
Tnonclassical = [0.0,2.7536,2.7536+1.06965j,2.88978+1.06965j]
Tnonclassical = [0.0,3.7428,3.7428+0.1842j,3.9446+0.1842j,3.9446,10.0]#,3.8446,10.0]#3.7921,3.7921+0.1842j,3.8446 + 0.1846j]
Tnonclassical = [0.0,4.87,4.87+1j,5+1j]#4.9594,4.9594+1.6701j]
Tnonclassical = [0.0,3.624,3.624+0.6102j]#,3.7889+0.4102j]#4.709,4.709+1.103j,5.408+1.103j,5.408+1.548j,6.354+1.548j]
Tnonclassical = [0.0,4.265,4.265+0.5374j,4.9922+0.5374j,4.9922,10.0]#6.4295+1.0877j]#,4.702,10.0] # 4.4482,4.4482+0.4913j
Tnonclassical = [0.0,2.6052,2.6052-0.2282j,3.0285-0.2282j,3.0285,10.0]
Tnonclassical = [0.0,4.123,4.123+1.485j,4.353+1.485j,4.353,10.0]
tcoord = 3.115+1.4562j
tlb = (real(tcoord)-0.05) + 1j*(imag(tcoord)-0.05)
trb = (real(tcoord)+0.05) + 1j*(imag(tcoord)-0.05)
tlt = (real(tcoord)-0.05) + 1j*(imag(tcoord)+0.05)
trt = (real(tcoord)+0.05) + 1j*(imag(tcoord)+0.05)
    
Tnonclassical = [0.0,real(tlb),tlb,tlt,trt,trb,real(trb),10.0]
Tnonclassical = [0.0,2.6645,2.6645+0.08j,3.3+0.08j,3.3+1j,3.6+1j,3.6+0.08j,3.6,10.0]#3.0+0.9j,3.6+0.9j,3.6]#3.3,10.0]#3.2+0.7564j,]#2.9059,10.0]

Tnonclassicalreal = real(Tnonclassical)
Tnonclassicalimag = imag(Tnonclassical)

Tnonclassical = [0.0,2.84,2.84+0.1j,3.0+0.1j,3.0,10.0]
Tnonclassical = [0.0,2.73,2.73+1.23j,2.81+1.23j,2.81,10.0]
Tnonclassical = [0.0,2.63,2.63+1.33j,2.83+1.33j,2.83,10.0]
Tnonclassical = [0.0,2.33,2.33+1.96j,2.4027+1.96j,2.4027]
Tnonclassical = [0.0,2.6113,2.6113+1.5189j,]#2.7130+1.5189j,2.7130]
Tnonclassical = [0.0,2.59,2.59+1.49j,2.79+1.49j,2.79]#2.84]

#3.3336+0.08j,3.3336+0.8874j,3.5191+0.8874j,3.5191,10.0]#,3.0093+0.8874j,]#3.2242+0.8874j]
#[0.0,3.2,3.2+0.05j,3.3+0.05j,3.3,10.0]
#[0.0,3.2,3.2-0.05j,3.3-0.05j,3.3+0.05j,3.2+0.05j,3.2-0.05j,3.3-0.05j,3.3,10.0]
#[0.0,4.24,4.24+0.11j,4.24+0.2040j,4.13+0.2040j,4.13+0.11j,4.24+0.11j,4.24,10.0]

#[0.0,real(singpoint)+0.05,singpoint+0.05-0.05j,singpoint+0.05+0.05j,singpoint-0.05+0.05j,singpoint-0.05-0.05j,\
#singpoint+0.05-0.05j,real(singpoint)+0.05,10.0]
#[0.0,1.65,1.65+0.9j,1.87+0.9j,1.87,4.0]
#[0.0,2.15,2.15+1.051j,2.25+1.051j,2.25+0.957j,2.15+0.957j,2.15, 4.0]
#[0.0,2.0,2.0+0.96j,2.17+0.96j,2.17+0.87j,2.0+0.87j,2.0,4.0]#[0.0,1.88,1.88+0.2j,2.5+0.2j,2.5,4.0]
Ttarget = Tusual#Tnonclassical

op = open('/home/vijay/Codes/Pickle_files/Time_contour_for_point_{}.pkl'.format(point),'wb')
pickle.dump(Ttarget,op)
op.close()

trajdatatunnel = contour_integration(xinit,yinit,pxinit,pyinit,gammaf,gammai,trajdata,S20,Ttarget,x0,px0,y0,py0)


#---------------------------------------------------------------------------------------------

pxarrinv = 1/np.array(pxarr)
pyarrinv = 1/np.array(pyarr)
#print('pole',trajdatatunnel[1],trajdatatunnel[1]*(trajdatatunnel[0]-1j*pi/2)**4)
xarr = array(xarr) #- pi*1j/2#trajdata[0]


n=0.5 #x doesn't quite follow t^0.5, it goes as t^0.55 
# plt.plot(abs(tiarr),(abs(pxarrinv))*(abs(tiarr[5])**n)/abs(pxarrinv[5]),color='g',marker='x')#,s=2)#-abs(pxarrinv[len(pxarrinv)-1])
# plt.plot(abs(tiarr),(abs(tiarr))**n,color='b')#,s=2)
# #plt.plot(abs(tiarr),2*(abs(tiarr))**0.5,color='m')
# plt.show()

if(plot==1):

    xre= np.arange(-20,20,0.05)
    xim = np.arange(-20,20,0.05)
    #print(xarr)
    xreal,ximag = np.meshgrid(xre,xim)
    # for xc in range(len(yarr)):
        # plt.scatter(real(yarr[xc]),imag(yarr[xc]))
    # Complex_plotter.plotcomplex(potential(0.0,xreal+1j*ximag),1,1,xre[0],xre[len(xre)-1],xim[0],xim[len(xim)-1])
    # plt.show()
    # for xc in range(len(xarr)):
        # plt.scatter(real(xarr[xc]),imag(xarr[xc]))
    # Complex_plotter.plotcomplex(potential(xreal+1j*ximag,0.0),1,1,xre[0],xre[len(xre)-1],xim[0],xim[len(xim)-1])
    # plt.show()

    plt.figure(1)
    plt.suptitle('Trajectory for initial conditions: ({},{})  Coordinate index: {}'.format(pointarr[0][0],pointarr[0][1],coordin))
    for xc in range(len(xarr)):
        plt.scatter(real(xarr[xc]),real(yarr[xc]),color=colorarr[xc])
    Complex_plotter.plotcomplex(potential(xreal,ximag),1,1,xre[0],xre[len(xre)-1],xim[0],xim[len(xim)-1])
    plt.show()
    #plt.show()
#fig = plt.figure(2)
#fig.figsize = fig_size
#ax = fig.gca(projection='3d')
#ax.set_xlim3d(1.4,1.8)
#ax.set_ylim3d(-0.2,0.2)
#xfinbar[timereal<1.4]=nan
#xfinbar[timereal>1.8]=nan
# for i in range(5):
    # magnitudeplot = np.minimum(abs(potential(xreal+1j*ximag,i*1j))/10.0,1.0)
    # phaseplot = (angle(-potential(xreal+1j*ximag,i*1j))+np.pi)/2/np.pi
    # #phaseplot-= min(min(x) for x in phaseplot)
    # #print(phaseplot.shape,magnitudeplot.shape)
    # z_data_rgb = rgb((phaseplot))
    # datap=z_data_rgb*magnitudeplot[:,:,np.newaxis]
    # datap[...,3]=1.
    # surf = ax.plot_surface(xreal,ximag,i*ones_like(ximag),rstride=10, cstride=10,facecolors=datap,shade=False)#,antialiased=True)#, alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)

# for xc in range(len(xarr)):
    # plt.scatter(real(xarr[xc]),imag(xarr[xc]),imag(yarr[xc]))
#plt.show()


####-------------------------------------
# fig = plt.figure(2)
# ax = fig.gca(projection='3d')

# for xc in range(len(xarr)):
    # ax.scatter(real(xarr[xc]),imag(xarr[xc]),real(yarr[xc]))
    
# plt.show()
####-------------------------------------

