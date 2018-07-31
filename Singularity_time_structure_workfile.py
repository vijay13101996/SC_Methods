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
import scipy.ndimage
from Parameter_file_Singularity_time_structure import *
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
import Complex_plotter
#from Singularity_time_structure_2D import Ttarget
#from mayavi import *

#-----------------------------------------------------
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
#--------------------------------------------------------------------------------

#print(xinit,yinit)
print(pointarr)

#print(tre,tim)

xfinbar = zeros((len(tre),len(tim)),dtype=complex)
yfinbar = zeros((len(tre),len(tim)),dtype=complex)
pxfinbar = zeros((len(tre),len(tim)),dtype=complex)
pyfinbar = zeros((len(tre),len(tim)),dtype=complex)

xbar = x0*ones((len(tre),len(tim)),dtype=complex)
ybar = y0*ones((len(tre),len(tim)),dtype=complex)
pxbar = px0*ones((len(tre),len(tim)),dtype=complex)
pybar = py0*ones((len(tre),len(tim)),dtype=complex)

Sfinbar = zeros((len(tre),len(tim)),dtype=complex)

mpp = zeros((len(tre),len(tim),2,2),dtype=complex)
mqp = zeros((len(tre),len(tim),2,2),dtype=complex)
mpq = zeros((len(tre),len(tim),2,2),dtype=complex)
mqq = zeros((len(tre),len(tim),2,2),dtype=complex)

gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]
S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
data = zeros((len(tre),len(tim),20),dtype = complex)
j0 = zeros((len(tre),len(tim)),dtype=complex)
j0x = zeros((len(tre),len(tim)),dtype=complex)

pfi = zeros((len(tre),len(tim),2),dtype=complex)
pin = zeros((len(tre),len(tim),2),dtype=complex)
Sin= zeros((len(tre),len(tim),2),dtype=complex)
Sfi= zeros((len(tre),len(tim),2),dtype=complex)

dsigma= zeros((len(tre),len(tim),2),dtype=complex)
dsigmax= zeros((len(tre),len(tim)),dtype=complex)
dsigmay= zeros((len(tre),len(tim)),dtype=complex)
dsigmamag= zeros((len(tre),len(tim)),dtype=complex)

timereal,timeimag = meshgrid(tre,tim)

Time = zeros((len(tre),len(tim)),dtype=complex)
t1 =(1.7600502483227527+0.2454058275120668j)
t1 = 3.616+0.558j#1.025-0.0437j#3.595+0.558j #1.025-0.0437j
t1 = 1.26-1.35j
for tr in range(len(tre)):
    for ti in range(len(tim)):
        Time[tr][ti] = tre[tr] + 1j*tim[ti]
        #print(Time[tr][ti])

for i in range(len(pointarr)):
    for j in [2]:        
        op = open('/home/vijay/Codes/Pickle_files/Singularities_double_slit_pointarr_{}_branch_{}.pkl'.format(i,j),'rb')
        data = pickle.load(op)
        op.close()




        xfinbar,pxfinbar,yfinbar,pyfinbar,Sfinbar = data[:,:,:5].transpose((2,0,1))
        mpp=data[:,:,5:9].reshape(data.shape[:2]+(2,2))
        mpq=data[:,:,9:13].reshape(data.shape[:2]+(2,2))
        mqp=data[:,:,13:17].reshape(data.shape[:2]+(2,2))
        mqq=data[:,:,17:].reshape(data.shape[:2]+(2,2))


        DET=lambda M:M[...,0,0]*M[...,1,1]-M[...,0,1]*M[...,1,0]

        xsi = 2*np.matmul(gammaf,mqq) + 2*np.matmul(np.matmul(gammaf,mqp),S20) -1j*mpq - 1j*np.matmul(mpp,S20)
        xsixx = xsi[...,0,0]
        xsixy = xsi[...,0,1]
        xsiyx = xsi[...,1,0]
        xsiyy = xsi[...,1,1]
        j0 = 0.25/(np.linalg.det(gammaf))*DET(xsi)

        j0xx = 0.25/(np.linalg.det(gammaf))*xsixx
        j0xy = 0.25/(np.linalg.det(gammaf))*xsixy
        j0yx = 0.25/(np.linalg.det(gammaf))*xsiyx
        j0yy = 0.25/(np.linalg.det(gammaf))*xsiyy

        

            
        J = -(0.25/np.linalg.det(gammaf))*abs(j0)**2 + 0j
        Integrand = -(-0.25)**2*(4/(pi**2*np.linalg.det(gammaf)))**0.75*J/(j0**0.5)

        xfin = real((2*gammaxf*xfinbar - 1j*pxfinbar))/(2*(gammaxf))
        yfin = real((2*gammayf*yfinbar - 1j*pyfinbar))/(2*(gammayf))
        pxfin = -imag((2*gammaxf*xfinbar - 1j*pxfinbar)) 
        pyfin = -imag((2*gammayf*yfinbar - 1j*pyfinbar))
        
        ExponentS = 1j*Sfinbar
        Exponentx = -gammaxf*xfinbar**2 + gammaxf*xfin**2 +\
                         1j*pxfinbar*xfin - 1j*pxfin*xfinbar
        Exponenty = -gammayf*yfinbar**2 + gammayf*yfin**2 +\
                         1j*pyfinbar*yfin - 1j*pyfin*yfinbar
        sigma = Exponentx + Exponenty + ExponentS
        Joverlap = Integrand*exp(sigma)

        # for tr in range(len(tre)):
            # for ti in range(len(tim)):
                # pfi[tr][ti] = [pxfinbar[tr][ti],pyfinbar[tr][ti]]
                # pin[tr][ti] = [pxbar[tr][ti],pybar[tr][ti]]
                # Sin[tr][ti] = 1j*np.matmul(pfi[tr][ti],xsi[tr][ti])
                # Sfi[tr][ti] = 1j*array(pin[tr][ti])
                # dsigma[tr][ti] = Sfi[tr][ti]-Sin[tr][ti]
                # dsigmax[tr][ti] = dsigma[tr][ti][0]
                # dsigmay[tr][ti] = dsigma[tr][ti][1]
                # dsigmamag[tr][ti] = dsigmax[tr][ti]**2 + dsigmay[tr][ti]**2
                 
         #print(dxsi,pfi,pin,1j*np.matmul(pfi,dxsi))

##        Sin = 1j*np.matmul(xsi,pfi)
##        Sfi = 1j*array(pin)
##
##         #print(Sfi,Sin)
##
##        dsigma = Sfi-Sin
##
##        dsigmax = dsigma[0]
##        dsigmay = dsigma[1]
##        dsigmamag = dsigmax**2 + dsigmay**2
##

##op = open('Singularities_double_slit_corrected_point_{}_J0.pkl'.format(point),'wb')
##pickle.dump(j0,op)
##op.close()
##
##op = open('Singularities_double_slit_corrected_point_{}_J0xx.pkl'.format(point),'wb')
##pickle.dump(j0xx,op)
##op.close()
##
##op = open('Singularities_double_slit_corrected_point_{}_J0xy.pkl'.format(point),'wb')
##pickle.dump(j0xy,op)
##op.close()
##
##op = open('Singularities_double_slit_corrected_point_{}_J0yx.pkl'.format(point),'wb')
##pickle.dump(j0yx,op)
##op.close()
##
##op = open('Singularities_double_slit_corrected_point_{}_J0yy.pkl'.format(point),'wb')
##pickle.dump(j0yy,op)
##op.close()


        pxfinbar = matrix.transpose(pxfinbar)
        pyfinbar = matrix.transpose(pyfinbar)

        yfinbar = matrix.transpose(yfinbar)
        xfinbar = matrix.transpose(xfinbar)
        Sfinbar = matrix.transpose(Sfinbar)
        Joverlap = matrix.transpose(Joverlap)
        j0 = matrix.transpose(j0)
        j0xx = matrix.transpose(j0xx)
        j0xy = matrix.transpose(j0xy)
        j0yx = matrix.transpose(j0yx)
        j0yy = matrix.transpose(j0yy)
        dsigmamag = matrix.transpose(dsigmamag)
        print(pointarr[i][0],pointarr[i][1])

        Time = matrix.transpose(Time)
        
        op = open('/home/vijay/Codes/Pickle_files/Time_contour_for_point_{}.pkl'.format(point),'rb')
        Ttarget = pickle.load(op)
        op.close()
    
        Timecontour = Ttarget
        
        
        

        ax=plt.subplots()
        plt.suptitle('Singularity time structure for initial conditions: ({},{})  Coordinate index: {}'.format(pointarr[0][0],pointarr[0][1],coordin))
        
        ax = Complex_plotter.plotcomplex(pyfinbar,1,10,tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1])
        for tc in range(len(Timecontour)):
            plt.scatter(real(Timecontour[tc]),imag(Timecontour[tc]))
            
        plt.plot(real(Timecontour),imag(Timecontour))
        plt.show()
        # #fig = plt.figure(2)
        # #fig.figsize = fig_size
        # ax = fig.gca(projection='3d')
        # ax.set_xlim3d(1.4,1.8)
        # ax.set_ylim3d(-0.2,0.2)
        # xfinbar[timereal<1.4]=nan
        # xfinbar[timereal>1.8]=nan
        # magnitudeplot = np.minimum(abs(pxfinbar)/10.0,1.0)
        # phaseplot = (angle(-pxfinbar)+np.pi)/2/np.pi
        # #phaseplot-= min(min(x) for x in phaseplot)
        # #print(phaseplot.shape,magnitudeplot.shape)
        # z_data_rgb = rgb((phaseplot))
        # datap=z_data_rgb*magnitudeplot[:,:,np.newaxis]
        # datap[...,3]=1.
        #surf = ax.plot_surface(timereal,timeimag,real(xfinbar),rstride=5, cstride=5,facecolors=datap,shade=False)#,antialiased=True)#, alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
        #surf = ax.plot_surface(Time, Y, Fric_map, alpha = 1, rstride=1, cstride=1, cmap=cm.autumn,linewidth=0.5, antialiased=True)
        
        #plt.show()
        # for tax in range(len(singularity_times)):
            # plt.scatter(real(singularity_times[tax]),imag(singularity_times[tax]))
        #plt.show()
        #Complex_plotter.plotcomplex(Joverlap,1,10,tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1])
        #Complex_plotter.plotcomplex(dsigmamag,1,1000,tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1])


#plt.figure(2)        
#plt.show()
count=0
total = 0
for tr in range(len(tre)):
    for ti in range(len(tim)):
        total+=1
        
        
        #print(Time[tr][ti],xfinbar[tr][ti])
        if (isnan(pxfinbar[tr][ti]) == False):
            #plt.scatter(tre[tr],tim[ti])
            count+=1
            
print(count,total)
#plt.show()
        
        
        


