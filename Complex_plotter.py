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
#from Singularity_time_structure_2D import V0,w,E,m,gammaxi,gammayi,gammaxf,gammayf,dt,tre,tim,treal,timag,xinit,yinit,pxinit,pyinit
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
from pylab import *

#----------------------------------------------------------------

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



#[c('red'), c('yellow'), 0.333, c('yellow'), c('green'), c('cyan'), 0.5, c('cyan'),
    # c('blue'), c('magenta'), 0.833, c('magenta'), c('red')])

def click(event):
   """If the left mouse button is pressed: draw a little square. """
   tb = get_current_fig_manager().toolbar
   if event.button==1 and event.inaxes and tb.mode == '':
       x,y = event.xdata,event.ydata
       plot([x],[y],'rs')
       draw()
#----------------------------------------------------------------

def plotcomplex(function,height,scale,xrl,xrr,xil,xir):

    #norm = mcolors.Normalize(0,height)

    magnitudeplot = np.minimum(abs(function)/scale,1.0)
    #magnitudeplot = abs(function)
    #print(magnitudeplot)
    phaseplot = (angle(-function)+np.pi)/2/np.pi
    #phaseplot-= min(min(x) for x in phaseplot)
    #print(phaseplot.shape,magnitudeplot.shape)
    z_data_rgb = rgb((phaseplot))
    #print('max,min',max(max(x) for x in phaseplot),min(min(x) for x in phaseplot))
    #print(z_data_rgb.shape)
    intensity = (magnitudeplot)
    #print(intensity.shape)
    #ls = mcolors.LightSource(0,90)

    plt.figure(1)

    

    #datap=ls.blend_hsv(z_data_rgb, intensity)#*magnitudeplot[:,:,np.newaxis]
    datap=z_data_rgb*magnitudeplot[:,:,np.newaxis]
    datap[...,3]=1.
    #datap=ls.blend_hsv(z_data_rgb, intensity)*np.minimum(magnitudeplot[:,:,np.newaxis]/scale,1)
    #datap=ls.blend_hsv(z_data_rgb, intensity)*mag[:,:,np.newaxis]
    #print (datap.shape,datap[:,:,0].max())
    
    ax = plt.imshow(datap,origin = 'lower',extent=[xrl,xrr,xil,xir],interpolation='none')
    gca().set_autoscale_on(False)
    connect('button_press_event',click)

    T1 = [0.0,1.0+0.1j,1.2+0.1j,1.2-0.1j,3.0] # Encloses one caustic point, the left most one, but no singularities
    T2 = [0.0,1.0+0.4j,3.0]#Encloses both singularity 
    T3 = [0.0,1.0-0.2j,3.0]#Encloses no caustic point, no singularity
    T4 = [0.0,1.0+0.4j,1.2+0.4j,1.2-0.1j,3.0-0.1j,3.0]#Encloses left singularity
    T5 = [0.0,1.0+0.1j,1.2+0.1j,1.2+0.3j,3.0+0.3j,3.0]#Encloses right singularity
    T6 = [0.0,1.0+0.1j,1.2+0.1j,1.2+0.3j,1+0.3j,1.0+0.1j,1.4,3.0]#Goes around left singularity
    T7 = [0.0,1.0+0.1j,1.2+0.1j,1.2-0.1j,1.8-0.1j,1.8+0.3j,1.2+0.3j,1.2-0.1j,3.0]#Goes around right singularity
##    plt.plot(real(T1),imag(T1),color='r')
##    plt.plot(real(T2),imag(T2),color='b')
##    plt.plot(real(T3),imag(T3),color='g')
##    plt.plot(real(T4),imag(T4),color='m')
##    plt.plot(real(T5),imag(T5),color='y')
##    plt.plot(real(T6),imag(T6),color='c')
##    plt.plot(real(T7),imag(T7),color='k')
    
    #plt.savefig('J0.png')
    #return fig
    #plt.show()
    return ax


#plotcomplex(function,5.0)
