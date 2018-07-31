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


#----------------------------------------------------------------
##dt=0.05##
##tre = arange(-2.0,2.0 + 0.1*dt,dt)
##tim = arange(-2.0,2.0+0.1*dt,dt)
##
##treal,timag = meshgrid(tre,tim)
##
##t = treal + 1j*timag
##function=t
##print(function)

def plotcomplex(function,height):

    norm = mcolors.Normalize(0,height)

    magnitudeplot = abs(function)
    print(magnitudeplot)
    phaseplot = (angle(-function)+np.pi)/2/np.pi
    #phaseplot-= min(min(x) for x in phaseplot)
    print(phaseplot.shape,magnitudeplot.shape)
    z_data_rgb = rgb((phaseplot))
    print('max,min',max(max(x) for x in phaseplot),min(min(x) for x in phaseplot))
    print(z_data_rgb.shape)
    intensity = norm(magnitudeplot)
    print(intensity.shape)
    ls = mcolors.LightSource()

    datap=ls.blend_hsv(z_data_rgb, intensity)*np.minimum(magnitudeplot[:,:,np.newaxis],1)
    #datap=ls.blend_hsv(z_data_rgb, intensity)*mag[:,:,np.newaxis]
    print (datap.shape,datap[:,:,0].max())
    plt.imshow(datap,origin = 'lower',extent=[tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1]],interpolation='none')
    #plt.savefig('J0.png')
    plt.show()


#plotcomplex(function,5.0)
