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
from Parameter_file_Singularity_time_structure_1D import *
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
import Complex_plotter

#-----------------------------------------------------------
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
#-----------------------------------------------------------
     
qfinbar = zeros((len(tre),len(tim)),dtype=complex)
pfinbar = zeros((len(tre),len(tim)),dtype=complex)

timereal,timeimag = meshgrid(tre,tim)
Time = zeros((len(tre),len(tim)),dtype=complex)

#---------------------------------------------------------------------------------------
         
op = open('/home/vijay/Codes/Pickle_files/Singularities_{}_point_{}.pkl'.format(potkey,1),'rb')
data = pickle.load(op)
op.close()

qfinbar,pfinbar = data[:,:,:2].transpose((2,0,1))

qfinbar = matrix.transpose(qfinbar)
pfinbar = matrix.transpose(pfinbar)

op = open('/home/vijay/Codes/Pickle_files/Time_contour_1D_for_point_1.pkl','rb')
Timecontour = pickle.load(op)
op.close()

for tc in range(len(Timecontour)):
            plt.scatter(real(Timecontour[tc]),imag(Timecontour[tc]))

Complex_plotter.plotcomplex(qfinbar,1,10,tre[0],tre[len(tre)-1],tim[0],tim[len(tim)-1])
plt.show()
        
        
