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
import random
import warnings
from matplotlib import pyplot as plt
import Singularity_time_finder
import time
from Parameter_file_Singularity_time_structure import *


start_time = time.time()
xtarget_list = [(0.0-1.5667183672681677j)]#,(0.0+1.5667183672681677j),(0.0-4.7067885713578885j),(0.0+4.7067885713578885j) ]

singpoints = [[] for i in range(len(xtarget_list))]
count=0

for xtrgt in xtarget_list:
    Singularity_time_finder.set_target(xtrgt)
    singpoints[count] = Singularity_time_finder.singtime_finder(xinit,yinit,pxinit,pyinit)
    print('Singularities',count,singpoints[count])
    count+=1
 
print("Time taken:",time.time()-start_time) 
    
op = open('Singularities_times_for_point_{}.pkl'.format(point),'wb')
pickle.dump(singpoints[count-1],op)
op.close()

