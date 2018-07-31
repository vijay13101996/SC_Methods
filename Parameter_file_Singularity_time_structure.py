
import numpy as np
from numpy import *
#import trial_code
import matplotlib.backends.backend_tkagg as backend
from matplotlib.figure import Figure
import pickle, pprint
from matplotlib import pyplot as plt

V0 = 16.0
w = 4.0
E = 1.0
m = 1.0
Vc=0.0#0.8776
print('Vc',Vc)

tPrint=0
tPrintSpace=0.01

gammaxi = 0.5
gammayi = 0.5
gammaxf = 0.5
gammayf = 0.5

div = 400
dt=0.02
tre = np.arange(0.0,4+0.1*dt,dt)
tim = np.arange(-2.,2.+0.1*dt,dt)

treal,timag = np.meshgrid(tre,tim)

x0 = -5.0
y0 = 0.0
px0 = 3.0
py0 = 1.0
##xinit = -5.0
##yinit = 0.0
##pxinit = 3.0
##pyinit = 1.0

lamda = -0.11
ws=1.0
wu=1.1
A=0.01
B=1

potkey = 'gaspard'#'diffraction'#'reaction_coordinate'#'

#-----------------------------------------
##(-3.4-1.6j) -0.4j     # Magnitude >50
##(-3.4-1.6j) 0.4j
##(-7-0.8j) (-0.8+1.6j)
##(-7-0.8j) (0.8-1.6j)
##(-4.2-0.4j) (-0.4+0.4j)
##(-4.2-0.4j) (0.4-0.4j)
##(-4.6+0j) -1.2j
##(-4.6+0j) 1.2j

    ### This is another set - magnitude between 20 and 50
##2655 (-4.6-1.2j) 0j     # ok
##3070 (-7-0.8j) (0.8-2j) # ok, seems alright
##3555 (-5-0.8j) 0j       # Reflecting, but ok
##3655 (-4.6-0.8j) 0j     # Borderline
##4060 (-7-0.4j) (0.4-2j) # Seems alright, but further check required : Grazing through central hump, after reflection.
##4555 (-5-0.4j) 0j       # ok
##4655 (-4.6-0.4j) 0j     # ok
##5527 (-5+0j) (-1.2+0.8j)# Seems good, but further check required: A grazing tangential collission. Needs to be thought about.
##5555 (-5+0j) 0j         # ok
##5583 (-5+0j) (1.2-0.8j) # Reflection: Above
##5680 (-4.6+0j) (1.2-2j) # Problem only after time 3.0
##7436 (-5.4+0.8j) (-0.8+0.4j) # Probably bad : A deep collision
##7474 (-5.4+0.8j) (0.8-0.4j) # Reflection

# To remove: [4060,5527,5583,7436,7474]

#----------------------------------------- Obtained from N=10...
 



point0 = [x0,y0]   # temp2 - [x0,y0-1j], temp1 = [x0+0.5,y0-1j], temp3 = [x0,y0], temp4 = [x0+0.5,y0]

#(-4.35 + 0.27j,0.0)     (-5.898 - 1.517j,0.0) point 6,7 respectively. They are apparently forming an order 1 caustic and are symmetrically positioned in a submanifold

# -3.785 +0.507j  Point 8,Seems to be hitting a singularity

#-6.199+0.600j Point 9, Seems to be in the stokes' divergence region

#(-6.8-1.4j) (-1+0.2j)  Point 10, These points cause divergence in the final wave-function - Hitting the weird singularity time.

#(-4.6-0.8j) (1.8-1.2j) Point 11, Same as point 10 - Hitting a singularity

#(-6.2-1.8j) (0.8-1.8j) Point 12, Same as above - Hitting a weird singularity

# (-5.4-1.8j) (1.6-1.8j) Point 13 Same as above - Close to the weird singularity, but not quite head-on

# (-4.8+1j) (0.6+1.2j) Point 14 This point is currently removed by the filter inside the integration step.

#(-5.4-1.4j) (0.4-0.4j) Point 15. Same as above

point10 = [-6.8-1.4j,-1+0.2j]
point11 = [-4.6-0.8j,1.8-1.2j]
point12 = [(-6.2-1.8j),(0.8-1.8j)]


point13 = [(-6-1.8j) ,(0.2-0.2j)]

point14 = [(-5.2-1.6j), (-1.6+1.4j)]
point15 = [(-4-1j),(0.2-0.2j)]
point16 = [(-3.8-0.6j), (-0.2+1.4j)]
point17= [(-4.6-0.2j), (-1.6+0.8j)]
point18 = [(-6.8-1.2j), (0.4-0.2j)]
point19 = [(-5-0.6j), (-1.4+1.2j)]
point20 = [(-5.6+0.6j),(0.4-0.4j)]

point21 = [(-5.2+1.6j),0.4j]
point22 = [(-4.6+1.8j),(0.8-0.6j)]
point23 = [-4.55+0.86j,-0.2+0.2j]   #This is a point where j0 diverges, for no apparent reason, while the Joverlap vanishes. This is simply baffling,

#On one hand, it isn't clear why there is a singularity time occuring due to an incidence with the potential wall. This is a confusing feature here, and apparently,
#no singularities in the potential seems to be nearby the path of the trajectory. 

#This point can be used to show how y component distorts the riemann sheet.

point24 = [-6.509-0.895j,-0.2+0.2j]

point25 = [-5.0+0.2j,0.0]
point26=  [(-4.2-0.4j),(0.4-0.4j)]

point27 = [-6.452+0.56j,0.8j]

point28 = [(-3.4-1.6j), 0.4j]
point29 = [(-5-2j), (-0.2+0.2j)]
point30 = [(-4.2+0.4j) ,(-0.8-2j)]
point31 =[(-5.6+1.2j),(0.6-0.4j)]
point32 =[(-5-0.6j),(-1.4+1.2j)]

point2655 =  [(-4.6-1.2j) ,0j] # Stokes
point3070 =  [(-7-0.8j) ,(0.8-2j)]
point3555 =  [(-5-0.8j), 0j]#Might be good
point3655 =  [(-4.6-0.8j) ,0j]#
point4060 =  [(-7-0.4j),(0.4-2j)]#Stokes
point4555 =  [(-5-0.4j), 0j]#Might be good
point4655 =  [(-4.6-0.4j), 0j]#
point5527 =  [(-5+0j),(-1.2+0.8j)]#Might be good
point5555 =  [(-5+0j),0j]#Good
point5583 =  [(-5+0j),(1.2-0.8j)]#Good
point5680 =  [(-4.6+0j),(1.2-2j)]#Might be good   [5555,3555,3655,4555,4655,5680]


point33 = [-4.737+0.175j,-1.2j]
point34 = [-4.529 +0.982j,0.0]
point35 = [0.0,0.8-1.25j]
point36 = [-4.779-1.724j,-0.8+0.8j]
point37 = [-5.0,(0.0-1j)*8/8]  #This initial leads to the trajectory going through the slit
point38 = [-5.0,(2.0-4j)*8/16]
point39 = [-5.0,-3.88-0.54j]
point40 = [-5.0,-0.97-0.458j]
point41 = [-5.0,-0.18+0.96j]
point42 = [-5.0,(-0.557+0.557j)*8/8]#(-0.6+0.55j)*0/8]#-0.603+0.52j]#-0.640+0.440j]#(-0.557+0.557j)*8/8]#
#Region to target: (-0.6,-0.5)x(0.55j,0.65j)
point41 = [-5.0,-0.18+0.96j] # Transmitted region, 3 bounces...
point42 = [-5.0,(-0.557+0.557j)*8/8] # Transmitted region, 2 bounces...
point43 = [-5.0,(-0.956+0.494j)*8/8] # Transmitted region, one bounce.
point44 =[-5.0,-0.418+0.835j]
point45  =[-5.0,(-0.424906+1.20407j)*8/8]#-0.458+1.255j] #On the borders of the Reflecting region....
point46 = [-5.0,-0.416878+1.19604j] #On the borders of the Transmitted region with 3 bounces...
point47 = [-5.0,-0.910+0.4609j] # On single bounce borders
point48 = [-5.0,-0.608+0.4609j] # On double bounce borders

c1 = 1.0
c2 = 1 - c1
point49 = c1*np.array(point47) + c2*np.array(point48)

point50 = [-5.0,-0.9114+0.4732j]
point51 = [-5.0,-0.3396 +1.0611j]
point52 = [-5.0,-0.5469 +0.5946j]
point53 = c1*np.array(point51) + c2*np.array(point52)
point54 = [-5.0,-0.4278+0.8727j]#-0.854+0.734j]#-0.4278+0.8727j]#
point55 = [-5.0,-1.133 + 0.897j]#-1.1582+0.8118j]#-1.133 + 0.897j
point56 = [-5.0,-1.186+0.738j]#-1.0503 +0.6992j]-1.186+0.738j
point57 = [-5.0,-1.177+0.6041j]#-1.1927+0.6666j]#-1.009+0.5454j]
point58 = [-5.0,-0.6391+0.6801j]
point59 = [-5.0,-0.5520+1.1257j]
point60 = [-5.0,(-0.4192 +0.9166j)*0/8]#-0.5989+1.2447j] -0.4192 +0.9166j
point61 = [-5.0,-0.6056 +1.1068j]#-0.5468+1.1222j] -0.6056 +1.1068j
point62 = [-5.0,-0.9218+0.8177j]
point63 = [-5.0,-0.864+0.958j]

point64 = [-5.0,-0.760+0.966j]
point65 = [-5.0,-0.736+0.966j]
point66 = [-5.0,-0.749+0.948j]

xc = np.array([0.0,0.0213,0.0243,0.0260,0.02669,0.02815,0.02962,0.0305,0.0511,0.1053,0.1333,0.1660,0.1846,0.1893,0.2080,0.1800,0.2033,0.1660])
yc = np.array([0.0,0.1446,0.2482,0.2900,0.31850,0.35760,0.40500,0.4572,0.4619,0.5319,0.6018,0.6112,0.6765,0.7278,0.7698,0.9191,1.1198,1.1804])*-1

coordin = len(xc)-2
point67 = [-5.0,xc[coordin]+1j*yc[coordin]]
point68 = [-5.0,0.1994-0.7498j]
point69 = [-5.0,0.3567-0.4250j]
point70 = [-5.0,0.02373-0.2194j]

xc = np.array([0.19-1.03j,0.15-1.05j,0.1-1.08j,0.11-1.12j,0.14-1.15j,0.17-1.18j])
xc2 = np.array([0.71-0.30j,0.725-0.396j,0.727-0.431j,0.738-0.555j,0.617-0.584j,0.574-0.659j,0.601-0.737j,0.648-0.741j])
xc3 = np.array([0.71-0.30j,0.725-0.396j,0.727-0.431j,0.738-0.555j,0.852-0.600j,0.957-0.659j,1.108-0.709j,1.154-0.672j,1.158-0.595j,1.150-0.524j])
coordin = 7
point71 = [-5.0,xc2[coordin]]
point72 = [-5.0,0.718-0.666j]
point73 = [-5.0,0.499-0.799j]
point74 = [-0.2+0j,-0.2+0j]
 
point75 =  [-5.0,-2.18-0.67j] 
point76 = [-5.0,-1.32+1.06j]
point77 = [-5.0,0.49+2.42j]
point78 = [-5.0,2.5+3.43j]

coordin = 3
xc4 = np.array([0.50+0.25j,0.39+0.37j,0.39+0.91j,0.39+1.04,0.42+1.22j,0.51+1.34j,0.66+1.41j,1.08+0.62j])
point79 = [-5.0,0]
xinit = point79[0]#(-6.1+1.3j,0), (-5.5+0.7j,0), (-6.5 + 0.3j,0),(-6.42 + 1.48j,0), (-5.823 + 1.411j,-0.8),(-6.547 -1.365j,1.0)point 1,2,3,4 respectively
yinit = point79[1]
pxinit = -1j*(2*gammaxi*x0 + 1j*px0 - 2*gammaxi*xinit)
pyinit = -1j*(2*gammayi*y0 + 1j*py0 - 2*gammayi*yinit)

print('px',pxinit,'py',pyinit)
point = 79
#pointarr = [point2655,point3070,point3555,point3655,point4060,point4555,point4655,point5527,point5555,point5583,point5680]
pointarr = [point79]
##op = open('pointarr_N_20.pkl','rb')
##pointarr = pickle.load(op)
##op.close()
#
####print(len(pointarr))
##k = 978
##xinit = pointarr[k][0]#(-6.1+1.3j,0), (-5.5+0.7j,0), (-6.5 + 0.3j,0),(-6.42 + 1.48j,0), (-5.823 + 1.411j,-0.8),(-6.547 -1.365j,1.0)point 1,2,3,4 respectively
##yinit = pointarr[k][1]
##pxinit = -1j*(2*gammaxi*x0 + 1j*px0 - 2*gammaxi*xinit)
##pyinit = -1j*(2*gammayi*y0 + 1j*py0 - 2*gammayi*yinit)
##
##
##pointarr = [pointarr[k]]


#(-4.4+0.4j) (0.4-0.6j)  - If at all, something exhibits 'dynamic singularity',
#this should. Check up

if(potkey=='double_slit'):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4
    def potential(x,y):
                return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0) + Vc)*E/(np.cosh(x)**2)

    def dxpotential(x,y):
                return -2*E*(V0 + Vc - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*np.sinh(x)/np.cosh(x)**3

    def dypotential(x,y):
                return E*(-1.0*m*w**2*y + m**2*w**4*y**3/(4*V0))/np.cosh(x)**2

    def ddpotential1(x,y):
                return E*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)*(16*V0 + 16*Vc - 8.0*m*w**2*y**2 + m**2*w**4*y**4/V0)/(8*np.cosh(x)**2)

    def ddpotential2(x,y):
                return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

    def ddpotential3(x,y):
                return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

    def ddpotential4(x,y):
                return -E*m*w**2*(1 - 3*m*w**2*y**2/(4*V0))/np.cosh(x)**2

# def potential(x,y):
            # return (A*x**2 + B*y**2 + 1/(x-y)**2)
    
# def dxpotential(x,y):
            # return (2*A*x - 2/(x - y)**3)

# def dypotential(x,y):
            # return (2*B*y + 2/(x - y)**3)

# def ddpotential1(x,y):
            # return (2*(A + 3/(x - y)**4))
            
# def ddpotential2(x,y):
            # return (-6/(x - y)**4)
            
# def ddpotential3(x,y):
            # return (-6/(x - y)**4)
            
# def ddpotential4(x,y):
            # return (2*(B + 3/(x - y)**4))
            
if(potkey=='henon_heiles'):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4  
    def potential(x,y):
                return (0.5*(x**2+y**2) + lamda*(x**2*y - y**3/3))
                
    def dxpotential(x,y):
                return (2*lamda*x*y + 1.0*x)
                
    def dypotential(x,y):
                return (lamda*(x**2 - y**2) + 1.0*y)

    def ddpotential1(x,y):
                return (2*lamda*y + 1.0)
                
    def ddpotential2(x,y):
                return (2*lamda*x)
                
    def ddpotential3(x,y):
                return (2*lamda*x)

    def ddpotential4(x,y):
                return (-2*lamda*y + 1.0)

if(potkey=='Davis_Heller'):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4  
    def potential(x,y):
                return (0.5*ws**2*x**2 + 0.5*wu**2*y**2 + lamda*y**2*x)
                
    def dxpotential(x,y):
                return lamda*y**2 + 1.0*ws**2*x
                
    def dypotential(x,y):
                return 2*lamda*x*y + 1.0*wu**2*y

    def ddpotential1(x,y):
                return 1.0*ws**2
                
    def ddpotential2(x,y):
                return 2*lamda*y
                
    def ddpotential3(x,y):
                return 2*lamda*y

    def ddpotential4(x,y):
                return 2*lamda*x + 1.0*wu**2

if(potkey=='reaction_coordinate'):

        global w,E,Vc,D,alpha
        
        alpha = 1.0
        D = 1.0
        w = 4.0
        E = 1.0
        Vc = 8.0

        def potential(x,y):
            return E*(D*(1 - np.exp(-alpha*y))**2 + Vc)/np.cosh(x)**2
        
        def dxpotential(x,y):
            return -2*E*(D*(1 - np.exp(-alpha*y))**2 + Vc)*np.sinh(x)/np.cosh(x)**3
        
        def dypotential(x,y):
            return 2*D*E*alpha*(1 - np.exp(-alpha*y))*np.exp(-alpha*y)/np.cosh(x)**2
        
        def ddpotential1(x,y):
            return 2*E*(D*(1 - np.exp(-alpha*y))**2 + Vc)*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2
        
        def ddpotential2(x,y):
            return -4*D*E*alpha*(1 - np.exp(-alpha*y))*np.exp(-alpha*y)*np.sinh(x)/np.cosh(x)**3
        
        def ddpotential3(x,y):
            return -4*D*E*alpha*(1 - np.exp(-alpha*y))*np.exp(-alpha*y)*np.sinh(x)/np.cosh(x)**3
        
        def ddpotential4(x,y):
            return 2*D*E*alpha**2*(-1 + 2*np.exp(-alpha*y))*np.exp(-alpha*y)/np.cosh(x)**2
    # global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4 
    # def potential(x,y):
                # return ((Vc+w*y**2)*E/np.cosh(x)**2)
                
    # def dxpotential(x,y):
                # return (-2*E*(Vc + w*y**2)*np.sinh(x)/np.cosh(x)**3)
                
    # def dypotential(x,y):
                # return (2*E*w*y/np.cosh(x)**2)

    # def ddpotential1(x,y):
                # return (2*E*(Vc + w*y**2)*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2)
                
    # def ddpotential2(x,y):
                # return (-4*E*w*y*np.sinh(x)/np.cosh(x)**3)
                
    # def ddpotential3(x,y):
                # return (-4*E*w*y*np.sinh(x)/np.cosh(x)**3)

    # def ddpotential4(x,y):
                # return (2*E*w/np.cosh(x)**2)
                
if(potkey=='diffraction'):
    global potential,dxpotential,dypotential,ddpotential1,ddpotential2,ddpotential3,ddpotential4 
    def potential(x,y):
                return (V0/(np.cosh(x)**2*np.cosh(y)**2))
                
    def dxpotential(x,y):
                return (-2*V0*np.sinh(x)/(np.cosh(x)**3*np.cosh(y)**2))
                
    def dypotential(x,y):
                return (-2*V0*np.sinh(y)/(np.cosh(x)**2*np.cosh(y)**3))

    def ddpotential1(x,y):
                return (2*V0*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/(np.cosh(x)**2*np.cosh(y)**2))
                
    def ddpotential2(x,y):
                return (4*V0*np.sinh(x)*np.sinh(y)/(np.cosh(x)**3*np.cosh(y)**3))
                
    def ddpotential3(x,y):
                return (4*V0*np.sinh(x)*np.sinh(y)/(np.cosh(x)**3*np.cosh(y)**3))
                
    def ddpotential4(x,y):
                return (2*V0*(3*np.sinh(y)**2/np.cosh(y)**2 - 1)/(np.cosh(x)**2*np.cosh(y)**2))
                
if(potkey=='gaspard'):
    
        global V0,Vc,edge,x1,y1,x2,y2,x3,y3
        V0 = 64.0
        Vc=0.0
        edge = 8.0
        x1 = 0.0
        y1 = -0.5*edge
        x2 = 0.5*edge
        y2 = (0.5*3**0.5 -0.5)*edge
        x3 = -0.5*edge
        y3 = (0.5*3**0.5 -0.5)*edge
        
        
        def potential(x,y):
                    return V0/(cosh(x - x3)**2*cosh(y - y3)**2) + V0/(cosh(x - x2)**2*cosh(y - y2)**2) + V0/(cosh(x - x1)**2*cosh(y - y1)**2)
                    
        def dxpotential(x,y):
                    return -2*V0*sinh(x - x1)/(cosh(x - x1)**3*cosh(y - y1)**2) - 2*V0*sinh(x - x2)/(cosh(x - x2)**3*cosh(y - y2)**2) - 2*V0*sinh(x - x3)/(cosh(x - x3)**3*cosh(y - y3)**2)
                    
        def dypotential(x,y):
                    return -2*V0*sinh(y - y1)/(cosh(x - x1)**2*cosh(y - y1)**3) - 2*V0*sinh(y - y2)/(cosh(x - x2)**2*cosh(y - y2)**3) - 2*V0*sinh(y - y3)/(cosh(x - x3)**2*cosh(y - y3)**3)

        def ddpotential1(x,y):
                    return 2*V0*(3*sinh(x - x1)**2/(cosh(x - x1)**4*cosh(y - y1)**2) + 3*sinh(x - x2)**2/(cosh(x - x2)**4*cosh(y - y2)**2) + 3*sinh(x - x3)**2/(cosh(x - x3)**4*cosh(y - y3)**2) - 1/(cosh(x - x3)**2*cosh(y - y3)**2) - 1/(cosh(x - x2)**2*cosh(y - y2)**2) - 1/(cosh(x - x1)**2*cosh(y - y1)**2))
                    
        def ddpotential2(x,y):
                    return 4*V0*(sinh(x - x1)*sinh(y - y1)/(cosh(x - x1)**3*cosh(y - y1)**3) + sinh(x - x2)*sinh(y - y2)/(cosh(x - x2)**3*cosh(y - y2)**3) + sinh(x - x3)*sinh(y - y3)/(cosh(x - x3)**3*cosh(y - y3)**3))

                    
        def ddpotential3(x,y):
                    return 4*V0*(sinh(x - x1)*sinh(y - y1)/(cosh(x - x1)**3*cosh(y - y1)**3) + sinh(x - x2)*sinh(y - y2)/(cosh(x - x2)**3*cosh(y - y2)**3) + sinh(x - x3)*sinh(y - y3)/(cosh(x - x3)**3*cosh(y - y3)**3))

                    
        def ddpotential4(x,y):
                    return 2*V0*(3*sinh(y - y1)**2/(cosh(x - x1)**2*cosh(y - y1)**4) + 3*sinh(y - y2)**2/(cosh(x - x2)**2*cosh(y - y2)**4) + 3*sinh(y - y3)**2/(cosh(x - x3)**2*cosh(y - y3)**4) - 1/(cosh(x - x3)**2*cosh(y - y3)**2) - 1/(cosh(x - x2)**2*cosh(y - y2)**2) - 1/(cosh(x - x1)**2*cosh(y - y1)**2))
                


                
                







                
    




                



# def potential(x,y):
    # return (V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*E/(np.cosh(x)**2) +0j#0.5*k*x**2 + 0.5*k*y**2 #V0*E/(np.cosh(x)**2*np.cosh(y)**2)#*(exp(-x**2/s**2)
# #(V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0)))*##

# def dxpotential(x,y):
    # return -2*E*(V0 - 0.5*m*w**2*y**2 + m**2*w**4*y**4/(16*V0))*np.sinh(x)/np.cosh(x)**3

# def dypotential(x,y):
    # return E*(-1.0*m*w**2*y + m**2*w**4*y**3/(4*V0))/np.cosh(x)**2

# def ddpotential1(x,y):
    # return E*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)*(16*V0 - 8.0*m*w**2*y**2 + m**2*w**4*y**4/V0)/(8*np.cosh(x)**2)

# def ddpotential2(x,y):
    # return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

# def ddpotential3(x,y):
    # return 2*E*m*w**2*y*(1 - m*w**2*y**2/(4*V0))*np.sinh(x)/np.cosh(x)**3

# def ddpotential4(x,y):
    # return -E*m*w**2*(1 - 3*m*w**2*y**2/(4*V0))/np.cosh(x)**2

# def potential(x,y):
    # return (V0)*E/(np.cosh(x)**2)

# def dxpotential(x,y):
    # return -2*E*V0*np.sinh(x)/np.cosh(x)**3

# def dypotential(x,y):
    # return 0

# def ddpotential1(x,y):
    # return 2*E*V0*(3*np.sinh(x)**2/np.cosh(x)**2 - 1)/np.cosh(x)**2

# def ddpotential2(x,y):
    # return 0

# def ddpotential3(x,y):
    # return 0

# def ddpotential4(x,y):
    # return 0
    


# divergearr

##(-6.4-2j) (-0.4+1j)
##(-6.2-1.8j) (0.8-1.8j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-6.2-1.8j) (0.8-1.8j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-6.2-1.8j) (0.8-1.8j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-4.6-0.8j) (1.8-1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-4.6-0.8j) (1.8-1.2j)
##(-4+0.2j) (-0.2+1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)
##(-4+0.2j) (-0.2+1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-6.8-1.4j) (-1+0.2j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)
##(-5.4-1.8j) (1.6-1.8j)
##(-5.4-1.4j) (0.4-0.4j)
##(-4.6-0.8j) (1.8-1.2j)



