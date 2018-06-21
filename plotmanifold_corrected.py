import numpy as np
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)


gammaxi = 0.5
gammayi = 0.5
gammaxf = 0.5
gammayf = 0.5

S20= ([2*gammaxi*1j,0],[0,2*gammayi*1j])
gammaf = [[gammaxf,0],[0,gammayf]]
gammai = [[gammaxi,0],[0,gammayi]]


sel=None

time=6
largeSpan=80
smallSpan=20
#spins=[10,10,10,10]
spins=[7,7,7,7]
dims=[0,2,3]
plotType=['slides','3d','cloud'][1]
coords=['q0','qt'][0]
col=lambda obj:obj.J0**-0.0*abs(obj.J0)**0*np.exp(obj.sigma)*3e0
#col=lambda obj:np.exp(obj.sigma.real)
#col=lambda obj:obj.J0*5e-1
#col=lambda obj:obj.magnitude*np.exp(1j*obj.phaseangle)
#sel=lambda obj:(obj.sigma.real<0) & (obj.xbar.real>-6)
#sel=lambda obj:(obj.magnitude.real<10)
#sel=lambda obj:(abs(obj.J0)<2e0)
#J0=lambda obj:(obj.mqq+1j*obj.mqp-1j*(obj.mpq+1j*obj.mpp))
DET=lambda mat:mat[...,0,0]*mat[...,1,1]-mat[...,0,1]*mat[...,1,0]
#col=lambda obj:DET(J0(obj))/2
#sel=lambda obj:(abs(obj.J0).max((-1,-2))<3e0)

if coords=='q0':
  axes=[lambda obj:obj.xbar.real,lambda obj:obj.xbar.imag,
        lambda obj:obj.ybar.real,lambda obj:obj.ybar.imag]
  plotCmd='u 1:2:3:4 w pm3d lc rgb variable notitle'
elif coords=='qt':
  axes=[lambda obj:obj.xfinbar.real,lambda obj:obj.xfinbar.imag,
        lambda obj:obj.yfinbar.real,lambda obj:obj.yfinbar.imag]
  plotCmd='u 1:2:3:4 w p pt 5 ps 1 lc rgb variable notitle'

mmap=True
if 0:
  keys=['xGrid','yGrid','wf',\
        'xbar','pxbar','ybar','pybar',\
        'xfinbar','pxfinbar','yfinbar','pyfinbar',\
        'Sfinbar','gammafin',\
        'mpp','mpq','mqp','mqq',\
        'Turningpoints','Squantum',\
        'xfin','yfin','pxfin','pyfin',\
        'magnitude','phaseangle','sigma','J0']
  path='/home/wkoch/midorivijay/Codes/'
  #path='/home/vijay/Codes/'
  filePattern='FINCO_double_slit_eckart_{}_{}_{}_{}_at_time_{}'
  loadType='all'
else:
  keys=['wf','mpp','mpq','mqp','mqq','Squantum','xfin','yfin','pxfin','pyfin','magnitude','phaseangle','sigma','J0']
  path='/home/wkoch/midorivijay/Codes'
  path='/home/vijay/Codes'
  filePattern='FINCO_double_slit_eckart_modified_{}_{}_{}_{}_at_time_{}'
  manifoldFilePattern='FINCO_further_modified_ketmanifold_{}_{}_{}_{}'
  manifoldKeys=['xbar','pxbar','ybar','pybar']
  loadType='subset'

class dObj():
  pass

def convert(fileName,path,keys=[]):
  import subprocess
  converter=subprocess.check_call(['python3','-c','''
import pickle,sys,os,io
import numpy as np
with open(os.path.join('{1}','{0}.pkl'),"rb") as inFile:
  try:
    data=[]
    while (1):
      data.append(pickle.load(inFile))
      #print(data[-1],file=sys.stderr)
      #print(data[-1].shape,file=sys.stderr)
      #np.save(sys.stdout.buffer,data[-1])
  except EOFError:
    pass
keys=[{2}]
if len(keys)>0:
  for key,arr in zip(keys,data):
    np.save('{0}_'+key+'.npy',arr,allow_pickle=False)
else:
  np.savez_compressed('{0}.npz',*tuple(data))
          '''.format(fileName,path,','.join("'"+k+"'" for k in keys))])

def loadNPY(fileName,keys,thisData,mmap=False):
  if mmap:
    for ind,key in enumerate(keys):
      setattr(thisData,key,np.load(fileName+'_'+key+'.npy',mmap_mode='r'))
  else:
    dataFile=open(fileName+'.npz','r')
    # with  as dataFile:
    thisDataNpz=np.load(dataFile)#,mmap_mode='r')
    #print thisDataNpz.keys(),len(keys)
    for ind,key in enumerate(keys):
      #print ind,key,thisDataNpz[ind].shape,thisDataNpz[ind][0,0]
      setattr(thisData,key,thisDataNpz['arr_{}'.format(ind)])
      # setattr(thisData,key,pickle.loads(dataString))
      # setattr(thisData,key,pickle.load(dataFile))
    thisDataNpz.close()

def openPickle(fileName,path,keys,thisData,mmap=False):
  try:
    loadNPY(fileName,keys,thisData,mmap)
  except IOError:
    convert(fileName,path,keys=[[],keys][mmap])
    loadNPY(fileName,keys,thisData,mmap)


def myreset(spinPar=spins,wxapp=None,term=None):
  import os,pickle
  data={}
  extra={}
  ukeys={}
  #print >>plot.gnuplot,'reset'
  np.set_printoptions(edgeitems=240,linewidth=200)
  
  for x1 in range(0,3):
    for x2 in range(x1+1,4):
      for time in [0,2,3,6]:#range(7):
        print time,x1,x2
        counts=[smallSpan]*4
        counts[x1]=largeSpan
        counts[x2]=largeSpan
        thisData=dObj()
        fileName=filePattern.format(*tuple(counts+[time]))
        if loadType=='all':
          openPickle(fileName,path,keys,thisData,mmap)
        elif loadType=='subset':
          openPickle(fileName,path,keys,thisData,mmap)
          if time==0:
            fileName=manifoldFilePattern.format(*tuple(counts+[time]))
            openPickle(fileName,path,manifoldKeys,thisData,mmap)
          else:
            for ind,key in enumerate(manifoldKeys):
              #print ind,key,thisDataNpz[ind].shape,thisDataNpz[ind][0,0]
              setattr(thisData,key,getattr(data['{}_{}_{}'.format(0,x1,x2)][0],key))
            
        data['{}_{}_{}'.format(time,x1,x2)]=[thisData]
        
  return data,extra,ukeys

remapDims=[1,0,2,3]
def subSlice(data,key,time,dims,indices,noSlice=False):
  data=getattr(data['{}_{}_{}'.format(time,dims[0],dims[1])][0],key)
  counts=[smallSpan]*4
  ind=list(indices)
  for d in dims:
    counts[remapDims[d]]=largeSpan
    if not noSlice:
      ind[remapDims[d]]=slice(None)
  data=data.reshape(tuple(counts)+data.shape[1:])
  return data[ind]

class dataObj(object):
  def __init__(self,data,time,dims,indices,noSlice=False):
    self.data=data
    self.time=time
    self.dims=dims
    self.indices=indices
    self.noSlice=noSlice

  def __getattr__(self,name):
    return subSlice(self.data,name,self.time,self.dims,self.indices,self.noSlice)

#MATMUL = lambda mat1,mat2:np.mat[...,0,0]*mat[...,1,1]-mat[...,0,1]*mat[...,1,0]
def mycomputerank(data):
  for entry in data:
    dxidq0=2*matmul(gammaf,(entry.mqq+matmul(entry.mqp,S20))) -1j*(entry.mpq+matmul(entry.mpp,S20))
    
    realJ0xx = [[real(dxidq0[:,0,0]),imag(dxidq0[:,0,0])],[-imag(dxidq0[:,0,0]),real(dxidq0[:,0,0])]]
    realJ0xy = [[real(dxidq0[:,0,1]),imag(dxidq0[:,0,1])],[-imag(dxidq0[:,0,1]),real(dxidq0[:,0,1])]]
    realJ0yx = [[real(dxidq0[:,1,0]),imag(dxidq0[:,1,0])],[-imag(dxidq0[:,1,0]),real(dxidq0[:,1,0])]]
    realJ0yy = [[real(dxidq0[:,1,1]),imag(dxidq0[:,1,1])],[-imag(dxidq0[:,1,1]),real(dxidq0[:,1,1])]]

    realJ0 = [[realJ0xx,realJ0xy],[realJ0yx,realJ0yy]]

    rank = np.linalg.matrix_rank(realJ0)
    entry.rank=rank
    
def myplot(data,extra,ukeys,spinPar=spins,\
           wxapp=None,allstores=None,SetRange=None,gplId=None):
  global gnuplot
  if 0:
    import nogui
    plot.gnuplot,plot.gplProc=plot.initGnuplot()#plot.gplProc)
    nogui.gnuplot=plot.gnuplot
  #gnuplot=plot.gnuplot
  plot.emptyStore()
  #plot.storePlot('x**2')
  print >>plot.gnuplot,'''
#  reset
  set autoscale fix;set xyplane 0
set pm3d depthorder
unset colorbox
# set xrange [-7:-3]
# set yrange [-2:2]
# set zrange [-2:2]
#set hidden3d
'''
  print >>gnuplot,'set style fill solid border'
  dims.sort()
  print >>plot.gnuplot,'set xlabel "{}";set ylabel "{}";set zlabel "{}"'\
    .format(*(['x_r','x_i','y_r','y_i'][d] for d in dims))
  if plotType=='slides':
    for x1 in range(4):
      for x2 in range(x1+1,4):
        thisDat=dataObj(data,time,(x1,x2),spinPar)
        pDat=col(thisDat)
        x=np.array(np.meshgrid(np.arange(pDat.shape[0]),np.arange(pDat.shape[1]))).T
        if x1==0:
          x[...,0]+=(x2-1)*pDat.shape[0]
        else:
          x[...,1]+=pDat.shape[1]
          if x1==2:
            x[...,0]+=(x2-1)*pDat.shape[0]
          else:
            x[...,0]+=(x2-2)*pDat.shape[0]
        plot.d2plot(pDat,x,pmax=0.1,rgb=1,noImage=0)
        pDat=subSlice(data,'magnitude',time,(x1,x2),spinPar)\
              *np.exp(1j*subSlice(data,'phaseangle',time,(x1,x2),s))
  elif plotType=='3d':
    for x1,x2 in [[0,1],[0,2],[1,2]]:
      # for spinPar[1] in range (5,8):
        thisDat=dataObj(data,time,(dims[x1],dims[x2]),spinPar)
        pos=tuple(axes[d](thisDat) for d in dims)
        color=col(thisDat)

        colorRGB=plot.getRGB(color)*np.maximum(
          np.minimum(1.,np.abs(color[...,np.newaxis])),0)*255
        if sel is not None:
          selector=sel(thisDat)
          colorRGB[~selector]=255-(255-colorRGB[~selector])*0.25
        color=plot.reduceRGB(colorRGB)

        allDat=np.concatenate(tuple(k[...,np.newaxis] for k in pos+(color,)),color.ndim)
        plot.storePlot(plotCmd,allDat)
  elif plotType=='cloud':
    allDat=[]
    allFlex=list(spinPar)
    for d in dims:
      allFlex[remapDims[d]]=slice(None)
    #allFlex=[slice(None)]*4
    for x1 in range(4):
      for x2 in range(x1+1,4):
        thisDat=dataObj(data,time,(x1,x2),[spinPar,allFlex][1],True)
        #print thisDat.mpp.shape
    # for x1,x2 in [[0,1],[0,2],[1,2]]:
    #     thisDat=dataObj(data,time,(dims[x1],dims[x2]),[spinPar,allFlex][1])
        pos=tuple(axes[d](thisDat) for d in dims)
        color=col(thisDat)
        colorRGB=plot.getRGB(color)*np.maximum(
          np.minimum(1.,np.abs(color[...,np.newaxis])),0)*255
        color=plot.reduceRGB(colorRGB)
        if sel is not None:
          selector=sel(thisDat)
          pos=tuple(a[selector] for a in pos)
          color=color[selector]

        allDat.append(np.concatenate(tuple(k[...,np.newaxis] 
                                           for k in pos+(color,)),1))
    allDat=np.r_[tuple(allDat)]
    plot.storePlot('u 1:2:3:4 w p pt 5 ps 1 lc rgb variable notitle',allDat)
        

  if 0:
    z=np.linspace(-2,2,40)[:,np.newaxis]+1j*np.linspace(-3,3,40)[np.newaxis,:]
    V0=16
    omega=4
    mass=1.
    pot=np.cosh(z.real)**-2*(V0-mass*omega**2/2*z.imag**2+(mass*omega**2)**2*z.imag**4/16/V0)
    plot.storePlot('u 1:2:3 w l',np.dstack((z.real,z.imag,pot-5)))
  plot.runPlots(disp=True,d2=plotType!='slides',viafile=False)
  print 'done plotting'
  return plot.storedPlots

def standAlone():
  global plot,data,store,gnuplot
  import plot
  gnuplot,_=plot.initGnuplot()
  plot.gnuplot=gnuplot
  execfile('plotmanifold.py');data=myreset()
  execfile('plotmanifold.py');store=myplot(*data)
