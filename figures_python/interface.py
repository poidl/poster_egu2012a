from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np


vec=[64,58,66,61,60,62,63];

ifc=np.nan*np.ones((np.size(vec),70))
ifc0=np.nan*np.ones(np.size(vec))
Qm=np.nan*np.ones(np.size(vec))
E=np.nan*np.ones(np.size(vec))
cnt=0
for ii in vec:     
    ddir='../run'+str(ii)+'/';
    sname=ddir+'saves/save0.00e00.070.019.nc';
    gname=ddir+'data/grd.nc'
    gname2=ddir+'saves/grid.70.19.nc'

    ff=NF(sname,'r')
    e=ff.variables['e'][:]
    u=ff.variables['u'][:]
    thick=ff.variables['h'][:]
    ff.close()
    
    ff=NF(gname2,'r')
    dyh=ff.variables['dyh'][:]
    ff.close()
    
    u=u[:,:,:,:-1]+0.5*np.diff(u,axis=3) #regrid
    dyh=np.tile(dyh,(u.shape[0],u.shape[1],1,1))
    U=thick*u*dyh
    U=np.sum(U,axis=2)    
    ind=33.   
    U=U[:,:,ind]/1e6
    Qm[cnt]=np.mean(U[-60:,0],axis=0)
    Q1=Qm[cnt]
    Q2=-Q1
    
    beta=0.00075
    r0=1027.6
    r2=1030.3
    S0=36.5
    S1=36.5
    S2=(1./beta)*(r2/r0-1.)+S0
    E[cnt]=-2.*(S1*Q1+S2*Q2)/(S1+S2)
    
    ifc[cnt,:]=np.mean(e[-60:,1,9,:],axis=0)
    ifc0[cnt]=e[0,1,9,-1]
    cnt+=1

ff=NF(gname,'r')
D=ff.variables['D'][:]
ff.close()
D=-D[9,:]
sill_i=np.where(D==-284.)[0]
narrow_i=np.where(D==-880.)[0]

for ii in range(10):
    plt.close() 

plt.plot(ifc[0,:],'b--',label=r'$D_R=$'+str(ifc0[0])+' m, Q='+str(np.round(Qm[0],decimals=3))+' Sv., E='+str(np.round(E[0],decimals=3))+' Sv.')    
for ii in np.arange(1,ifc.shape[0]):
    plt.plot(ifc[ii,:],label=r'$D_R=$'+str(ifc0[ii])+' m, Q='+str(np.round(Qm[ii],decimals=3))+' Sv., E='+str(np.round(E[ii],decimals=3))+' Sv.')
plt.axvline(narrow_i,ymin=0,ymax=1,color='r')    
plt.legend(loc=4)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
llines = leg.get_lines()
plt.setp(ltext, fontsize='small')
plt.setp(llines, linewidth=2.0) 
plt.grid(True,color='grey')
plt.xlim([sill_i,60.])
plt.ylim([-300.,0.])
plt.xlabel('x [dx=2km]')
plt.ylabel('Depth [m]')

plt.annotate('Sill', [23.,-284.], xytext=[23,-325], xycoords='data',
         arrowprops=dict(facecolor='black', shrink=0.05),textcoords='data',horizontalalignment='center')
plt.annotate('Contraction', [narrow_i,-284.], xytext=[narrow_i,-325], xycoords='data',
         arrowprops=dict(facecolor='black'),textcoords='data',horizontalalignment='center')

plt.text(19.,0.,'6)',fontsize=17)

plt.savefig('../figures/interface.pdf')
    
    
    
    