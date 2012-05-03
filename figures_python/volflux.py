from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import os as os


Q=np.nan*np.ones((4,121,2))
Qm=np.nan*np.ones((4,2))

cnt=0
ddir=os.getcwd();
for ii in np.arange(3,7,1): # run53-run56    
    ddir='../run5'+str(ii)+'/';
    sname=ddir+'saves/save0.00e00.070.019.nc';
    gname=ddir+'saves/grid.70.19.nc'


    ff=NF(sname,'r')
    u=ff.variables['u'][:]
    thick=ff.variables['h'][:]
    time=ff.variables['Time'][:]
    ff.close()
    
    ff=NF(gname,'r')
    dyh=ff.variables['dyh'][:]
    ff.close()
    
    u=u[:,:,:,:-1]+0.5*np.diff(u,axis=3) #regrid
    dyh=np.tile(dyh,(u.shape[0],u.shape[1],1,1))
    
    U=thick*u*dyh
    U=np.sum(U,axis=2)
    
    ind=33.
    
    U=U[:,:,ind]/1e6
    Q[cnt,...]=U
    Qm[cnt,:]=np.mean(U[-60:-1,:],axis=0)
    cnt+=1

for ii in range(10):
    plt.close()

labels=['1028.5','1029.6','1030.8','1031.9']
cols=['r','g','b','k']

for ii in range(Q.shape[0]):  
    plt.plot(time,Q[ii,:,0],color=cols[ii],label=r"$\rho_2=$"+labels[ii])
    plt.plot(time,-Q[ii,:,1],color=cols[ii],linestyle='--')
#    plt.axhline(Qm[ii,0],xmin=0,xmax=1,color=cols[ii])
#    plt.axhline(-Qm[ii,1],xmin=0,xmax=1,color=cols[ii],linestyle='--')
plt.legend(loc=4)
plt.grid(True,color='grey')
plt.ylim([0.6,1.7])
#plt.ylim(0.,0.8)
#ax = plt.axes()
#p = ax.imshow(np.squeeze(u[tndx,layer,...]),interpolation='nearest',origin='lower',extent=[(lonq[0]-0.5*dy),(lonq[-1]+0.5*dy),lath[0]-0.5*dz,lath[-1]+0.5*dz],aspect=fac)
#c=ax.contour(lonh,lath,maskp,[0.1])
#plt.colorbar(p).set_label(r'$u\,[ms^{-1}]$',fontsize=15)
plt.xlabel('Time [d]'); 
plt.ylabel('Layer volume transp. [Sv]');
plt.text(-5,1.7,'5)',fontsize=17)
plt.savefig('../figures/volflux.pdf')