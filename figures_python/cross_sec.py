from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import os as os

ddir=os.getcwd();
gname1=ddir+'/../run53/data/grd.nc'
print gname1
gname=ddir+'/../run53/saves/grid.70.19.nc'

ff=NF(gname1,'r')
xh=ff.variables['xh'][:]
yh=ff.variables['yh'][:]
D=ff.variables['D'][:]
ff.close()


ff=NF(gname,'r')
wet=ff.variables['wet'][:]
dxh=ff.variables['dxh'][:]
dyh=ff.variables['dyh'][:]
ff.close()

a=np.sum( (D*wet)*dyh, axis=0)

Dmax=np.max(D,axis=0)

W=np.sum(wet,axis=0)*dyh[0]

sill_i=np.argmin(Dmax)
narrow_i=np.argmin(W)
narrow_is=np.where(W==W[narrow_i])[0] # multiple minima

Dn=880.
narrow_target_depth_i=np.where(Dmax==Dn)


for ii in range(10):
    plt.close()
    
fig=plt.figure()
ax1=plt.subplot(311)
topo=wet*D
topo[wet==0]=np.nan
plt.contourf(np.arange(xh.shape[0]),yh/1e3,topo,20)
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_target_depth_i, ymin=0, ymax=1,color='r')

plt.xlim([0,70])
plt.ylabel('y [km]')
plt.grid(True, color='gray')
xticklabels = ax1.get_xticklabels()
plt.setp(xticklabels, visible=False)
plt.colorbar(orientation='horizontal',shrink=0.7,aspect=100)
pos= ax1.get_position().get_points()
ax1.set_position([pos[0,0],pos[0,1],pos[1,0]-pos[0,0],pos[1,1]-pos[0,1]+0.07])
plt.figtext(0.46,0.63,'Depth [m]')
plt.text(-10,17,'1a)',fontsize=15)


ax2=plt.subplot(312)
plt.plot(a*1e-7)
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_target_depth_i, ymin=0, ymax=1,color='r')

plt.ylim([0.,2])
plt.xlim([0,70])
plt.ylabel('Area of \n Cross-sec.\n [1e7 m^2]', multialignment='center')
xticklabels = ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)
plt.grid(True, color='gray')
plt.text(-10,2.,'b)',fontsize=15)

plt.subplot(313)
plt.plot(-Dmax,'.')
plt.plot(-Dmax)
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_target_depth_i, ymin=0, ymax=1,color='r')
plt.ylabel('Depth [m]')
plt.xlabel('x [dx=2km]')
plt.grid(True, color='gray')
plt.text(-10,-200,'c)',fontsize=15)

print "min(Dmax): "+str(np.min(Dmax))
print "W[sill_i]: "+str(W[sill_i])
print "W[narrow_target_depth_i]: "+str(W[narrow_target_depth_i])
print ""+str()
print ""+str()
print ""+str()
plt.savefig('../figures/channel_topo.pdf')