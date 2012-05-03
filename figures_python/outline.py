from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import os as os


ddir=os.getcwd();
sname=ddir+'/../run54/saves/save0.00e00.070.019.nc';
gname=ddir+'/../run54/saves/grid.70.19.nc'

ff=NF(sname,'r')
e=ff.variables['e'][:]
u=ff.variables['u'][:]
time=ff.variables['Time'][:]
ff.close()
v=u


x=np.arange(e.shape[3])
x=np.r_[[x]*e.shape[1]] 

for ii in range(10):
    plt.close()
    
fig=plt.figure(figsize=(8,3))
ax=fig.add_subplot(111)
pc=ax.pcolormesh(x,e[40,:,9,:],v[40,:,9,:],vmin=-2.,vmax=2.,cmap='RdGy')
for ii in range(e.shape[1]):
    myline=plt.Line2D(x[ii,:],e[40,ii,9,:],color='k')    
    ax.add_line(myline)    
cb=plt.colorbar(pc)
cb.set_label('u [m/s]')
plt.xlabel('x [dx=2km]')
plt.ylabel('Depth [m]')
plt.xlim([0,69])
plt.ylim([-1700,100])

pos= ax.get_position().get_points()
ax.set_position([pos[0,0],pos[0,1]+0.05,pos[1,0]-pos[0,0],pos[1,1]-pos[0,1]])


verts = [
    (0., -1700.), # left, bottom
    (0., 0.), # left, top
    (8.,0.), # right, top
    (8., -1700.), # right, bottom
    (0., -1700.), # ignored
    ]
codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]
path = Path(verts, codes)
patch = patches.PathPatch(path, facecolor='blue',alpha=0.2, lw=0)
ax.add_patch(patch)

verts = [
    (61., -1700.), # left, bottom
    (61., 0.), # left, top
    (69.,0.), # right, top
    (69., -1700.), # right, bottom
    (61., -1700.), # ignored
    ]
codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]
path = Path(verts, codes)
patch = patches.PathPatch(path, facecolor='blue',alpha=0.2, lw=0)
ax.add_patch(patch)


plt.text(-11,0.26,'3)',fontsize=17)


plt.savefig('../figures/outline.pdf')




    