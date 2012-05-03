from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np
import mod_data as md
import mod_triang as tr
import os as os
import inspect as inspect

rundir= os.path.dirname(inspect.getfile(inspect.currentframe()))
grdname=rundir+'/../data/grd.nc';


#set up basin
domx=750e3 #distance Cam. sill-Gibraltar ~40km
domx+=-10.723*(15+2*180)
domy=170e3
#domy+=2e3
depth=2e3

north=tr.make_coast2()
bottom=tr.make_bottom()

#res_c=15
#res_b=180.
#nx=res_c+2.*(res_b-2.);
nx=float(north.shape[0])
ny=np.round((domy/domx)*nx)
#
xh=domx*(np.arange(nx)/nx-0.5)
yh=domy*(np.arange(ny)/ny)
yh=yh-yh[np.floor(yh.shape[0]/2.)]
#

narrow_i=np.argmin(north)
north[narrow_i+1]=north[narrow_i] # dirty

bot=np.nan*np.ones((ny,nx))
for jj in np.arange(bot.shape[0]):
    for ii in np.arange(bot.shape[1]):
        dist=north[ii]-np.abs(yh[jj])
#+0.5*np.diff(yh)[0]
        if dist<0.:
            bot[jj,ii]=0.
        else:
            bot[jj,ii]=bottom[ii]*dist/north[ii]
            

sill_i=np.argmin(bottom)

plt.close()
plt.close()
plt.close()

plt.close()
plt.close()
plt.close()
plt.figure()
plt.plot(north)
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_i, ymin=0, ymax=1,color='k')
plt.figure()
plt.plot(bottom,'.')
plt.axvline(sill_i, ymin=0, ymax=1,color='k')
plt.axvline(narrow_i, ymin=0, ymax=1,color='k')
plt.figure()
#bot[bot==0.]=np.nan
plt.imshow(bot,interpolation='nearest')

print "Depth at Sill:    " +str(bottom[sill_i])
print "Depth at Narrows: " +str(bottom[narrow_i])
W=2*north
print "Width at Sill:    " +str(W[sill_i])
print "Width at Narrows: " +str(W[narrow_i])

print "Distance Sill-Narrow: " +str((domx/nx)*(narrow_i-sill_i))

print "dx (domx/nx): "+str(domx/nx)
print "dy (domy/ny): "+str(domy/ny)

#plt.contourf(xh,yh,bot)
#plt.colorbar()
 

if 0: #large domain
    xs=140.
    xh=xh[xs:]
    bot=bot[:,xs:]
if 0: #small
    xs=150.; xe=xs+90.
    xh=xh[xs:xe]
    bot=bot[:,xs:xe]
if 1: #channel
    xs=160.; xe=xs+70.
    w=10.
    ys=43.-w; ye=43.+w-1;
    xh=xh[xs:xe]
    yh=yh[ys:ye]
    bot=bot[ys:ye,xs:xe]   

bot[bot==0.]=0.1
bot[0,:]=0.1
bot[-1,:]=0.1 # dirty 

md.create_grd(grdname,yh,xh,'x')

grid=NF(grdname,'a')
grid.variables['xh'][:]=xh
grid.variables['yh'][:]=yh
grid.variables['D'][:]=bot
grid.close()
#
print ""
print "lenlat: " + str(np.diff(yh)[0]*yh.shape[0])
print "lenlon: " + str(np.diff(xh)[0]*xh.shape[0])
print ""
print "lowlat: " + str(min(yh))
print "westlon: " + str(min(xh))
#
print "size: " + str(bot.shape)



plt.figure()
bot[bot==0.1]=np.nan
#plt.imshow(bot,interpolation='nearest')
plt.contourf(xh,yh,bot)
plt.colorbar()

