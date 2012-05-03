from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import numpy as np

gname='/home/stefan/arbeit/him/run81/saves/grid.213.85.nc';
ff=NF(gname,'r')
wet=ff.variables['wet'][:]
ff.close()

tndx=599

for ii in range(10):
    plt.close()
plt.figure(figsize=(18,10))
n=3
m=3
cnt=1
###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run81/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=0
ifc=0

ff=NF(sname,'r')
e=np.mean(ff.variables['e'][tndx:,ifc,...],axis=0)
ff.close()
e[wet==0]=np.nan
emin=-0.08
emax=0.2;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = emin, vmax = emax, clip = False)
g=plt.contourf(e,np.arange(emin,emax,0.01),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
plt.colorbar(g)

plt.ylabel('y [dy=2km]')
plt.text(-35,83,'8a)',fontsize=12)

###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run81/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=1.8;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.25),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
    
plt.quiver(yg,xg,utmp,vtmp,scale=15)
plt.colorbar(g)
plt.text(-35,83,'b)',fontsize=12)

###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run81/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=1
ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=0.5;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.05),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
utmp[np.abs(utmp)>0.3]=np.nan
vtmp[np.abs(vtmp)>0.3]=np.nan
plt.quiver(yg,xg,utmp,vtmp,scale=8)
plt.colorbar(g)

plt.text(-35,83,'c)',fontsize=12)
###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run83/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=0
ifc=0

ff=NF(sname,'r')
e=np.mean(ff.variables['e'][tndx:,ifc,...],axis=0)
ff.close()
e[wet==0]=np.nan
emin=-0.08
emax=0.2;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = emin, vmax = emax, clip = False)
g=plt.contourf(e,np.arange(emin,emax,0.01),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
plt.colorbar(g)

plt.ylabel('y [dy=2km]')
plt.text(-35,83,'d)',fontsize=12)

###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run83/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=1.8;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.25),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
    
plt.quiver(yg,xg,utmp,vtmp,scale=15)
plt.colorbar(g)
plt.text(-35,83,'e)',fontsize=12)

###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run83/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=1
ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=0.5;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.05),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
utmp[np.abs(utmp)>0.3]=np.nan
vtmp[np.abs(vtmp)>0.3]=np.nan
plt.quiver(yg,xg,utmp,vtmp,scale=8)
plt.colorbar(g)

plt.text(-35,83,'f)',fontsize=12)
###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run82/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=0
ifc=0

ff=NF(sname,'r')
e=np.mean(ff.variables['e'][tndx:,ifc,...],axis=0)
ff.close()
e[wet==0]=np.nan
emin=-0.07
emax=0.07;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = emin, vmax = emax, clip = False)
g=plt.contourf(e,np.arange(emin,emax,0.005),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
plt.colorbar(g)
plt.ylabel('y [dy=2km]')
plt.text(-35,83,'g)',fontsize=12)

###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run82/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=0.9;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.05),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
plt.quiver(yg,xg,utmp,vtmp,scale=15)
plt.colorbar(g)
plt.text(-35,83,'h)',fontsize=12)
###################################################
plt.subplot(m,n,cnt)
cnt+=1
ddir='/home/stefan/arbeit/him/run82/saves/';
sname=ddir+'save0.00e00.213.085.nc';
gname=ddir+'grid.213.85.nc';

lay=1
ff=NF(sname,'r')
u=np.mean(ff.variables['u'][tndx:,lay,...],axis=0)
v=np.mean(ff.variables['v'][tndx:,lay,...],axis=0)
ff.close()

u=u[:,:-1]+0.5*np.diff(u,axis=1)
v=v[:-1,:]+0.5*np.diff(v,axis=0)

xg=np.arange(e.shape[0])
yg=np.arange(e.shape[1])
isub=5
jsub=5
utmp=u[::isub,::jsub]
vtmp=v[::isub,::jsub]
xg=xg[::isub]
yg=yg[::jsub]
u[wet==0]=np.nan
v[wet==0]=np.nan
ke=np.sqrt(u**2+v**2)
kemin=0.
kemax=0.1;
palette = matplotlib.cm.jet
norm1=matplotlib.colors.Normalize(vmin = kemin, vmax = kemax, clip = False)
g=plt.contourf(ke,np.arange(kemin,kemax,0.01),\
    cmap=palette,\
    norm = norm1, extend='both'
    )
utmp[np.abs(utmp)>0.3]=np.nan
vtmp[np.abs(vtmp)>0.3]=np.nan
plt.quiver(yg,xg,utmp,vtmp,scale=2)
plt.colorbar(g)
plt.text(-35,83,'i)',fontsize=12)


#subplots_adjust(left=0.08, bottom=0.1, right=0.885, top=0.9)

plt.figtext(0.03,0.8,'EXP1',fontsize=18)
plt.figtext(0.03,0.5,'EXP2',fontsize=18)
plt.figtext(0.03,0.2,'EXP3',fontsize=18)

plt.figtext(0.16,0.95,'Free Surface [m]',fontsize=18)
plt.figtext(0.41,0.95,r'$\sqrt{u^2+v^2}$'+' Layer 1 [m/s]',fontsize=18)
plt.figtext(0.68,0.95,r'$\sqrt{u^2+v^2}$'+' Layer 2 [m/s]',fontsize=18)


plt.savefig('../figures/all.pdf')