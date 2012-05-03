from netCDF4 import Dataset as NF
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import os as os

def h1(H1,D):
    return H1*(D-0.5*H1)/(D-H1)
    
def h2(H1,D):
    return 0.5*(D-H1)

def f(x1,x2):
    return 1.-x1/x2

def g(x1,x2):
    return 2.-x1/x2    

def U2s(gint,H1s,Ds):
    h1s=h1(H1s,Ds);
    h2s=h2(H1s,Ds);    
    t1=(1./h1s)*(Ds**2/H1s**2)*f(H1s,Ds)**4/g(H1s,Ds)**2
    t2=(1./h2s)
    return -np.sqrt(gint/(t1+t2))    

def U1n(gint,H1n,Dn):
    h1n=h1(H1n,Dn);
    h2n=h2(H1n,Dn);    
    t1=(1./h2n)*(H1n**2/Dn**2)*g(H1n,Dn)**2/f(H1n,Dn)**4
    t2=(1./h1n)
    return np.sqrt(gint/(t1+t2))

def U2n(gint,H1n,Dn):
    U1n_=U1n(gint,H1n,Dn)
    bla= U1n_**2*H1n**2*g(H1n,Dn)**2/(Dn**2*f(H1n,Dn)**4)    
    return -np.sqrt(bla)
    
def U1s(gint,H1s,Ds):
    U2s_=U2s(gint,H1s,Ds)
    bla= U2s_**2*Ds**2*f(H1s,Ds)**4/(H1s**2*g(H1s,Ds)**2)    
    return np.sqrt(bla)   
    
def mom(gint,H1s,H1n,Ds,Dn): 
    U1s_=U1s(gint,H1s,Ds)
    U2s_=U2s(gint,H1s,Ds)
    U1n_=U1n(gint,H1n,Dn)
    U2n_=U2n(gint,H1n,Dn)
    left=U1s_**2-U2s_**2+2.*gint*H1s
    right=U1n_**2-U2n_**2+2.*gint*H1n
    return (left-right)
    
def two(gint,H1s,H1n,Ds,Dn,Ws,Wn):
    U2s_=U2s(gint,H1s,Ds)
    left=U1n(gint,H1n,Dn)
    right=-U2s_*Ds*Ws*f(H1s,Ds)**2/(H1n*Wn*g(H1n,Dn))
    return (left-right)

def fun(ar,gint,Ds,Dn,Ws,Wn):
    H1s=ar[0]; H1n=ar[1];
    mom_=mom(gint,H1s,H1n,Ds,Dn)
    two_=two(gint,H1s,H1n,Ds,Dn,Ws,Wn)
    return (mom_,two_)

gint=9.81*1.5*1e-3 #any value; gint cancels out

Ds=284.
Dn=880.
#Ws=22.3e3
#Wn=13.8e3
Ws=22e3 #HIM 2km resolution
Wn=14e3 #HIM 2km resolution

[H1s,H1n]=opt.fsolve(fun, [200.,100.], (gint,Ds,Dn,Ws,Wn))

# end (Bryden and Kinder)
########

if 1: #
###########################################
# analytic 

    A1s=Ws*H1s*(Ds-H1s/2.)/Ds
    A2s=0.5*Ws*Ds-A1s

    grav=9.81
    beta=0.00075
    r0=1027.6
    S0=36.5
    
    S1=36.5
    r1=r0*(1.+beta*(S1-S0))   
    
    S2=np.linspace(36.5,42.08,20.)
    E=np.nan*S2
    r2=np.nan*S2
    Q=np.nan*S2
    for ii in np.arange(S2.shape[0]):

        r2[ii]=r0*(1.+beta*(S2[ii]-S0))
        gint=grav*(r2[ii]-r1)/r2[ii]        
       
        u1s=U1s(gint,H1s,Ds)
        u1n=U1n(gint,H1n,Dn)
        u2s=U2s(gint,H1s,Ds)
        u2n=U2n(gint,H1n,Dn)
        
        Q1=u1s*A1s
        Q2=u2s*A2s
        Q[ii]=Q1
        
        E[ii]=-2*(S1*Q1+S2[ii]*Q2)/(S1+S2[ii])


###########################################
# him
S2_him=S2[4::5]
E_him=np.nan*S2_him
sd_him=np.nan*S2_him
Q_him=np.nan*S2_him


cnt=0
#for ii in np.arange(5,9,1): # run45-run48
    #ddir='/home/stefan/arbeit/him/run4'+str(ii)+'/';

for ii in np.arange(3,7,1): # run53-run56    
    ddir='../run5'+str(ii)+'/';
    sname=ddir+'saves/save0.00e00.070.019.nc';
    print sname
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
    U=U[:,:,ind]
    
    Q1=np.mean(U[-60:-1,0])
    sd=np.std(U[-60:-1,0])
    Q2=-Q1
    Q_him[cnt]=Q1
    sd_him[cnt]=sd
    E_him[cnt]=-2.*(S1*Q1+S2_him[cnt]*Q2)/(S1+S2_him[cnt])
    cnt+=1
###########################################

for ii in range(10):
    plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hold(True)
#pt=11
p1=ax1.plot(S2,E*1e-6,'k')
p1=ax1.plot(S2,E*1e-6,'k.')
p1=ax1.plot(S2_him,E_him*1e-6,'bo',markersize=10)
#plt.errorbar(S2_him,E_him*1e-6, xerr=0., yerr=sd_him*1e-6,fmt=None)
#ax1.plot(S2[pt],E[pt]*1e-6,'xk')
ax1.set_xlabel(ur'$S_2$ (\u2030)')
ax1.set_ylabel('E [Sv]')
xtick=ax1.get_xticks()
rtick=r0*(1.+beta*(xtick-S0))
ax1.grid(True,color='grey')
ax2 = ax1.twinx()
ax2.plot(S2,Q*1e-6,'r')
ax2.plot(S2,Q*1e-6,'r.')
ax2.plot(S2_him,Q_him*1e-6,'yo',markersize=10)
plt.errorbar(S2_him,Q_him*1e-6, xerr=0., yerr=sd_him*1e-6,fmt=None,capsize=20,ecolor='k')
#ax2.plot(S2[pt],Q[pt]*1e-6,'xk')
#ax2.grid(True,color='r')
ax2.set_ylabel('$Q_1 (=-Q_2)$ [Sv]', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax3=ax1.twiny()
ax3.set_xticks(xtick)
ax3.set_xlim([xtick[0],xtick[-1]])
ax3.set_xticklabels(np.round(rtick-1e3, decimals=2))
ax3.set_xlabel(r'$\sigma_{\theta}$ $\rm [kg/m^3]$')
plt.text(35.2,0.26,'4)',fontsize=17)
plt.savefig('../figures/him_bryKin.pdf')




    