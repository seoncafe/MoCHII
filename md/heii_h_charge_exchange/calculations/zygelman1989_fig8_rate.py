import numpy as np, math
# Zygelman et al. 1989 Fig. 8 nonradiative points, read from the 400-dpi scan
# (view pixels: x frame 103 = 1 eV, 817 = 100 eV; y 230 = 1e-2, 377 = 1e-3, 522 = 1e-4 a0^2)
px=[(347,572),(406,390),(440,341),(451,327),(485,297),(542,256),(603,220),(674,183)]
E=np.array([10**((x-103)/357.0) for x,y in px]); s=np.array([10**(-2-(y-230)/146.0) for x,y in px])
a0sq=2.80028520e-17
for e,v in zip(E,s): print('E %6.2f eV  sigma %.2e a0^2 = %.2e cm^2'%(e,v,v*a0sq))
lE=np.log(E); ls=np.log(s*a0sq)
slope_hi=(ls[-1]-ls[-3])/(lE[-1]-lE[-3]); slope_lo=(ls[1]-ls[0])/(lE[1]-lE[0])
def sig(e):
    le=np.log(e)
    out=np.interp(le,lE,ls)
    out=np.where(le>lE[-1], ls[-1]+slope_hi*(le-lE[-1]), out)
    out=np.where(le<lE[0], ls[0]+slope_lo*(le-lE[0]), out)   # steep drop below 4.8 eV
    return np.exp(out)
u=1.66053907e-24; eV=1.602176634e-12; kB=8.617333262e-5
mu=(4.002602-5.4858e-4)*(1.007276+5.4858e-4)/(4.002602+1.007276)*u
def rate(T, cap=None):
    kT=kB*T; e=np.exp(np.linspace(np.log(0.05),np.log(80*kT),20000)); sg=sig(e)
    if cap: sg=np.where(e>cap, sig(np.array([cap]))[0]*(e/cap)**0, sg) if cap=='flat' else sg
    return math.sqrt(8/(math.pi*mu))*(kT*eV)**-1.5*np.trapz(e*sg*np.exp(-e/kT)*e,np.log(e))*eV**2
def kf96(T):
    t=min(T,1e5)/1e4; return 7.47e-15*t**2.06*(1+9.93*math.exp(-3.89*t))
print('slopes: low %.1f high %.2f'%(slope_lo,slope_hi))
print('   T     Fig8 singlet   Fig8 x 1/4     KF96     KF96/Fig8')
for T in (6e3,1e4,2e4,5e4,1e5):
    r=rate(T); print('%6.0f  %.3e   %.3e   %.3e   %.3f'%(T,r,r/4,kf96(T),kf96(T)/r))
