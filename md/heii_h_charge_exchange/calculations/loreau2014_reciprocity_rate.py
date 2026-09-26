import numpy as np, math
# Loreau, Ryabchenko & Vaeck (2014) Table 2, process (1): H+ + He(1s2) -> H(1s) + He+(1s)
a=[-5.167,-5.200,-0.3265,-0.3786,-0.04399,-0.05318,-0.06957,0.01072]
b=[None,0.1131,1.279,0.4867,0.1944,0.1305,0.03501,-0.04530]
w=0.7800
def sig_rev(x):   # x in eV/u ; returns cm^2 ; log = log10? check both
    L=np.log10(x)
    s=a[0]+sum(a[n]*np.cos(n*w*L)+b[n]*np.sin(n*w*L) for n in range(1,8))
    return s
mp=1.007276; mHe=4.002602; me=5.48580e-4
mHep=mHe-me; mH=mp+me
mu_r=mp*mHe/(mp+mHe); mu_f=mH*mHep/(mH+mHep)
DE=24.587387-13.598434
# lab energy per u of the proton projectile -> CM: E_cm = x*mp*mHe/(mp+mHe)
def x_of_Ecm(Ecm): return Ecm*(mp+mHe)/(mp*mHe)
print('threshold eV/u', x_of_Ecm(DE), ' data start 15 eV/u -> E_fwd', 15/x_of_Ecm(1.0)-DE)
# check log base: print sigma at some energies with log10 and ln interpretation
for x in (15,30,100,1000,2e4):
    L=sig_rev(x)
    print(x, 'log10-> %.3e'%(10**L), ' ln-> %.3e'%(math.exp(L)))

u=1.66053907e-24; kB_eV=8.617333262e-5; eV=1.602176634e-12
def srev_cm2(x): return 10**sig_rev(x)*1e-16
def k_fwd(T, low='fit'):
    kT=kB_eV*T
    Ef=np.exp(np.linspace(np.log(1e-4), np.log(60*kT), 40000))   # forward CM energy [eV]
    Er=Ef+DE; x=x_of_Ecm(Er)
    s=srev_cm2(x)
    if low=='zero15': s=np.where(x<15,0.0,s)
    if low=='flat15': s=np.where(x<15,srev_cm2(15.0),s)
    sf=0.25*(mu_r*Er)/(mu_f*Ef)*s
    integ=np.trapz(Ef*sf*np.exp(-Ef/kT)*Ef, np.log(Ef))       # int E sigma e^{-E/kT} dE
    return math.sqrt(8/(math.pi*mu_f*u))*(kT*eV)**-1.5*integ*eV**2
def kf96(T):
    t=min(T,1e5)/1e4; return 7.47e-15*t**2.06*(1+9.93*math.exp(-3.89*t))
print('\n   T      k_fwd(Loreau fit) [zero<15, flat<15]      KF96       KF96/4    radiative   ratio KF96/Loreau')
for T in (3e3,6e3,1e4,2e4,5e4,1e5):
    kl=k_fwd(T); kz=k_fwd(T,'zero15'); kf=k_fwd(T,'flat15')
    print('%7.0f  %.3e [%.3e, %.3e]  %.3e  %.3e  %.3e   %.0f'%(T,kl,kz,kf,kf96(T),kf96(T)/4,1.25e-15*(T/300)**0.25,kf96(T)/kl))
