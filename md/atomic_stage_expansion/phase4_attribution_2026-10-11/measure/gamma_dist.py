import sys,h5py,numpy as np
W='/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/'+sys.argv[1]
X=sys.argv[2]
def L(x,k,f='rates'):
    return h5py.File(f'{W}/{x}/bb1e5_{f}.h5','r')[k+'/data'][:]
V=L('d','LeafSize')**3
e0=L(X,'emis_s_3','emis')[:,5]; e1=L('d','emis_s_3','emis')[:,5]
w=V*0.5*(e0+e1); m=w>0; w=w[m]/w[m].sum()
G0=L(X,'Gamma_s_stages')[m]; G1=L('d','Gamma_s_stages')[m]
J0=L(X,'J_nu')[m]; J1=L('d','J_nu')[m]; E=L('d','E_bin'); dE=L('d','dE_bin')
for j in [1,2]:
    r=np.log(G1[:,j]/G0[:,j]); o=np.argsort(r); c=np.cumsum(w[o])
    print(f'dln Gamma stage {j}: wmean {np.sum(w*r):+.3e}  wpercentiles 5/25/50/75/95:', ' '.join(f'{r[o][np.searchsorted(c,p)]:+.3e}' for p in [.05,.25,.5,.75,.95]))
# crude Gamma from J: sigma ~ E^-2.5 above 23.34 (shape only) 
for thr in [23.34,34.86]:
    s=(E>thr)
    sig=np.where(s,(E/thr)**-2.5,0)
    g0=(J0*sig*dE/E).sum(1); g1=(J1*sig*dE/E).sum(1)
    r=np.log(g1/g0); o=np.argsort(r); c=np.cumsum(w[o])
    print(f'crude J-based rate above {thr}: wmean {np.sum(w*r):+.3e}, median {r[o][np.searchsorted(c,.5)]:+.3e}')
