import h5py, numpy as np, sys
R='/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase4_20261011/runs'
def ld(x):
    f=h5py.File(f'{R}/{x}/bb1e5/bb1e5_p4_rates.h5')
    return {k:f[k+'/data'][:] for k in ['T_e','n_e','x_s_stages','Gamma_s_stages','LeafSize','x_HII','x_HeIII','x_o_stages','Gamma_HI']}
A=ld('before'); 
for y in ['work','noise']:
    B=ld(y)
    V=A['LeafSize']**3
    # sulfur stage weights: S2+ index 2 (0=S0)
    for nm,i in [('S+',1),('S2+',2),('S3+',3)]:
        wa=np.sum(V*A['n_e']*A['x_s_stages'][:,i]); wb=np.sum(V*B['n_e']*B['x_s_stages'][:,i])
        print(y,nm,'sum V ne x', f'{wb/wa-1:+.4e}')
    w=V*A['n_e']*A['x_s_stages'][:,1]
    for j,nm in [(1,'Gamma S+'),(2,'Gamma S2+')]:
        print(y,nm,'x(S+)ne-weighted', f"{np.sum(w*B['Gamma_s_stages'][:,j])/np.sum(w*A['Gamma_s_stages'][:,j])-1:+.4e}")
    w2=V*A['n_e']*A['x_s_stages'][:,2]
    print(y,'Te S2+-wt',f"{np.sum(w2*B['T_e'])/np.sum(w2*A['T_e'])-1:+.3e}",'ne S2+-wt',f"{np.sum(w2*B['n_e'])/np.sum(w2*A['n_e'])-1:+.3e}")
