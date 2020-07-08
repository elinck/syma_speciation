import allel
import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt

# build dd
fs=moments.Spectrum.from_file("syma_unlinked_sfs.txt")
fs=fs.fold()
ns=fs.sample_sizes
projections = [12,16]
fs.mask[1,:]=True #mask singletons in all populations
fs.mask[:,1]=True

def IM(params, ns):
    nToro,nMega,tSplit,m12,m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nToro,nMega], tSplit, m = np.array([[0, m12], [m21, 0]]))
    return fs

upper_bound = [10,10,10,10,10]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4]

for i in range(20):
        poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(5)]
        poptg=moments.Inference.optimize_log(poptg, fs, IM,lower_bound=lower_bound,
                                             upper_bound=upper_bound,verbose=True,
                                             maxiter=3)
        model=IM(poptg, ns)
        ll_model=moments.Inference.ll_multinom(model,fs)
        theta = moments.Inference.optimal_sfs_scaling(model, fs)
        L=9.89e8*(50521/78882912.0)
        Nref=theta/(4*2.3e-9*L)
        nToro=poptg[0]*Nref
        nMega=poptg[1]*Nref
        tSplit=poptg[2]*2*Nref
        m12=(poptg[3]/(2*Nref))*nToro
        m21=(poptg[4]/(2*Nref))*nMega
        out1=[nToro,nMega,tSplit,m12,m21,ll_model,theta]
        out1="\t".join(map(str,out1))+"\n"
        f=open("IM_realparams.txt",'a')
        f.write(out1)
        f.close()
        out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
        f=open("IM_modelparams.txt",'a')
        f.write(out2)
        f.close()