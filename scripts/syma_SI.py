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

def SI(params, ns):
    s,nToro,nMega,tSplit = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    nu1_func = lambda t: s*(nToro/s)**(t/tSplit)
    nu2_func = lambda t: (1-s)*(nMega/(1-s))**(t/tSplit)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    fs.integrate(nu_func, tSplit, dt_fac=0.01)
    return  fs

upper_bound = [1,10,10,10]
lower_bound = [1e-4,1e-4,1e-4,1e-4]

for i in range(20):
        poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(4)]
        poptg=moments.Inference.optimize_log(poptg, fs, SI,lower_bound=lower_bound,
                                             upper_bound=upper_bound,verbose=True,
                                             maxiter=3)
        model=SI(poptg, ns)
        ll_model=moments.Inference.ll_multinom(model,fs)
        theta = moments.Inference.optimal_sfs_scaling(model, fs)
        L=9.89e8*(50521/78882912.0)
        Nref=theta/(4*2.3e-9*L)
        nToro=poptg[1]*Nref
        nMega=poptg[2]*Nref
        tSplit=poptg[3]*2*Nref
        out1=[nToro,nMega,tSplit,ll_model,theta]
        out1="\t".join(map(str,out1))+"\n"
        f=open("SIg_realparams.txt",'a')
        f.write(out1)
        f.close()
        out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
        f=open("SIg_modelparams.txt",'a')
        f.write(out2)
        f.close()
