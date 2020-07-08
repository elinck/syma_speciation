import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt
import pandas as pd

# load sfs
fs=moments.Spectrum.from_file("syma_unlinked_sfs.txt")
fs=fs.fold()
ns=fs.sample_sizes
projections = [15,15]
fs.mask[1,:]=True #mask singletons in all populations
fs.mask[:,1]=True
fs=fs.sample()
ns=fs.sample_sizes

def SC(params, ns):
    nToro,nMega,t1,t2,m12,m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nToro,nMega], t1, m = np.array([[0, 0], [0, 0]]))
    fs.integrate([nToro,nMega], t2, m = np.array([[0, m12], [m21, 0]]))
    return  fs

upper_bound = [10,10,10,10,10,10]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]

top_params=pd.read_csv("SC_modelparams.txt",sep="\t",header=None)
top_params=top_params.sort_values(6,ascending=False)
poptg=np.array(top_params)[0,:]
poptg=poptg[range(6)]
poptg=moments.Inference.optimize(poptg, fs, SC, lower_bound=lower_bound, upper_bound=upper_bound, verbose=True, maxiter=1)

model=SC(poptg, ns)
ll_model=moments.Inference.ll_multinom(model,fs)
theta = moments.Inference.optimal_sfs_scaling(model, fs)
L=9.89e8*(50521/78882912.0)
Nref=theta/(4*2.3e-9*L)
nToro=poptg[0]*Nref
nMega=poptg[1]*Nref
t1=poptg[2]*2*Nref
t2=poptg[3]*2*Nref
m12=poptg[4]/(2*Nref)*nToro
m21=poptg[5]/(2*Nref)*nMega
out1=[nToro,nMega,t1,t2,m12,m21,ll_model,theta]
out1="\t".join(map(str,out1))+"\n"
f=open("SC_realparams_boots.txt",'a')
f.write(out1)
f.close()
out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
f=open("SC_modelparams_boots.txt",'a')
f.write(out2)
f.close()

#can run this script in parallel with eg
#parallel -N0 -j 10 python SC_bootstraps.py ::: {1..10}
