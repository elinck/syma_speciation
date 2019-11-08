import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt
import pandas as pd
os.chdir("/media/burke/bigMac/ethan/moments_revision/")

# set print options
np.set_printoptions(precision=3)

# build dd
vcf_filename = "syma.6x.moments.unlinked.recode.vcf"
popinfo_filename = "moments.pops.txt"
pop_ids=["mega", "toro"]
projections = [16,12]
dd = Misc.make_data_dict_vcf(vcf_filename, popinfo_filename)
fs = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
fs.mask[1,:]=True #mask singletons in all populations
fs.mask[:,1]=True
ns=fs.sample_sizes

def IM(params, ns):
    nMeg0,nTor0,nMeg1,nTor1,T0,T1,mMT0,mMT1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nMeg0,nTor1], T0, m = np.array([[0, mMT0], [mMT0, 0]]))
    fs.integrate([nMeg0,nTor1], T1, m = np.array([[0, mMT1], [mMT1, 0]]))
    return  fs

upper_bound = [10,10,10,10,2,2,2,2]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]

for i in range(10):
        poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(8)]
        poptg=moments.Inference.optimize(poptg, fs, IM,lower_bound=lower_bound,
                                             upper_bound=upper_bound,verbose=True,
                                             maxiter=1)
        model=IM(poptg, ns)
        ll_model=moments.Inference.ll_multinom(model,fs)
        theta = moments.Inference.optimal_sfs_scaling(model, fs)
        L=9.89e8*(4492/78882912.0)
        Nref=theta/(4*2.3e-9*L)
        nMeg0=poptg[0]*Nref
        nTor0=poptg[1]*Nref
        nMeg1=poptg[2]*Nref
        nTor1=poptg[3]*Nref
        T0=poptg[4]*2*Nref
        T1=poptg[5]*2*Nref
        mMT0=(poptg[6]/(2*Nref))
        mMT1=(poptg[7]/(2*Nref))
        out1=[nMeg0,nTor0,nMeg1,nTor1,T0,T1,mMT0,mMT1,ll_model,theta]
        out1="\t".join(map(str,out1))+"\n"
        f=open("/media/burke/bigMac/ethan/moments_revision/IM_realparams.txt",'a')
        f.write(out1)
        f.close()
        out2="\t".join(map(str,poptg))+"\t"+str(ll_model)+"\t"+str(theta)+"\n"
        f=open("/media/burke/bigMac/ethan/moments_revision/IM_modelparams.txt",'a')
        f.write(out2)
        f.close()
