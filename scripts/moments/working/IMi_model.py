# load libraries
import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys

# import data
infile = "/media/burke/bigMac/ethan/moments_inference/syma.moments.revised.txt"
pop_ids=["mega", "toro"]
projections = [16,12]
dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)

# print info about sfs
print "projection", projections
print "sample sizes", data.sample_sizes
sfs_sum = np.around(data.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

# define model
params=[1,1,1,1,1,1,1,1,1,1,1]

def IMi(params, ns):
    nu1_0,nu2_0,nu1,nu2,T0,T1,m12,m21,m12i,m21i,P = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T1)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1,nu2], T0, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1,nu2], T0, dt_fac=0.01, m=np.array([[0, m12i], [m21i, 0]]))
    fsi.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12i], [m21i, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs

func=IMi
upper_bound = [200,200,200,200,200,200,200,200,200,200,1]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)

# run initial optimization
poptg = moments.Inference.optimize_log(params, data, func,lower_bound=lower_bound,
                                upper_bound=upper_bound,verbose=False, maxiter=10)

# run 25 optimizations with less pertubation
for i in range(25):
        poptg=moments.Misc.perturb_params(poptg, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
        poptg=moments.Inference.optimize_log(poptg, data, func,lower_bound=lower_bound,
                                        upper_bound=upper_bound,verbose=False, maxiter=10)
        sample_sizes = ns
        model=func(poptg, sample_sizes)
        ll_model=moments.Inference.ll_multinom(model, data)
        print('Maximum log composite likelihood: {0}'.format(ll_model))
        theta = moments.Inference.optimal_sfs_scaling(model, data)
        print('Optimal value of theta: {0}'.format(theta))
        L=9.89e8*(19067.0/78882912.0)
        Nref=theta/(4*2.3e-9*L)
        nu1_0=poptg[0]*Nref
        nu2_0=poptg[1]*Nref
        nu1=poptg[2]*Nref
        nu2=poptg[3]*Nref
        T0=poptg[4]*2*Nref
        T1=poptg[5]*2*Nref
        m12=(poptg[6]/(2*Nref)*nu1)
        m21=(poptg[7]/(2*Nref)*nu2)
        m12i=(poptg[8]/(2*Nref)*nu1)
        m21i=(poptg[9]/(2*Nref)*nu2)
        fst=data.Fst()
        realparams=[Nref,nu1_0,nu2_0,nu1,nu2,T0,T1,m12,m21,m12i,m21i,fst]
        roundedparams=['%.3f' % e for e in realparams]
        print(['Nref','nu1_0','nu2_0','nu1','nu2','T0','T1','m12','m21','m12i','m21i','Fst'])
        print(roundedparams)
        IMimod=open('IMi_modelparams.txt','a')
        IMimod.write(str(poptg[0])+'\t'+
            str(poptg[1])+'\t'+
            str(poptg[2])+'\t'+
            str(poptg[3])+'\t'+
            str(poptg[4])+'\t'+
            str(poptg[5])+'\t'+
            str(poptg[6])+'\t'+
            str(poptg[7])+'\t'+
            str(poptg[8])+'\t'+
            str(poptg[9])+'\t'+
            str(theta)+'\t'+
            str(fst)+'\t'+
            str(ll_model)+'\t'+
            str(L)+'\n')
        IMimod.close()
        IMireal=open('IMi_realparams.txt','a')
        IMireal.write(str(nu1_0)+'\t'+
            str(nu2_0)+'\t'+
            str(nu1)+'\t'+
            str(nu2)+'\t'+
            str(T0)+'\t'+
            str(T1)+'\t'+
            str(m12)+'\t'+
            str(m21)+'\t'+
            str(m12i)+'\t'+
            str(m21i)+'\t'+
            str(theta)+'\t'+
            str(fst)+'\t'+
            str(ll_model)+'\t'+
            str(L)+'\n')
        IMireal.close()
