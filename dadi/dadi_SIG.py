import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions

# Usage: python dadi_Run_Optimizations.py

# import data
snps = "/media/burke/bigMac/ethan/dadi/syma_dadi_file.tsv"
dd = dadi.Misc.make_data_dict(snps)
pop_ids=["1", "2"]
proj = [12,16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

# print info about sfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

# define model
def SIG(params, ns, pts):
    s, nu1, nu2, Tdiv = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_func = lambda t: s * (nu1/s)**(t/Tdiv)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/Tdiv)
    phi=dadi.Integration.two_pops(phi,xx,Tdiv,nu1=nu1_func,nu2=nu2_func,m12=0,m21=0)
    model_fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return model_fs

# define optimization parameters
pts = [40,50,60]
p_labels = "s, nu1, nu2, Tdiv"
upper = [1,20,20,15]
lower = [0.01,0.01,0.01,0.1]
reps = [10,20,50]
maxiters = [5,10,20]
folds = [3,2,1]

# run 5 optimizations
for i in range(1,6):
    prefix = "V5_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "SIG", SIG, 3, 4,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)