import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions

# Usage: python dadi_Run_Optimizations.py

# import data
snps = "/media/burke/bigMac/ethan/dadi/dadi_input_wgs.tsv"
dd = dadi.Misc.make_data_dict(snps)
pop_ids=["1", "2"]
proj = [18,18]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

# print info about sfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

# define model
def IIM(params, ns, pts):
    nu1, nu2, m12, m21, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi=dadi.Integration.two_pops(phi,xx,T1,nu1,nu2,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,T2,nu1,nu2,m12=0,m21=0)
    model_fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return model_fs

# define optimization parameters
pts = [40,50,60]
p_labels = "nu1, nu2, m12, m21, T1, T2"
upper = [20,20,10,10,10,15]
lower = [0.01,0.01,0.01,0.01,0.1,0.5]
reps = [10,20,50]
maxiters = [5,10,20]
folds = [3,2,1]

# run 5 optimizations
for i in range(1,6):
    prefix = "V5_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "IIM", IIM, 3, 6,  param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)