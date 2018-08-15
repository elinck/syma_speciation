import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions_GOF

# load data
snps = "/media/burke/bigMac/ethan/dadi/syma_dadi_file.tsv"
dd = dadi.Misc.make_data_dict(snps)
pop_ids=["1", "2"]
proj = [12,16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

# define model
def SC(params, ns, pts):
    nu1, nu2, m12, m21, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi=dadi.Integration.two_pops(phi,xx,T1,nu1,nu2,m12=0,m21=0)
    phi=dadi.Integration.two_pops(phi,xx,T2,nu1,nu2,m12,m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

# define optimized params
pts = [40,50,60]
emp_params = [1.0864,1.7594,9.8827,4.4261,9.8656,0.978]

# fit the model using these parameters and return the folded model SFS (scaled by theta).
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, pts, "Empirical", "SC", SC, emp_params)

# simulation params
sims = 100
p_num = 6
rounds = 3
reps = [10,20,50]
maxiters = [5,10,20]
folds = [3,2,1]

# Execute the optimization routine for each of the simulated SFS.
Optimize_Functions_GOF.Perform_Sims(sims, scaled_fs, pts, "SC", SC, rounds, p_num, reps=reps, maxiters=maxiters, folds=folds)
