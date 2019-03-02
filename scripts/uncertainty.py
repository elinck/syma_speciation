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

# build dd
infile = "/media/burke/bigMac/ethan/moments_inference/syma.moments.revised.txt"
pop_ids=["mega", "toro"]
projections = [16,12]
dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)

# load w/ no theta to calculate uncertainty
with open('IM_modelparams.txt','r') as f:
        IMmod=[x.strip().split('\t') for x in f]
        IMmod=np.array(IMmod).astype(np.float)

IMbest = np.where(IMmod == np.max(IMmod[:,10]))
IMbestpm = IMmod[IMbest[0]]
IMbestpm = IMbestpm.flatten()
IMbestpm = IMbestpm[0:8]

# set up model
def IM(params, ns):
    nu1_0,nu2_0,nu1,nu2,T0,T1,m12,m21 = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T1)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_0,nu2_0], T0, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return  fs

funcIM=IM
all_boot_IM=Misc.bootstrap(dd,pop_ids,projections,polarized=False)
uncertIM=moments.Godambe.GIM_uncert(funcIM,all_boot_IM,IMbestpm,data)

# reload with theta
with open('IM_modelparams.txt','r') as f:
        IMmod=[x.strip().split('\t') for x in f]
        IMmod=np.array(IMmod).astype(np.float)

IMbest = np.where(IMmod == np.max(IMmod[:,10]))
IMbestpm = IMmod[IMbest[0]]
IMbestpm = IMbestpm.flatten()
IMbestpm = IMbestpm[0:9]

# calc real values
L=7e8*(19067.0/78882912.0)
Nref=9.132e+02/(4*2.3e-9*L)
IMnu1_0=IMbestpm[0]*Nref
IMnu2_0=IMbestpm[1]*Nref
IMnu1=IMbestpm[2]*Nref
IMnu2=IMbestpm[3]*Nref
IMT0=IMbestpm[4]*2*Nref
IMT1=IMbestpm[4]*2*Nref
IMm12=(IMbestpm[5]/(2*Nref)*IMnu1)
IMm21=(IMbestpm[6]/(2*Nref)*IMnu2)
IM_uncert_nu1_0=uncertIM[0]*Nref
IM_uncert_nu2_0=uncertIM[1]*Nref
IM_uncert_nu1=uncertIM[2]*Nref
IM_uncert_nu2=uncertIM[3]*Nref
IM_uncert_T0=uncertIM[4]*2*Nref
IM_uncert_T1=uncertIM[4]*2*Nref
IM_uncert_m12=(uncertIM[5]/(2*Nref)*IMnu1)
IM_uncert_m21=(uncertIM[6]/(2*Nref)*IMnu2)

# write IM stuff
IMbest=open('IM_real_bestparams.txt','a')
IMbest.write(str(IMnu1_0)+'\t'+str(IM_uncert_nu1_0)+'\n'+
    str(IMnu2_0)+'\t'+str(IM_uncert_nu2_0)+'\n'+
    str(IMnu1)+'\t'+str(IM_uncert_nu1)+'\n'+
    str(IMnu2)+'\t'+str(IM_uncert_nu2)+'\n'+
    str(IMT0)+'\t'+str(IM_uncert_T0)+'\n'+
    str(IMT1)+'\t'+str(IM_uncert_T1)+'\n'+
    str(IMm12)+'\t'+str(IM_uncert_m12)+'\n'+
    str(IMm21)+'\t'+str(IM_uncert_m21)+'\n')
IMbest.close()

# open SC
with open('SC_modelparams.txt','r') as f:
        SCmod=[x.strip().split('\t') for x in f]
        SCmod=np.array(SCmod).astype(np.float)

SCbest = np.where(SCmod == np.max(SCmod[:,10]))
SCbestpm = SCmod[SCbest[0]]
SCbestpm = SCbestpm.flatten()
SCbestpm = SCbestpm[0:8]

# set up next model
def SC(params, ns):
    nu1_0,nu2_0,nu1,nu2,T0,T1,m12,m21 = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T1)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_0,nu2_0], T0, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return  fs

funcSC=SC
all_boot_SC=Misc.bootstrap(dd,pop_ids,projections,polarized=False)
uncertSC=moments.Godambe.GIM_uncert(funcSC,all_boot_SC,SCbestpm,data)

# reload w/ theta
with open('SC_modelparams.txt','r') as f:
        SCmod=[x.strip().split('\t') for x in f]
        SCmod=np.array(SCmod).astype(np.float)

SCbest = np.where(SCmod == np.max(SCmod[:,10]))
SCbestpm = SCmod[SCbest[0]]
SCbestpm = SCbestpm.flatten()
SCbestpm = SCbestpm[0:9]

# calc real values
L=7e8*(19067.0/78882912.0)
Nref=9.132e+02/(4*2.3e-9*L)
SCnu1_0=SCbestpm[0]*Nref
SCnu2_0=SCbestpm[1]*Nref
SCnu1=SCbestpm[2]*Nref
SCnu2=SCbestpm[3]*Nref
SCT0=SCbestpm[4]*2*Nref
SCT1=SCbestpm[5]*2*Nref
SCm12=(SCbestpm[6]/(2*Nref)*SCnu1)
SCm21=(SCbestpm[7]/(2*Nref)*SCnu2)
SC_uncert_nu1_0=uncertSC[0]*Nref
SC_uncert_nu2_0=uncertSC[1]*Nref
SC_uncert_nu1=uncertSC[2]*Nref
SC_uncert_nu2=uncertSC[3]*Nref
SC_uncert_T0=uncertSC[4]*2*Nref
SC_uncert_T1=uncertSC[5]*2*Nref
SC_uncert_m12=(uncertSC[6]/(2*Nref)*SCnu1)
SC_uncert_m21=(uncertSC[7]/(2*Nref)*SCnu2)

# write SC stuff
SCbest=open('SC_real_bestparams.txt','a')
SCbest.write(str(SCnu1_0)+'\t'+str(SC_uncert_nu1_0)+'\n'+
    str(SCnu2_0)+'\t'+str(SC_uncert_nu2_0)+'\n'+
    str(SCnu1)+'\t'+str(SC_uncert_nu1)+'\n'+
    str(SCnu2)+'\t'+str(SC_uncert_nu2)+'\n'+
    str(SCT0)+'\t'+str(SC_uncert_T0)+'\n'+
    str(SCT1)+'\t'+str(SC_uncert_T1)+'\n'+
    str(SCm12)+'\t'+str(SC_uncert_m12)+'\n'+
    str(SCm21)+'\t'+str(SC_uncert_m21)+'\n')
SCbest.close()

# open SC w/ no theta for LRT
with open('SC_modelparams.txt','r') as f:
        SCmod=[x.strip().split('\t') for x in f]
        SCmod=np.array(SCmod).astype(np.float)

SCbest = np.where(SCmod == np.max(SCmod[:,10]))
SCbestpm = SCmod[SCbest[0]]
SCbestpm = SCbestpm.flatten()
SCbestpm = SCbestpm[0:8]

# likelihood ratio test
adj = moments.Godambe.LRT_adjust(funcIM, all_boot_IM, SCbestpm, data, nested_indices=[4], multinom=False)

#adj = -0.00014546845535846927
D=adj*2*(-538.724045378618--566.0220678660236) #0.02381559338597699
weights = (0, 1)
p = moments.Godambe.sum_chi2_ppf(D, weights) #p = 0.8773550502494083; fail to reject

# write remainder of bestparam tables
with open('SI_modelparams.txt','r') as f:
        SImod=[x.strip().split('\t') for x in f]
        SImod=np.array(SImod).astype(np.float)

SIbest = np.where(SImod == np.max(SImod[:,8]))
SIbestpm = SImod[SIbest[0]]
SIbestpm = SIbestpm.flatten()
SIbestpm = SIbestpm[0:6]

def SI(params, ns):
    nu1_0,nu2_0,nu1,nu2,T0,T1 = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T1)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_0,nu2_0], T0, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return  fs

funcSI=SI
all_boot_SI=Misc.bootstrap(dd,pop_ids,projections,polarized=False)
uncertSI=moments.Godambe.GIM_uncert(funcSI,all_boot_SI,SIbestpm,data)

# reload w/ theta
with open('SI_modelparams.txt','r') as f:
        SImod=[x.strip().split('\t') for x in f]
        SImod=np.array(SImod).astype(np.float)

SIbest = np.where(SImod == np.max(SImod[:,8]))
SIbestpm = SImod[SIbest[0]]
SIbestpm = SIbestpm.flatten()
SIbestpm = SIbestpm[0:6]

# calc real values
L=7e8*(19067.0/78882912.0)
Nref=9.132e+02/(4*2.3e-9*L)
SInu1_0=SIbestpm[0]*Nref
SInu2_0=SIbestpm[1]*Nref
SInu1=SIbestpm[2]*Nref
SInu2=SIbestpm[3]*Nref
SIT0=SIbestpm[4]*2*Nref
SIT1=SIbestpm[5]*2*Nref
SI_uncert_nu1_0=uncertSI[0]*Nref
SI_uncert_nu2_0=uncertSI[1]*Nref
SI_uncert_nu1=uncertSI[2]*Nref
SI_uncert_nu2=uncertSI[3]*Nref
SI_uncert_T0=uncertSI[4]*2*Nref
SI_uncert_T1=uncertSI[5]*2*Nref

# write SC stuff
SIbest=open('SI_real_bestparams.txt','a')
SIbest.write(str(SInu1_0)+'\t'+str(SI_uncert_nu1_0)+'\n'+
    str(SInu2_0)+'\t'+str(SI_uncert_nu2_0)+'\n'+
    str(SInu1)+'\t'+str(SI_uncert_nu1)+'\n'+
    str(SInu2)+'\t'+str(SI_uncert_nu2)+'\n'+
    str(SIT0)+'\t'+str(SI_uncert_T0)+'\n'+
    str(SIT1)+'\t'+str(SI_uncert_T1)+'\n')
SIbest.close()

# load w/ no theta to calculate uncertainty
with open('IIM_modelparams.txt','r') as f:
        IIMmod=[x.strip().split('\t') for x in f]
        IIMmod=np.array(IIMmod).astype(np.float)

IIMbest = np.where(IIMmod == np.max(IIMmod[:,10]))
IIMbestpm = IIMmod[IIMbest[0]]
IIMbestpm = IIMbestpm.flatten()
IIMbestpm = IIMbestpm[0:8]

# set up model
def IIM(params, ns):
    nu1_0,nu2_0,nu1,nu2,T0,T1,m12,m21 = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T1)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_0,nu2_0], T0, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T1, dt_fac=0.01, m=np.array([[0, 0], [0, 0]]))
    return  fs

funcIIM=IIM
all_boot_IIM=Misc.bootstrap(dd,pop_ids,projections,polarized=False)
uncertIIM=moments.Godambe.GIM_uncert(funcIIM,all_boot_IIM,IIMbestpm,data)

# reload with theta
with open('IIM_modelparams.txt','r') as f:
        IIMmod=[x.strip().split('\t') for x in f]
        IIMmod=np.array(IIMmod).astype(np.float)

IIMbest = np.where(IIMmod == np.max(IIMmod[:,10]))
IIMbestpm = IIMmod[IIMbest[0]]
IIMbestpm = IIMbestpm.flatten()
IIMbestpm = IIMbestpm[0:9]

# calc real values
L=7e8*(19067.0/78882912.0)
Nref=9.132e+02/(4*2.3e-9*L)
IIMnu1_0=IIMbestpm[0]*Nref
IIMnu2_0=IIMbestpm[1]*Nref
IIMnu1=IIMbestpm[2]*Nref
IIMnu2=IIMbestpm[3]*Nref
IIMT0=IIMbestpm[4]*2*Nref
IIMT1=IIMbestpm[4]*2*Nref
IIMm12=(IIMbestpm[5]/(2*Nref)*IIMnu1)
IIMm21=(IIMbestpm[6]/(2*Nref)*IIMnu2)
IIM_uncert_nu1_0=uncertIIM[0]*Nref
IIM_uncert_nu2_0=uncertIIM[1]*Nref
IIM_uncert_nu1=uncertIIM[2]*Nref
IIM_uncert_nu2=uncertIIM[3]*Nref
IIM_uncert_T0=uncertIIM[4]*2*Nref
IIM_uncert_T1=uncertIIM[4]*2*Nref
IIM_uncert_m12=(uncertIIM[5]/(2*Nref)*IIMnu1)
IIM_uncert_m21=(uncertIIM[6]/(2*Nref)*IIMnu2)

# write IM stuff
IIMbest=open('IIM_real_bestparams.txt','a')
IIMbest.write(str(IIMnu1_0)+'\t'+str(IIM_uncert_nu1_0)+'\n'+
    str(IIMnu2_0)+'\t'+str(IIM_uncert_nu2_0)+'\n'+
    str(IIMnu1)+'\t'+str(IIM_uncert_nu1)+'\n'+
    str(IIMnu2)+'\t'+str(IIM_uncert_nu2)+'\n'+
    str(IIMT0)+'\t'+str(IIM_uncert_T0)+'\n'+
    str(IIMT1)+'\t'+str(IIM_uncert_T1)+'\n'+
    str(IIMm12)+'\t'+str(IIM_uncert_m12)+'\n'+
    str(IIMm21)+'\t'+str(IIM_uncert_m21)+'\n')
IIMbest.close()
