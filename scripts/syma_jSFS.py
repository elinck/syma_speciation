import allel
import moments
import numpy as np

# load vcf and extract folded sfs; write to file
vcf_file = "syma.moments.bi.recode.vcf"
f = allel.read_vcf(vcf_file, fields='*')
gt = f['calldata/GT']
gt = allel.GenotypeArray(gt)
toro=gt[:, np.r_[0,1,2,7,12,13,17,18]]
mega=gt[:, np.r_[3,4,5,6,8,9,10,14,15,16]]
gn_t = toro.to_n_alt(fill=-1)
gn_m = mega.to_n_alt(fill=-1)
t_unlinked = allel.locate_unlinked(gn_t)
m_unlinked = allel.locate_unlinked(gn_m)
pass_linkage = np.logical_and(t_unlinked, m_unlinked)
final = gt.compress(pass_linkage, axis=0)
toro2=final[:, np.r_[0,1,2,7,12,13,17,18]]
mega2=final[:, np.r_[3,4,5,6,8,9,10,14,15,16]]
toro_ac = toro2.count_alleles()
mega_ac = mega2.count_alleles()
fsfs = allel.joint_sfs_folded(toro_ac, mega_ac)
m_fsfs = moments.Spectrum(fsfs)
m_fsfs.to_file("syma_unlinked_sfs.txt")
