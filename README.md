# Analysis of *Syma* spp. genotypes, phenotypes for Linck et al. 2019

Repository for analysis of *Syma torotoro* and *Syma megarhyncha* genotype and phenotype data. A digital notebook recording commands for processing and analyzing data is found in `pipeline.md.`  

### Scripts  

(in `scripts/` subdirectory)  

`IM_moments.py:` Model for demographic inference with [moments](https://bitbucket.org/simongravel/moments)  

`IM_moments_bootstraps.py:` Model for demographic inference with [moments](https://bitbucket.org/simongravel/moments)  

`assembly.R:` Pipeline for assembling sequencing reads against *Halcyon senegaloides* draft genome.  

`d_tests.py` Script for performing *D*-tests.  

`genotype_analysis.R:` Analysis of genotypic data.

`morphology_analysis.R:` Analysis of morphological data.  

`msa.py:` Commands to assembly and align mtDNA genomes.

`oneliners.sh:` Various bash oneliners.  

`plotting.R:` Make figures for manuscript.  

`ref_nuDNA.sh:` Generate alignment for SVDquartets.

`smcpp_run.sh:` Run `SMC++` for bootstrapped data.   

`unified_genotyper_all.sh:` Command to run [UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) on all samples.  

`vocalization_analysis.R:` Analysis of bioacoustic data.  

### Data   

(in `data/` subdirectory)  

`Syma_*:` Shape files for plotting range maps.

`chr_labels.csv:` Chromosome coordinates for manhattan plots.  

`coverage.txt:` Per-sample coverage.  

`d_d3*csv:` *D*-test results.   

`morphology.csv:` Morphological data from specimens.  

`params.csv:` Bioacoustic parameter data.  

`percent_mapped.txt:` Read mapping results.  

`syma.svd.tre:` SVDquartets lineage tree.

`syma_ND2.tree:` ND2 tree.

`syma_nd2_final.xml:` BEAST .xml input file.

`syma_nuc.fasta:` Fasta alignment of nuclear DNA data.

`syma_nuc.nexus:` Nexus alignment of nuclear DNA data for SVDquartets.

`syma_nuc_d3.fasta` Fasta alignment of nuclear DNA data for *D*-tests.

`syma_spp_calls.csv:` Processed bioacoustic data.  

`syma_spp_master.csv:` Specimen, labwork data.  

`syma_spp_morphology.csv:` Morphology data before processing.  

`syma_spp_pcs.csv:` PCA results.  

`syma_sptree.nexus:` Nexus alignment of nuclear DNA data for SVDquartets species tree analysis.  

`syma_sptree_concordant.tre:` Species tree from SVDquartets.

`window_correlations.csv:` Dataframe for plotting correlations between windowed summary statistics.  

`window_stats.csv:` Raw data from sliding window analyses (big).   

`window_stats_chr.csv:` Processed data from sliding window analyses for plotting (big).  
