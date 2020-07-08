# Analysis of *Syma* spp. genotypes, phenotypes for Linck et al. *in prep*

Repository for analysis of *Syma torotoro* and *Syma megarhyncha* genotype and phenotype data. A digital notebook recording commands for processing and analyzing data is found in `pipeline.md.`  

### Scripts  

(in `scripts/` subdirectory)  

`assembly.R:` Pipeline for assembling sequencing reads against *Halcyon senegalensis* draft genome.  

`genotype_analysis.R:` Analysis of genotypic data and demographic inference output.

`morphology_analysis.R:` Analysis of morphological data.  

`msa.py:` Commands to assembly and align mtDNA genomes.

`oneliners.sh:` Various bash oneliners.  

`plotting.R:` Make figures for manuscript.  

`syma_*.py`: Python scripts to generate joint site frequency spectra and fit demographic models in *moments*

`unified_genotyper_all.sh:` Command to run [UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) on all samples.  

`vocalization_analysis.R:` Analysis of bioacoustic data.  

### Data   

(in `data/` subdirectory)  

`Syma_*:` Shape files for plotting range maps.

`coverage.txt:` Per-sample coverage.  

`meg_sizes.csv:` Population size estimates for *S. megarhyncha* inferred from *moments*.   

`migration_rate.csv:` Inferred rates of gene flow from *moments*.  

`morphology.csv:` Morphological data from specimens.  

`params.csv:` Bioacoustic parameter data.  

`percent_mapped.txt:` Read mapping results.  

`syma_ND2.tree:` ND2 tree.

`syma_nd2_final.xml:` BEAST .xml input file.

`syma_spp_calls.csv:` Processed bioacoustic data.  

`syma_spp_master.csv:` Specimen, labwork data.  

`syma_spp_morphology.csv:` Morphology data before processing.  

`syma_spp_morphology_log.csv:` Log-transformed morphology data.   

`syma_spp_pcs.csv:` PCA results.  

`syma_spp_pcs_old.csv:` PCA results without modern samples.   

`syma_spp_pcs_wgs.csv:` PCA results without hyRAD data.   

`tor_sizes.csv:` Population size estimates for *S. torotoro* inferred from *moments*.   
