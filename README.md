# Analysis of *Syma* spp. genotypes, phenotypes for Linck et al. 2019

Repository for analysis of *Syma torotoro* and *Syma megarhyncha* genotype and phenotype data. A complete digital notebook recording commands for processing and analyzing data is found in `pipeline.md.`  
  
### Supplemental Material  
  
All supplemental figures, tables, and captions are found in the `supplement/` subdirectory.  
    
### Scripts  
  
(in `scripts/` subdirectory)  
     
`*_model.py:` Models for demographic inference with [moments](https://bitbucket.org/simongravel/moments)  
  
`abbababa.R:` Script for performing ABBA-BABA tests.  
  
`assembly.R:` Pipeline for assembling sequencing reads against *Halcyon senegaloides* draft genome.  
  
`genotype_analysis.R:` Analysis of song data.  
    
`morphology_analysis.R:` Analysis of morphological data.  
  
`mtDNA.sh:` Commands to assembly mtDNA genomes.  
  
`oneliners.sh:` Various bash oneliners.  
      
`plotting.R:` Make figures for manuscript.  
  
`smcpp_data.sh:` Prepare data for [SMC++](https://github.com/popgenmethods/smcpp)  
  
`smcpp_run.sh:` Run `SMC++` for bootstrapped data.   
  
`smcpp_boot.R:` Generate list of files for bootstrapping for `SMC++`.  
  
`treemix.R:` Run and plot [treemix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).  
  
`uncertainty.py:` Quantify uncertainty from `moments` analysis; perform likelihood ratio test of top models.  
   
`unified_genotyper_all.sh:` Command to run [UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) on all samples.  
  
`vocalization_analysis.R:` Analysis of bioacoustic data.  
  
### Data   
  
(in `data/` subdirectory)  
     
`chr_labels.csv:` Chromosome coordinates for manhattan plots.  
  
`coverage.txt:` Per-sample coverage.  
  
`d_stat_results.csv:` Individual *D*-statistics from ABBA-BABA tests.  
  
`d_summary.csv:` *D*-statistics summarized.  
  
`moments_plotting.csv:` Data to plot `moments` inference results.   
  
`morphology.csv:` Morphological data from specimens.  
  
`params.csv:` Bioacoustic parameter data.  
  
`percent_mapped.txt:` Read mapping results.  
  
`pi.csv:` Windowed π calculation for *megarhyncha* chromosome 5 and simulated data.  
  
`smcpp_df.csv:` Dataframe for plotting `SMC++` trajectories.  
  
`sweeps.csv:` Selective sweep analysis data.  
  
`syma_spp_calls.csv:` Processed bioacoustic data.  
  
`syma_spp_master.csv:` Specimen, labwork data.  
  
`syma_spp_morphology.csv:` Morphology data before processing.  
  
`syma_spp_pcs.csv:` PCA results.  
  
`window_correlations.csv:` Dataframe for plotting correlations between windowed summary statistics.  
  
`window_stats.csv:` Raw data from sliding window analyses (big).   
  
`window_stats_chr.csv:` Processed data from sliding window analyses for plotting (big).  




      
