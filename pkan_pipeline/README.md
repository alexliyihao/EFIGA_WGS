# Pkan pipeline

The pipeline running on Columbia Neurology HPC cluster.

Pkan project has 120 samples, but for data re-useability, some larger-scale alternative is used rather then small one (e.g. GATK GenomicsDBImport rather than GATK CombineGVCFs)

The pipeline order should be as follows:
1. pkan_bwa_mem (N_task = 120)
2. pkan_processing (N_task = 120, depends on the output of 1)
3. pkan_statistics (N_task = 120, which depends on an intermediate output of 2)
4. pkan_genotyping (N_task = 25, depends on the output of 2)
5. pkan_merging (N_task = 1, depends on the output of 4)
6. pkan_VariantRecalibrator_indel and pkan_VariantRecalibrator_snp (N_task = 1 for each, can run simultaneously, depends on the output of 5)
7. pkan_VQSR (N_task = 1, depends on the output of 6, both output necessary)
8. pkan_CalculateGenotypePosteriors (N_task = 1)

If the pipeline are written into multiple scripts usually because some steps requires way different hardware requirement from its preceding and following steps. pkan_VariantRecalibrator_indel and pkan_VariantRecalibrator_snp are using different settings.
