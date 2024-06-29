# EFIGA FBS Family WGS

This is a repo for WGS pipeline on *EFIGA FBS FAMILY* dataset on both SLURM and MMVerge pipeline.

The pipeline is trying to reproduce [NIAGADS VCPA pipeline](https://bitbucket.org/NIAGADS/vcpa-pipeline/src/master/) on GATK 4 environment while following the [GATK Best Practices Workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). Some of the steps and benchmarks also from Qi Yu et.al's [tutorial on GATK 4](https://hpc.nih.gov/training/gatk_tutorial/index.html).

## Dataset
 - The local HPC pipeline runs 120 samples for *pkan* project on Columbia Neurology HPC
 - TODO: The AWS cloud pipeline is trying to run ~4000 samples on MMVerge.

## References
 - The pipeline is working on GRCh38_full_analysis_set_plus_decoy_hla.fa
 - Other resource please check the references.txt

## Authorship
The pipeline is implemented and maintained by Yihao Li(yl4326@cumc.columbia.edu) and the pipeline designs are supervised by Dr.Badri Vardarajan from the Gerturde H. Sergievsky Center, Department of Neurology, Columbia University Irving Medical Center. Please email Yihao if you need any help or can help with anything:)
