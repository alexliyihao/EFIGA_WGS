# Memverge EFIGA deployment

This folder includes the testing codes deploying the pkan_pipeline on Memverge

TODO:

We are having writing I/O issue with [samtools](https://www.htslib.org/) and [bamutils](https://genome.sph.umich.edu/wiki/BamUtil) writing to AWS S3 bucket.

For memverge colleagues:

1. Mounting setting:

- s3://\<bucket\>/P4-41393:/data/P4-41393
- s3://\<bucket\>/hg38_reference:/data/reference

2. You can skip the following:

- individual_step_cleaned.sh
- individual_step_with_output.sh

and working on:
- Step_1_alignment.sh is split into step 1_1 and 1_2, the output of 1_1 is piped into 1_2
- Step 2_1 is extracted from step_2
