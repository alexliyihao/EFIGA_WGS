# Pipeline stage 1: alignment

## Input:

Fastq files, the pipeline uses Illumina.Inc's naming format ${PREFIX}_S*_L00*_R1_001.fastq.gz and ${PREFIX}_S*_L00*_R2_001.fastq.gz

## Running:
This stage runs the following:
 - bwa mem
 - samtools sort (-n)
 - samblaster --addMateTags
 - samtools sort (no -n anymore)

## Output
  - bam file

## Notes:
  - samblaster only take *.sam* format as inflow and output, so the cache is massive. The overall memory is approximately 20-25GB at its peak.
  - All of these programs more or less can benefited from multi-threading, modify the *THREADS* environment from mmfloat_stage_1_submit.sh if you want to switch to 4 or 16, etc.
