# Pipeline stage 2: BQSR

## Input:

BAM files from Phase 1

## Running:
This stage runs the following:
 - bamutil dedup_lowmem
 - sambamba index
 - gatk BaseRecalibrator
 - gatk ApplyBQSR

## Output
  - BQSRed bam file

## Notes:
 - Multi-threading: As of the time this repository is prepared. GATK's Spark implementation is still in it's Beta. bamutil is single threaded from it's implementation. So we only used 2nd threads for JAVA's garbage collection. If your instance is memory-emphasized which will be perfect.
 - bamutil dedup_lowmem and sambamba index take very little resources so it start from small instances and supposed to "float up" from Memverge's WaveRaider options.
 - The output for human WGS can be up to 200 GBs, prepare enough volumes.
