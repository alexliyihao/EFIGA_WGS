# Pipeline stage 3: Haplotypecaller

## Input:

BAM files from Phase 2

## Running:
This stage runs the following:
 - gatk HaplotypeCaller
 - samtools view (converting bam to cram)

## Output
  - cram file
  - g.vcf file

## Notes:
 - This stage was a part of stage 2, for gatk HaplotypeCaller is requiring more 20+GBs memory while this type of OnSpot machine is not very often captured in us-east-1. Thus this pipeline often got OOM crash or unwanted *NoAvailableInstance* interruption. So we make it an individual part.
 - Multi-threading: As of the time this repository is prepared. GATK's Spark implementation is still in it's Beta. bamutil is single threaded from it's implementation. So we only used 2nd threads for JAVA's garbage collection. If your instance is memory-emphasized which will be perfect
