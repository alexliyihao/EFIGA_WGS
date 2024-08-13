#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an memverge environmental variable,
#please use P4-41393 (no quotes) here
#THREADS need to be used as an memverge environment variable,
#please use the minimal number of cpus defined here, I used 4.

INPUT_PATH="/data"

REFERENCE_PATH="$INPUT_PATH/genome"
REF="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_DBSNP="${REFERENCE_PATH}/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_KNOWN="${REFERENCE_PATH}/Homo_sapiens_assembly38.known_indels.vcf.gz"
REF_GOLD="${REFERENCE_PATH}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
BED="${REFERENCE_PATH}/hg38_main_chr.bed"
if [ -s $REF ] &&
  [ -s $REF_DBSNP ] &&
  [ -s $REF_KNOWN ] &&
  [ -s $REF_GOLD ] &&
  [ -s $BED ]; then
  echo "Initializing: confirming that all the reference are properly mounted."
else
  echo "Initializing: not all the files are properly mounted"
fi

#PREFIX need to be used as an environment variable
echo "Initializing: the sample prefix from environment variable is ${PREFIX}"

#THREADS need to be used as an environment variable
echo "Initializing: the number of threads from environment variable is ${THREADS}"

WORKING_DIR=${INPUT_PATH}/${PREFIX}
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

echo "Setting: The list of fastq.gz detected"
echo "R1:"
echo $(ls -1l ${WORKING_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz)
echo "R2"
echo $(ls -1l ${WORKING_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz)

PROGRESS_LOG=${WORKING_DIR}/${PREFIX}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG

#-------------------------------- Section 2: Running --------------------------------
start=`date +%s`

echo "Running: checking if alignment finished..."
if (grep -q "mate tages finished: 0" $PROGRESS_LOG); then
  echo "Running: last mate tages finished, going into next step";
  echo "mate tages finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last mate tages not finished, restart";
  echo "Running: Reading Group: @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"
  cat ${WORKING_DIR}/${PREFIX}_bwa_mem_test.sam | \
    $SAMTOOLS sort \
      -@${THREAD} \
      -l 0 \
      -n \
      -T "sorted_step1" \
      -O sam | \
      $SAMBLASTER \
        --addMateTags \
        --ignoreUnmated | \
        $SAMTOOLS sort \
          -@${THREAD} \
          -l 5 \
          -T "sorted_step2" \
          -O bam \
          -o ${WORKING_DIR}/${PREFIX}_with_mated_tags_test.bam
fi

timeAlignment=`date +%s`
runtime=$((timeAlignment-start))
echo "Running: time cost for mate tages ${runtime}s"
