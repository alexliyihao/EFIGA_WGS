#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an memverge environmental variable,
#please use P4-41393 (no quotes) here
#THREADS need to be used as an memverge environment variable,
#please use the minimal number of cpus defined here, I used 2(not 4, all the step 2 related are actually single thread).


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

PROGRESS_LOG=${WORKING_DIR}/${PREFIX}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG

#-------------------------------- Section 2: Running --------------------------------
start=`date +%s`

echo "Running: checking if dedup_lowmem finished..."
if (grep -q "BamUtil dedup_lowmem finished: 0" $PROGRESS_LOG); then
  echo "Running: last bamutil finished, going into next step";
  echo "BamUtil dedup_lowmem finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last bamUtil dedup_lowmem not finished, restart";
  rm -f ${WORKING_DIR}/${PREFIX}_dedup_lowmem.log
  rm -f ${WORKING_DIR}/${PREFIX}_dup_marked.bam
  $BAMUTIL dedup_lowmem \
    --in ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam \
    --out ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
    --log ${WORKING_DIR}/${PREFIX}_dedup_lowmem.log \
    --force \
    --excludeFlags 0xB00
  exit_code=$?
  echo "BamUtil finished with exit code $exit_code"
  echo "BamUtil dedup_lowmem finished: $exit_code" >> $PROGRESS_LOG
fi

timeDedup=`date +%s`
runtime=$((timeDedup-start))
echo "Running: time cost for dedup_lowmem ${runtime}s"
