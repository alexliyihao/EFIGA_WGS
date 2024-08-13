#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an memverge environmental variable,
#please use P4-41393 (no quotes) here
#THREADS need to be used as an memverge environment variable,
#please use the minimal number of cpus defined here, I used 4.

INPUT_PATH="/data"

REFERENCE_PATH="$INPUT_PATH/reference"
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
if (grep -q "bwa mem finished: 0" $PROGRESS_LOG); then
  echo "Running: last bwa mem finished, going into next step";
  echo "bwa mem finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last bwa mem not finished, restart";
  echo "Running: Reading Group: @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"
  $BWA mem \
    -t ${THREADS} \
    -v 1 \
    -K 100000000 \
    -Y \
    -R "@RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}" \
    $REF \
    <(cat ${WORKING_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz) \
    <(cat ${WORKING_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz) \
    > ${WORKING_DIR}/${PREFIX}_bwa_mem_test.sam
  exit_code=$?
  echo "bwa mem finished with exit code $exit_code"
  echo "bwa mem finished: $exit_code" >> $PROGRESS_LOG
fi

timeAlignment=`date +%s`
runtime=$((timeAlignment-start))
echo "Running: time cost for bwa mem step ${runtime}s"
