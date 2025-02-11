#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an memverge environmental variable,
#please use P4-41393 (no quotes) here
#THREADS need to be used as an memverge environment variable,
#please use the minimal number of cpus defined here, I used 2(not 4, all the step 2 related are actually single thread).


INPUT_PATH="/data"

echo "Initializing: adding reference paths"
REFERENCE_PATH="$INPUT_PATH/reference"
REF="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_INDEX="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
REF_DBSNP="${REFERENCE_PATH}/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_KNOWN="${REFERENCE_PATH}/Homo_sapiens_assembly38.known_indels.vcf.gz"
REF_GOLD="${REFERENCE_PATH}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
BED="${REFERENCE_PATH}/hg38_main_chr.bed"
if [ -s $REF ] &&
  [ -s $REF_INDEX ] &&
  [ -s $REF_DBSNP ] &&
  [ -s $REF_KNOWN ] &&
  [ -s $REF_GOLD ] &&
  [ -s $BED ]; then
  echo "Initializing: confirming that all the reference are properly mounted."
else
  echo "Initializing: WARNING - not all the files are properly mounted"
fi

echo "Initializing: adding AWS CLI credentials"
CRED_PATH="${INPUT_PATH}/creds"
CRED_CREDENTIALS="${CRED_PATH}/credentials"
readarray -t CRED_INFO < $CRED_CREDENTIALS
aws configure set aws_access_key_id ${CRED_INFO[0]}
aws configure set aws_secret_access_key ${CRED_INFO[1]}
aws configure set region ${CRED_INFO[2]}
aws configure set output ${CRED_INFO[3]}
echo "Initializing: region set to ${CRED_INFO[2]}"

#PREFIX need to be used as an environment variable
echo "Initializing: the sample prefix from environment variable is ${PREFIX}"

#THREADS need to be used as an environment variable
echo "Initializing: the number of threads from environment variable is ${THREADS}"

WORKING_DIR=${INPUT_PATH}/${PREFIX}
echo "Initializing: the S3 directory is set to $WORKING_DIR"

#TEMP_DIR="/tmp/data/${PREFIX}"
#echo $(ls /tmp/data/)
TEMP_DIR="/efiga/data/${PREFIX}"
echo "Initializing: the file in ${TEMP_DIR} detected:"
echo $(ls $TEMP_DIR)
echo "Initializing: the temporary working directory is set to $TEMP_DIR"
cd $TEMP_DIR

PROGRESS_LOG=${WORKING_DIR}/${PREFIX}_Progress.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG

OVERSEER_LOG=${CRED_PATH}/overseer_single_thread.txt
echo "Setting: the finishing of this pipeline will be report to ${OVERSEER_LOG}"

#-------------------------------- Section 2: Running --------------------------------
start=`date +%s`

if (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG) &&
   (grep -q "${PREFIX}_dup_marked.bam deleted with exit code 0" $PROGRESS_LOG); then
  echo "Running: the pipeline has already been finished correctly, exit"
  exit 0
fi

if ([ ! -f "${TEMP_DIR}/${PREFIX}_with_mated_tags.bam" ] ||
  !(grep -q "Alignment steps finished with exit code 0" $PROGRESS_LOG)) &&
  !(grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: no input file ${PREFIX}_with_mated_tags.bam detected while no previous running finished BamUtil dedup_lowmem, the pipeline terminated";
  exit 255
fi

echo "Running: checking if dedup_lowmem finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: last bamutil dedup_lowmem finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "last bamutil dedup_lowmem finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: last bamUtil dedup_lowmem not finished, restart - $(date '+%Y-%m-%d %X')";
  rm -f ${TEMP_DIR}/${PREFIX}_dedup_lowmem.log
  rm -f ${TEMP_DIR}/${PREFIX}_dup_marked.bam
  $BAMUTIL dedup_lowmem \
    --in ${TEMP_DIR}/${PREFIX}_with_mated_tags.bam \
    --out ${TEMP_DIR}/${PREFIX}_dup_marked.bam \
    --log ${TEMP_DIR}/${PREFIX}_dedup_lowmem.log \
    --force \
    --excludeFlags 0xB00
  exit_code=$?
  echo "BamUtil dedup_lowmem finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "BamUtil dedup_lowmem finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeDedup=`date +%s`
runtime=$((timeDedup-start))
echo "Running: time spent for bamutil dedup_lowmem: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running dedup_lowmem";
   echo "The pipeline terminated due to non-zero exit running dedup_lowmem"  >> $PROGRESS_LOG;
   exit 1;
fi

# once ${PREFIX}_dup_marked.bam is generated, ${PREFIX}_with_mated_tags.bam is safe to be deleted
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  [ -f "${TEMP_DIR}/${PREFIX}_with_mated_tags.bam" ]; then
  echo "Running: removing ${PREFIX}_with_mated_tags.bam"
  rm -f ${TEMP_DIR}/${PREFIX}_with_mated_tags.bam
  exit_code=$?
  echo "${PREFIX}_with_mated_tags.bam deleted with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "${PREFIX}_with_mated_tags.bam deleted with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
else
  echo "${TEMP_DIR}/${PREFIX}_with_mated_tags.bam not removed when BamUtil dedup_lowmem finished - $(date '+%Y-%m-%d %X')"
  echo "${TEMP_DIR}/${PREFIX}_with_mated_tags.bam not removed when BamUtil dedup_lowmem finished - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit removing ${PREFIX}_with_mated_tags.bam";
   echo "The pipeline terminated due to non-zero exit removing ${PREFIX}_with_mated_tags.bam"  >> $PROGRESS_LOG;
   exit 2;
fi


echo "Running: checking if Sambamba index finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: Sambamba index finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "Sambamba index finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: last sambamba index not finished, restart"
  rm -f ${TEMP_DIR}/${PREFIX}_dup_marked.bam.bai
  echo "Running: Sambamba Index"
  $SAMBAMBA index \
    -t ${THREADS} \
    ${TEMP_DIR}/${PREFIX}_dup_marked.bam
  exit_code=$?
  echo "Sambamba Index finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "Sambamba Index finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeSambambaIndex=`date +%s`
runtime=$((timeSambambaIndex-timeDedup))
echo "Running: time spent for sambamba index: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running Sambamba index";
   echo "The pipeline terminated due to non-zero exit running Sambamba index" >> $PROGRESS_LOG;
   exit 3;
fi

echo "Running: checking if BaseRecalibrator finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK BaseRecalibrator finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: BaseRecalibrator finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "BaseRecalibrator finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: last BaseRecalibrator not finished, restart"
  rm -f ${TEMP_DIR}/${PREFIX}_recal_data.table
  echo "Running: GATK BaseRecalibrator..."
  $GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms30G -Xmx120G -XX:ParallelGCThreads=2" \
  	BaseRecalibrator \
  	--reference $REF \
  	--input ${TEMP_DIR}/${PREFIX}_dup_marked.bam \
  	--use-original-qualities \
  	--known-sites $REF_DBSNP \
  	--known-sites $REF_KNOWN \
  	--known-sites $REF_GOLD \
  	--output ${TEMP_DIR}/${PREFIX}_recal_data.table \
  	> ${TEMP_DIR}/${PREFIX}_BaseRecalibrator.log
  exit_code=$?
  echo "GATK BaseRecalibrator finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "GATK BaseRecalibrator finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeBaseRecalibrator=`date +%s`
runtime=$((timeBaseRecalibrator-timeSambambaIndex))
echo "Running: time spent for GATK BaseRecalibrator: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running GATK BaseRecalibrator";
   echo "The pipeline terminated due to non-zero exit running GATK BaseRecalibrator" >> $PROGRESS_LOG;
   exit 4;
fi

echo "Running: checking if GATK ApplyBQSR finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK BaseRecalibrator finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: GATK ApplyBQSR finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "GATK ApplyBQSR finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: last ApplyBQSR not finished, restart"
  rm -f ${TEMP_DIR}/${PREFIX}_BQSR.bam
  rm -f ${TEMP_DIR}/${PREFIX}_BQSR.bai
  echo "Running: GATK ApplyBQSR..."
  $GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms30G -Xmx120G -XX:ParallelGCThreads=2" \
    ApplyBQSR \
    --input ${TEMP_DIR}/${PREFIX}_dup_marked.bam \
    --reference $REF \
    --bqsr-recal-file ${TEMP_DIR}/${PREFIX}_recal_data.table \
    --output ${TEMP_DIR}/${PREFIX}_BQSR.bam \
    > ${TEMP_DIR}/${PREFIX}_ApplyBQSR.log
  exit_code=$?
  echo "GATK ApplyBQSR finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "GATK ApplyBQSR finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeBQSR=`date +%s`
runtime=$((timeBQSR-timeBaseRecalibrator))
echo "Running: time spent for GATK ApplyBQSR: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running GATK ApplyBQSR";
   echo "The pipeline terminated due to non-zero exit running GATK ApplyBQSR" >> $PROGRESS_LOG;
   exit 5;
fi

# only when GATK ApplyBQSR is finished in TEMP directory and the file is in TEMP_DIR
if ([ -s ${TEMP_DIR}/${PREFIX}_dup_marked.bam.bai ]) &&
  ([ -s ${TEMP_DIR}/${PREFIX}_dup_marked.bam ]); then
  echo "Running: detected ${PREFIX}_dup_marked.bam and ${PREFIX}_dup_marked.bam.bai after ApplyBQSR"
  echo "Running: checking if last ApplyBQSR finished successfully..."
  # if the ApplyBQSR finished with exit code 0, we have a result saved somewhere
  # (maybe in juicefs, maybe already transferred to S3 bucket)
  if (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG); then
    echo "Running: last ApplyBQSR finished successfully, removing ${PREFIX}_dup_marked.bam"
    rm -f ${TEMP_DIR}/${PREFIX}_dup_marked.bam.bai
    rm -f ${TEMP_DIR}/${PREFIX}_dup_marked.bam
    exit_code=$?
    echo "${PREFIX}_dup_marked.bam deleted with exit code $exit_code - $(date '+%Y-%m-%d %X')"
    echo "${PREFIX}_dup_marked.bam deleted with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  else
    echo "${PREFIX}_dup_marked.bam was not removed when GATK ApplyBQSR finished"
    echo "${PREFIX}_dup_marked.bam was not removed when GATK ApplyBQSR finished" >> $PROGRESS_LOG
  fi
else
  echo "Running: No ${PREFIX}_dup_marked.bam and ${PREFIX}_dup_marked.bam.bai detected, skip the cleaning - $(date '+%Y-%m-%d %X')"
  echo "No ${PREFIX}_dup_marked.bam and ${PREFIX}_dup_marked.bam.bai detected, skip the cleaning - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit removing ${PREFIX}_dup_marked.bam";
   echo "The pipeline terminated due to non-zero exit removing ${PREFIX}_dup_marked.bam" >> $PROGRESS_LOG;
   exit 6;
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
