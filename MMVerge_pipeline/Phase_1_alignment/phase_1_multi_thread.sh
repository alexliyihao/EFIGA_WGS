#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an memverge environmental variable,
#please use P4-41393 (no quotes) here
#THREADS need to be used as an memverge environment variable,
#please use the minimal number of cpus defined here, I used 4.


INPUT_PATH="/data"

REFERENCE_PATH="$INPUT_PATH/reference"
REF="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_INDEX="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
if [ -s $REF ] && [ -s $REF_INDEX ]; then
  echo "Initializing: confirming that the reference is properly mounted."
else
  echo "Initializing: the reference is not properly mounted"
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

#PREFIX needs to be used as an environment variable
echo "Initializing: the sample prefix from environment variable is ${PREFIX}"

#THREADS need to be used as an environment variable
echo "Initializing: the number of threads from environment variable is ${THREADS}"

WORKING_DIR=${INPUT_PATH}/${PREFIX}
echo "Initializing: the S3 directory is set to $WORKING_DIR"

#TEMP_DIR="/tmp/data/${PREFIX}"
#echo $(ls /tmp/data/)
TEMP_DIR="/efiga/data/${PREFIX}"
echo $(ls /efiga/data/)
mkdir -p $TEMP_DIR
echo $(ls $TEMP_DIR)
rm -rf $TEMP_DIR/*
echo "Initializing: the temporary working directory is set to $TEMP_DIR"
#cd $TEMP_DIR

echo "Setting: The list of fastq.gz detected"
echo "R1:"
echo $(ls ${WORKING_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz)
echo "R2"
echo $(ls ${WORKING_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz)

PROGRESS_LOG=${WORKING_DIR}/${PREFIX}_Progress.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG

OVERSEER_LOG=${CRED_PATH}/overseer_multi_thread.txt
echo "Setting: the finishing of this pipeline will be report to ${OVERSEER_LOG}"
#-------------------------------- Section 2: Running --------------------------------
start=`date +%s`


echo "Running: checking if alignment finished..."
if (grep -q "Alignment steps finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: last Alignment steps finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "Running: last Alignment steps finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: last Alignment stepes not finished, restart - $(date '+%Y-%m-%d %X')";
  echo "Running: copy the source file to temp directory"
  # Maybe AWS CLI can provide a more reliable HTTPS traffic, especially for large files
  aws s3 cp s3://<your_bucket_here>/${PREFIX}/ ${TEMP_DIR}/ \
    --only-show-errors \
    --recursive \
    --exclude "*" \
    --include "${PREFIX}*.fastq.gz"
  #cp ${WORKING_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz ${TEMP_DIR}/
  #cp ${WORKING_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz ${TEMP_DIR}/
  echo "Running: Reading Group: @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"
  echo "Running: Starting alignment pipeline..."
  $BWA mem \
    -t ${THREADS} \
    -v 1 \
    -K 100000000 \
    -Y \
    -R "@RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}" \
    $REF \
    <(cat ${TEMP_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz) \
    <(cat ${TEMP_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz) | \
    $SAMTOOLS sort \
      -@${THREADS} \
      -l 0 \
      -n \
      -T "${TEMP_DIR}/sorted_step1" \
      -O sam | \
      $SAMBLASTER \
        --addMateTags \
        --ignoreUnmated | \
        $SAMTOOLS sort \
          -@${THREADS} \
          -l 5 \
          -T "${TEMP_DIR}/sorted_step2" \
          -O bam \
          -o ${TEMP_DIR}/${PREFIX}_with_mated_tags.bam
  exit_code=$?
  echo "Alignment steps finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "Alignment steps finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running alignment";
   echo "The pipeline terminated due to non-zero exit running alignment"  >> $PROGRESS_LOG;
   exit 1;
fi

echo "Running: checking if pipeline output is properly generated..."
if [ -s ${TEMP_DIR}/${PREFIX}_with_mated_tags.bam ] &&
  (grep -q "Alignment steps finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: confirmed that the alignment pipeline finished with exit code 0"
  echo "Running: confirmed that the output file ${PREFIX}_with_mated_tags.bam exist and are not empty"
  echo "Running: cleaning the input fastq.gz files"
  # Please note that the .bam file is intentionally left there, only remove the input
  rm -f ${TEMP_DIR}/${PREFIX}*.fastq.gz
  exit_code=$?
  echo "All the fastq.gz input removed with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "All the fastq.gz input removed with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  echo -e "${PREFIX}\tfinished correctly\t$(date '+%Y-%m-%d %X')" >> $OVERSEER_LOG
else
  echo "Not all the output are correctly generated, the pipeline exited with fastq.gz kept at ${TEMP_DIR}"
  echo "Not all the output are correctly generated, the pipeline exited with fastq.gz kept at ${TEMP_DIR}" >> $PROGRESS_LOG
  echo -e "${PREFIX}\tdid not finish correctly\t$(date '+%Y-%m-%d %X')" >> $OVERSEER_LOG
fi

# (only for testing purpose) Maybe AWS CLI can provide a more reliable HTTPS traffic, especially for large files
#aws s3 cp ${TEMP_DIR}/${PREFIX}_with_mated_tags.bam s3://<your_bucket_here>/${PREFIX}/${PREFIX}_with_mated_tags_jfs_cli.bam \
#  --only-show-errors

timeAlignment=`date +%s`
runtime=$((timeAlignment-start))
echo "Running: time cost for Alignment stepes ${runtime}s"
