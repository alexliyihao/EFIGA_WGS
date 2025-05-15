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
if (grep -q "${PREFIX}\tfinished correctly" $OVERSEER_LOG) &&
   (grep -q "All the output copied to S3 bucket with exit code 0" $PROGRESS_LOG) &&
   (grep -q "All the intermediate output removed with exit code 0" $PROGRESS_LOG); then
  echo "Running: the pipeline has already been finished correctly, exit"
  exit 0
fi

if (!(grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG)); then
  echo "Running: GATK ApplyBQSR not finished, the pipeline terminated"
  exit 255
fi


echo "Running: checking if GATK HaplotypeCaller finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK BaseRecalibrator finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK HaplotypeCaller finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: last HaplotypeCaller finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "last HaplotypeCaller finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: GATK HaplotypeCaller not finished, restart"
  rm -f ${TEMP_DIR}/${PREFIX}.g.vcf.gz
  rm -f ${TEMP_DIR}/${PREFIX}.g.vcf.gz.tbi
  echo "Running: GATK HaplotypeCaller..."
  $GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms30G -Xmx120G -XX:ParallelGCThreads=2"  \
  	HaplotypeCaller \
  	--reference $REF \
  	--input ${TEMP_DIR}/${PREFIX}_BQSR.bam \
  	--dbsnp $REF_DBSNP \
  	--min-pruning 2 \
  	--standard-min-confidence-threshold-for-calling 30 \
  	--emit-ref-confidence GVCF \
  	--pcr-indel-model NONE \
  	--verbosity INFO \
  	--output ${TEMP_DIR}/${PREFIX}.g.vcf.gz \
  	> ${TEMP_DIR}/${PREFIX}_HaplotypeCaller.log
    exit_code=$?
  echo "GATK HaplotypeCaller finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "GATK HaplotypeCaller finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeHaplotypeCaller=`date +%s`
runtime=$((timeHaplotypeCaller-timeBQSR))
echo "Running: time spent for GATK HaplotypeCaller: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running GATK HaplotypeCaller";
   echo "The pipeline terminated due to non-zero exit running GATK HaplotypeCaller" >> $PROGRESS_LOG;
   exit 1;
fi

echo "Running: checking if BAM convert to CRAM finished..."
if (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK BaseRecalibrator finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK HaplotypeCaller finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "BQSR.bam converted to cram file finished with exit code 0" $PROGRESS_LOG); then
  echo "Running: BAM convert to CRAM finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')";
  echo "BAM convert to CRAM finished in previous attempts, going into next step - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  exit_code=0
else
  echo "Running: ${TEMP_DIR}/${PREFIX}_BQSR.bam not converted to cram, restart"
  rm -f ${TEMP_DIR}/${PREFIX}_BQSR.cram
  echo "Running: converting ${PREFIX}_BQSR.bam to cram"
  $SAMTOOLS view \
    -@{THREADS} \
    --reference ${REF} \
    --cram \
    -o ${TEMP_DIR}/${PREFIX}_BQSR.cram \
    ${TEMP_DIR}/${PREFIX}_BQSR.bam
  exit_code=$?
  echo "BQSR.bam converted to cram file finished with exit code $exit_code - $(date '+%Y-%m-%d %X')"
  echo "BQSR.bam converted to cram file finished with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
fi

timeSamtoolsToCram=`date +%s`
runtime=$((timeSamtoolsToCram-timeHaplotypeCaller))
echo "Running: time spent for SAMTOOLS converting BAM to CRAM: ${runtime}s"

if [ "$exit_code" -ne "0" ]; then
   echo "Running: the pipeline terminated due to non-zero exit running Samtools generating CRAM files";
   echo "The pipeline terminated due to non-zero exit running Samtools generating CRAM files" >> $PROGRESS_LOG;
   exit 2;
fi

echo "Running: checking if all pipeline output are properly generated..."
if [ -s ${TEMP_DIR}/${PREFIX}.g.vcf.gz ] &&
  [ -s ${TEMP_DIR}/${PREFIX}.g.vcf.gz.tbi ] &&
  [ -s ${TEMP_DIR}/${PREFIX}_BQSR.cram ] &&
  (grep -q "BamUtil dedup_lowmem finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "Sambamba Index finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK BaseRecalibrator finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK ApplyBQSR finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "GATK HaplotypeCaller finished with exit code 0" $PROGRESS_LOG) &&
  (grep -q "BQSR.bam converted to cram file finished with exit code 0" $PROGRESS_LOG) ; then
  echo "Running: confirmed that the pipeline finished with exit code 0"
  echo "Running: confirmed that the output files exist and are not empty"
  echo "Running: copying all the output to S3 bucket"
  echo $(ls $TEMP_DIR)
  #AWS CLI copy result back
  #cp for mounted AWS S3 bucket cannot take continuous ingress streams more than ~89-95GB in experimeents
  aws s3 cp ${TEMP_DIR}/ \
    s3://<your_bucket_here>/${PREFIX}/ \
    --recursive \
    --exclude "*" \
    --include "${PREFIX}.g.vcf.gz" \
    --include "${PREFIX}.g.vcf.gz.tbi" \
    --include "${PREFIX}_BQSR.cram"
  exit_code_aws_cli=$?
  echo "All the output copied to S3 bucket with exit code $exit_code_aws_cli - $(date '+%Y-%m-%d %X')"
  echo "All the output copied to S3 bucket with exit code $exit_code_aws_cli - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
  echo $(ls ${WORKING_DIR})
  if [ $exit_code_aws_cli -eq 0 ] &&
    [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz ] &&
    [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz.tbi ] &&
    [ -s ${WORKING_DIR}/${PREFIX}_BQSR.cram ]; then
    echo "Running: Confirm that the outputs for ${PREFIX} are correctly transferred to S3 bucket"
    echo "Running: cleaning temp working directory"
    rm -rf ${TEMP_DIR}
    exit_code=$?
    echo "All the intermediate output removed with exit code $exit_code - $(date '+%Y-%m-%d %X')"
    echo "All the intermediate output removed with exit code $exit_code - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
    echo -e "${PREFIX}\tfinished correctly\t$(date '+%Y-%m-%d %X')" >> $OVERSEER_LOG
  else
    echo "Not all result was correctly transferred via AWS CLI s3 cp operation. The intermediate output for ${PREFIX} was not removed - $(date '+%Y-%m-%d %X')"
    echo "Not all result was correctly transferred via AWS CLI s3 cp operation. The intermediate output for ${PREFIX} was not removed - $(date '+%Y-%m-%d %X')" >> $PROGRESS_LOG
    echo -e "${PREFIX}\tdid not finish correctly\t$(date '+%Y-%m-%d %X')" >> $OVERSEER_LOG
  fi
else
  echo "Not all the output are correctly generated, the pipeline exited with intermediate result kept at ${TEMP_DIR}"
  echo "Not all the output are correctly generated, the pipeline exited with intermediate result kept at ${TEMP_DIR}" >> $PROGRESS_LOG
  echo -e "${PREFIX}\tdid not finish correctly\t$(date '+%Y-%m-%d %X')" >> $OVERSEER_LOG
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
