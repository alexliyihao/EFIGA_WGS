#!/bin/bash

#-------------------------------- Section 1: Initializing and setting --------------------------------
#PREFIX need to be used as an environment variable
#THREADS need to be used as an environment variable

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
if (grep -q "Alignment stepes finished: 0" $PROGRESS_LOG); then
  echo "Running: last Alignment stepes finished, going into next step";
  echo "Alignment stepes finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last Alignment stepes not finished, restart";
  echo "Running: Reading Group: @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"
  $BWA mem \
    -t ${THREADS} \
    -v 1 \
    -K 100000000 \
    -Y \
    -R "@RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}" \
    $REF \
    <(cat ${WORKING_DIR}/${PREFIX}_S*_L00*_R1_001.fastq.gz) \
    <(cat ${WORKING_DIR}/${PREFIX}_S*_L00*_R2_001.fastq.gz) | \
    $SAMTOOLS sort \
      -@${THREADS} \
      -l 0 \
      -n \
      -m 8G \
      -T ${WORKING_DIR}/sorted_step1 \
      -O sam | \
      $SAMBLASTER \
        --addMateTags \
        --ignoreUnmated | \
        $SAMTOOLS sort \
          -@${THREADS} \
          -l 5 \
          -m 8G \
          -T ${WORKING_DIR}/sorted_step2 \
          -O bam \
          -o ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam
  exit_code=$?
  echo "Alignment stepes finished with exit code $exit_code"
  echo "Alignment stepes finished: $exit_code" >> $PROGRESS_LOG
fi

timeAlignment=`date +%s`
runtime=$((timeAlignment-start))
echo "Running: time cost for Alignment stepes ${runtime}s"

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

if [ $exit_code -eq 0 ] &&
   [ -s ${WORKING_DIR}/${PREFIX}_dup_marked.bam ] &&
   (grep -q "BamUtil dedup_lowmem finished: 0" $PROGRESS_LOG) ; then
  echo "removing ${PREFIX}_with_mated_tags.bam"
  rm -f ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam
  exit_code=$?
  echo "with_mated_tags.bam deleted with exit code $exit_code"
  echo "with_mated_tags.bam deleted: $exit_code" >> $PROGRESS_LOG
else
  echo "${WORKING_DIR}/${PREFIX}_with_mated_tags.bam not removed when BamUtil dedup_lowmem finished, probably BamUtil dedup_lowmem did not successfully finished"
  echo "${WORKING_DIR}/${PREFIX}_with_mated_tags.bam not removed when BamUtil dedup_lowmem finished, probably BamUtil dedup_lowmem did not successfully finished" >> $PROGRESS_LOG
fi

timeDedup=`date +%s`
runtime=$((timeDedup-timeAlignment))
echo "Running: time cost for bamutil dedup_lowmem: ${runtime}s"

echo "Running: checking if Sambamba index finished..."
if [ -f ${WORKING_DIR}/${PREFIX}_dup_marked.bam.bai ] &&
   (grep -q "Sambamba Index finished: 0" $PROGRESS_LOG); then
  echo "Running: Sambamba index finished, going into next step";
  echo "Sambamba Index finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last sambamba index not finished, restart"
  rm -f ${WORKING_DIR}/${PREFIX}_dup_marked.bam.bai
  echo "Running: Sambamba Index"
  $SAMBAMBA index \
    -t ${THREADS} \
    ${WORKING_DIR}/${PREFIX}_dup_marked.bam
  exit_code=$?
  echo "Sambamba Index finished with exit code $exit_code"
  echo "Sambamba Index finished: $exit_code" >> $PROGRESS_LOG
fi

timeSambambaIndex=`date +%s`
runtime=$((timeSambambaIndex-timeDedup))
echo "Running: time cost for sambamba index: ${runtime}s"

echo "Running: checking if BaseRecalibrator is finished"
if [ -f ${WORKING_DIR}/${PREFIX}_recal_data.table ] &&
   ((grep -q "SUCCESS" ${WORKING_DIR}/${PREFIX}_BaseRecalibrator.log) ||
   (grep -q "GATK BaseRecalibrator finished: 0" $PROGRESS_LOG)); then
  echo "Running: last BaseRecalibrator finished, going into next step";
  echo "GATK BaseRecalibrator finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last BaseRecalibrator not finished, restart"
  rm -f ${WORKING_DIR}/${PREFIX}_recal_data.table
  echo "Running: GATK BaseRecalibrator..."
  $GATK --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
  	BaseRecalibrator \
  	--reference $REF \
  	--input ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
  	--use-original-qualities \
  	--known-sites $REF_DBSNP \
  	--known-sites $REF_KNOWN \
  	--known-sites $REF_GOLD \
  	--output ${WORKING_DIR}/${PREFIX}_recal_data.table \
  	> ${WORKING_DIR}/${PREFIX}_BaseRecalibrator.log
  exit_code=$?
  echo "GATK BaseRecalibrator finished with exit code $exit_code"
  echo "GATK BaseRecalibrator finished: $exit_code" >> $PROGRESS_LOG
fi

timeBaseRecalibrator=`date +%s`
runtime=$((timeBaseRecalibrator-timeSambambaIndex))
echo "Running: time cost for GATK BaseRecalibrator: ${runtime}s"

echo "Running: checking if GATK ApplyBQSR is finished"
if ([ -s ${WORKING_DIR}/${PREFIX}_BQSR.bam ] ||
   [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz.tbi ]) &&
   (grep -q "GATK ApplyBQSR finished: 0" $PROGRESS_LOG) &&
   (! grep -q "error" ${WORKING_DIR}/${PREFIX}_ApplyBQSR.log); then
  echo "Running: last ApplyBQSR finished, going into next step";
  echo "GATK ApplyBQSR finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last ApplyBQSR not finished, restart"
  rm -f ${WORKING_DIR}/${PREFIX}_BQSR.bam
  rm -f ${WORKING_DIR}/${PREFIX}_BQSR.bai
  echo "Running: GATK ApplyBQSR..."
  $GATK --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
    ApplyBQSR \
    --input ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
    --reference $REF \
    --bqsr-recal-file ${WORKING_DIR}/${PREFIX}_recal_data.table \
    --output ${WORKING_DIR}/${PREFIX}_BQSR.bam \
    > ${WORKING_DIR}/${PREFIX}_ApplyBQSR.log
  exit_code=$?
  echo "GATK ApplyBQSR finished with exit code $exit_code"
  echo "GATK ApplyBQSR finished: $exit_code" >> $PROGRESS_LOG
fi

if #(grep -q "Depth of coverage computing finished: 0" $PROGRESS_LOG) &&
   #(grep -q "Sambamba view MQ30 finished: 0" $PROGRESS_LOG) &&
   (grep -q "GATK ApplyBQSR finished: 0" $PROGRESS_LOG) &&
   #[ -s ${WORKING_DIR}/${PREFIX}_MQ30.txt ] &&
   #[ -s ${WORKING_DIR}/${PREFIX}_depth.txt ] &&
   #[ -s ${WORKING_DIR}/${PREFIX}_depth_of_coverage.txt ] &&
   ([ -s ${WORKING_DIR}/${PREFIX}_BQSR.bam ] || [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz ]); then
  echo "removing ${PREFIX}_dup_marked.bam"
  rm -f ${WORKING_DIR}/${PREFIX}_dup_marked.bam
  rm -f ${WORKING_DIR}/${PREFIX}_dup_marked.bam.bai
  echo "${PREFIX}_dup_marked.bam removed"
  #echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished"
  #echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished" >> $PROGRESS_LOG
  echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam removed when GATK ApplyBQSR finished"
  echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam removed when GATK ApplyBQSR finished" >> $PROGRESS_LOG
else
  #echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably Depth of coverage computing not finished or GATK ApplyBQSR not successfully finished"
  #echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably Depth of coverage computing not finished or GATK ApplyBQSR not successfully finished" >> $PROGRESS_LOG
  echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably GATK ApplyBQSR not successfully finished"
  echo "${WORKING_DIR}/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably GATK ApplyBQSR not successfully finished" >> $PROGRESS_LOG
fi

timeBQSR=`date +%s`
runtime=$((timeBQSR-timeBaseRecalibrator))
echo "Running: time cost for GATK ApplyBQSR: ${runtime}s"

echo "Running: checking if GATK HaplotypeCaller is finished"
if [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz.tbi ] &&
   (grep -q "GATK HaplotypeCaller finished: 0" $PROGRESS_LOG) &&
   (! grep -q "error" ${WORKING_DIR}/${PREFIX}_HaplotypeCaller.log); then
  echo "Running: last HaplotypeCaller finished, going into next step";
  echo "GATK HaplotypeCaller finished: 0" >> $PROGRESS_LOG
else
  echo "Running: GATK HaplotypeCaller not finished, restart"
  rm -f ${WORKING_DIR}/${PREFIX}.g.vcf.gz
  rm -f ${WORKING_DIR}/${PREFIX}.g.vcf.gz.tbi
  echo "Running: GATK HaplotypeCaller..."
  $GATK --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
  	HaplotypeCaller \
  	--reference $REF \
  	--input ${WORKING_DIR}/${PREFIX}_BQSR.bam \
  	--dbsnp $REF_DBSNP \
  	--min-pruning 2 \
  	--standard-min-confidence-threshold-for-calling 30 \
  	--emit-ref-confidence GVCF \
  	--pcr-indel-model NONE \
  	--verbosity INFO \
  	--output ${WORKING_DIR}/${PREFIX}.g.vcf.gz \
  	> ${WORKING_DIR}/${PREFIX}_HaplotypeCaller.log
    exit_code=$?
  echo "GATK HaplotypeCaller finished with exit code $exit_code"
  echo "GATK HaplotypeCaller finished: $exit_code" >> $PROGRESS_LOG
fi

timeHaplotypeCaller=`date +%s`
runtime=$((timeHaplotypeCaller-timeBQSR))
echo "Running: time cost for GATK HaplotypeCaller: ${runtime}s"

echo "Running: checking if BAM convert to CRAM is finished"
if (grep -q "BQSR.bam converted to cram file: 0" $PROGRESS_LOG) &&
   [ -s ${PREFIX}_BQSR.cram ]; then
  echo "Running: BAM convert to CRAM finished, going into next step";
  echo "BQSR.bam converted to cram file: 0" >> $PROGRESS_LOG
else
  echo "Running: ${WORKING_DIR}/${PREFIX}_BQSR.bam not converted to cram, restart"
  rm -f ${WORKING_DIR}/${PREFIX}_BQSR.cram
  echo "Running: converting ${PREFIX}_BQSR.bam to cram"
  $SAMTOOLS view \
    -@{THREADS} \
    --reference ${REF} \
    --cram \
    -o ${WORKING_DIR}/${PREFIX}_BQSR.cram \
    ${WORKING_DIR}/${PREFIX}_BQSR.bam
  exit_code=$?
  echo "BQSR.bam converted to cram file with exit code $exit_code"
  echo "BQSR.bam converted to cram file: $exit_code" >> $PROGRESS_LOG
  if [ $exit_code -eq 0 ] &&
     [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz ] &&
     [ -s ${WORKING_DIR}/${PREFIX}.g.vcf.gz.tbi ] &&
     (! grep -q "error" ${WORKING_DIR}/${PREFIX}_HaplotypeCaller.log) &&
     (grep -q "GATK HaplotypeCaller finished: 0" $PROGRESS_LOG) &&
     (grep -q "BQSR.bam converted to cram file: 0" $PROGRESS_LOG); then
    echo "removing ${PREFIX}_BQSR.bam"
    rm -f ${WORKING_DIR}/${PREFIX}_BQSR.bai
    rm -f ${WORKING_DIR}/${PREFIX}_BQSR.bam
    exit_code=$?
    echo "BQSR.bam deleted with exit code $exit_code"
    echo "BQSR.bam deleted: $exit_code" >> $PROGRESS_LOG
  else
    echo "${WORKING_DIR}/${PREFIX}_BQSR.bam not removed when GATK ApplyBQSR finished, probably GATK HaplotypeCaller did not successfully finished"
    echo "${WORKING_DIR}/${PREFIX}_BQSR.bam not removed when GATK ApplyBQSR finished, probably GATK HaplotypeCaller did not successfully finished" >> $PROGRESS_LOG
  fi
fi


echo "Running: checking if GATK VariantEval is finished"
if [ -s ${WORKING_DIR}/${PREFIX}_variant_eval.table ] &&
  (grep -q "GATK VariantEval finished: 0" $PROGRESS_LOG) &&
  (! grep -q "error" ${WORKING_DIR}/${PREFIX}_VariantEval.log); then
  echo "Running: last HaplotypeCaller finished, going into next step";
  echo "GATK VariantEval finished: 0" >> $PROGRESS_LOG
else
  echo "Running: GATK VariantEval not finished, restart"
  rm -f ${WORKING_DIR}/${PREFIX}_variant_eval.table
  echo "Running: GATK VariantEval"
  $GATK --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
  	VariantEval \
  	--reference $REF \
  	--dbsnp $REF_DBSNP \
  	--gold-standard $REF_GOLD \
  	--eval-module CountVariants \
  	--eval-module VariantSummary \
  	--eval ${WORKING_DIR}/${PREFIX}.g.vcf.gz \
  	--merge-evals \
  	--output ${WORKING_DIR}/${PREFIX}_variant_eval.table \
  	> ${WORKING_DIR}/${PREFIX}_VariantEval.log
    exit_code=$?
  echo "GATK VariantEval finished with exit code $exit_code"
  echo "GATK VariantEval finished: $exit_code" >> $PROGRESS_LOG
fi

timeVariantEval=`date +%s`
runtime=$((timeVariantEval-timeHaplotypeCaller))
echo "Running: time cost for GATK VariantEval: ${runtime}s"

#cp ${WORKING_DIR}/${PREFIX}.g.vcf.gz ${OUTPUT_PATH}/${PREFIX}.g.vcf.gz
#cp ${WORKING_DIR}/${PREFIX}_BQSR.cram ${OUTPUT_PATH}/${PREFIX}_BQSR.cram
#cp ${WORKING_DIR}/${PREFIX}_variant_eval.table ${OUTPUT_PATH}/${PREFIX}_variant_eval.table

#echo "Finishing: Removing all the output"
#rm -f $WORKING_DIR
end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
