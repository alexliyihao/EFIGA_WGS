source /usr/local/bin/_activate_current_env.sh

echo "Setting: the prefix is ${PREFIX}"
INPUT_PATH="/home/$username/data/${PREFIX}"
OUTPUT_PATH=$INPUTPATH
WGS_PATH="/home/$username/wgs"
REFERENCE_PATH="/home/$username/data/genome"
REF="${REFERENCE_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_DBSNP="${REFERENCE_PATH}/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_KNOWN="${REFERENCE_PATH}/Homo_sapiens_assembly38.known_indels.vcf.gz"
REF_GOLD="${REFERENCE_PATH}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
BED="${REFERENCE_PATH}/hg38_main_chr.bed"
THREADS=8 # aws r5/r6/r7 only have 8 cores for 64 GB memory

start=`date +%s`
readarray -t INPUTFILES < ${DATA_INFLOW}

WORKING_DIR=${WGS_PATH}/${PREFIX}
echo "Setting: the working directory is set to $WORKING_DIR"
mkdir -p WORKING_DIR

set -o pipefail
PROGRESS_LOG=${WORKING_DIR}/${PREFIX}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG


echo "Running: Aligning with Reading Group: @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"
$BWA mem \
  -t ${THREAD} \
  -v 1 \
  -K 100000000 \
  -Y \
  -R "@RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}" \
  $REF \
  <(cat ${INPUT_PATH}/${PREFIX}_S*_L00*_R1_001.fastq.gz) \
  <(cat ${INPUT_PATH}/${PREFIX}_S*_L00*_R2_001.fastq.gz) | \
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
        -o ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam
exit_code=$?
echo "Alignment stepes finished with exit code $exit_code"
echo "Alignment stepes finished: $exit_code" >> $PROGRESS_LOG
timeAlignment=`date +%s`
runtime=$((timeAlignment-start))
echo "Running: time cost for Alignment stepes ${runtime}s"


echo "Running: BAMUTIL dedup_lowmem..."
$BAMUTIL dedup_lowmem \
  --in ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam \
  --out ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
  --log ${WORKING_DIR}/${PREFIX}_dedup_lowmem.log \
  --force \
  --excludeFlags 0xB00
exit_code=$?
echo "BamUtil finished with exit code $exit_code"
echo "BamUtil dedup_lowmem finished: $exit_code" >> $PROGRESS_LOG
timeDedup=`date +%s`
runtime=$((timeDedup-timeAlignment))
echo "Running: time cost for bamutil dedup_lowmem: ${runtime}s"


echo "Running: Removing the ${PREFIX}_with_mated_tags.bam"
rm ${WORKING_DIR}/${PREFIX}_with_mated_tags.bam
echo "Running: ${PREFIX}_with_mated_tags.bam removed"


echo "Running: Sambamba Index..."
$SAMBAMBA index \
  -t ${THREADS} \
  ${WORKING_DIR}/${PREFIX}_dup_marked.bam
exit_code=$?
echo "Sambamba Index finished with exit code $exit_code"
echo "Sambamba Index finished: $exit_code" >> $PROGRESS_LOG
timeSambambaIndex=`date +%s`
runtime=$((timeSambambaIndex-timeDedup))
echo "Running: time cost for sambamba index: ${runtime}s"


echo "Running: GATK BaseRecalibrator..."
$GATK --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
	BaseRecalibrator \
	--reference $REF \
	--input ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
	--use-original-qualities \
	--known-sites $REF_DBSNP \
	--known-sites $REF_KNOWN \
	--known-sites $REF_GOLD \
	--output ${WORKING_DIR}/${PREFIX}_recal_data.table
exit_code=$?
echo "GATK BaseRecalibrator finished with exit code $exit_code"
echo "GATK BaseRecalibrator finished: $exit_code" >> $PROGRESS_LOG
timeBaseRecalibrator=`date +%s`
runtime=$((timeBaseRecalibrator-timeSambambaIndex))
echo "Running: time cost for GATK BaseRecalibrator: ${runtime}s"


echo "Running: GATK ApplyBQSR..."
$GATK --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
  ApplyBQSR \
  --input ${WORKING_DIR}/${PREFIX}_dup_marked.bam \
  --reference $REF \
  --bqsr-recal-file ${WORKING_DIR}/${PREFIX}_recal_data.table \
  --output ${WORKING_DIR}/${PREFIX}_BQSR.bam
exit_code=$?
echo "GATK ApplyBQSR finished with exit code $exit_code"
echo "GATK ApplyBQSR finished: $exit_code" >> $PROGRESS_LOG
timeBQSR=`date +%s`
runtime=$((timeBQSR-timeBaseRecalibrator))
echo "Running: time cost for GATK ApplyBQSR: ${runtime}s"


echo "removing ${PREFIX}_dup_marked.bam"
rm ${WORKING_DIR}/${PREFIX}_dup_marked.bam
rm ${WORKING_DIR}/${PREFIX}_dup_marked.bam.bai
echo "${PREFIX}_dup_marked.bam removed"


echo "Running: GATK HaplotypeCaller..."
$GATK --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
	HaplotypeCaller \
	--reference $REF \
	--input ${WORKING_DIR}/${PREFIX}_BQSR.bam \
	--dbsnp $REF_DBSNP \
	--min-pruning 2 \
	--standard-min-confidence-threshold-for-calling 30 \
	--emit-ref-confidence GVCF \
	--pcr-indel-model NONE \
	--verbosity INFO \
	--output ${WORKING_DIR}/${PREFIX}.g.vcf.gz
exit_code=$?
echo "GATK HaplotypeCaller finished with exit code $exit_code"
echo "GATK HaplotypeCaller finished: $exit_code" >> $PROGRESS_LOG
timeHaplotypeCaller=`date +%s`
runtime=$((timeHaplotypeCaller-timeBQSR))
echo "Running: time cost for GATK HaplotypeCaller: ${runtime}s"


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


echo "removing ${PREFIX}_BQSR.bam"
rm ${WORKING_DIR}/${PREFIX}_BQSR.bai
rm ${WORKING_DIR}/${PREFIX}_BQSR.bam
exit_code=$?
echo "BQSR.bam deleted with exit code $exit_code"
echo "BQSR.bam deleted: $exit_code" >> $PROGRESS_LOG


echo "Running: GATK VariantEval"
$GATK --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
	VariantEval \
	--reference $REF \
	--dbsnp $REF_DBSNP \
	--gold-standard $REF_GOLD \
	--eval-module CountVariants \
	--eval-module VariantSummary \
	--eval ${WORKING_DIR}/${PREFIX}.g.vcf.gz \
	--merge-evals \
	--output ${WORKING_DIR}/${PREFIX}_variant_eval.table \
exit_code=$?
echo "GATK VariantEval finished with exit code $exit_code"
echo "GATK VariantEval finished: $exit_code" >> $PROGRESS_LOG

timeVariantEval=`date +%s`
runtime=$((timeVariantEval-timeHaplotypeCaller))
echo "Running: time cost for GATK VariantEval: ${runtime}s"

cp ${WORKING_DIR}/${PREFIX}.g.vcf.gz ${OUTPUT_PATH}/${PREFIX}.g.vcf.gz
cp ${WORKING_DIR}/${PREFIX}_BQSR.cram ${OUTPUT_PATH}/${PREFIX}_BQSR.cram
cp ${WORKING_DIR}/${PREFIX}_variant_eval.table ${OUTPUT_PATH}/${PREFIX}_variant_eval.table
cp $PROGRESS_LOG ${OUTPUT_PATH}/${PREFIX}_Progress_Records.log

if [ -s ${OUTPUT_PATH}/${PREFIX}.g.vcf.gz ] &&
  [ -s ${OUTPUT_PATH}/${PREFIX}_BQSR.cram ] &&
  [ -s ${OUTPUT_PATH}/${PREFIX}_variant_eval.table ]; then
    rm ${INPUTPATH}/${PREFIX}*.fastq.gz
    echo "${PREFIX}*.fastq.gz removed" >> ${OUTPUT_PATH}/${PREFIX}_Progress_Records.log

echo "Finishing: Removing WORKING_DIR ${WORKING_DIR}"
rm $WORKING_DIR
end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
