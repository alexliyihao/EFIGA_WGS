#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_processing_tidy_up

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

#array
#SBATCH --array 3-120

# Memory requirement
#SBATCH --mem=60G
#SBATCH --exclude=node[50-70]

# Time requirement
#SBATCH -t 5-00:00:00

# Change directory to the current
#SBATCH --chdir=./

# Specify the outerr file
#SBATCH --output=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/SLURM_logs/%x/%j.out
#SBATCH --error=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/SLURM_logs/%x/%j.err

# Partition
#SBATCH --partition=GEN

## Section 2: Transfering Data & Path Settings ----------------------------------------------------

INPUTPATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs"
REF="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_DBSNP="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_KNOWN="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Homo_sapiens_assembly38.known_indels.vcf.gz"
REF_GOLD="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
BED="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/data_inflow/hg38_main_chr.bed"
DATA_INFLOW="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/data_inflow/data_inflow.txt"
BAMUTIL="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/external_software/bam"
SAMBAMBA="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/external_software/sambamba_v0.6.6"
THREADS=2


## Section 3: Deploying the Program ----------------------------------------------------

module load GATK/4.2.1.0
module load SAMTOOLS/1.18

## Section 4:Executing the Program ----------------------------------------------------
start=`date +%s`

echo "Setting: the task ID is ${SLURM_ARRAY_TASK_ID}"

readarray -t INPUTFILES < ${DATA_INFLOW}
PREFIX=${INPUTFILES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Setting: the prefix is ${PREFIX}"

WORKING_DIR=${INPUTPATH}/${PREFIX}
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

rm $WORKING_DIR/*.config
rm $WORKING_DIR/*.so

PROGRESS_LOG=$WORKING_DIR/${PREFIX}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG

echo "Running: checking if dedup_lowmem finished..."
if (grep -q "Dedup_LowMem complete!" $WORKING_DIR/${PREFIX}_dedup_lowmem.log) ||
   (grep -q "BamUtil dedup_lowmem finished: 0" $PROGRESS_LOG) ||
   ([ ! -f $WORKING_DIR/${PREFIX}_with_mated_tags.bam ]); then
  echo "Running: last bamutil finished, going into next step";
  echo "BamUtil dedup_lowmem finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last bamUtil dedup_lowmem not finished, restart";
  rm $WORKING_DIR/${PREFIX}_dedup_lowmem.log
  rm $WORKING_DIR/${PREFIX}_dup_marked.bam
  $BAMUTIL dedup_lowmem \
    --in $WORKING_DIR/${PREFIX}_with_mated_tags.bam \
    --out $WORKING_DIR/${PREFIX}_dup_marked.bam \
    --log $WORKING_DIR/${PREFIX}_dedup_lowmem.log \
    --force \
    --excludeFlags 0xB00
  exit_code=$?
  echo "BamUtil finished with exit code $exit_code"
  echo "BamUtil dedup_lowmem finished: $exit_code" >> $PROGRESS_LOG
fi
rm $WORKING_DIR/${PREFIX}_with_mated_tags.bam

timeDedup=`date +%s`
runtime=$((timeDedup-start))
echo "Running: time cost for bamutil dedup_lowmem: ${runtime}s"

echo "Running: checking if Sambamba index finished..."
if [ -f $WORKING_DIR/${PREFIX}_dup_marked.bam.bai ] &&
   (grep -q "Sambamba Index finished: 0" $PROGRESS_LOG); then
  echo "Running: Sambamba index finished, going into next step";
  echo "Sambamba Index finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last sambamba index not finished, restart"
  rm $WORKING_DIR/${PREFIX}_dup_marked.bam.bai
  $SAMBAMBA index \
    -t ${THREADS} \
    $WORKING_DIR/${PREFIX}_dup_marked.bam
  exit_code=$?
  echo "Sambamba Index finished with exit code $exit_code"
  echo "Sambamba Index finished: $exit_code" >> $PROGRESS_LOG
fi

timeSambambaIndex=`date +%s`
runtime=$((timeSambambaIndex-timeDedup))
echo "Running: time cost for sambamba index: ${runtime}s"

echo "Running: checking if BaseRecalibrator is finished"
if [ -f $WORKING_DIR/${PREFIX}_recal_data.table ] &&
   ((grep -q "SUCCESS" $WORKING_DIR/${PREFIX}_BaseRecalibrator.log) ||
   (grep -q "GATK BaseRecalibrator finished: 0" $PROGRESS_LOG)); then
  echo "Running: last BaseRecalibrator finished, going into next step";
  echo "GATK BaseRecalibrator finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last BaseRecalibrator not finished, restart"
  rm $WORKING_DIR/${PREFIX}_recal_data.table
  echo "Running: GATK BaseRecalibrator..."
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
  	BaseRecalibrator \
  	--reference $REF \
  	--input $WORKING_DIR/${PREFIX}_dup_marked.bam \
  	--use-original-qualities \
  	--known-sites $REF_DBSNP \
  	--known-sites $REF_KNOWN \
  	--known-sites $REF_GOLD \
  	--output $WORKING_DIR/${PREFIX}_recal_data.table \
  	> $WORKING_DIR/${PREFIX}_BaseRecalibrator.log
  exit_code=$?
  echo "GATK BaseRecalibrator finished with exit code $exit_code"
  echo "GATK BaseRecalibrator finished: $exit_code" >> $PROGRESS_LOG
fi

timeBaseRecalibrator=`date +%s`
runtime=$((timeBaseRecalibrator-timeSambambaIndex))
echo "Running: time cost for GATK BaseRecalibrator: ${runtime}s"

echo "Running: checking if GATK ApplyBQSR is finished"
if ([ -s $WORKING_DIR/${PREFIX}_BQSR.bam ] ||
   [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz.tbi ]) &&
   (grep -q "GATK ApplyBQSR finished: 0" $PROGRESS_LOG) &&
   (! grep -q "error" $WORKING_DIR/${PREFIX}_ApplyBQSR.log); then
  echo "Running: last ApplyBQSR finished, going into next step";
  echo "GATK ApplyBQSR finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last ApplyBQSR not finished, restart"
  rm $WORKING_DIR/${PREFIX}_BQSR.bam
  rm $WORKING_DIR/${PREFIX}_BQSR.bai
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2" \
	ApplyBQSR \
	--input $WORKING_DIR/${PREFIX}_dup_marked.bam \
	--reference $REF \
 	--bqsr-recal-file $WORKING_DIR/${PREFIX}_recal_data.table \
	--output $WORKING_DIR/${PREFIX}_BQSR.bam \
	> $WORKING_DIR/${PREFIX}_ApplyBQSR.log
  exit_code=$?
  echo "GATK ApplyBQSR finished with exit code $exit_code"
  echo "GATK ApplyBQSR finished: $exit_code" >> $PROGRESS_LOG
fi

if ((grep -q "Depth of coverage computing finished: 0" $PROGRESS_LOG) &&
   (grep -q "Sambamba view MQ30 finished: 0" $PROGRESS_LOG) &&
   (grep -q "GATK ApplyBQSR finished: 0" $PROGRESS_LOG) &&
   [ -s $WORKING_DIR/${PREFIX}_MQ30.txt ] &&
   [ -s $WORKING_DIR/${PREFIX}_depth.txt ] &&
   [ -s $WORKING_DIR/${PREFIX}_depth_of_coverage.txt ] &&
   ([ -s $WORKING_DIR/${PREFIX}_BQSR.bam ] || [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz ])); then
  rm $WORKING_DIR/${PREFIX}_dup_marked.bam
  rm $WORKING_DIR/${PREFIX}_dup_marked.bam.bai
  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished"
  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished" >> $PROGRESS_LOG
else
  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably Depth of coverage computing not finished or GATK ApplyBQSR not successfully finished"
  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam not removed when GATK ApplyBQSR finished, probably Depth of coverage computing not finished or GATK ApplyBQSR not successfully finished" >> $PROGRESS_LOG
fi

timeBQSR=`date +%s`
runtime=$((timeBQSR-timeBaseRecalibrator))
echo "Running: time cost for GATK ApplyBQSR: ${runtime}s"

echo "Running: checking if GATK HaplotypeCaller is finished"
if [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz.tbi ] &&
   (grep -q "GATK HaplotypeCaller finished: 0" $PROGRESS_LOG) &&
   (! grep -q "error" $WORKING_DIR/${PREFIX}_HaplotypeCaller.log); then
  echo "Running: last HaplotypeCaller finished, going into next step";
  echo "GATK HaplotypeCaller finished: 0" >> $PROGRESS_LOG
else
  echo "Running: GATK HaplotypeCaller not finished, restart"
  rm $WORKING_DIR/${PREFIX}.g.vcf.gz
  rm $WORKING_DIR/${PREFIX}.g.vcf.gz.tbi
  # -newQual is now by default since 4.1.0.0 https://gatk.broadinstitute.org/hc/en-us/articles/4404604697243-HaplotypeCaller
  # --genotyping_mode DISCOVERY is removed from the manual but still by default https://gatk.broadinstitute.org/hc/en-us/community/posts/360077648352-GATK-4-2-0-0-Haplotype-caller-genotyping-mode-DISCOVERY-Option
  # --pcr-indel-model NONE for Badri confirmed that the dataset is PCR-free
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
	HaplotypeCaller \
	--reference $REF \
	--input $WORKING_DIR/${PREFIX}_BQSR.bam \
	--dbsnp $REF_DBSNP \
	--min-pruning 2 \
	--standard-min-confidence-threshold-for-calling 30 \
	--emit-ref-confidence GVCF \
	--pcr-indel-model NONE \
	--verbosity INFO \
	--output $WORKING_DIR/${PREFIX}.g.vcf.gz \
	> $WORKING_DIR/${PREFIX}_HaplotypeCaller.log
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
  echo "Running: $WORKING_DIR/${PREFIX}_BQSR.bam not converted to cram, restart"
  rm $WORKING_DIR/${PREFIX}_BQSR.cram
  samtools view \
    -@{THREADS} \
    --reference ${REF} \
    --cram \
    -o $WORKING_DIR/${PREFIX}_BQSR.cram \
    $WORKING_DIR/${PREFIX}_BQSR.bam
  exit_code=$?
  echo "BQSR.bam converted to cram file with exit code $exit_code"
  echo "BQSR.bam converted to cram file: $exit_code" >> $PROGRESS_LOG
  if [ $exit_code -eq 0 ] &&
     [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz ] &&
     [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz.tbi ] &&
     (! grep -q "error" $WORKING_DIR/${PREFIX}_HaplotypeCaller.log) &&
     (grep -q "GATK HaplotypeCaller finished: 0" $PROGRESS_LOG) &&
     (grep -q "BQSR.bam converted to cram file: 0" $PROGRESS_LOG); then
    rm $WORKING_DIR/${PREFIX}_BQSR.bai
    rm $WORKING_DIR/${PREFIX}_BQSR.bam
    exit_code=$?
    echo "BQSR.bam deleted with exit code $exit_code"
    echo "BQSR.bam deleted: $exit_code" >> $PROGRESS_LOG
  else
    echo "$WORKING_DIR/${PREFIX}_BQSR.bam not removed when GATK ApplyBQSR finished, probably GATK HaplotypeCaller did not successfully finished"
    echo "$WORKING_DIR/${PREFIX}_BQSR.bam not removed when GATK ApplyBQSR finished, probably GATK HaplotypeCaller did not successfully finished" >> $PROGRESS_LOG
  fi
fi


echo "Running: checking if GATK VariantEval is finished"
if [ -s $WORKING_DIR/${PREFIX}_variant_eval.table ] &&
  (grep -q "GATK VariantEval finished: 0" $PROGRESS_LOG) &&
  (! grep -q "error" $WORKING_DIR/${PREFIX}_VariantEval.log); then
  echo "Running: last HaplotypeCaller finished, going into next step";
  echo "GATK VariantEval finished: 0" >> $PROGRESS_LOG
else
  echo "Running: GATK VariantEval not finished, restart"
  rm $WORKING_DIR/${PREFIX}_variant_eval.table
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms50G -Xmx50G -XX:ParallelGCThreads=2"  \
	VariantEval \
	--reference $REF \
	--dbsnp $REF_DBSNP \
	--gold-standard $REF_GOLD \
	--eval-module CountVariants \
	--eval-module VariantSummary \
	--eval $WORKING_DIR/${PREFIX}.g.vcf.gz \
	--merge-evals \
	--output $WORKING_DIR/${PREFIX}_variant_eval.table \
	> $WORKING_DIR/${PREFIX}_VariantEval.log
  exit_code=$?
  echo "GATK VariantEval finished with exit code $exit_code"
  echo "GATK VariantEval finished: $exit_code" >> $PROGRESS_LOG
fi

timeVariantEval=`date +%s`
runtime=$((timeVariantEval-timeHaplotypeCaller))
echo "Running: time cost for GATK VariantEval: ${runtime}s"

rm $WORKING_DIR/*.config
rm $WORKING_DIR/*.so

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
