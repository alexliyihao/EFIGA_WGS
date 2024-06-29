#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_ApplyVQSR

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

# Memory requirement
#SBATCH --mem=32G
#SBATCH --exclude node[50-70]

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

META_PATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs/meta"
REF_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome"
REF="${REF_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
THREADS=2

## Section 3: Deploying the Program ----------------------------------------------------

module load GATK/4.2.1.0

## Section 4:Executing the Program ----------------------------------------------------
start=`date +%s`

WORKING_DIR=$META_PATH/merged
mkdir -p $WORKING_DIR
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

PROGRESS_LOG=$WORKING_DIR/merge_VQSR_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"

echo "Running: checking if GATK ApplyVQSR_indel_step is finished"
if [ -s $WORKING_DIR/pkan_INDEL_recalibrated.vcf.gz ] &&
   (grep -q "GATK ApplyVQSR_indel_step finished: 0" $PROGRESS_LOG); then
  echo "Running: last ApplyVQSR_indel_step finished, going into next step";
  echo "GATK ApplyVQSR_indel_step finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last ApplyVQSR_indel_step not finished, restart"
  rm $WORKING_DIR/pkan_INDEL_recalibrated.vcf.gz
  echo "Running: GATK ApplyVQSR_indel_step..."
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xmx20g -Xmx20g -XX:ParallelGCThreads=2" \
    ApplyVQSR \
    --variant $WORKING_DIR/pkan_merged.vcf.gz \
    --output $WORKING_DIR/pkan_INDEL_recalibrated.vcf.gz \
    --reference $REF \
    --recal-file $WORKING_DIR/pkan_merged_recalibrate_INDEL.recal \
    --tranches-file $WORKING_DIR/pkan_merged_recalibrate_INDEL.tranches \
    --mode INDEL \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
    -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
    -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 \
    -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    > $WORKING_DIR/pkan_ApplyVQSR_indel_step.log
  exit_code=$?
  echo "GATK ApplyVQSR_indel_step finished with exit code $exit_code"
  echo "GATK ApplyVQSR_indel_step finished: $exit_code" >> $PROGRESS_LOG
fi

timeApplyVQSR_indel=`date +%s`
runtime=$((timeApplyVQSR_indel-start))
echo "Running: time cost for GATK ApplyVQSR indel mode: ${runtime}s"


echo "Running: checking if GATK ApplyVQSR_snp_step is finished"
if [ -s $WORKING_DIR/pkan_SNP_recalibrated.vcf.gz ] &&
   (grep -q "GATK ApplyVQSR_snp_step finished: 0" $PROGRESS_LOG); then
  echo "Running: last ApplyVQSR_snp_step finished, going into next step";
  echo "GATK ApplyVQSR_snp_step finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last ApplyVQSR_snp_step not finished, restart"
  rm $WORKING_DIR/pkan_SNP_recalibrated.vcf.gz
  echo "Running: GATK ApplyVQSR_snp_step..."
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xmx20g -Xmx20g -XX:ParallelGCThreads=2" \
    ApplyVQSR \
    --variant $WORKING_DIR/pkan_INDEL_recalibrated.vcf.gz \
    --output $WORKING_DIR/pkan_SNP_recalibrated.vcf.gz \
    --reference $REF \
    --recal-file $WORKING_DIR/pkan_merged_recalibrate_SNP.recal \
    --tranches-file $WORKING_DIR/pkan_merged_recalibrate_SNP.tranches \
    --mode SNP \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
    -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
    -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 \
    -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    > $WORKING_DIR/pkan_ApplyVQSR_snp_step.log
  exit_code=$?
  echo "GATK ApplyVQSR_snp_step finished with exit code $exit_code"
  echo "GATK ApplyVQSR_snp_step finished: $exit_code" >> $PROGRESS_LOG
fi

timeApplyVQSR_snp=`date +%s`
runtime=$((timeApplyVQSR_snp-timeApplyVQSR_indel))
echo "Running: time cost for GATK ApplyVQSR snp mode: ${runtime}s"

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
