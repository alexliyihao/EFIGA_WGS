#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_merging

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

# Memory requirement
#SBATCH --mem=500G
#SBATCH --exclusive
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

INPUTPATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs"
META_PATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs_meta"
PICARD="/mnt/vast/hpc/mgmt/apps/Picard3/picard.jar"
REF_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome"
REF="${REF_PATH}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
THREADS=2


## Section 3: Deploying the Program ----------------------------------------------------

module load Picard

## Section 4:Executing the Program ----------------------------------------------------
start=`date +%s`

WORKING_DIR=$META_PATH/merged
mkdir -p $WORKING_DIR
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

PROGRESS_LOG=$WORKING_DIR/merge_VQSR_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"

#from https://hpc.nih.gov/training/gatk_tutorial/genotype-gvcfs.html#optimized-script-4
# For GATK 4 no longer take multi vcfs this should be done by the side
echo "Running: checking if Picard GatherVcfs is finished"
if [ -s $WORKING_DIR/pkan_merged.vcf.gz ] &&
   (grep -q "Picard GatherVcfs finished: 0" $PROGRESS_LOG); then
  echo "Running: last Picard GatherVcfs finished, going into next step";
  echo "Picard GatherVcfs finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last Picard GatherVcfs not finished, restart"
  echo "Running: Picard GatherVcfs..."
  java -Xms450G -Xmx450G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$WORKING_DIR \
    -jar $PICARD \
    GatherVcfs \
    --INPUT $META_PATH/chr1/pkan_chr1.vcf.gz \
    --INPUT $META_PATH/chr2/pkan_chr2.vcf.gz \
    --INPUT $META_PATH/chr3/pkan_chr3.vcf.gz \
    --INPUT $META_PATH/chr4/pkan_chr4.vcf.gz \
    --INPUT $META_PATH/chr5/pkan_chr5.vcf.gz \
    --INPUT $META_PATH/chr6/pkan_chr6.vcf.gz \
    --INPUT $META_PATH/chr7/pkan_chr7.vcf.gz \
    --INPUT $META_PATH/chr8/pkan_chr8.vcf.gz \
    --INPUT $META_PATH/chr9/pkan_chr9.vcf.gz \
    --INPUT $META_PATH/chr10/pkan_chr10.vcf.gz \
    --INPUT $META_PATH/chr11/pkan_chr11.vcf.gz \
    --INPUT $META_PATH/chr12/pkan_chr12.vcf.gz \
    --INPUT $META_PATH/chr13/pkan_chr13.vcf.gz \
    --INPUT $META_PATH/chr14/pkan_chr14.vcf.gz \
    --INPUT $META_PATH/chr15/pkan_chr15.vcf.gz \
    --INPUT $META_PATH/chr16/pkan_chr16.vcf.gz \
    --INPUT $META_PATH/chr17/pkan_chr17.vcf.gz \
    --INPUT $META_PATH/chr18/pkan_chr18.vcf.gz \
    --INPUT $META_PATH/chr19/pkan_chr19.vcf.gz \
    --INPUT $META_PATH/chr20/pkan_chr20.vcf.gz \
    --INPUT $META_PATH/chr21/pkan_chr21.vcf.gz \
    --INPUT $META_PATH/chr22/pkan_chr22.vcf.gz \
    --INPUT $META_PATH/chrX/pkan_chrX.vcf.gz \
    --INPUT $META_PATH/chrY/pkan_chrY.vcf.gz \
    --INPUT $META_PATH/chrM/pkan_chrM.vcf.gz \
    --OUTPUT $WORKING_DIR/pkan_merged.vcf.gz \
    --TMP_DIR $WORKING_DIR \
    --REFERENCE_SEQUENCE $REF \
    > $WORKING_DIR/Picard_GatherVcfs.log
  exit_code=$?
  echo "Picard GatherVcfs finished with exit code $exit_code"
  echo "Picard GatherVcfs finished: $exit_code" >> $PROGRESS_LOG
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time cost for Picard GatherVcfs: ${runtime}s"
