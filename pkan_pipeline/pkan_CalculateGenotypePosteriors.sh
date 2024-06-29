#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_CalculateGenotypePosteriors

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

# Memory requirement
#SBATCH --mem=8G
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

echo "Running: checking if CalculateGenotypePosteriors is finished"
if [ -s $WORKING_DIR/pkan_SNP_GenotypePosteriors.vcf.gz ] &&
   (grep -q "GATK ApplyVQSR_snp_step finished: 0" $PROGRESS_LOG); then
  echo "Running: last ApplyVQSR_snp_step finished, going into next step";
  echo "GATK ApplyVQSR_snp_step finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last ApplyVQSR_snp_step not finished, restart"
  rm $WORKING_DIR/pkan_SNP_recalibrated.vcf.gz
  echo "Running: GATK ApplyVQSR_snp_step..."
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xmx5g -Xmx5g -XX:ParallelGCThreads=2" \
    CalculateGenotypePosteriors \
    --variant $WORKING_DIR/pkan_SNP_recalibrated.vcf.gz \
    --output $WORKING_DIR/pkan_SNP_GenotypePosteriors.vcf.gz \
    > $WORKING_DIR/pkan_CalculateGenotypePosteriors.log
  exit_code=$?
  echo "GATK ApplyVQSR_snp_step finished with exit code $exit_code"
  echo "GATK ApplyVQSR_snp_step finished: $exit_code" >> $PROGRESS_LOG
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time cost for GATK CalculateGenotypePosteriors snp mode: ${runtime}s"
