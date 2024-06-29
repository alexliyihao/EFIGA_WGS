#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_bwa_mem

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=8

#array
#SBATCH --array 1-120

# Memory requirement
#SBATCH --mem=50G

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

INPUTPATH="/mnt/vast/hpc/mayeux_lab/data/EFIGA_FBS_FAMILY_WGS/Nov_2023/pkan"
OUTPUTPATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs"
REF="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
DATA_INFLOW="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/EFIGA_FBS_FAMILY_WGS/data_inflow/data_inflow.txt"
THREAD=8
## Section 3: Deploying the Program ----------------------------------------------------

module load BWA/0.7.17
module load SAMTOOLS/1.18
module load Samblaster/0.1.26

## Section 4:Executing the Program ----------------------------------------------------

start=`date +%s`
readarray -t INPUTFILES < ${DATA_INFLOW}
PREFIX=${INPUTFILES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Setting: the task ID is ${SLURM_ARRAY_TASK_ID}"
echo "Setting: the prefix is ${PREFIX}"
mkdir -p ${OUTPUTPATH}/${PREFIX}
cd ${OUTPUTPATH}/${PREFIX}

echo "Running: bwa-mem with @RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}"

bwa mem \
  -t ${THREAD} \
  -v 1 \
  -K 100000000 \
  -Y \
  -R "@RG\tID:${PREFIX}\tPU:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA\tLB:${PREFIX}" \
  $REF \
  <(cat ${INPUTPATH}/${PREFIX}_S*_L00*_R1_001.fastq.gz) \
  <(cat ${INPUTPATH}/${PREFIX}_S*_L00*_R2_001.fastq.gz) | \
  samtools sort \
    -@${THREAD} \
    -l 0 \
    -n \
    -T "sorted_step1" \
    -O sam | \
    samblaster \
      --addMateTags \
      --ignoreUnmated | \
      samtools sort \
        -@${THREAD} \
        -l 5 \
        -T "sorted_step2" \
        -O bam \
        -o $OUTPUTPATH/${PREFIX}/${PREFIX}_with_mated_tags.bam

end=`date +%s`
runtime=$((end-start))
echo "time duration: ${runtime}s"
