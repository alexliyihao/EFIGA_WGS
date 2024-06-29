#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_genotyping

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

#array
#SBATCH --array 1-25

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

INPUTPATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs"
META_PATH="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs_meta"
REF="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_DBSNP="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_KNOWN="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Homo_sapiens_assembly38.known_indels.vcf.gz"
REF_GOLD="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/GRCh38_reference_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
GVCF_MAP="/mnt/vast/hpc/bvardarajan_lab/NotBackedUp/pkan/wgs_meta/cohort_list.map"

## Section 3: Deploying the Program ----------------------------------------------------

module load GATK/4.2.1.0

## Section 4:Executing the Program ----------------------------------------------------
start=`date +%s`

echo "Setting: the task ID is ${SLURM_ARRAY_TASK_ID}"

CHRS=($(seq 1 1 22))
CHRS+=('X' 'Y' 'M')

CHR=${CHRS[$SLURM_ARRAY_TASK_ID - 1]}
echo "Setting: working on Chromosome $CHR"

WORKING_DIR=${META_PATH}/chr${CHR}
mkdir -p $WORKING_DIR
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

PROGRESS_LOG=$WORKING_DIR/chr${CHR}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"
echo "Progress Records: $(date "+%Y-%m-%d %X")" >> $PROGRESS_LOG

echo "Running: checking if GenomicsDBImport is finished"
if [ -s $WORKING_DIR/pkan_chr${CHR}.gdb ] &&
   (grep -q "GATK GenomicsDBImport finished: 0" $PROGRESS_LOG); then
  echo "Running: last GenomicsDBImport finished, going into next step";
  echo "GATK GenomicsDBImport finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last GenomicsDBImport not finished, restart"
  #given --overwrite-existing-genomicsdb-workspace flag
  #rm $WORKING_DIR/pkan_chr${CHR}.gdb
  echo "Running: GATK GenomicsDBImport..."
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms6G -Xmx6G -XX:ParallelGCThreads=2" \
    GenomicsDBImport \
    --genomicsdb-workspace-path $WORKING_DIR/pkan_chr${CHR}.gdb \
    --intervals chr${CHR} \
    --sample-name-map $GVCF_MAP \
    --batch-size 0 \
    --overwrite-existing-genomicsdb-workspace \
    --reference ${REF} \
    --validate-sample-name-map \
    --add-output-vcf-command-line \
    > $WORKING_DIR/pkan_chr${CHR}_GenomicsDBImport.log
  exit_code=$?
  echo "GATK GenomicsDBImport finished with exit code $exit_code"
  echo "GATK GenomicsDBImport finished: $exit_code" >> $PROGRESS_LOG
fi

timeGenomicsDBImport=`date +%s`
runtime=$((timeGenomicsDBImport-start))
echo "Running: time cost for GATK GenomicsDBImport: ${runtime}s"

echo "Running: checking if GATK GenotypeGVCFs is finished"
if [ -s $WORKING_DIR/pkan_chr${CHR}.vcf.gz ] &&
   (grep -q "GATK GenotypeGVCFs finished: 0" $PROGRESS_LOG); then
  echo "Running: last GenotypeGVCFs finished, going into next step";
  echo "GATK GenotypeGVCFs finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last GenotypeGVCFs not finished, restart"
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms6G -Xmx6G -XX:ParallelGCThreads=2" \
    GenotypeGVCFs \
    --variant gendb://$WORKING_DIR/pkan_chr${CHR}.gdb \
    --output $WORKING_DIR/pkan_chr${CHR}.vcf.gz \
    --reference $REF \
    --dbsnp $REF_DBSNP \
    --intervals chr${CHR} \
    --add-output-vcf-command-line \
    > $WORKING_DIR/${PREFIX}_GenotypeGVCFs.log
  exit_code=$?
  echo "GATK GenotypeGVCFs finished with exit code $exit_code"
  echo "GATK GenotypeGVCFs finished: $exit_code" >> $PROGRESS_LOG
fi

timeGenotypeGVCFs=`date +%s`
runtime=$((timeGenotypeGVCFs-timeGenomicsDBImport))
echo "Running: time cost for GATK GenotypeGVCFs: ${runtime}s"

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
