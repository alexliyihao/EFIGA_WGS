#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_VariantRecalibrator_snp

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=2

# Memory requirement
#SBATCH --mem=500G
#SABTCH --exclusive
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
REF_DBSNP="${REF_PATH}/Homo_sapiens_assembly38.dbsnp138.vcf"
REF_HIGH="${REF_PATH}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
REF_OMNI="${REF_PATH}/1000G_omni2.5.hg38.vcf.gz"
REF_HAPMAP="${REF_PATH}/hapmap_3.3.hg38.vcf.gz"
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

echo "Running: checking if VariantRecalibrator_mode_snp is finished"
if [ -s $WORKING_DIR/pkan_merged_recalibrate_SNP.recal ] &&
   [ -s $WORKING_DIR/pkan_merged_recalibrate_SNP.tranches ] &&
   [ -s $WORKING_DIR/pkan_merged_ecalibrate_SNP_plots.R ] &&
   (grep -q "GATK VariantRecalibrator_mode_snp finished: 0" $PROGRESS_LOG); then
  echo "Running: last VariantRecalibrator_mode_snp finished, going into next step";
  echo "GATK VariantRecalibrator_mode_snp finished: 0" >> $PROGRESS_LOG
else
  echo "Running: last VariantRecalibrator_mode_snp not finished, restart"
  rm $WORKING_DIR/pkan_merged_recalibrate_SNP.recal
  rm $WORKING_DIR/pkan_merged_recalibrate_SNP.tranches
  rm $WORKING_DIR/pkan_merged_ecalibrate_SNP_plots.R
  echo "Running: GATK VariantRecalibrator_mode_snp..."
  # Compared with GATK3 taking multi-vcfs as input, we gather the input first and the input is via
  # The original code gives --allPoly (--trust-all-polymorphic now) but didn't give any values, the default value doesn't change so copy their code...
  gatk --java-options "-Djava.io.tmpdir=$WORKING_DIR -Xms440G -Xmx440G -XX:ParallelGCThreads=2" \
  	VariantRecalibrator \
    --variant $WORKING_DIR/pkan_merged.vcf.gz \
    --reference $REF \
    --output $WORKING_DIR/pkan_merged_recalibrate_SNP.recal \
    --tranches-file $WORKING_DIR/pkan_merged_recalibrate_SNP.tranches \
    --rscript-file $WORKING_DIR/pkan_merged_ecalibrate_SNP_plots.R \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
    -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
    -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 \
    -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7:$REF_DBSNP \
    -resource:hapmap,known=false,training=true,truth=true,prior=15:$REF_HAPMAP \
    -resource:omni,known=false,training=true,truth=true,prior=12:$REF_OMNI \
    -resource:1000G,known=false,training=true,truth=false,prior=10:$REF_HIGH \
    --use-annotation QD \
    --use-annotation FS \
    --use-annotation DP \
    --use-annotation SOR \
    --use-annotation MQ \
    --use-annotation ReadPosRankSum \
    --use-annotation MQRankSum \
    --use-annotation InbreedingCoeff \
    -mode SNP \
    --trust-all-polymorphic \
    -tranche 100.0 \
    -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.7 \
    -tranche 99.5 \
    -tranche 99.3 \
    -tranche 99.0 \
    -tranche 98.5 \
    -tranche 98.0 \
    -tranche 97.0 \
    -tranche 95.0 \
    -tranche 90.0 \
  	> $WORKING_DIR/pkan_VariantRecalibrator_snp.log
  exit_code=$?
  echo "GATK VariantRecalibrator_mode_snp finished with exit code $exit_code"
  echo "GATK VariantRecalibrator_mode_snp finished: $exit_code" >> $PROGRESS_LOG
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time cost for GATK VariantRecalibrator snp mode: ${runtime}s"
