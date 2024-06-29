#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=pkan_statistics

# Multicore
#SBATCH --nodes=1
#SBATCH --ntasks=8

#array
#SBATCH --array 1-120

# Memory requirement
#SBATCH --mem=20G

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
THREADS=8


## Section 3: Deploying the Program ----------------------------------------------------

module load SAMTOOLS/1.18

function get_read_size(){
 local BAM="$1"
 [[ -z "$BAM" ]] && exit 1
 echo -n $(samtools view "$BAM" | head -1 | awk '{print length($10);}')
}

## Section 4:Executing the Program ----------------------------------------------------

start=`date +%s`

echo "Setting: the task ID is ${SLURM_ARRAY_TASK_ID}"

readarray -t INPUTFILES < ${DATA_INFLOW}
PREFIX=${INPUTFILES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Setting: the prefix is ${PREFIX}"

WORKING_DIR=${INPUTPATH}/${PREFIX}
echo "Setting: the working directory is set to $WORKING_DIR"
cd $WORKING_DIR

PROGRESS_LOG=$WORKING_DIR/${PREFIX}_Progress_Records.log
echo "Setting: the progress report is written to $PROGRESS_LOG"

echo "Running: checking if Sambamba index finished..."
if ([ -f $WORKING_DIR/${PREFIX}_dup_marked.bam.bai ] &&
   (grep -q "Sambamba Index finished: 0" $PROGRESS_LOG)); then
  echo "Running: Sambamba index finished, start the statistics gathering"

  if (grep -q "Sambamba view MQ30 finished: 0" $PROGRESS_LOG); then
    echo "Sambamba view MQ30 finished, going into next step"
    echo "Sambamba view MQ30 finished: 0" >> $PROGRESS_LOG
  else
    echo "Running: last sambamba view MQ>30 not finished, restart"
    $SAMBAMBA view \
    	-t ${THREADS} \
    	-c \
    	-F ' mapping_quality>=30 ' \
    	$WORKING_DIR/${PREFIX}_dup_marked.bam \
    	> $WORKING_DIR/${PREFIX}_MQ30.txt
    exit_code=$?
    echo "Sambamba view MQ30 finished with exit code $exit_code"
    echo "Sambamba view MQ30 finished: $exit_code" >> $PROGRESS_LOG
  fi
  timeMq=`date +%s`
  runtime=$((timeMq-start))
  echo "Running: time cost for sambamba viewing MQ30: ${runtime}s"


  if (grep -q "Sambamba depth region finished: 0" $PROGRESS_LOG); then
    echo "Sambamba depth region finished, going into next step"
    echo "Sambamba depth region finished: 0" >> $PROGRESS_LOG
  else
    echo "Running: last sambamba depth region not finished, restart"
    $SAMBAMBA depth region \
      -t ${THREADS} \
      --combined \
      -L $BED \
      -T 5 -T 10 -T 20 -T 30 -T 40 -T 50 \
      $WORKING_DIR/${PREFIX}_dup_marked.bam \
      > $WORKING_DIR/${PREFIX}_depth.txt
    exit_code=$?
    echo "Sambamba depth region finished with exit code $exit_code"
    echo "Sambamba depth region finished: $exit_code" >> $PROGRESS_LOG
  fi
  timeDepthRegion=`date +%s`
  runtime=$((timeDepthRegion-timeMq))
  echo "Running: time cost for depth region: ${runtime}s"


  if (grep -q "Depth of coverage computing finished: 0" $PROGRESS_LOG); then
    echo "Depth of coverage computing finished, going into next step"
    echo "Depth of coverage computing finished: 0" >> $PROGRESS_LOG
  else
    echo "Running: last Depth of coverage computing not finished, restart"
    DOC_INFILE=$WORKING_DIR/${PREFIX}_depth.txt
    DOC_OUTFILE=$WORKING_DIR/${PREFIX}_depth_of_coverage.txt
    if ([ -f $DOC_OUTFILE ] && (grep -Eq "^RSIZE\t([0-9]+)" $DOC_OUTFILE)); then
      RSIZE=$(grep -E "^RSIZE\t([0-9]+)" $DOC_OUTFILE | sed -nr 's/^RSIZE\t([0-9]+)/\1/p')
      exit_code=$?
      echo "Running: obtained the RSIZE=${RSIZE} from a exist but unfinished ${PREFIX}_depth_of_coverage.txt"
    else
      RSIZE=$(get_read_size $WORKING_DIR/${PREFIX}_dup_marked.bam)
      exit_code=$?
      echo "Running: obtained the RSIZE=${RSIZE} from samtools"
    fi
    echo "DepthOfCoverage stats:" > $DOC_OUTFILE
    echo -e "RSIZE\t$RSIZE" >> $DOC_OUTFILE
    ATGC_refSize=2745186602 #hg38 - w/o N's
    echo -e "ATGC_refSize\t$ATGC_refSize" >> $DOC_OUTFILE
    #chrom chromStart      chromEnd        readCount       meanCoverage     5x	10x	20x	30x	40x	50x
    # $1   $2              $3              $4              $5              $6       $7      $8      $9      $10     $11
    P5=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize '{SUM += $6 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P5\t$P5" >> $DOC_OUTFILE
    P10=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize  '{SUM += $7 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P10\t$P10" >> $DOC_OUTFILE
    P20=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize '{SUM += $8 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P20\t$P20" >> $DOC_OUTFILE
    P30=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize '{SUM += $9 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P30\t$P30" >> $DOC_OUTFILE
    P40=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize '{SUM += $10 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P40\t$P40" >> $DOC_OUTFILE
    P50=$(grep "^chr" $DOC_INFILE | awk -v GSIZE=$ATGC_refSize '{SUM += $11 * $3} END{ printf "%0.2f", SUM/GSIZE }')
    echo -e "P50\t$P50" >> $DOC_OUTFILE
    # average coverage: chromosome read count * read size summed, divided by number of non-N nucleotides
    AVGCOV=$(grep "^chr" $DOC_INFILE | awk -v RSIZE=$RSIZE -v GSIZE=$ATGC_refSize  '{SUM += $4 * RSIZE} END{ printf "%0.2f", SUM/GSIZE }')
    echo "Depth of coverage computing finished with exit code $exit_code"
    echo -e "AVG_COV\t$AVGCOV" >> $DOC_OUTFILE
    echo "Depth of coverage computing finished: $exit_code" >> $PROGRESS_LOG
  fi
  timeDepthCoverage=`date +%s`
  runtime=$((timeDepthCoverage-timeDepthRegion))
  echo "Running: time cost for Depth of coverage computing: ${runtime}s"

#if ((grep -q "Depth of coverage computing finished: 0" $PROGRESS_LOG) &&
#   (grep -q "Sambamba view MQ30 finished: 0" $PROGRESS_LOG) &&
#   (grep -q "GATK ApplyBQSR finished: 0" $PROGRESS_LOG) &&
#   [ -s $WORKING_DIR/${PREFIX}_MQ30.txt ] &&
#   [ -s $WORKING_DIR/${PREFIX}_depth.txt ] &&
#   [ -s $WORKING_DIR/${PREFIX}_depth_of_coverage.txt ] &&
#   ([ -s $WORKING_DIR/${PREFIX}_BQSR.bam ] || [ -s $WORKING_DIR/${PREFIX}.g.vcf.gz ])); then
#  rm $WORKING_DIR/${PREFIX}_dup_marked.bam
#  rm $WORKING_DIR/${PREFIX}_dup_marked.bam.bai
#  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished"
#  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam removed when both Depth of coverage computing and GATK ApplyBQSR finished" >> $PROGRESS_LOG
#else
#  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam not removed when Depth of coverage computing finished, GATK ApplyBQSR probably not finished or Depth of coverage computing not successfully finished"
#  echo "$WORKING_DIR/${PREFIX}_dup_marked.bam not removed when Depth of coverage computing finished, GATK ApplyBQSR probably not finished or Depth of coverage computing not successfully finished" >> $PROGRESS_LOG
#fi

else
  echo "Running: Sambamba index not finished, stops the statistics gathering"
fi

end=`date +%s`
runtime=$((end-start))
echo "Finishing: time duration: ${runtime}s"
