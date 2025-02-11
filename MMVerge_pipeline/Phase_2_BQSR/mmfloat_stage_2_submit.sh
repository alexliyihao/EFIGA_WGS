#!/bin/bash
PREFIX=$1
JFS=$2
JFS_AS=$3
ACCESSKEY=$4
SECRETKEY=$5
IMAGE=$6

mmfloat submit \
  -n efiga_stage_2_genotyping_BQSR \
  -j phase_2_double_thread.sh \
  -i ${IMAGE} \
  --tag 0.3 \
  -c 2 \
  -m 8 \
  --env PREFIX=${PREFIX} \
  --env THREADS=2 \
  --dataVolume ${JFS}:${JFS_AS} \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/${PREFIX}:/data/${PREFIX} \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/hg38_reference:/data/reference \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/creds:/data/creds \
  --vmPolicy [spotOnly=true,retryLimit=20,retryInterval=900s] \
  --migratePolicy [cpu.upperBoundRatio=99,cpu.lowerBoundRatio=1,cpu.upperBoundDuration=600s,cpu.lowerBoundDuration=600s,cpu.limit=8,cpu.lowerLimit=2,mem.upperBoundRatio=80,mem.lowerBoundRatio=10,mem.upperBoundDuration=1s,mem.lowerBoundDuration=600s,mem.limit=64,mem.lowerLimit=8,stepAuto=true,evadeOOM=true]
