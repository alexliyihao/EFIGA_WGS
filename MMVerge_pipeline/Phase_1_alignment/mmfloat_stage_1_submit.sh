#!/bin/bash
PREFIX=$1
JFS=$2
JFS_AS=$3
ACCESSKEY=$4
SECRETKEY=$5
IMAGE=$6

mmfloat submit \
  -n efiga_stage_1_alignment \
  -j phase_1_multi_thread.sh \
  -i ${IMAGE} \
  --tag 0.3 \
  -c 8:8 \
  -m 25:32 \
  --env PREFIX=${PREFIX} \
  --env THREADS=8 \
  --dataVolume ${JFS}:${JFS_AS} \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/${PREFIX}:/data/${PREFIX} \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/hg38_reference:/data/reference \
  --dataVolume [accessKey=${ACCESSKEY},secret=${SECRETKEY},endpoint=s3.us-east-1.amazonaws.com,mode=rw]s3://<your_bucket_here>/creds:/data/creds \
  --vmPolicy [spotOnly=true,retryLimit=20,retryInterval=900s]
