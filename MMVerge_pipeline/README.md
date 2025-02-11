
# Memverge EFIGA deployment

This folder includes the codes deploying the pipeline on Memverge.

## Cloud Computing Infrastructures

 - Cloud Services:
  See *EFIGA_Memverge_juiceFS_approach.png*, we used a JuiceFS file system as a supplement for S3 bucket, for AWS S3 mount to EC2 is unable to handle >95GBs continuous ingress (writing from EC2 & ingress to S3) mount simultaneously. An AWS EBS is fine as well but it's a bit more expensive.

 - S3 bucket:
  See *bucket_file_structure.txt*

## Files

### Pipeline Scripts
For each phase,  we have one folder with three files:
  - *phase_\<number\>_**.sh*: the actual script running
  - *mmfloat_stage_\<number\>_submit.sh*: the script submit running script to mmfloat
  - *mmfloat_stage_\<number\>_batch_submission.sh*: a wrapper submit a batch of jobs from submission script to a specific OpCenter

### Supporting Files:
  All formats are provided, confidential info are neutralized
  - *opcenter_info.txt*: some metadata used in job submission. The submission should run as follows:
  **sh ./mmfloat_stage_\<number\>_batch_submission.sh opcenter_info.txt**
  I didn't left the start and end outside for I use that as my own notes, you can do it yourself as variable $2 and $3
  - credentials: see *bucket_file_structure.txt*, for the pipeline' login to AWS CLI downloading input fastqs and upload final crams and gvcfs
  - data_inflow_efiga.txt: see *mmfloat_stage_\<number\>_batch_submission.sh*, providing an overall orders for file uploading and jobs submissions
