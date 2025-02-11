#!/bin/bash
CONFIGURE=$1
OPCENTER=$(sed "1q;d" $CONFIGURE)
ROLE=$(sed "2q;d" $CONFIGURE)
CODE=$(sed "3q;d" $CONFIGURE)
ACCESSKEY=$(sed "4q;d" $CONFIGURE)
SECRETKEY=$(sed "5q;d" $CONFIGURE)
JFS=$(sed "6q;d" $CONFIGURE)
JFS_AS=$(sed "7q;d" $CONFIGURE)
IMAGE_NAME=$(sed "8q;d" $CONFIGURE)
INFLOW="data_inflow_efiga.txt"

mmfloat login \
  -a ${OPCENTER} \
  -u ${ROLE} \
  -p ${CODE}

for i in {1001..1800}
do
  echo "submitting job ${i}"
  sh mmfloat_stage_3_submit.sh $(sed "${i}q;d" ${INFLOW}) ${JFS} ${JFS_AS} ${ACCESSKEY} ${SECRETKEY} ${IMAGE_NAME}
done
