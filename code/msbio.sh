#!/bin/bash
conda activate singularity

cd /media/ssd/sdb1/data/ljh/software/msbio_v3.5.20250101/

# singularity build --force msbio2.sif docker://docker.io/metadocker8/msbio2:latest
# singularity build --force msbio2.sif docker://dockerproxy.net/metadocker8/msbio2:latest
# singularity build --force msbio2.sif docker://hub.c.163.com/library/metadocker8/msbio2:latest
# singularity build --force msbio2.sif docker://hub.c.163.com/library/msbio2:latest

if [ ! -d "data" ]; then
  mkdir data
fi
chmod a+rwx data

bin/sms.sh -u -o /data/output_single_id_txt /data/example/single_list_id.txt
