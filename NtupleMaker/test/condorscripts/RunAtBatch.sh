#!/bin/sh
#
#
#
export X509_USER_PROXY=/afs/cern.ch/user/n/nbinnorj/myProxy
cd ${JOBWORKDIR}
#
#
#
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
#
#
#
SAMPLE_NAME=${1}
FILE_NUM=${2}
mkdir -p TMP_WORKDIR_${SAMPLE_NAME}_${FILE_NUM}
cd TMP_WORKDIR_${SAMPLE_NAME}_${FILE_NUM}
#
#
#
source ${JOBWORKDIR}/RunStep1.sh ${SAMPLE_NAME} ${FILE_NUM}

