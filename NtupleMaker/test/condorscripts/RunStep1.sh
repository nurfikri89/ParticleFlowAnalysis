#!/bin/sh

SAMPLE_NAME=${1}
FILE_NUM=${2}
EOS_DIR=/eos/user/n/nbinnorj/RECOPFStudiesV2/

#
# Config file
#
CFG_TEMPLATE_FILE_NAME=${JOBWORKDIR}/Template_CMSRunConfig_Ntuples_${SAMPLE_NAME}.py
CFG_FILE_NAME=CMSRunConfig_Ntuples_${SAMPLE_NAME}_Run${FILE_NUM}.py
cp ${CFG_TEMPLATE_FILE_NAME} ./${CFG_FILE_NAME}

#
# Input file
#
OUTPUT_FILE_NAME=Ntuple_${SAMPLE_NAME}_Run${FILE_NUM}.root
TMP_OUTPUT_FILE_PATH=${TMPDIR}/${OUTPUT_FILE_NAME}

#
# Output file
#
DEST_OUTPUT_FILE_PATH=${EOS_DIR}/${OUTPUT_FILE_NAME}

#
#
#
NEVENTS=2500
# use 2500 events/lumi
declare -i LUMIBLOCK=${FILE_NUM}*2500/${NEVENTS}+1
declare -i FIRSTEVENT=${FILE_NUM}*${NEVENTS}
echo "Event" ${FIRSTEVENT} "lumi" ${LUMIBLOCK}

echo "process.maxEvents.input=cms.untracked.int32(${NEVENTS})" >> ${CFG_FILE_NAME}
echo "process.TFileService.fileName=cms.string('${TMP_OUTPUT_FILE_PATH}');" >> ${CFG_FILE_NAME}
echo "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper;randHelper=RandomNumberServiceHelper(process.RandomNumberGeneratorService);randHelper.populate()" >> ${CFG_FILE_NAME}
echo "process.source.firstLuminosityBlock=cms.untracked.uint32(${LUMIBLOCK});" >> ${CFG_FILE_NAME}
echo "process.source.firstEvent=cms.untracked.uint32(${FIRSTEVENT});" >> ${CFG_FILE_NAME}
echo "" >> ${CFG_FILE_NAME}
echo "" >> ${CFG_FILE_NAME}

echo "START"
date

cmsRun ${CFG_FILE_NAME}

#
# Copy config file to eos
#
cp -fv ${CFG_FILE_NAME} ${EOS_DIR}/

#
#
#
echo ""
echo "Check file size"
du -cskh ${TMP_OUTPUT_FILE_PATH}


#
# Transfer file to eos
#
echo ""
echo "Transfer from TMPDIR: ${TMP_OUTPUT_FILE_PATH}"
echo "To eos: ${DEST_OUTPUT_FILE_PATH}"
xrdcp -fv ${TMP_OUTPUT_FILE_PATH} root://eosuser.cern.ch/${DEST_OUTPUT_FILE_PATH}


#
# Delete file in TMPDIR
#
echo ""
echo "Delete file in TMPDIR"
rm -rvf ${TMP_OUTPUT_FILE_PATH}

date
echo "END"

