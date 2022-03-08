#!/bin/bash
#PBS -P gw26
#PBS -q biodev
#PBS -l walltime=24:00:00
#PBS -l mem=64GB
#PBS -l ncpus=1
#PBS -m e
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/
#PBS -l storage=scratch/te53+gdata/te53

		
TAR_OUT=/g/data/te53/MitoImpute.tar.gz
DIR_TO_TAR=/g/data/te53/MitoImpute/
OUTPUT_DIR=/g/data/te53/

ls -lh ${TAR_OUT}
ls -lh ${OUTPUT_DIR}

echo
echo "DIRECTORY TO COMPRESS:	${DIR_TO_TAR}"	
echo "OUTPUT TAR FILE:			${TAR_OUT}"

if [ ! -s ${TAR_OUT} ] && [ -d ${DIR_TO_TAR} ]
then
	echo "COMPRESSING DIRECTORIES FOR ${CHR} OF ${POP}"

	tar -zcvf ${TAR_OUT} -C ${OUTPUT_DIR} ${POP} && rm -rf ${DIR_TO_TAR}
	#echo tar -cvf ~/Desktop/SANDBOX/test_tar.tar -C ~/Desktop/SANDBOX/ TreeTime_testTar && rm -R ~/Desktop/SANDBOX/TreeTime_testTar/
elif [ -s ${TAR_OUT} ] && [ -d ${DIR_TO_TAR} ]
then
	echo "WARNING: COMPRESSED DIRECTORY AND UNCOMPRESSED STILL EXISTING"
	echo "CONTINUE COMPRESSING DIRECTORIES FOR ${CHR} OF ${POP}"
	tar -zcvf ${TAR_OUT} -C ${OUTPUT_DIR} ${POP} && rm -rf ${DIR_TO_TAR}
else
	echo "COMPRESSED DIRECTORY ALREADY FOUND FOR ${CHR} OF ${POP}	...	PASSING"
	echo "FOUND AT:	${TAR_OUT}"
fi








# END!
