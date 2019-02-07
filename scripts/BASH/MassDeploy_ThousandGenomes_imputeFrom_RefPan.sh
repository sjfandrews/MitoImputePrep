#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt
MCMC_list=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/MCMC_list.txt
MCMC="1"
BURN="0"
KHAP="500"
NE="20000"

REFpanel=$1

echo
echo "USAGE:"
echo "sh ~/GitCode/MitoImputePrep/scripts/BASH/MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh <REF_PANEL>"

#MCMC 1,0	5,1	10,3	20,6	30,10

for i in `cat ${STRAND_LIST}`; do
	qsub -v REFpanel=${REFpanel},MtPlatforms=${i},mcmc=${MCMC},burn=${BURN} ~/GitCode/MitoImputePrep/scripts/BASH/ThousandGenomes_imputeFrom_RefPan.sh
	echo "JOB FOR ${i} SUBMITTED FOR ${REFpanel} WITH PARAMETERS:	MCMC LENGTH = ${MCMC} (${BURN} BURN-INS),	K_HAP = ${KHAP},	NE = ${NE}"
done
