#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt
MCMC_list=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/MCMC_list.txt
MCMC="1"
BURN="0"

REFpanel=$1

#MCMC 1,0	5,1	10,3	20,6	30,10

for i in `cat ${STRAND_LIST}`; do
	qsub -v REFpanel=${REFpanel},MtPlatforms=${i},mcmc=${MCMC},burn=${BURN} ~/GitCode/MitoImputePrep/scripts/BASH/ThousandGenomes_imputeFrom_RefPan_v6.sh
	echo "JOB FOR ${i} SUBMITTED WITH MCMC LENGTH = ${MCMC} AND ${BURN} BURN-INS"
done
