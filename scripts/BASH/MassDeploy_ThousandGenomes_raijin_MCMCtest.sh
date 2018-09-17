#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt
MCMC_list=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/MCMC_list.txt
MCMC="1"

for i in `cat ${STRAND_LIST}`; do
	qsub -v MtPlatforms=${i},mcmc=${MCMC} ~/GitCode/MitoImputePrep/scripts/BASH/ThousandGenomes_raijin.sh
	echo "JOB FOR ${i} SUBMITTED"
done
