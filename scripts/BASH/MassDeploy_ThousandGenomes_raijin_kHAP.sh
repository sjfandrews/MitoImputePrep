#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt
MCMC_list=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/MCMC_list.txt
MCMC="1"
BURN="0"
KHAP="250"

#-k_hap 100, 250, 500, 1000, 2500, 5000, 10000, 20000, 30000

for i in `cat ${STRAND_LIST}`; do
	qsub -v MtPlatforms=${i},mcmc=${MCMC},burn=${BURN},khap=${KHAP} ~/GitCode/MitoImputePrep/scripts/BASH/ThousandGenomes_raijin_MCMCtest.sh
	echo "JOB FOR ${i} SUBMITTED WITH MCMC LENGTH = ${MCMC} (${BURN} BURN-INS) AND -k_hap = ${KHAP}" 
done
