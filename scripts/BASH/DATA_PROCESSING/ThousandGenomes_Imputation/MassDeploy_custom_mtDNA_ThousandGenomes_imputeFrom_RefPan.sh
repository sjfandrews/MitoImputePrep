#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/custom_mtDNA_platforms.txt
MCMC_list=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/MCMC_list.txt

REFpanel=$1
MCMC="1"
BURN="0"
KHAP="500"
NE="20000"

echo
echo "USAGE:"
echo "sh ~/GitCode/MitoImputePrep/scripts/BASH/MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh <REF_PANEL>"

#MCMC 1,0	5,1	10,3	20,6	30,10

for i in `cat ${STRAND_LIST}`; do
	qsub -v REFpanel=${REFpanel},MtPlatforms=${i},mcmc=${MCMC},burn=${BURN},khap=${KHAP},ne=${NE} -N ${i}_impute_SNPchip_1kGP ~/GitCode/MitoImputePrep/scripts/BASH/DATA_PROCESSING/ThousandGenomes_Imputation/ThousandGenomes_imputeFrom_RefPan.sh
	echo "JOB FOR ${i} SUBMITTED FOR ${REFpanel} WITH PARAMETERS:	MCMC LENGTH = ${MCMC} (${BURN} BURN-INS),	K_HAP = ${KHAP},	NE = ${NE}"
done
