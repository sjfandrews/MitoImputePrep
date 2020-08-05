#!/bin/bash

STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt
STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms_smallTest0.txt
#STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms_imputedOnly.txt
#STRAND_LIST=~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms_imputedOnly_smaller.txt
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
	#summary_file=/g/data1a/te53/MitoImpute/data/STRANDS/${i}/${REFpanel}/MCMC${MCMC}/chrMT_1kg_${i}_imputed_MCMC${MCMC}_SUMMARY.csv
	summary_file=/g/data1a/te53/MitoImpute/data/STRANDS/${i}/${REFpanel}/MCMC${MCMC}/chrMT_1kg_${i}_imputed_MCMC${MCMC}_SUMMARY.tsv
	
	#echo ${summary_file}
	
	if [ ! -s ${summary_file} ]
	then
		qsub -v nqstat_anu=nqstat_anu,REFpanel=${REFpanel},MtPlatforms=${i},mcmc=${MCMC},burn=${BURN},khap=${KHAP},ne=${NE} -N IMPUTE2_${i}_${REFpanel}_MCMC${MCMC}_kHAP${KHAP} ~/GitCode/MitoImputePrep/scripts/BASH/DATA_PROCESSING/ThousandGenomes_Imputation/ThousandGenomes_imputeFrom_RefPan.sh
		echo "JOB FOR ${i} SUBMITTED FOR ${REFpanel} WITH PARAMETERS:	MCMC LENGTH = ${MCMC} (${BURN} BURN-INS),	K_HAP = ${KHAP},	NE = ${NE}"
	fi
done
