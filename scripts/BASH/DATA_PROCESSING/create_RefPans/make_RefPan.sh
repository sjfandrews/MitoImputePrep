#!/bin/bash

# INFO:
#THIS SCRIPT CREATES VCFs AND OTHER FILE TYPES FOR A REFERENCE PANEL
#REFpanel=ReferencePanel_v2
#REFpanel=ReferencePanel_v1

# RUN COMMANDS

# SPECIFY CURRENT MASTER ALIGNMENT
MASTER_ALN=$1 # INPUT THE CURRENT MASTER ALIGNMENT FOR THE REFERENCE PANEL
MASTER_ALN=`basename $MASTER_ALN`
echo ${MASTER_ALN}

if [ ${MASTER_ALN} == "McInerney_Master_Alignment_July18_2018.fasta" ]
then
	REFpanel=ReferencePanel_v1
	echo
	echo "MASTER ALIGNMENT VERSION 1 FOUND (18 JULY, 2018)"
	echo ${MASTER_ALN}
	echo
else
	REFpanel=ReferencePanel_vX
	echo
	echo "MASTER ALIGNMENT NOT FOUND... PLEASE CHECK YOUR FILE"
	echo "KILLING THE SCRIPT HERE"
	echo
	exit 1
fi

# SPECIFY THE MINOR ALLELE FREQUENCY
MAF_IN=$2
MAF_PC=`echo "${MAF_IN} * 100" | bc`
MAF_PC=`printf "%8.2f\n" "${MAF_PC}"`
echo
echo "MINOR ALLELE FREQUENCY: ${MAF_IN} (${MAF_PC}% )"
echo

WORKING_VERSION=${REFpanel}_${MAF_IN}
echo ${WORKING_VERSION}

#CURRENT=McInerney_Master_Alignment_July18_2018.fasta
MT_DIR=/Volumes/TimMcInerney/MitoImpute/data/
ALN=${MT_DIR}FASTA/masters/${CURRENT}
ALN_DIR=`dirname $ALN`/
ALN_BASE=`basename $ALN .fasta`

# CONVERT MISSING OR NON-N AMBIGUOUS CHARACTER STATES TO N AMBIGUOUS CHARACTER STATE
ALN_AMB=${MT_DIR}FASTA/ambiguous2missing/${ALN_BASE}"_ambig2missing.fasta"
if [ -f ${ALN_AMB} ]
then
	echo
	echo "${ALN_AMB} FOUND ... PASSING STEP"
else
	echo echo "${ALN_AMB} NOT FOUND ... CONVERTING MISSING OR NON-N AMBIGUOUS CHARACTER STATES TO N AMBIGUOUS CHARACTER STATE"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py -i ${ALN} -o ${ALN_AMB} -v
fi

## CONVERT THE CURRENT MASTER ALIGNMENT FASTA TO VCF FORMAT
#VCF_AMB=${MT_DIR}VCF/`basename ${ALN_AMB} .fasta`.vcf.gz
#if [ -f ${VCF_AMB} ]
#then
#	echo
#	echo "${VCF_AMB} FOUND ... PASSING STEP"
#else
#	echo 
#	echo "${VCF_AMB} NOT FOUND ... CONVERTING THE CURRENT MASTER ALIGNMENT FASTA TO VCF FORMAT"
#	python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${ALN_AMB} -o ${VCF_AMB} -v
#fi
#
## RUN THE CURRENT MASTER ALIGNMENT VCF THROUGH BCFTOOLS TO MAKE SURE IT CONFORMS TO STANDARDS
#REF_VCF=${MT_DIR}VCF/${REFpanel}.vcf.gz
#if [ -f ${REF_VCF} ]
#then
#	echo
#	echo "${REF_VCF} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${REF_VCF} NOT FOUND ... RUNNING THE CURRENT MASTER ALIGNMENT VCF THROUGH BCFTOOLS TO MAKE SURE IT CONFORMS TO STANDARDS"
#	bcftools view -Oz -o ${REF_VCF} ${VCF_AMB}
#	bcftools index ${REF_VCF}
#fi
#
## REMOVE LOW QUALITY SEQUENCES
#HQ_FILE=/Volumes/TimMcInerney/MitoImpute/metadata/`basename ${ALN_AMB} .fasta`"_highQual.txt"
#VCF_HQ=${MT_DIR}VCF/`basename ${REF_VCF} .vcf.gz`"_highQual.vcf.gz"
#if [ -f ${VCF_HQ} ]
#then
#	echo
#	echo "${VCF_HQ} FOUND ... PASSING STEP"
#else
#	echo "${REF_VCF} NOT FOUND ... REMOVING LOW QUALITY SEQUENCES"
#	Rscript ~/GitCode/MitoImputePrep/scripts/R/removeLowQuality_cmdline.R ${ALN_AMB} ${HQ_FILE}
#	bcftools view -S ${HQ_FILE} -Oz -o ${VCF_HQ} ${REF_VCF}
#	bcftools index ${VCF_HQ}
#fi
#
## APPLY SITE FILTRATION CRITERIA
#VCF_FILT=${MT_DIR}VCF/`basename ${VCF_HQ} .vcf.gz`"_MAF${MAF_IN}_filtered.vcf.gz"
#MAF_L=${MAF_IN}
#MAF_U=`echo 1.0 - ${MAF_L} | bc`
#if [ -f ${VCF_FILT} ]
#then
#	echo
#	echo "${VCF_FILT} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${VCF_FILT} NOT FOUND ... APPLYING SITE FILTRATION CRITERIA"
#	#vt decompose ${VCF_HQ} | bcftools +fill-tags | bcftools view -i 'ALT!="-"' | bcftools view -q 0.01 -Q 0.99 | bcftools view -Oz -o ${VCF_FILT}
#	vt decompose ${VCF_HQ} | bcftools +fill-tags | bcftools view -i 'ALT!="-"' | bcftools view -q ${MAF_L} -Q ${MAF_U} | bcftools view -Oz -o ${VCF_FILT}
#	bcftools index ${VCF_FILT}
#fi
#
## EXTRACT SAMPLE LIST
#SAMPS=/Volumes/TimMcInerney/MitoImpute/metadata/`basename ${VCF_FILT} .vcf.gz`"_sampleList.txt"
#if [ -f ${SAMPS} ]
#then
#	echo
#	echo "${SAMPS} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${SAMPS} NOT FOUND ... EXTRACTING SAMPLE LIST"
#	bcftools query -l ${VCF_FILT} > ${SAMPS}
#fi
#
## ADD SEX LABEL
#SEX=/Volumes/TimMcInerney/MitoImpute/metadata/`basename ${SAMPS} .txt`"_sex.txt"
#if [ -f ${SEX} ]
#then
#	echo
#	echo "${SEX} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${SEX} NOT FOUND ... ADDING SEX LABEL"
#	Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${SAMPS} ${SEX}
#fi
#
## CONVERT TO OXFORD FORMAT
#OXF=/Volumes/TimMcInerney/MitoImpute/data/OXFORD/${WORKING_VERSION}
#if [ -f ${OXF}.hap.gz ] && [ -f ${OXF}.legend.gz ] && [ -f ${OXF}.samples ]
#then
#	echo
#	echo "${OXF}.hap.gz AND ${OXF}.legend.gz AND ${OXF}.samples FOUND ... PASSING STEP"
#else
#	echo
#	echo "${OXF}.hap.gz OR ${OXF}.legend.gz OR ${OXF}.samples NOT FOUND ... CONVERTING TO OXFORD HAP, LEGEND, AND SAMPLES FORMAT"
#	bcftools1.4.1 convert --haplegendsample ${OXF} --sex ${SEX} ${VCF_FILT}
#fi
#
## CONVERT TO PLINK FORMAT
#PLK=/Volumes/TimMcInerney/MitoImpute/data/PLINK/${WORKING_VERSION}
#if [ -f ${PLK}.map ] && [ -f ${PLK}.ped ]
#then
#	echo
#	echo "${OXF}.map AND ${OXF}.ped FOUND ... PASSING STEP"
#else
#	echo
#	echo "${OXF}.map OR ${OXF}.ped NOT FOUND ... CONVERTING TO PLINK FORMAT"
#	plink1.9 --vcf ${VCF_FILT} --recode --double-id --keep-allele-order --out ${PLK}
#fi
#
## CONVERT TO GEN and SAMPLE FORMAT
#GEN=/Volumes/TimMcInerney/MitoImpute/data/OXFORD/${WORKING_VERSION}
#if [ -f ${GEN}.gen.gz ]
#then
#	echo
#	echo "${GEN}.gen.gz FOUND ... PASSING STEP"
#else
#	echo
#	echo "${GEN}.gen.gz NOT FOUND ... CONVERTING TO GEN FORMAT"
#	bcftools convert --gensample ${GEN} --sex ${SEX} ${VCF_FILT}
#fi
#
## MAKE RECOMBINATION MAP AND STRAND FILES
#MAP=/Volumes/TimMcInerney/MitoImpute/data/REF_PANEL/${WORKING_VERSION}_MtMap.txt
#STRAND=/Volumes/TimMcInerney/MitoImpute/data/REF_PANEL/${WORKING_VERSION}_MtStrand.txt
#if [ -f ${MAP} ] && [ -f ${STRAND} ]
#then
#	echo
#	echo "${MAP} AND ${STRAND} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${MAP} AND ${STRAND} NOT FOUND ... MAKING RECOMBINATION MAP AND STRAND FILES"
#	Rscript ~/GitCode/MitoImputePrep/scripts/R/mt_recombination_map.R ${VCF_FILT} ${MAP} ${STRAND}
#fi
#
## COPY RELEVANT FILES TO GitHub DIRECTORY
#GIT_DIR=~/GitCode/MitoImputePrep/DerivedData/${WORKING_VERSION}/
#if [ -d ${GIT_DIR} ]
#then
#	echo
#	echo "${GIT_DIR} EXISTS ... PASSING"
#else
#	echo
#	echo "${GIT_DIR} NOT FOUND ... CREATING DIRECTORY"
#	mkdir ${GIT_DIR}
#fi
#
#cp ${GEN}.gen.gz ${GIT_DIR}
#cp ${OXF}.hap.gz ${GIT_DIR}
#cp ${OXF}.legend.gz ${GIT_DIR}
#cp ${OXF}.samples ${GIT_DIR}
#cp ${MAP} ${GIT_DIR}
#cp ${STRAND} ${GIT_DIR}
#cp ${VCF_FILT} ${GIT_DIR}
#cp ${VCF_FILT}.csi ${GIT_DIR}
#cp ${SEX} ${GIT_DIR}
#
## END