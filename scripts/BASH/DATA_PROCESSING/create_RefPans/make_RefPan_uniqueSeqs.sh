#!/bin/bash

# INFO:
#THIS SCRIPT CREATES VCFs AND OTHER FILE TYPES FOR A REFERENCE PANEL
#REFpanel=ReferencePanel_v2
#REFpanel=ReferencePanel_v1

# RUN COMMANDS

# SPECIFY HAPLOGREP VERSION
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.25/haplogrep-2.1.25.jar

# SPECIFY CURRENT MASTER ALIGNMENT

CURRENT="McInerney_Master_Alignment_July18_2018_duplicatesRemoved.fasta"

MASTER_ALN=$1 # INPUT THE CURRENT MASTER ALIGNMENT FOR THE REFERENCE PANEL
MASTER_ALN=`basename $MASTER_ALN`
echo ${MASTER_ALN}

if [ ${MASTER_ALN} == "McInerney_Master_Alignment_July18_2018_duplicatesRemoved.fasta" ]
then
	REFpanel=ReferencePanel_v1-unique
	echo
	echo "MASTER ALIGNMENT VERSION 1 [UNIQUE SEQUENCES ONLY] FOUND (18 JULY, 2018)"
	echo ${MASTER_ALN}
	echo
else
	REFpanel=ReferencePanel_vX-unique
	echo
	echo "MASTER ALIGNMENT NOT FOUND... PLEASE CHECK YOUR FILE"
	echo "DEFAULTING TO THE CURRENT REFERENCE PANEL:	${CURRENT}"
	#echo "KILLING THE SCRIPT HERE"
	echo
	#exit 1
fi

# SPECIFY THE MINOR ALLELE FREQUENCY
MAF_IN=$2
MAF_PC=`echo "${MAF_IN} * 100" | bc`
MAF_PC=`printf "%1.2f\n" "${MAF_PC}"`
echo
echo "MINOR ALLELE FREQUENCY: ${MAF_IN} ( ${MAF_PC}% )"
echo

WORKING_VERSION=${REFpanel}_${MAF_IN}
echo "WORKING VERSION:	${WORKING_VERSION}"

#CURRENT=McInerney_Master_Alignment_July18_2018.fasta
MT_DIR=/Volumes/TimMcInerney/MitoImpute/data/
ALN=${MT_DIR}FASTA/masters/${MASTER_ALN}
ALN_DIR=`dirname $ALN`/
ALN_BASE=`basename $ALN .fasta`

# CHECK FOR OR CREATE DIRECTORIES
echo

MASTER_VCF_DIR=${MT_DIR}VCF/Masters/`basename ${MASTER_ALN} .fasta`/
echo "SAVING MASTER ALIGNMENT VCF FILES TO:								${MASTER_VCF_DIR}"

REF_PAN_VCF_DIR=${MT_DIR}VCF/ReferencePanels/${REFpanel}/
REF_PAN_VCF_DIR_FILT=${REF_PAN_VCF_DIR}${WORKING_VERSION}/
echo "SAVING REFERNECE PANEL VCF FILES TO:								${REF_PAN_VCF_DIR} AND ${REF_PAN_VCF_DIR_FILT}"

REF_PAN_SAMP_DIR=/Volumes/TimMcInerney/MitoImpute/metadata/ReferencePanels/${REFpanel}/
REF_PAN_SAMP_DIR_FILT=${REF_PAN_SAMP_DIR}${WORKING_VERSION}/
echo "SAVING REFERNECE PANEL SAMPLE INFORMATION FILES TO:						${REF_PAN_SAMP_DIR_FILT}"

REF_PAN_OXFORD_DIR=${MT_DIR}OXFORD/ReferencePanels/${REFpanel}/
REF_PAN_OXFORD_DIR_FILT=${REF_PAN_OXFORD_DIR}${WORKING_VERSION}/
echo "SAVING REFERNECE PANEL OXFORD FILES (.hap.gz, .legend.gz, .samples, .gen.gz) TO:		${REF_PAN_OXFORD_DIR_FILT}"

REF_PAN_PLINK_DIR=${MT_DIR}PLINK/ReferencePanels/${REFpanel}/
REF_PAN_PLINK_DIR_FILT=${REF_PAN_PLINK_DIR}${WORKING_VERSION}/
echo "SAVING REFERNECE PANEL PLINK FILES (.ped, .map) TO:						${REF_PAN_PLINK_DIR_FILT}"

REF_PAN_RECOMB_DIR=${MT_DIR}REF_PANEL/ReferencePanels/${REFpanel}/
REF_PAN_RECOMB_DIR_FILT=${REF_PAN_RECOMB_DIR}${WORKING_VERSION}/
echo "SAVING REFERNECE PANEL RECOMBINATION MAP FILES TO:						${REF_PAN_RECOMB_DIR_FILT}"

if [ ! -d ${MASTER_VCF_DIR} ]
then
	echo "${MASTER_VCF_DIR} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${MASTER_VCF_DIR}
fi

if [ ! -d ${REF_PAN_VCF_DIR_FILT} ]
then
	echo "${REF_PAN_VCF_DIR} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${REF_PAN_VCF_DIR_FILT}
fi

if [ ! -d ${REF_PAN_SAMP_DIR_FILT} ]
then
	echo "${REF_PAN_SAMP_DIR_FILT} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${REF_PAN_SAMP_DIR_FILT}
fi

if [ ! -d ${REF_PAN_OXFORD_DIR_FILT} ]
then
	echo "${REF_PAN_OXFORD_DIR_FILT} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${REF_PAN_OXFORD_DIR_FILT}
fi

if [ ! -d ${REF_PAN_PLINK_DIR_FILT} ]
then
	echo "${REF_PAN_PLINK_DIR_FILT} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${REF_PAN_PLINK_DIR_FILT}
fi

if [ ! -d ${REF_PAN_RECOMB_DIR_FILT} ]
then
	echo "${REF_PAN_RECOMB_DIR_FILT} NOT FOUND	...	CREATING DIRECTORY"
	mkdir -p ${REF_PAN_RECOMB_DIR_FILT}
fi
exit
# CONVERT MISSING OR NON-N AMBIGUOUS CHARACTER STATES TO N AMBIGUOUS CHARACTER STATE
ALN_AMB=${MT_DIR}FASTA/ambiguous2missing/${ALN_BASE}"_ambig2missing.fasta"
ALN_AMB_GP=${MT_DIR}FASTA/ambiguous2missing/${ALN_BASE}"_ambigANDgap2missing.fasta"
if [ -f ${ALN_AMB} ]
then
	echo
	echo "${ALN_AMB} FOUND ... PASSING STEP"
else
	echo
	echo "${ALN_AMB} NOT FOUND ... CONVERTING MISSING OR NON-N AMBIGUOUS CHARACTER STATES TO N AMBIGUOUS CHARACTER STATE"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py -i ${ALN} -o ${ALN_AMB} -v
fi

if [ -f ${ALN_AMB_GP} ]
then
	echo
	echo "${ALN_AMB_GP} FOUND ... PASSING STEP"
else
	echo
	echo "${ALN_AMB_GP} NOT FOUND ... CONVERTING MISSING OR NON-N AMBIGUOUS AND GAP CHARACTER STATES TO N AMBIGUOUS CHARACTER STATE"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py -i ${ALN} -o ${ALN_AMB_GP} -g -v
fi



# ASSIGN HAPLOGROUPINGS VIA HAPLOGREP
ALN_AMB_HG=${MT_DIR}FASTA/ambiguous2missing/${ALN_BASE}"_ambig2missing_HaploGrep.txt"
ALN_AMB_GP_HG=${MT_DIR}FASTA/ambiguous2missing/${ALN_BASE}"_ambigANDgap2missing_HaploGrep.txt"

if [ -f ${ALN_AMB_HG} ]
then
	echo
	echo "${ALN_AMB_HG} FOUND ... PASSING STEP"
else
	echo
	echo "${ALN_AMB_HG} NOT FOUND ... ASSIGNING HAPLOGROUPS VIA HAPLOGREP"
	java -jar ${HAPLOGREP} --in ${ALN_AMB} --format fasta --extend-report --out ${ALN_AMB_HG} # assign haplogreps
fi

if [ -f ${ALN_AMB_GP_HG} ]
then
	echo
	echo "${ALN_AMB_GP_HG} FOUND ... PASSING STEP"
else
	echo
	echo "${ALN_AMB_GP_HG} NOT FOUND ... ASSIGNING HAPLOGROUPS VIA HAPLOGREP"
	java -jar ${HAPLOGREP} --in ${ALN_AMB_GP} --format fasta --extend-report --out ${ALN_AMB_GP_HG} # assign haplogreps
fi

# CONVERT THE CURRENT MASTER ALIGNMENT FASTA TO VCF FORMAT
VCF_AMB=${MASTER_VCF_DIR}`basename ${ALN_AMB} .fasta`.vcf.gz
VCF_AMB_GP=${MASTER_VCF_DIR}`basename ${ALN_AMB_GP} .fasta`.vcf.gz

if [ -f ${VCF_AMB} ]
then
	echo
	echo "${VCF_AMB} FOUND ... PASSING STEP"
else
	echo 
	echo "${VCF_AMB} NOT FOUND ... CONVERTING THE CURRENT MASTER ALIGNMENT FASTA TO VCF FORMAT"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${ALN_AMB} -o ${VCF_AMB} -v
fi

if [ -f ${VCF_AMB_GP} ]
then
	echo
	echo "${VCF_AMB_GP} FOUND ... PASSING STEP"
else
	echo 
	echo "${VCF_AMB_GP} NOT FOUND ... CONVERTING THE CURRENT MASTER ALIGNMENT FASTA TO VCF FORMAT"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${ALN_AMB_GP} -o ${VCF_AMB_GP} -v
fi

# RUN THE CURRENT MASTER ALIGNMENT VCF THROUGH BCFTOOLS TO MAKE SURE IT CONFORMS TO STANDARDS
REF_VCF_TMP=${REF_PAN_VCF_DIR}${REFpanel}_tmp.vcf.gz
REF_VCF=${REF_PAN_VCF_DIR}${REFpanel}.vcf.gz
if [ -f ${REF_VCF} ]
then
	echo
	echo "${REF_VCF} FOUND ... PASSING STEP"
else
	echo
	echo "${REF_VCF} NOT FOUND ... RUNNING THE CURRENT MASTER ALIGNMENT VCF THROUGH BCFTOOLS TO MAKE SURE IT CONFORMS TO STANDARDS"
	bcftools view -Oz -o ${REF_VCF_TMP} ${VCF_AMB}
	bcftools +fill-tags ${REF_VCF_TMP} -Oz -o ${REF_VCF}
	bcftools index ${REF_VCF}
	rm ${REF_VCF_TMP}
fi

## RUN HAPLOGREP ON THE REFERENCE PANEL VCF
#REF_VCF_HG=${REF_PAN_VCF_DIR}${REFpanel}_HaploGrep.txt
#if [ -f ${REF_VCF_HG} ]
#then
#	echo
#	echo "${REF_VCF_HG} FOUND ... PASSING STEP"
#else
#	echo
#	echo "${REF_VCF_HG} NOT FOUND ... ASSIGNING HAPLOGROUPS VIA HAPLOGREP"
#	java -jar ${HAPLOGREP} --in ${REF_VCF} --format vcf --extend-report --out ${REF_VCF_HG} # assign haplogreps
#fi

# REMOVE LOW QUALITY SEQUENCES
HQ_FILE=/Volumes/TimMcInerney/MitoImpute/metadata/`basename ${ALN_AMB} .fasta`"_highQual.txt"
VCF_HQ=${REF_PAN_VCF_DIR}`basename ${REF_VCF} .vcf.gz`"_highQual.vcf.gz"
if [ -f ${VCF_HQ} ]
then
	echo
	echo "${VCF_HQ} FOUND ... PASSING STEP"
else
	echo "${REF_VCF} NOT FOUND ... REMOVING LOW QUALITY SEQUENCES"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R ${ALN_AMB} ${HQ_FILE}
	bcftools view -S ${HQ_FILE} -Oz -o ${VCF_HQ} ${REF_VCF}
	bcftools index ${VCF_HQ}
fi

# APPLY SITE FILTRATION CRITERIA
VCF_FILT=${REF_PAN_VCF_DIR_FILT}`basename ${VCF_HQ} .vcf.gz`"_MAF${MAF_IN}_filtered.vcf.gz"
MAF_L=${MAF_IN}
MAF_U=`echo 1.0 - ${MAF_L} | bc`

echo
echo "LOWER MAF:	$MAF_L"
echo "UPPER MAF:	$MAF_U"
echo

if [ -f ${VCF_FILT} ]
then
	echo
	echo "${VCF_FILT} FOUND ... PASSING STEP"
else
	echo
	echo "${VCF_FILT} NOT FOUND ... APPLYING SITE FILTRATION CRITERIA"
	#vt decompose ${VCF_HQ} | bcftools +fill-tags | bcftools view -i 'ALT!="-"' | bcftools view -q 0.01 -Q 0.99 | bcftools view -Oz -o ${VCF_FILT}
	vt decompose ${VCF_HQ} | bcftools +fill-tags | bcftools view -i 'ALT!="-"' | bcftools view -q ${MAF_L} -Q ${MAF_U} | bcftools view -Oz -o ${VCF_FILT}
	bcftools index ${VCF_FILT}
fi

# EXTRACT SAMPLE LIST
SAMPS=${REF_PAN_SAMP_DIR_FILT}`basename ${VCF_FILT} .vcf.gz`"_sampleList.txt"
if [ -f ${SAMPS} ]
then
	echo
	echo "${SAMPS} FOUND ... PASSING STEP"
else
	echo
	echo "${SAMPS} NOT FOUND ... EXTRACTING SAMPLE LIST"
	bcftools query -l ${VCF_FILT} > ${SAMPS}
fi

# ADD SEX LABEL
SEX=${REF_PAN_SAMP_DIR_FILT}`basename ${SAMPS} .txt`"_sex.txt"
if [ -f ${SEX} ]
then
	echo
	echo "${SEX} FOUND ... PASSING STEP"
else
	echo
	echo "${SEX} NOT FOUND ... ADDING SEX LABEL"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${SAMPS} ${SEX}
fi

# CONVERT TO OXFORD FORMAT
OXF=${REF_PAN_OXFORD_DIR_FILT}${WORKING_VERSION}
if [ -f ${OXF}.hap.gz ] && [ -f ${OXF}.legend.gz ] && [ -f ${OXF}.samples ]
then
	echo
	echo "${OXF}.hap.gz AND ${OXF}.legend.gz AND ${OXF}.samples FOUND ... PASSING STEP"
else
	echo
	echo "${OXF}.hap.gz OR ${OXF}.legend.gz OR ${OXF}.samples NOT FOUND ... CONVERTING TO OXFORD HAP, LEGEND, AND SAMPLES FORMAT"
	bcftools1.4.1 convert --haplegendsample ${OXF} --sex ${SEX} ${VCF_FILT}
fi

# CONVERT TO PLINK FORMAT
PLK=${REF_PAN_PLINK_DIR_FILT}${WORKING_VERSION}
if [ -f ${PLK}.map ] && [ -f ${PLK}.ped ]
then
	echo
	echo "${OXF}.map AND ${OXF}.ped FOUND ... PASSING STEP"
else
	echo
	echo "${OXF}.map OR ${OXF}.ped NOT FOUND ... CONVERTING TO PLINK FORMAT"
	plink1.9 --vcf ${VCF_FILT} --recode --double-id --keep-allele-order --out ${PLK}
fi

# CONVERT TO GEN and SAMPLE FORMAT
GEN=${REF_PAN_OXFORD_DIR_FILT}${WORKING_VERSION}
if [ -f ${GEN}.gen.gz ]
then
	echo
	echo "${GEN}.gen.gz FOUND ... PASSING STEP"
else
	echo
	echo "${GEN}.gen.gz NOT FOUND ... CONVERTING TO GEN FORMAT"
	bcftools convert --gensample ${GEN} --sex ${SEX} ${VCF_FILT}
fi

# MAKE RECOMBINATION MAP AND STRAND FILES
MAP=${REF_PAN_RECOMB_DIR_FILT}${WORKING_VERSION}_MtMap.txt
STRAND=${REF_PAN_RECOMB_DIR_FILT}${WORKING_VERSION}_MtStrand.txt
if [ -f ${MAP} ] && [ -f ${STRAND} ]
then
	echo
	echo "${MAP} AND ${STRAND} FOUND ... PASSING STEP"
else
	echo
	echo "${MAP} AND ${STRAND} NOT FOUND ... MAKING RECOMBINATION MAP AND STRAND FILES"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/mt_recombination_map.R ${VCF_FILT} ${MAP} ${STRAND}
fi

# COPY RELEVANT FILES TO GitHub DIRECTORY
GIT_DIR=~/GitCode/MitoImputePrep/DerivedData/${WORKING_VERSION}/
if [ -d ${GIT_DIR} ]
then
	echo
	echo "${GIT_DIR} EXISTS ... PASSING"
else
	echo
	echo "${GIT_DIR} NOT FOUND ... CREATING DIRECTORY"
	mkdir ${GIT_DIR}
fi

cp ${GEN}.gen.gz ${GIT_DIR}
cp ${OXF}.hap.gz ${GIT_DIR}
cp ${OXF}.legend.gz ${GIT_DIR}
cp ${OXF}.samples ${GIT_DIR}
cp ${MAP} ${GIT_DIR}
cp ${STRAND} ${GIT_DIR}
cp ${VCF_FILT} ${GIT_DIR}
cp ${VCF_FILT}.csi ${GIT_DIR}
cp ${SEX} ${GIT_DIR}

# END