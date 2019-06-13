#!/bin/bash

# ADNI RESEQUENCED VCF:		${WK_DIR}MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# ADNI GENOTYPED VCF:		${WK_DIR}MitoImpute/data/ADNI/Timpute/ADNI/mito_snps_rcrs_ed.vcf

# ADNI FILE (BELOW) HAD TO BE MANUALLY FIXED TO CONTAIN THE LINE "##contig=<ID=26,length=16569>"
# ${WK_DIR}MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# RENAMED TO: ${WK_DIR}MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214_fixed.vcf

# MAKE SURE DIRECTORIES TO M DRIVE WORK
if [ -d "/Volumes/TimMcInerney/" ]
then
	WK_DIR="/Volumes/TimMcInerney/"
else
	WK_DIR="/Volumes/MHS/"
fi

# SPECIFY REFERENCE PANEL
#REFpanel="ReferencePanel_v6"
REFpanel=$1
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
mcmc=1
burn=0

# SET THE ORIGINAL VCF FILES
TYP

echo
echo "REFERENCE PANEL:	${REFpanel}"
echo "HAPLOGREP:		${HAPLOGREP}"
echo "MCMC LENGTH:		${mcmc}"
echo "MCMC BURN-IN:		${burn}"
echo
echo "WGS VCF FILE:		${WGS_VCF}"
echo "GENOTYPED VCF:	${TYP_VCF}"
echo "IMPUTED VCF:		${IMP_VCF}"
echo

# SPECIFY REFERENCE FASTA FILES
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 

# ASSIGN THE PREFIX
ADNI_PREFIX=ADNI_12GO.Merged.chrM

# MAKE A VCF FOR THE PURPOSES OF SAMPLE ID EXTRACTION
echo
echo "MAKING A VCF FOR THE PURPOSES OF SAMPLE ID EXTRACTION"

BFILE=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/PLINK/ORIGINALS/${ADNI_PREFIX}
OUT_VCF=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/VCF/${ADNI_PREFIX}
ADNI_12GO_VCF=${OUT_VCF}.vcf.gz

plink1.9 --bfile ${BFILE} --recode vcf --out ${OUT_VCF}

bcftools norm --check-ref s -f ${REF26} -m + ${OUT_VCF}.vcf -Oz -o ${ADNI_12GO_VCF}
bcftools index ${ADNI_12GO_VCF}
bcftools index -t ${ADNI_12GO_VCF}


# EXTRACT SAMPLE IDS 
echo 
echo "EXTRACT SAMPLE IDS"

TYP_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/INFO/${ADNI_PREFIX}_samples.txt
TYP_SEX_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/INFO/${ADNI_PREFIX}_samples_sex.txt


bcftools query -l ${ADNI_12GO_VCF} > ${TYP_SAMPLE}
Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${TYP_SAMPLE} ${TYP_SEX_SAMPLE}



# GENERATE PLINK FILES (PED AND MAP FILES!)
echo
echo "GENERATE PLINK FILES (PED AND MAP FILES!)"

TYP_PLINK=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/PLINK/${ADNI_PREFIX}

plink1.9 --bfile ${BFILE} --recode --double-id --keep-allele-order --out ${TYP_PLINK}

# GENERATE GEN FILES
echo 
echo "GENERATE GEN AND SAMPLE FILES"
TYP_GEN_OUT=${WK_DIR}MitoImpute/data/ADNI_12GO/GENOTYPED/HAP_LEGEND_GEN/${ADNI_PREFIX}

bcftools convert --gensample ${TYP_GEN_OUT} ${ADNI_12GO_VCF} --sex ${TYP_SEX_SAMPLE}

# FIX SEX ID COLUMN IN .samples FILE
echo
echo "FIXING SEX ID COLUMN IN .samples FILE"
Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/FixSamplesFile_raijin.R ${TYP_GEN_OUT}.samples

# RUN IMPUTE2
map=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt
hap=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
leg=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
gen=${TYP_GEN_OUT}.gen.gz
sam=${TYP_GEN_OUT}.samples
out_dir=${WK_DIR}MitoImpute/data/ADNI_12GO/IMPUTED/${REFpanel}/IMPUTE2/
impute2_out=${out_dir}MitoImpute_${REFpanel}_imputed

if [ ! -d ${out_dir} ]
then
	mkdir -p ${out_dir}
fi

if [ -f ${impute2_out} ]
then
	echo "${impute2_out} FOUND! ... PASSING"
else
	echo "${impute2_out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${map} -h ${hap} -l ${leg} -g ${gen} -sample_g ${sam} -int 1 16569 -Ne 20000 -o ${impute2_out} -iter ${mcmc} -burnin ${burn}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
#impute2_out=${WK_DIR}MitoImpute/data/ADNI_12GO/IMPUTED/IMPUTE2/${REFpanel}/MitoImpute_ADNI_imputed
impute2_out_fixed=${impute2_out}_ChromFixed
awk '{{$1 = "26"; print}}' ${impute2_out} > ${impute2_out_fixed}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
impute2_gen=${impute2_out}_ChromFixed
impute2_sam=${impute2_out}_samples
#impute2_out=${impute2_out}

plink1.9 --gen ${impute2_gen} --sample ${impute2_sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${impute2_out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
impute2_gen=${impute2_out}_ChromFixed
impute2_sam=${impute2_out}_samples
#impute2_out=${impute2_out}

plink1.9 --gen ${impute2_gen} --sample ${impute2_sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${impute2_out}

# MAKE SURE IMPUTED VCF HAS REFERENCE ALLELES THE SAME AS rCRS AND CONVERT TO MULTI-ALLELIC FOR THIS PURPOSE
echo
echo "MAKING SURE IMPUTED VCF HAS REFERENCE ALLELES THE SAME AS rCRS, AND CONVERT TO MULTI-ALLELIC FOR THIS PURPOSE"
IMP_VCF_rCRS=${impute2_out}_rCRS.vcf
IMP_VCF=${impute2_out}_FINAL.vcf.gz
bcftools norm --check-ref s -f ${REF26} -m + ${impute2_out}.vcf -Oz -o ${IMP_VCF_rCRS}

# DECOMPOSE
vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
bcftools index ${IMP_VCF}

# PERFORM HAPLOGREP2 HAPLOGROUP ASSIGNMENT
HG2_FILE=${impute2_out}_FINAL_HaploGrep.txt
if [ -f ${HG2_FILE} ]
then
	echo
	echo "${HG2_FILE} FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${HG2_FILE} NOT FOUND ... RECODING TO plink VCF FILE"
	java -jar ${HAPLOGREP} --in ${impute2_out}_FINAL.vcf.gz --format vcf --chip --out ${HG2_FILE} # assign haplogreps
fi

# PERFORM Hi-MC HAPLOGROUP ASSIGNMENT
HiMC_FILE=${impute2_out}_FINAL_HiMC.csv
PED_FILE=${impute2_out}.ped
MAP_FILE=${impute2_out}.map
if [ -f ${HiMC_FILE} ]
then
	echo
	echo "${HiMC_FILE} FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${HiMC_FILE} NOT FOUND ... RECODING TO plink VCF FILE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/HiMC/HiMC_haplogroup_assignment.R ${PED_FILE} ${MAP_FILE} ${HiMC_FILE}  # assign haplogreps
fi

