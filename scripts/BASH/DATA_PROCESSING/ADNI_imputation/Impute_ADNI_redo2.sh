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
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.25/haplogrep-2.1.25.jar
mcmc=1
burn=0

# SET THE ORIGINAL VCF FILES
WGS_VCF=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed.vcf.gz # Whole genome resequence (ADNI 3)
TYP_VCF=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed.vcf.gz # Genotyped data (ADNI 1)
#IMP_VCF=${WK_DIR}

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


# SUBSET EACH VCF FILE SO THEY ONLY CONTAIN THE n=258 FOUND IN BOTH ADNI 1 and ADNI 3
echo
echo "SUBSET EACH VCF FILE SO THEY ONLY CONTAIN THE n=258 FOUND IN BOTH ADNI 1 and ADNI 3"
TYP_n258=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258.vcf.gz
WGS_n258=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258.vcf.gz

if [ ! -s ${TYP_n258} ]
then
	echo
	echo "${TYP_n258} NOT FOUND	...	RUNNING"
	bcftools view -S ~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.txt ${TYP_VCF} | bcftools +fill-tags -Oz -o ${TYP_n258}
	bcftools index ${TYP_n258}
fi

if [ ! -s ${WGS_n258} ]
then
	echo
	echo "${WGS_n258} NOT FOUND	...	RUNNING"
	bcftools view -S ~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH_reseq.txt ${WGS_VCF} | bcftools +fill-tags -Oz -o ${WGS_n258}
	bcftools index ${WGS_n258}
fi

# FILTER VCF TO CONTAIN ONLY SNVs 
echo
echo "FILTER VCF TO CONTAIN ONLY SNVs"
TYP_n258_tmp=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_tmp.vcf.gz
WGS_n258_tmp=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_tmp.vcf.gz
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 

#bcftools annotate -x INFO,^FORMAT/GT ${WGS_n258} | bcftools norm --check-ref s -f ${REF26} -m - | bcftools view -V indels,mnps | bcftools +fill-tags -Oz -o ${WGS_n258_biallelic}
#bcftools annotate -x INFO,^FORMAT/GT ${TYP_n258} | bcftools norm --check-ref s -f ${REFMT} -m - | bcftools view -V indels,mnps | bcftools +fill-tags -Oz -o ${TYP_n258_biallelic}

if [ ! -s ${TYP_n258_tmp} ] && [ ! -s ${WGS_n258_tmp} ]
then
	echo
	echo "${TYP_n258_tmp} OR ${WGS_n258_tmp} NOT FOUND	...	RUNNING"
	bcftools annotate -x INFO,^FORMAT/GT ${WGS_n258} | bcftools norm --check-ref s -f ${REF26} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${WGS_n258_tmp}
	bcftools annotate -x INFO,^FORMAT/GT ${TYP_n258} | bcftools norm --check-ref s -f ${REFMT} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${TYP_n258_tmp}

	bcftools index ${WGS_n258_tmp}
	bcftools index ${TYP_n258_tmp}
fi

# DECOMPOSE VCF FILES AND SPLIT MULTIALLELIC RECORDS INTO BIALLELIC 
echo
echo "DECOMPOSE VCF FILES AND SPLIT MULTIALLELIC RECORDS INTO BIALLELIC"

TYP_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_biallelic.vcf.gz
WGS_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic.vcf.gz


if [ ! -s ${TYP_n258_biallelic} ] && [ ! -s ${WGS_n258_biallelic} ]
then
	echo
	echo "${TYP_n258_biallelic} OR ${WGS_n258_biallelic} NOT FOUND	...	RUNNING"
	vt decompose ${WGS_n258_tmp} | bcftools +fill-tags -Oz -o ${WGS_n258_biallelic}
	vt decompose ${TYP_n258_tmp} | bcftools +fill-tags -Oz -o ${TYP_n258_biallelic}
	
	bcftools index ${WGS_n258_biallelic}
	bcftools index ${TYP_n258_biallelic}
fi


# RELABEL WGS TO ADNI1 SAMPLE IDs
samples_csv=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv
WGS_relab=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled.vcf

if [ ! -s ${WGS_relab}.gz ]
then
	echo
	echo "${WGS_relab}.gz NOT FOUND	...	RUNNING"
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/fix_vcf_names.py -i ${TYP_n258_biallelic} -o ${WGS_relab} -c ${samples_csv} -v
	# TIM YOU HAVE TO MAKE THIS COMMAND LINE >:(

	# GZIP USING BCFTOOLS TO MAKE SURE VCF CONFORMS TO STANDARDS
	bcftools view ${WGS_relab} -Oz -o ${WGS_relab}.gz
fi


# EXTRACT SAMPLE IDS 
echo 
echo "EXTRACT SAMPLE IDS"

WGS_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/INFO/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled_sampleID.txt
WGS_SEX_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/INFO/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled_sampleID_sex.txt
TYP_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/INFO/mito_snps_rcrs_ed_n258_biallelic_sampleID.txt
TYP_SEX_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/INFO/mito_snps_rcrs_ed_n258_biallelic_sampleID_sex.txt

if [ ! -s ${TYP_SEX_SAMPLE} ] && [ ! -s ${WGS_SEX_SAMPLE} ]
then
	echo
	echo "${TYP_SEX_SAMPLE} OR ${WGS_SEX_SAMPLE} NOT FOUND	...	RUNNING"
	bcftools query -l ${WGS_relab} > ${WGS_SAMPLE}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${WGS_SAMPLE} ${WGS_SEX_SAMPLE}

	bcftools query -l ${TYP_n258_biallelic} > ${TYP_SAMPLE}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${TYP_SAMPLE} ${TYP_SEX_SAMPLE}
fi

# GENERATE PLINK FILES (PED AND MAP FILES!)
echo
echo "GENERATE PLINK FILES (PED AND MAP FILES!)"
WGS_PLINK=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/PLINK/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled
TYP_PLINK=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_biallelic

if [ ! -s ${TYP_PLINK}.ped ] && [ ! -s ${WGS_PLINK}.ped ]
then
	echo
	echo "${TYP_PLINK}.ped OR ${WGS_PLINK}.ped NOT FOUND	...	RUNNING"
	plink1.9 --vcf ${WGS_relab} --recode --double-id --keep-allele-order --out ${WGS_PLINK}
	plink1.9 --vcf ${TYP_n258_biallelic} --recode --double-id --keep-allele-order --out ${TYP_PLINK}
fi


# GENERATE GEN FILES
echo 
echo "GENERATE GEN AND SAMPLE FILES"
WGS_GEN_OUT=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/HAP_LEGEND_GEN/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled
TYP_GEN_OUT=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/HAP_LEGEND_GEN/mito_snps_rcrs_ed_n258_biallelic

if [ ! -s ${TYP_GEN_OUT}.gen.gz ] && [ ! -s ${WGS_GEN_OUT}.gen.gz ]
then
	echo
	echo "${TYP_SEX_SAMPLE}.gen.gz OR ${WGS_SEX_SAMPLE}.gen.gz NOT FOUND	...	RUNNING"
	bcftools convert --gensample ${WGS_GEN_OUT} ${WGS_relab} --sex ${WGS_SEX_SAMPLE}
	bcftools convert --gensample ${TYP_GEN_OUT} ${TYP_n258_biallelic} --sex ${TYP_SEX_SAMPLE}
fi


# FIX SEX ID COLUMN IN .samples FILE
echo
echo "FIXING SEX ID COLUMN IN .samples FILE"
Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/FixSamplesFile_raijin.R ${WGS_GEN_OUT}.samples
Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/FixSamplesFile_raijin.R ${TYP_GEN_OUT}.samples

# RUN IMPUTE2
map=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt
hap=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
leg=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
gen=${TYP_GEN_OUT}.gen.gz
sam=${TYP_GEN_OUT}.samples
out_dir=${WK_DIR}MitoImpute/data/ADNI_REDO/IMPUTED/${REFpanel}/IMPUTE2/
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
#impute2_out=${WK_DIR}MitoImpute/data/ADNI_REDO/IMPUTED/IMPUTE2/${REFpanel}/MitoImpute_ADNI_imputed
impute2_out_fixed=${impute2_out}_ChromFixed
if [ ! -s ${impute2_out_fixed} ]
then
	echo
	echo "${impute2_out_fixed} NOT FOUND	...	RUNNING"
	awk '{{$1 = "26"; print}}' ${impute2_out} > ${impute2_out_fixed}
fi

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
impute2_gen=${impute2_out}_ChromFixed
impute2_sam=${impute2_out}_samples
#impute2_out=${impute2_out}

if [ ! -s ${impute2_out}.ped ]
then
	echo
	echo "${impute2_out}.ped NOT FOUND	...	RUNNING"
	plink1.9 --gen ${impute2_gen} --sample ${impute2_sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${impute2_out}
fi


# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
impute2_gen=${impute2_out}_ChromFixed
impute2_sam=${impute2_out}_samples
#impute2_out=${impute2_out}

if [ ! -s ${impute2_out}.vcf ]
then
	echo
	echo "${impute2_out}.vcf NOT FOUND	...	RUNNING"
	plink1.9 --gen ${impute2_gen} --sample ${impute2_sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${impute2_out}
fi


# MAKE SURE IMPUTED VCF HAS REFERENCE ALLELES THE SAME AS rCRS AND CONVERT TO MULTI-ALLELIC FOR THIS PURPOSE
echo
echo "MAKING SURE IMPUTED VCF HAS REFERENCE ALLELES THE SAME AS rCRS, AND CONVERT TO MULTI-ALLELIC FOR THIS PURPOSE"
IMP_VCF_rCRS=${impute2_out}_rCRS.vcf
IMP_VCF=${impute2_out}_FINAL.vcf.gz

if [ ! -s ${IMP_VCF_rCRS} ]
then
	echo
	echo "${IMP_VCF_rCRS} NOT FOUND	...	RUNNING"
	bcftools norm --check-ref s -f ${REF26} -m + ${impute2_out}.vcf -Oz -o ${IMP_VCF_rCRS}
fi

# DECOMPOSE
if [ ! -s ${IMP_VCF} ]
then
	echo
	echo "${IMP_VCF} NOT FOUND	...	RUNNING"
	vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
	bcftools index ${IMP_VCF}
fi


# PERFORM HAPLOGREP2 HAPLOGROUP ASSIGNMENT
HG2_FILE=${impute2_out}_FINAL_HaploGrep.txt
if [ -f ${HG2_FILE} ]
then
	echo
	echo "${HG2_FILE} FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${HG2_FILE} NOT FOUND ... RECODING TO plink1.9 VCF FILE"
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
	echo "${HiMC_FILE} NOT FOUND ... RECODING TO plink1.9 VCF FILE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/HiMC/HiMC_haplogroup_assignment.R ${PED_FILE} ${MAP_FILE} ${HiMC_FILE}  # assign haplogreps
fi

## CALCULATE Matthew's Correlation Coefficient
WGS_VCF=${WGS_relab}.gz
TYP_VCF=${TYP_n258_biallelic}
#IMP_VCF=${impute2_out}.vcf
IMP_INFO=${impute2_out}_info
OUT_FILE=${impute2_out}
first_MCC_out=${OUT_FILE}_imputed_MCC.csv

if [ -f ${OUT_FILE}_imputed_MCC.csv ] & [ -f ${OUT_FILE}_typed_MCC.csv ]
then
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv FOUND ... PIPELINE COMPLETED"
else
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi


########################################################################


# CUTOFF BY IMPUTE2 INFO SCORE
impute2_out=${out_dir}MitoImpute_${REFpanel}_imputed


impute2_file=${impute2_out}
impute2_info_file=${impute2_file}_info
impute2_file_cutoff=${impute2_file}_cutoffRetained
impute2_info_file_cutoff=${impute2_info_file}_cutoffRetained
cutoff=0.3

if [ ! -s ${impute2_file_cutoff} ] && [ ! -s ${impute2_info_file_cutoff} ]
then
	echo
	echo "${impute2_file_cutoff} AND ${impute2_info_file_cutoff} NOT FOUND	...	CUTTING OFF AT INFO ≥ ${cutoff}"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/remove_impute2_cutoff.R ${impute2_file}
else
	echo
	echo "${impute2_file_cutoff} AND ${impute2_info_file_cutoff} FOUND	...	PASSING"
fi

# FIX CHROMOSOME NAMES

InFile=${impute2_file_cutoff}
OutFile=${InFile}_ChromFixed

if [ ! -s ${OutFile} ]
then
	echo
	echo "FIXING CHROMOSOME NAMES"
	awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}
fi

# CONVERT OXFORD TO PEDIGREE

out_prefix_cutoff=${impute2_file_cutoff}
gen=${out_prefix_cutoff}_ChromFixed
sam=${impute2_sam}
out=${out_prefix_cutoff}

if [ ! -s ${out_prefix_cutoff}.ped ] && [ ! -s ${out_prefix_cutoff}.map ]
then
	echo
	echo "CONVERTING OXFORD TO PEDIGREE"
	plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}
else
	echo 
	echo "PED AND MAP FILES FOUND	...	PASSING"
fi

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=${out_prefix_cutoff}_ChromFixed
sam=${impute2_sam}
out=${out_prefix_cutoff}

if [ ! -s ${out_prefix_cutoff}.vcf ] 
then
	echo
	echo "CONVERTING OXFORD TO VCF"
	plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

else
	echo 
	echo "VCF FILES FOUND	...	PASSING"
fi

# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=${out_prefix_cutoff}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

if [ -s ${imp_vcf}.vcf ] & [ -s ${norm_imp_vcf} ] & [ -s ${imp_fasta} ] & [ -s ${vcf_pos} ] & [ -s ${fixed_vcf} ] & [ -s ${final_vcf}.vcf.gz ]
then
	echo
	echo "THE FOLLOWING FILE WERE SUPPOSEDLY FOUND:"
	echo ${final_vcf}.txt
	echo ${imp_vcf}.vcf
	echo ${norm_imp_vcf}
	echo ${imp_fasta}
	echo ${vcf_pos}
	echo ${fixed_vcf}
	echo ${final_vcf}.vcf.gz
else
	echo
	echo "FIXING IMPUTED VCF	...	BRINGING UP TO BCFTOOLS STANDARDS AND ANNOTATING"
	echo "NORMALISED VCF SAVING TO:	${norm_imp_vcf}"
	bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
	bcftools index ${norm_imp_vcf} # index normalised vcf
	
	echo "VCF POSITION INFO SAVING TO:	${vcf_pos}" 
	bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/plink_sites_map.R ${vcf_pos} # add a column with the MT label
	perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
	
	echo "FASTA FILE SAVING TO:	${imp_fasta}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
	#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
	echo "FIXED VCF SAVING TO:	${fixed_vcf}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
	bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
	bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
	#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	echo "FINAL ANNOTATED HAPLOGREP COMPATIBLE VCF SAVING TO:	${final_vcf}.vcf.gz"
	bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	bcftools index ${final_vcf}.vcf.gz # index it
	#echo "HAPLOGREP FILE SAVING TO:	${final_vcf}.txt" 
	#java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi

if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo "HAPLOGREP FILE SAVING TO:	${final_vcf}.txt" 
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi


if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink1.9 VCF FILE"
	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink1.9 (haplogrep seems to love plink1.9 vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi



TYP_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_biallelic.vcf.gz
WGS_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic.vcf.gz

## CALCULATE Matthew's Correlation Coefficient
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 
WGS_VCF=${vcf_1kg}
TYP_VCF=${strand_dir}chrMT_1kg_${MtPlatforms}.vcf.gz
TYP_VCF_MULTIALLELIC=${strand_dir}chrMT_1kg_${MtPlatforms}_multiallelic.vcf.gz
TYP_VCF_DECOMPOSED=${strand_dir}chrMT_1kg_${MtPlatforms}_biallelic_decomposed.vcf.gz

TYP_VCF_DECOMPOSED=${TYP_n258_biallelic}

#IMP_VCF=${imp_ext}.vcf
IMP_INFO=${impute2_info_file_cutoff}
OUT_FILE=${out_prefix_cutoff}
IMP_VCF_rCRS=${out_prefix_cutoff}_rCRS.vcf
IMP_VCF=${out_prefix_cutoff}_FINAL.vcf.gz

if [ ! -s ${TYP_VCF_DECOMPOSED} ]
then
	echo
	echo "NORMALISE GENOTYPED VCF"
	# NORMALISE GENOTYPED VCF
	bcftools norm --check-ref s -f ${REFMT} -m + ${TYP_VCF} -Oz -o ${TYP_VCF_MULTIALLELIC}
	
	# DECOMPOSE GENOTYPED VCF
	vt decompose ${TYP_VCF_MULTIALLELIC} | bcftools +fill-tags -Oz -o ${TYP_VCF_DECOMPOSED}
	bcftools index ${TYP_VCF_DECOMPOSED}
else
	echo
	echo "NORMALISED GENOTYPED VCF FOUND ... MOVING ON"
fi


# NORMALISE IMPUTED VCF


if [ ! -s ${IMP_VCF_rCRS} ]
then
	echo
	echo "NORMALISE IMPUTED VCF"
	bcftools norm --check-ref s -f ${REF26} -m + ${imp_ext}.vcf -Oz -o ${IMP_VCF_rCRS}
else
	echo
	echo "NORMALISED IMPUTED VCF FOUND	...	PASSING"
fi

# DECOMPOSE IMPUTED VCF
if [ ! -s ${IMP_VCF} ]
then
	echo
	echo "DECOMPOSE VCF"
	vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
	bcftools index ${IMP_VCF}
else
	echo
	echo "DECOMPOSED VCF FOUND	... PASSING"
fi

mcc_imputed_cutoff=${OUT_FILE}_imputed_MCC.csv
mcc_typed_cutoff=${OUT_FILE}_typed_MCC.csv

WGS_VCF=${WGS_relab}.gz
TYP_VCF=${TYP_n258_biallelic}

if [ -s ${mcc_imputed_cutoff} ] & [ -f ${mcc_typed_cutoff} ]
then
	echo
	echo "${mcc_imputed_cutoff} AND ${mcc_typed_cutoff} FOUND ... PIPELINE COMPLETED"
	echo "ACTUALLY ... DO IT ANYWAY (DOUBLE CHECKING, REMOVE THIS LATER"
	#Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
else
	echo
	echo "${mcc_imputed_cutoff} AND ${mcc_typed_cutoff} NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi


# GENERATE HiMC HAPLOGROUPINGS
#imp_ext=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
#full_1kGP_pref=${mitoimpute_dir}data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes
#typed_1kGP_pref=${strand_dir}chrMT_1kg_${MtPlatforms}
#imputed_1kGP_pref=${imp_ext}
#imputed_cutoff_1kGP_pref=${imp_ext}_cutoffRetained

full_1kGP_pref=${WGS_PLINK}
typed_1kGP_pref=${TYP_PLINK}
#imp_ext=${out_prefix_cutoff}
imp_ext=${impute2_file}
imputed_1kGP_pref=${imp_ext}
imputed_cutoff_1kGP_pref=${imp_ext}_cutoffRetained

himc_hg_file=${imp_ext}_HiMC_haplogroups.csv


if [ -s ${himc_hg_file} ]
then
	echo
	echo "${himc_hg_file} FOUND!	...	PASSING"
	#Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/HiMC_haplogrouping.R ${full_1kGP_pref} ${typed_1kGP_pref} ${imputed_1kGP_pref} ${imputed_cutoff_1kGP_pref}
else
	echo
	echo "${himc_hg_file} NOT FOUND!	...	GENERATING HiMC HAPLOGROUPINGS"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/HiMC_haplogrouping.R ${full_1kGP_pref} ${typed_1kGP_pref} ${imputed_1kGP_pref} ${imputed_cutoff_1kGP_pref}
fi


# GENERATE HAPLOGREP HAPLOGROUPINGS
#imp_ext=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
#full_1kGP_hg=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/1000genomes_mtDNA_haplogrep.txt
#typed_1kGP_hg=${strand_dir}chrMT_1kg_${MtPlatforms}_diploid_haplogrep.txt
#imputed_1kGP_hg=${imp_ext}_haplogrep.txt
#imputed_cutoff_1kGP_hg=${imp_ext}_cutoffRetained_haplogrep.txt

imp_ext=${impute2_file}
full_1kGP_hg=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/INFO/adni_mito_genomes_180214_fixed_n258_biallelic_HaploGrep.txt
typed_1kGP_hg=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_haplogrep.txt
imputed_1kGP_hg=${out_dir}MitoImpute_${REFpanel}_imputed_FINAL_HaploGrep.txt
imputed_cutoff_1kGP_hg=${out_dir}MitoImpute_${REFpanel}_imputed_cutoffRetained_haplogrep.txt


haplogrep_hg_file=${imp_ext}_HaploGrep_haplogroups.tsv


if [ -s ${haplogrep_hg_file} ]
then
	echo
	echo "${haplogrep_hg_file} FOUND!	...	PASSING"
else
	echo
	echo "${haplogrep_hg_file} NOT FOUND!	...	GENERATING HiMC HAPLOGROUPINGS"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/HaploGrep_haplogrouping.R ${full_1kGP_hg} ${typed_1kGP_hg} ${imputed_1kGP_hg} ${imputed_cutoff_1kGP_hg}
fi



# SUMMARISE EVERYTHING!
final_summary_file=${out_dir}MitoImpute_${REFpanel}_imputed_SUMMARY.csv

mcmc_str="MCMC${mcmc}"
khap_str="kHAP${khap}"

MtPlatforms="ADNI_n258"
maf_str="MAF1%"
khap_str="kHAP500"
MTSnps=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/INFO/mito_snps_rcrs_ed_n258_SITE_LIST.txt

ls -lh ${final_summary_file}
echo
echo ${MtPlatforms}
echo ${mcmc_str}
echo ${maf_str}
echo ${khap_str}
echo ${himc_hg_file}
echo ${haplogrep_hg_file}
echo ${impute2_info_file}
echo ${impute2_info_file_cutoff}
echo ${MTSnps}
echo ${first_MCC_out}
echo ${mcc_imputed_cutoff}
echo


if [ -s ${final_summary_file} ]
then
	echo
	echo "${final_summary_file} FOUND	...	PASSING"
	echo "BUT YOU NEDE TO DO IT AGAIN!"
	#Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/summarise_HiMC_HaploGrep_perchip.R ${MtPlatforms} ${mcmc_str} ${maf_str} ${khap_str} ${himc_hg_file} ${haplogrep_hg_file} ${impute2_info_file} ${impute2_info_file_cutoff} ${MTSnps} ${first_MCC_out} ${mcc_imputed_cutoff}
else
	echo
	echo "${final_summary_file} NOT FOUND	...	SUMMARISING EVERYTHING"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/summarise_HiMC_HaploGrep_perchip.R ${MtPlatforms} ${mcmc_str} ${maf_str} ${khap_str} ${himc_hg_file} ${haplogrep_hg_file} ${impute2_info_file} ${impute2_info_file_cutoff} ${MTSnps} ${first_MCC_out} ${mcc_imputed_cutoff}
fi 

echo
echo "END!"
echo


########################################################################