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
REFpanel="ReferencePanel_v6"
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
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

bcftools view -S ~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.txt ${TYP_VCF} | bcftools +fill-tags -Oz -o ${TYP_n258}
bcftools index ${TYP_n258}

bcftools view -S ~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH_reseq.txt ${WGS_VCF} | bcftools +fill-tags -Oz -o ${WGS_n258}
bcftools index ${WGS_n258}

# FILTER VCF TO CONTAIN ONLY SNVs 
echo
echo "FILTER VCF TO CONTAIN ONLY SNVs"
TYP_n258_tmp=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_tmp.vcf.gz
WGS_n258_tmp=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_tmp.vcf.gz
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 

#bcftools annotate -x INFO,^FORMAT/GT ${WGS_n258} | bcftools norm --check-ref s -f ${REF26} -m - | bcftools view -V indels,mnps | bcftools +fill-tags -Oz -o ${WGS_n258_biallelic}
#bcftools annotate -x INFO,^FORMAT/GT ${TYP_n258} | bcftools norm --check-ref s -f ${REFMT} -m - | bcftools view -V indels,mnps | bcftools +fill-tags -Oz -o ${TYP_n258_biallelic}

bcftools annotate -x INFO,^FORMAT/GT ${WGS_n258} | bcftools norm --check-ref s -f ${REF26} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${WGS_n258_tmp}
bcftools annotate -x INFO,^FORMAT/GT ${TYP_n258} | bcftools norm --check-ref s -f ${REFMT} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${TYP_n258_tmp}

bcftools index ${WGS_n258_tmp}
bcftools index ${TYP_n258_tmp}

# DECOMPOSE VCF FILES AND SPLIT MULTIALLELIC RECORDS INTO BIALLELIC 
echo
echo "DECOMPOSE VCF FILES AND SPLIT MULTIALLELIC RECORDS INTO BIALLELIC"

TYP_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_biallelic.vcf.gz
WGS_n258_biallelic=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic.vcf.gz

vt decompose ${WGS_n258_tmp} | bcftools +fill-tags -Oz -o ${WGS_n258_biallelic}
vt decompose ${TYP_n258_tmp} | bcftools +fill-tags -Oz -o ${TYP_n258_biallelic}

bcftools index ${WGS_n258_biallelic}
bcftools index ${TYP_n258_biallelic}

# RELABEL WGS TO ADNI1 SAMPLE IDs
samples_csv=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv
WGS_relab=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled.vcf

python ~/GitCode/MitoImputePrep/scripts/PYTHON/fix_vcf_names.py -i ${TYP_n258_biallelic} -o ${WGS_relab} -c ${samples_csv} -v
# TIM YOU HAVE TO MAKE THIS COMMAND LINE >:(

# GZIP USING BCFTOOLS TO MAKE SURE VCF CONFORMS TO STANDARDS
bcftools view ${WGS_relab} -Oz -o ${WGS_relab}.gz

# EXTRACT SAMPLE IDS 
echo 
echo "EXTRACT SAMPLE IDS"

WGS_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/INFO/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled_sampleID.txt
WGS_SEX_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/INFO/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled_sampleID_sex.txt
TYP_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/INFO/mito_snps_rcrs_ed_n258_biallelic_sampleID.txt
TYP_SEX_SAMPLE=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/INFO/mito_snps_rcrs_ed_n258_biallelic_sampleID_sex.txt

bcftools query -l ${WGS_relab} > ${WGS_SAMPLE}
Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${WGS_SAMPLE} ${WGS_SEX_SAMPLE}

bcftools query -l ${TYP_n258_biallelic} > ${TYP_SAMPLE}
Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${TYP_SAMPLE} ${TYP_SEX_SAMPLE}

# GENERATE PLINK FILES (PED AND MAP FILES!)
echo
echo "GENERATE PLINK FILES (PED AND MAP FILES!)"
WGS_PLINK=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/PLINK/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled
TYP_PLINK=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/VCF/mito_snps_rcrs_ed_n258_biallelic
plink1.9 --vcf ${WGS_relab} --recode --double-id --keep-allele-order --out ${WGS_PLINK}
plink1.9 --vcf ${TYP_n258_biallelic} --recode --double-id --keep-allele-order --out ${TYP_PLINK}

# GENERATE GEN FILES
echo 
echo "GENERATE GEN AND SAMPLE FILES"
WGS_GEN_OUT=${WK_DIR}MitoImpute/data/ADNI_REDO/WGS/HAP_LEGEND_GEN/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled
TYP_GEN_OUT=${WK_DIR}MitoImpute/data/ADNI_REDO/GENOTYPED/HAP_LEGEND_GEN/mito_snps_rcrs_ed_n258_biallelic


bcftools convert --gensample ${WGS_GEN_OUT} ${WGS_relab} --sex ${WGS_SEX_SAMPLE}
bcftools convert --gensample ${TYP_GEN_OUT} ${TYP_n258_biallelic} --sex ${TYP_SEX_SAMPLE}

# FIX SEX ID COLUMN IN .samples FILE
echo
echo "FIXING SEX ID COLUMN IN .samples FILE"
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${WGS_GEN_OUT}.samples
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${TYP_GEN_OUT}.samples

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
	Rscript ~/GitCode/MitoImputePrep/scripts/R/HiMC_haplogroup_assignment.R ${PED_FILE} ${MAP_FILE} ${HiMC_FILE}  # assign haplogreps
fi

## CALCULATE Matthew's Correlation Coefficient
WGS_VCF=${WGS_relab}.gz
TYP_VCF=${TYP_n258_biallelic}
#IMP_VCF=${impute2_out}.vcf
IMP_INFO=${impute2_out}_info
OUT_FILE=${impute2_out}

if [ -f ${OUT_FILE}_imputed_MCC.csv ] & [ -f ${OUT_FILE}_typed_MCC.csv ]
then
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv FOUND ... PIPELINE COMPLETED"
else
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi

#
## GENERATE GEN SAMPLE
#echo
#echo "GENERATING GEN SAMPLE"
#sex=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH_SEX.txt 
#out=${WK_DIR}MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni
#
#bcftools convert --gensample ${out} ${vcf} --sex ${sex}
#Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${out}.samples
#
## GENERATE PLINK FILES
#echo
#echo "GENERATING PLINK FILES"
#out=${WK_DIR}MitoImpute/data/ADNI/Timpute/PLINK/mitoimpute_adni
#
#plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}
#
## RUN IMPUTE2
#echo
#echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
#m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
#h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
#l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
#g=${WK_DIR}MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni.gen.gz
#s=${WK_DIR}MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni.samples
#out=${WK_DIR}MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/mitoimpute_adni_imputed_MCMC${mcmc}
##g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
##s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
##out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/MCMC${mcmc}/chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
#
#if [ -d ${WK_DIR}MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/ ]
#then
#	echo "DIRECTORY FOUND"
#else
#	mkdir -p ${WK_DIR}MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/
#fi
#
#if [ -f ${out} ]
#then
#	echo "${out} FOUND! ... PASSING"
#else
#	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
#	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out} -iter ${mcmc} -burnin ${burn}
#fi
#
## FIX CHROMOSOME NAMES
#echo
#echo "FIXING CHROMOSOME NAMES"
#InFile=${WK_DIR}MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/mitoimpute_adni_imputed_MCMC${mcmc}
#OutFile=${InFile}_ChromFixed
#awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}
#
## CONVERT OXFORD TO PEDIGREE
#echo
#echo "CONVERTING OXFORD TO PEDIGREE"
#gen=${InFile}_ChromFixed
#sam=${InFile}_samples
#out=${InFile}
#
#plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}
#
## CONVERT OXFORD TO VCF
#echo
#echo "CONVERTING OXFORD TO VCF"
#gen=${InFile}_ChromFixed
#sam=${InFile}_samples
#out=${InFile}
#
#plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}
#
## CONVERT VCF TO FORMAT FOR HAPLOGREP2
#ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
#imp_ext=${InFile}
#imp_vcf=${imp_ext}.vcf
#norm_imp_vcf=${imp_ext}_norm.vcf.gz
#imp_fasta=${imp_ext}.fasta
#vcf_pos=${imp_ext}_norm_SNPpositions.txt
#fixed_vcf=${imp_ext}_fixed.vcf
#final_vcf=${imp_ext}_haplogrep
#
#bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
#bcftools index ${norm_imp_vcf} # index normalised vcf
#bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
#Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
#perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
#python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
##python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
#bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
#bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
##bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
#bcftools index ${final_vcf}.vcf.gz # index it
#java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
#
#if [ -f ${final_vcf}.txt ]
#then
#	echo
#	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
#else
#	echo
#	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
#	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
#	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
#fi
#