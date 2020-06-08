#!/bin/bash

# ADNI RESEQUENCED VCF:		/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# ADNI GENOTYPED VCF:		/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI/mito_snps_rcrs_ed.vcf

# ADNI FILE (BELOW) HAD TO BE MANUALLY FIXED TO CONTAIN THE LINE "##contig=<ID=26,length=16569>"
# /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# RENAMED TO: /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214_fixed.vcf

# SPECIFY REFERENCE PANEL
REFpanel="ReferencePanel_v5"
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
mcmc=1
burn=0

WGS_VCF=""
TYP_VCF=""
IMP_VCF=""

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

# GENERATE HAPLOGREP HAPLOGROUP ASSIGNMENTS FROM RESEQUENCED DATA
echo
echo "GENERATING HAPLOGREP HAPLOGROUP ASSIGNMENTS FROM RESEQUENCED DATA"

ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
reseq_ext=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214
reseq_vcf=${reseq_ext}_fixed.vcf
norm_vcf=${reseq_ext}_norm.vcf.gz
reseq_fasta=${reseq_ext}.fasta
vcf_pos=${reseq_ext}_norm_SNPpositions.txt
fixed_vcf=${reseq_ext}_fixed2.vcf
final_vcf=${reseq_ext}_haplogrep

bcftools annotate --set-id '.' -x INFO,^FORMAT/GT ${reseq_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
bcftools index ${norm_vcf} # index normalised vcf
bcftools query -f '%POS\n' ${norm_vcf} > ${vcf_pos} # extract genomic positions
Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_vcf} -o ${reseq_fasta} -v # convert to a fasta file
python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${reseq_fasta} -o ${fixed_vcf} -g -d -id -a -v # convert back to a vcf
bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools index ${final_vcf}.vcf.gz # index it
java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --out ${final_vcf}.txt # assign haplogreps
java -jar ${HAPLOGREP} --in ${norm_vcf} --format vcf --out ${final_vcf}2.txt # assign haplogreps

if [ -f ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --out ${final_vcf}.txt # assign haplogreps
fi

# DECOMPOSE RESEQUENCED DATA VCF
echo
echo "DECOMPOSING RESEQUENCED DATA VCF"

samp_file=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH_reseq.txt
reseq_ext=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214
decom_ext=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/DECOMPOSED/adni_mito_genomes_180214
orig_vcf=${reseq_ext}_fixed.vcf
norm_vcf=${decom_ext}_norm.vcf.gz
decom_vcf=${decom_ext}_norm_decomposed.vcf.gz
final_vcf=${decom_ext}_norm_decomposed_firstAlt.vcf.gz
plink_vcf=${decom_ext}_norm_decomposed_firstAlt
samps_adni=/Volumes/TimMcInerney/MitoImpute/metadata/SampleList_ADNI.txt
sex_adni=/Volumes/TimMcInerney/MitoImpute/metadata/SampleList_ADNI_sex.txt

bcftools annotate -x INFO,^FORMAT/GT ${orig_vcf} | bcftools norm -f ${ref_fasta_plink} -m - | bcftools view -V indels,mnps -S ${samp_file} | bcftools norm -m + | bcftools +fill-tags -Oz -o ${norm_vcf}
vt decompose ${norm_vcf} | bcftools +fill-tags -Oz -o ${decom_vcf}
#python3 ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt ${decom_vcf} | bcftools view -Oz -o ${vcf_1kg}
#bcftools index ${vcf_1kg}
plink1.9 --vcf ${decom_vcf} --recode --double-id --keep-allele-order --out ${plink_vcf}
bcftools query -l ${decom_vcf} > ${samps_adni}
Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${samps_adni} ${sex_adni}

# STRIP AWAY PROBLEMATIC COLUMNS
echo
echo "STRIPPING AWAY PROBLEMATIC COLUMNS"

orig_adni=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI/mito_snps_rcrs_ed.vcf
ref_fasta=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta

bcftools annotate -x INFO,^FORMAT/GT ${orig_adni} | bcftools norm --check-ref s -f ${ref_fasta} -m +any | bcftools +fill-tags -Oz -o ${orig_adni}.gz
bcftools index ${orig_adni}.gz

# SUBSET TO ONLY SAMPLES FOUND IN BOTH ADNI DATASETS
echo
echo "SUBSETTING TO ONLY SAMPLES FOUND IN BOTH ADNI DATASETS"
samp_file=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.txt
vcf=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/VCF/mitoimpute_adni/mitoimpute_adni.vcf.gz

bcftools view -S ${samp_file} ${orig_adni}.gz | bcftools +fill-tags -Oz -o ${vcf}
bcftools index ${vcf}

# CREATE DIPLOID VCF
echo
echo "GENERATING PLINK FILES (DIPLOID)"
geno_ext=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/VCF/mitoimpute_adni/mitoimpute_adni
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta
norm_vcf=${geno_ext}_norm.vcf.gz
vcf_pos=${geno_ext}_norm_SNPpositions.txt
geno_fasta=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/FASTA/mitoimpute_adni.fasta
fixed_vcf=${geno_ext}_fixed.vcf
diploid_vcf=${geno_ext}_diploid

#bcftools annotate --set-id '.' ${vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
#bcftools index ${norm_vcf} # index normalised vcf
bcftools query -f '%POS\n' ${vcf} > ${vcf_pos} # extract genomic positions
Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${vcf} -o ${geno_fasta} -v # convert to a fasta file
python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${geno_fasta} -o ${fixed_vcf} -g -d -id -a -v # convert back to a vcf
bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${diploid_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools index ${diploid_vcf}.vcf.gz # index it

java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf.gz --format vcf --chip --out ${diploid_vcf}.txt # assign haplogreps

if [ -f ${diploid_vcf}.txt ]
then
	echo
	echo "${diploid_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${diploid_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink1.9 --vcf ${diploid_vcf}.vcf.gz --recode vcf --out ${diploid_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --chip --out ${diploid_vcf}.txt # assign haplogreps
fi

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
sex=~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH_SEX.txt 
out=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/PLINK/mitoimpute_adni

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni.gen.gz
s=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/OXFORD/mitoimpute_adni.samples
out=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/mitoimpute_adni_imputed_MCMC${mcmc}
#g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
#s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
#out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/MCMC${mcmc}/chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}

if [ -d /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/ ]
then
	echo "DIRECTORY FOUND"
else
	mkdir -p /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/
fi

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out} -iter ${mcmc} -burnin ${burn}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC${mcmc}/mitoimpute_adni_imputed_MCMC${mcmc}
OutFile=${InFile}_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=${InFile}_ChromFixed
sam=${InFile}_samples
out=${InFile}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=${InFile}_ChromFixed
sam=${InFile}_samples
out=${InFile}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=${InFile}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
bcftools index ${norm_imp_vcf} # index normalised vcf
bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools index ${final_vcf}.vcf.gz # index it
java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps

if [ -f ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi
