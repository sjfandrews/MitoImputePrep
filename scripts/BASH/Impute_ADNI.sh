#!/bin/bash

# ADNI RESEQUENCED VCF:		/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# ADNI GENOTYPED VCF:		/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI/mito_snps_rcrs_ed.vcf

# ADNI FILE (BELOW) HAD TO BE MANUALLY FIXED TO CONTAIN THE LINE "##contig=<ID=26,length=16569>"
# /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214.vcf
# RENAMED TO: /Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214_fixed.vcf

# SPECIFY REFERENCE PANEL
REFpanel="ReferencePanel_v5"
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
echo
echo "REFERENCE PANEL:	${REFpanel}"
echo "HAPLOGREP:		${HAPLOGREP}"

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