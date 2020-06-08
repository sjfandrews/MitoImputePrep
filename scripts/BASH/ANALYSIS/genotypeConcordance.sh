#!/bin/bash

COMMAND=GenotypeConcordance
#SAMPLE_VCF=/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/IMPUTE2/MCMC1/mitoimpute_adni_imputed_MCMC1_norm.vcf.gz
SAMPLE_VCF=/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/VCF/mitoimpute_adni/mitoimpute_adni.vcf.gz
REF_VCF=/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/ADNI_reseq/adni_mito_genomes_180214_n258_subset_relabelled.vcf.gz
FASTA=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta
OUT_FILE=/Users/u5015730/Desktop/gatk_gc_test.grp

SAMPLE_COMP=ADNI_0074

/Applications/gatk-4.1.0.0/gatk ${COMMAND} \
	--CALL_VCF ${SAMPLE_VCF} \
	--TRUTH_VCF ${REF_VCF} \
	--OUTPUT ${OUT_FILE} \
	--IGNORE_FILTER_STATUS \
	--CALL_SAMPLE ${SAMPLE_COMP} \
	--TRUTH_SAMPLE ${SAMPLE_COMP} \
	--REFERENCE_SEQUENCE ${FASTA}