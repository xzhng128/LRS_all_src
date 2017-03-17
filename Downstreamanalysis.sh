#!bin/bash
set -e
set -u

date

INTERVALS_BED=IAD66395_209_Designed.bed
MPILEFILT=PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf.gz
GATKFILT=PHETE-FILTERED_final.GATK.break.vcf.recode.vcf.gz

#change names to sample names instead of barcodes. the python script has a help flag (-h)
#python rename.py -DICT dict.txt -VCF PHETE-FILTERED_final.GATK.break.vcf.recode.vcf
#python rename.py -DICT dict.txt -VCF PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf

if [ -f PHETE-FILTERED_final.GATK.break.vcf.recode.vcf ]
  then bgzip PHETE-FILTERED_final.GATK.break.vcf.recode.vcf
  tabix -p vcf PHETE-FILTERED_final.GATK.break.vcf.recode.vcf.gz
fi

if [ -f PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf ]
  then bgzip PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf
  tabix -p vcf PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf.gz
fi

if [ -f final.hapmap.recode.vcf ]
  then bgzip final.hapmap.recode.vcf
  tabix -p vcf final.hapmap.recode.vcf.gz
fi

#venn diagram
vcf-compare -a $MPILEFILT $GATKFILT final.hapmap.recode.vcf.gz > compare.txt

bgzip -d $GATKFILT
bgzip -d $MPILEFILT
bgzip -d final.hapmap.recode.vcf.gz



#figure 4
#version of snpeff was updated and didn't have this library in newer version of snpeff
snpEff -v -csvStats snpeffbed.txt AGPv3.27 -i bed $INTERVALS_BED -s summaryintervals > bed_targets.vcf
snpEff -v -csvStats snpeffsamtools.txt AGPv3.27 PHETE-FILTERED_allsamplesafterfiltering.dupremoved.mpileup.vcf -s summarysamtools > samtools_annotation.vcf
snpEff -v -csvStats snpeffGATK.txt AGPv3.27 PHETE-FILTERED_final.GATK.break.vcf.recode.vcf -s summaryGATK > GATK_annotation.vcf
snpEff -v -csvStats snpeffhapmap.txt AGPv3.27 final.hapmap.recode.vcf -s summaryhapmap > hapmap_annotation.vcf

#py python rename.py -DICT dict.txt -VCF PHETE-FILTERED_final.GATK.break.vcf.recode.vcf

#figure 5
vcftools --vcf PHETE-FILTERED_final.GATK.break.vcf.recode.vcf --plink --out GATKplink
plink --file GATKplink --cluster --mds-plot 2 --noweb


#FST

##get ready for tassel input
vcftools --vcf PHETE-FILTERED_final.GATK.break.vcf.recode.vcf --remove-indels --remove-filtered-all --maf 0.02 --recode --recode-INFO-all --out SNPs_only
bcftools annotate -x INFO,^FORMAT/GT SNPs_only.recode.vcf > converttassel.vcf
sed '/*/d' converttassel.vcf > converttasselsed.vcf
/Users/tjamann/bin/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile converttasselsed.vcf -outputFile sortedconverttasselsed.vcf -fileType VCF
date
#convert to text with TASSEL- I don't know how to do this with the command line. also apply default paramteters to remove monomorphic sites.
#run Jim's Fst script for r.



