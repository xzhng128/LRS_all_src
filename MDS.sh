
#need to do once: ulimit -n 3000
#plink version 1.07 
vcftools --vcf landrace_10_nopress.recode.vcf --plink --out chr10plink
plink --file chr10plink --mind 0.1 --cluster --mds-plot 2 --noweb --allow-no-sex