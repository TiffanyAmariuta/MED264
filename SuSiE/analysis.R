library(data.table)
library(susieR)

#Use the SuSiE tutorial website as a reference: https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html

#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 

library(data.table)
y <- fread("UKB_460K.cancer_BREAST.sumstats.gz", header = T)
#get SNPs in the window corresponding to that gene of interest. 
#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 
#1. Find the SNPs that are within 500kb +/- the gene body using the 1000G.EUR.10.bim file 
write.table(y, file = "BrCa_GWAS_ENSG00000078403_locus.txt", row.names = F, col.names = T, sep = "\t", quote = F)
system("gzip -f BrCa_GWAS_ENSG00000078403_locus.txt")

y <- fread("BrCa_GWAS_ENSG00000078403_locus.txt.gz", header = T)
susie_plot(sumstats$Z, y = "z")

R <- fread("LDmatrix_ENSG00000078403_locus.txt.gz", header = F) #Tiffany did this ahead of time using plink2R and computing LD using 1KG individuals for this locus
fitted_rss = susie_rss(z = sumstats$Z, R = as.matrix(R), n = 459324, L = 10,
                        estimate_residual_variance = TRUE)
#2. Plot the susie results and use the summary() function to understand the credible sets that were found 

                        
#3. What are the genomic coordinates of causal SNPs? 

