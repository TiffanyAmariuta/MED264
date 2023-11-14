library(data.table)
library(susieR)

#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 

library(data.table)
y <- fread("UKB_460K.cancer_BREAST.sumstats.gz", header = T)
#get SNPs in the window corresponding to that gene of interest. 
#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 
snps <- fread("1000G.EUR.10.bim", header = F)
w <- which(snps$V4 > 21524646 - 500000 & snps$V4 < 21743630 + 500000) #give a window around it
ww <- which(y$SNP %in% snps$V2[w]) #422
y <- y[ww,]
write.table(y, file = "BrCa_GWAS_ENSG00000078403_locus.txt", row.names = F, col.names = T, sep = "\t", quote = F)
system("gzip -f BrCa_GWAS_ENSG00000078403_locus.txt")

sumstats <- fread("BrCa_GWAS_ENSG00000078403_locus.txt.gz", header = T)
susie_plot(sumstats$Z, y = "z")

R <- fread("LDmatrix_ENSG00000078403_locus.txt.gz", header = F)
fitted_rss = susie_rss(z = sumstats$Z, R = as.matrix(R), n = 459324, L = 10,
                        estimate_residual_variance = TRUE)
susie_plot(fitted_rss, y="PIP")
summary(fitted_rss)$cs #there are 6 credible sets, each with one causal variant. 
                         
#3. What are the genomic coordinates of causal SNPs? 
pos <- snps$V4[match(sumstats$SNP[as.numeric(summary(fitted_rss)$cs$variable)], snps$V2)]
#chr10: 21849769 21989245 21983960 21856475 21884471 22062729
