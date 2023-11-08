library(data.table)
library(susieR)

#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 

library(data.table)
y <- fread("UKB_460K.cancer_BREAST.sumstats.gz", header = T)
#get SNPs in the window corresponding to that gene of interest. 
#Gene of interest: ENSG00000078403 Chromosome 10: 21,524,646-21,743,630 
snps <- fread("/expanse/lustre/projects/ddp412/tamariutabartell/1000Genomes/LDREF/1000G.EUR.10.bim", header = F)
w <- which(snps$V4 > 21524646 - 500000 & snps$V4 < 21743630 + 500000) #give a window around it
ww <- which(y$SNP %in% snps$V2[w]) #422
y <- y[ww,]
write.table(y, file = "BrCa_GWAS_ENSG00000078403_locus.txt", row.names = F, col.names = T, sep = "\t", quote = F)
system("gzip -f BrCa_GWAS_ENSG00000078403_locus.txt")

y <- fread("BrCa_GWAS_ENSG00000078403_locus.txt.gz", header = T)
z_scores <- y$Z
susie_plot(z_scores, y = "z", b=rep(0,length(z_scores))) #don't know the causal variants so nothing to put in "b" variable 

R <- fread("LDmatrix_ENSG00000078403_locus.txt.gz", header = F)
fitted_rss1 <- susie_rss(z = z_scores, n = 459324, R = R, L = 10,
                         estimate_residual_variance = TRUE)

