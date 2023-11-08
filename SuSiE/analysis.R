library(data.table)
library(susieR)
sumstats <- fread("chr10_LymCount_GWAS.txt.gz", header = T)
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)

R <- fread("LDmatrix_chr10locus.txt.gz", header = F)
fitted_rss1 <- susie_rss(bhat = sumstats$beta, shat = sumstats$se, n = 524867, R = R, var_y = 1), L = 10,
                         estimate_residual_variance = TRUE)

