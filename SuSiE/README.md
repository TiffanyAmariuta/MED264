# Finemapping exercises


Go to SuSiE tutorial: https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html

If on Expanse, launch an interactive session and load R:
```
module load gcc/9.2.0
module load r
```
### Launch R and install required libraries:
```
R
> install.packages("susieR") #could take up to 10 minutes. 
> library(susieR)
> set.seed(1)
> data(N3finemapping)
> attach(N3finemapping)
> #Y: phenotypes across 572 individuals. 2 phenotypes. 
> #X: genotype matrix (centered but not standardized) 
> #true_coef: true causal effect sizes for SNPs for each phenotype.
> dim(X) #574 x 1001 
> dim(Y) #572 x 2
> X_std <- scale(X) #can double check that the variance of each column is 1 and the mean is still 0. 
> b <- true_coef[,1] #true causal effect sizes for SNPs for first phenotype. 
> sumstats <- univariate_regression(X_std, Y[,1]) #performs a GWAS 
> z <- sumstats$betahat/sumstats$sebetahat
```

For the first in-class exercise, how many SNPs are significant at a threshold of p < 5x10^-8?
Hint: convert the z-score to a p-value using the pnorm() function.

Now, let's do some fine-mapping: 

```
> fitted <- susie(X_std, Y[,1],L = 10,verbose = TRUE)
> print(fitted$sets)
> which(b != 0)
```

Note that each of the three causal variants is reprenseted by a credible set identified by SuSiE. 

For the second in-class exercise, what happens if we change the parameter L (assumption on the number of causal variants)? Try the values 1:10. How do the credible sets match up with the true causal variants?   

For the third in-class exercise, explore the relationship between run time and size of locus (e.g. number of SNPs). Is the relationship quadratic? Linear? How is it affected by the number of causal signals present? Hint: use the following code in R to measure the time of running a piece of code:  

```
> start <- proc.time()
> run_code_here()
> proc.time() - start
```

For the fourth in-class exercise, how does SuSiE perform when the 3 causal variants are missing? How many credible sets? What is different than when all causal variants were present?
```
> w <- which(b != 0) #403 653 773
> X_miss <- X_std[,-w]
> fitted <- susie(X_miss, Y[,1],L = 10,verbose = TRUE)
```

Follow-up question: what could account for the difference in fine-mapping results using X_std vs X_miss? Hint: think about LD.   



