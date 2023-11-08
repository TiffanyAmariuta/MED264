traiti <- 76 #Breast cancer, index in the Sumstats2021_Update.txt.gz file 
ciswindow <- 2e6 #how many nucleotides we consider for co-regulation between genes. 
library(data.table) #install.packages("data.table")
#library(Hmisc) #install.packages("Hmisc") #for weighted variance function, but replaced by manual calculation

traits <- fread("Sumstats2021_Update.txt.gz", header= T) #full list of traits with h2 and sample size values
trait <- traits$Trait_Identifier[traiti] #match trait
N <- traits$N[traiti] #get sample size for GWAS

y <- fread("TissueGroups2.txt", header = T) #get tissue meta data
groups <- unique(y$GTEx_sameN_lowN) #tissues for analysis
tissues <- groups[-which(groups == "Remove")] #only keep tissues we want to analyze
small_tissues <- c(3,5,11,12,14,22,24,27:30,33:35,37:38) #with sample size < 320
normal_tissues <- c(1:length(tissues))[-small_tissues] #sample size = 320 

get_cov_alpha1alpha1_multitissue <- function(alpha_z,CoRegMat,nGWAS,weights1,weights2,weights3){ #function to run TCSC regression, TWAS chi-sq ~ co-regulation scores * N/G
y <- (alpha_z^2)/nGWAS #free intercept 
w1 <- 1/weights1 #heteroscedasticity of twas chisq
w2 <- 1/weights2 #co-regulation
w3 <- 1/weights3 #tissue redundancy
mod <- summary(glm(y ~ CoRegMat, weights = w1*w2*w3)) #multiple weights in weighted linear regression
cov_b1b1 <- coef(mod)[-1,1]  #remove intercept, get betas for each tissue
cov_b1b1 <- cov_b1b1[1:ncol(CoRegMat)] #only take value for each tissue
return(cov_b1b1)
}

load("TCSC_input_320_smalltoo_tissuecoreg_modjk_removemhc_PCbiascorr_removesmall_June21.RData") #co-regulation scores across tissues
gtex <- fread("gtexGeneV8.txt.gz", header = F, sep = "\t") #list of genes and their positions in genome

#### trait specific analysis #### 

#read transcript with it to decide what is protein coding and what isn't
alpha_z <- c()
for (i in 1:length(tissues)){
if(i %in% small_tissues){ # n < 320 sample size 
transcript_key <- fread(paste0("1KG_PredExp_Transcripts_v8_allEUR_blup/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1 #list of genes
keep <- fread(paste0("1KG_PredExp_Transcripts_v8_allEUR_blup/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1 #list of QC genes
z <- fread(paste0("Marginal_alphas_sumstats_1KG_v8_allEUR_blup/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1 #TWAS result gene-disease association
ww <- which(transcript_key %in% keep) 
transcript_key <- transcript_key[ww]  
z <- z[ww] 
m <- match(transcript_key,gtex$V7) #figure out if these genes are protein coding or not 
genetype <- gtex$V8[m]
ww <- which(genetype == "protein_coding") #only keep PC genes
z <- z[ww] 
alpha_z <- c(alpha_z,z)
#alpha_z <- c(alpha_z,z[which(genetype == "protein_coding")])
}else{ #same steps but for n = 320 tissues
transcript_key <- fread(paste0("1KG_PredExp_Transcripts_v8_320EUR/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("1KG_PredExp_Transcripts_v8_320EUR/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
z <- fread(paste0("Marginal_alphas_sumstats_1KG_v8_320EUR/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
ww <- which(transcript_key %in% keep) 
transcript_key <- transcript_key[ww] 
z <- z[ww] 
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
ww <- which(genetype == "protein_coding")
z <- z[ww] 
alpha_z <- c(alpha_z,z)
#alpha_z <- c(alpha_z,z[which(genetype == "protein_coding")])
}
}
alpha_z <- alpha_z[-remove_genes] #perform QC
alpha_z <- alpha_z[-remove_genes_from_smalltissues] #perform QC
if(length(w) > 0){alpha_z <- alpha_z[-w]}
w <- which(alpha_z^2 > max(80,0.001*N) | is.na(alpha_z)) #trait specific qc, remove outliers, if TWAS chi-sq is too high, we may not believe it. 
if(length(w)>0){
alpha_z <- alpha_z[-w]
transcripts <- transcripts[-w]
starts <- starts[-w]
chrs <- chrs[-w]
tissueassign <- tissueassign[-w]
X <- X[-w,]
}
N_tissuespecific <- as.numeric(table(tissueassign)) #how many genes are modeled in each tissue, indexed by tissue 
a_transcripts <- table(transcripts) #figure out how many times the same gene is modeled in different tissues 
weights3 <- sapply(1:length(transcripts), function(x) as.numeric(a_transcripts)[match(transcripts[x], names(a_transcripts))])  #how many tissues express each gene
tissueassign_tissues <- tissues[tissueassign] #indexed by gene, what tissue does each gene belong to 
y <- fread("TissueGroups2.txt", header = T)
groups <- unique(y$GTEx_sameN)
tissues <- groups[-which(groups == "Remove")] #23
tissueassign <- match(tissueassign_tissues,tissues)
weights2 <- sapply(1:nrow(X), function(x) sum(X[x,-tissueassign[x]])) + 1 #total co-regulation score weight
expected_true_cish2_genes <- N_tissuespecific #how many gene models go into the TCSC model per tissue, indexed by tissue 

#### set up jackknife for estimating standard error #### 
a <- cbind(1:length(starts),as.numeric(chrs),as.numeric(starts),unlist(sapply(1:length(N_tissuespecific), function(x) rep(x,N_tissuespecific[x]))))
a <- a[order(a[,2],a[,3],decreasing = F),]
chunks <- 200 #CHANGE THIS IF YOU LIKE
#assigns each gene to a chunk, each chunk is adjacent and of the same size across the genome
size_groups <- floor(length(starts)/chunks) #figure out how many genes need to be in each chunk
size_group5 <- length(starts) - (chunks-1)*size_groups
group_assignment <- c()
for (j in 1:chunks){
if(j == chunks){
group_assignment <- c(group_assignment, rep(j,size_group5))
}else{
group_assignment <- c(group_assignment, rep(j,size_groups))
}
}
a <- cbind(a, group_assignment) #add assignment for chunk
a <- a[order(a[,1], decreasing = F),] #sort by chunk number 
a <- cbind(a,transcripts) #add gene meta-data

#calculate heteroscedasticity weight for regression
mean_chisq <- mean(alpha_z^2, na.rm = T) #mean chi-squared
totcoreg <- rowSums(X) #co-regulation scores 
crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg)) #crude heritability estimate w/o regression. General tissue-wide signal
weights1 <- (1 + N*crude_h2_est*totcoreg)^2 #heteroscedasticity - proportional to co-regulation. 

variance_mat <- matrix(0,1+length(tissues),6) #estimate, jackknife SE, P
variance_mat[,1] <- c(tissues,"AllTissues")
variance_mat[1:length(tissues),2] <- get_cov_alpha1alpha1_multitissue(alpha_z,X,N,weights1,weights2,weights3) * expected_true_cish2_genes #run TCSC regression
variance_mat[(length(tissues)+1),2] <- sum(as.numeric(variance_mat[1:length(tissues),2])) #add up h2 across tissues

#estimate the SE for each tissue-specific h2 result
jk <- matrix(0,nrow = chunks,ncol = length(tissues))
jk_weights <- matrix(0,nrow = chunks,ncol = 1)
jk_sum <- matrix(0,nrow = chunks,ncol = 1)
for (chunk in 1:chunks){
print(chunk)
remove_genes <- which(a[,5] == chunk) #remove genes in each chunk at a time, then redo TCSC regression.
print(length(remove_genes))
tab <- table(a[remove_genes,4]) #freq of tissues
subtract_genes <- rep(0,length(tissues))
m <- match(1:length(tissues), names(tab))
w <- which(!is.na(m))
subtract_genes[w] <- as.numeric(tab)[m[w]]
N_tissuespecific_jk <- sapply(1:length(N_tissuespecific), function(x) N_tissuespecific[x] - subtract_genes[x])
alpha_z_jk <- alpha_z[-remove_genes]
mean_chisq <- mean(alpha_z_jk^2, na.rm = T)
X_jk <- X[-remove_genes,]
totcoreg <- rowSums(X_jk)
crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg)) #need crude h2 estimate to compute weights again
weights1_jk <- (1 + N*crude_h2_est*totcoreg)^2 #compute weights again 
weights2_jk <- weights2[-remove_genes] 
weights3_jk <- weights3[-remove_genes] #serves to upweight tissue-specific genes
jk[chunk,] <- get_cov_alpha1alpha1_multitissue(alpha_z_jk,X_jk,N,weights1_jk,weights2_jk,weights3_jk)* N_tissuespecific_jk 
jk_sum[chunk,1] <- sum(as.numeric(jk[chunk,]))
print(jk[chunk,1])
jk_weights[chunk,1] <- sum(1/(weights1_jk*weights2_jk*weights3_jk))
} 

write.table(jk, file = paste0("TCSC_",trait,"_jk.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
write.table(jk_weights, file = paste0("TCSC_",trait,"_jkweights.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

#calculate weighted variance which we need to interpret the jackknife results
calc_weighted_var <- function(jk_colx){
ws <- jk_weights[,1]
wf <- rep(1,chunks)
x <- jk_colx
xm <- weighted.mean(x, ws)
myvar <- sum(ws *(x-xm)^2)*(sum(ws)/(sum(ws)^2-sum(ws^2)))  # Variance
return(myvar)
}
variance_mat[1:ncol(jk),3] <- sapply(1:ncol(jk), function(x) sqrt(calc_weighted_var(jk[,x]))*sqrt(length(which(jk_weights[,1] != 0))))
variance_mat[nrow(variance_mat),3] <- sqrt(calc_weighted_var(jk_sum[,1]))*sqrt(length(which(jk_weights[,1] != 0)))

variance_mat[,4] <- pnorm(q = 0, mean = as.numeric(variance_mat[,2]), sd = as.numeric(variance_mat[,3])) #added tissues to first column therefore matrix nonnumeric now. 
variance_mat[1:ncol(jk),5] <- p.adjust(as.numeric(variance_mat[1:ncol(jk),4]), method = "fdr")
variance_mat[nrow(variance_mat),5] <- NA
variance_mat[,6] <- as.numeric(variance_mat[,2])/sum(as.numeric(variance_mat[1:ncol(jk),2]))

colnames(variance_mat) <- c("Tissue","Variance","JK_SE","P","FDRP","Variance_divideby_Total")
write.table(variance_mat, file = paste0("TCSC_",trait,"_",ciswindow,"_cislocus.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

## Analysis to do in class 

#1. Which tissue has the greatest BrCa heritability in the variance_mat matrix?

#2. Plot the distribution of heritability estimates for this tissue across each jackknife block using the jk matrix. 

#3. Which jackknife blocks have significantly lower heritability compared to a normal distribution characterized by mean and std across all jackknife blocks? 

#4. What are the TWAS z-scores of the genes that are in this block? You should use the "a" matrix and the alpha_z variable; in both of which the order of genes is the same.   
#Column 4 of the "a" matrix is the tissue index, Column 5 of the "a" matrix is the jackknife block. 

#5. What are the top two genes? 

#6. Are these genes specifically associated with BrCa in this tissue or do they have high TWAS associations in other tissues? 

#7. Is the gene that is specifically cis-heritable in Exposed Skin tissue co-regulated with other genes? The matrix X contains co-regulation scores and the variable weights2 contains the sum across the rows of matrix X. 


