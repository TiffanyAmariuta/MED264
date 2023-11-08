#For the following SNPs on chr10, positions: 21849769 21856475 21884471 21983960 21989245 22062729
#Find which IMPACT annotations out of the ones in this directory have cell-type-specific regulation potential. 
#Hypothesize about which variants might actually be causal


#IMPACT annotation key: 
#95: GSM798404	ESR1	BREAST	Tumour tissues	95	PRJNA147213_GSM798404_ESR1_peaks.narrowPeak	Tumor	1
#135: GSM951859	HSF1	BREAST	BPLER	135	PRJNA169333_GSM951859_HSF1_peaks.narrowPeak	Breast	1
#133: GSM951855	HSF1	BREAST	BPE	133	PRJNA169333_GSM951855_HSF1_peaks.narrowPeak	Breast	1
#30: GSM599839	HSF1	ADIPOCYTE	Adipocytes	30	PRJNA132733_GSM599839_HSF1_peaks.narrowPeak	Adipocytes	1
#41: GSM634630	PPARG	ADIPOCYTE	SGBS	41	PRJNA135749_GSM634630_PPARG_peaks.narrowPeak	Adipocytes	1
#214: GSM1139038	SOX2	LUNG	HCC95	214	PRJNA202445_GSM1139038_SOX2_peaks.narrowPeak	Lung	1
#460: GSM1010723	JUND	LUNG	A549	460	PRJNA63447_GSM1010723_JUND_peaks.narrowPeak	Lung	1

library(data.table)
positions <- c(21849769,21856475,21884471,21983960,21989245,22062729)
files <- list.files(pattern = "IMPACT_")
matrix_res <- matrix(0,length(positions), length(files))
#iterate through files and positions to find the IMPACT score from file i for the SNP at position j


rownames(matrix_res) <- paste0("SNP:",positions)
colnames(matrix_res) <- c("ESR1:BREAST", "HSF1:BREAST","HSF1:BREAST","HSF1:Adipocyte","PPARG:Adipocyte","SOX2:Lung","JUND:Lung")

