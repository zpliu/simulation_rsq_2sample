library(gsmr)
library(ggplot2)
library(foreach)
library(doMC)

registerDoMC(12)


## Calculate outlier
cal_pleio_outlier <- function (exposure_beta_file, outcome_beta_file, exposure_se_file, outcome_se_file, exposure_pval_file, outcome_pval_file, ldmat_file)
{
	exposureBETA = read.table(exposure_beta_file, sep="\t", header=T)
	outcomeBETA = read.table(exposure_beta_file, sep="\t", header=T)
	exposureSE = read.table(exposure_se_file, sep="\t", header=T)
	outcomeSE = read.table(exposure_se_file, sep="\t", header=T)
	exposurePVAL = read.table(exposure_pval_file, sep="\t", header=T)
	outcomePVAL = read.table(exposure_pval_file, sep="\t", header=T)
	ldrho = read.table(ldmat_file, sep="\t", header=T)

	sim_num = length(colnames(exposureBETA))
	snp_num = length(rownames(exposureBETA))
	# ldrho = matrix(0, nrow=snp_num, ncol=snp_num)
	# colnames(ldrho) = rownames(exposureBETA)
	# rownames(ldrho) = rownames(exposureBETA)
	snp_coeff_id = rownames(exposureBETA)

	n_ref = 5000    # Sample size of the reference sample for LD
	nsnps_thresh = 5   # the minimum number of instruments required for the GSMR analysis
	gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
	heidi_outlier_thresh = 0.01    # HEIDI-outlier threshold
	heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
	ld_r2_thresh = 0.1    # LD r2 threshold to remove SNPs in high LD
	ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments: default 0.05

	outlier = list()
	outlier <- foreach (i=1:sim_num, .packages='gsmr') %dopar% {
		bzx = exposureBETA[,i]    # SNP effects on the risk factor
		bzx_se = exposureSE[,i]    # standard errors of bzx
		bzx_pval = exposurePVAL[,i]   # p-values for bzx
		bzy = outcomeBETA[,i]    # SNP effects on the disease
		bzy_se = outcomeSE[,i]    # standard errors of bzy
		bzy_pval = outcomePVAL[,i]    # p-values for bzy
		# gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, heidi_outlier_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)    # GSMR analysis 
		heidi_results = heidi_outlier(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snp_coeff_id, n_ref, gwas_thresh, heidi_outlier_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)
		length(heidi_results$pleio_snps)
	}
	return(as.numeric(outlier))
}


## ------------------
##  Main function
## ------------------
beta = "1.0"
maf = '0.010.5'
snp_num = 500
sample_X = 50000
sample_Y = 20000
sample_ref = 5000
sim_num = 1000

rsqGX = 0.4
rsqXY = 0.05
rsqC = 0.2

pleio_ratio = 0.1
# pleio_snp_num = round(snp_num*as.numeric(pleio_ratio), digits=0)

rsqPG_list = c("0.03","0.04","0.05","0.06","0.07","0.08","0.09","0.1")
# rsqPG_list = c("0.1")
outlier_num = c()
# num_sd = c()

for (rsqPG in rsqPG_list) {
	message(rsqPG)
	# ratio = format(ratioList[j], nsmall=1)
	sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"Ref",sample_ref,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)

	outlier <- cal_pleio_outlier(
		exposure_beta_file = paste0(sim_para,"/BETA_XG.exp.xls.gz"),
		outcome_beta_file = paste0(sim_para,"/BETA_YG.out.xls.gz"),
		exposure_se_file = paste0(sim_para,"/SE_XG.exp.xls.gz"),
		outcome_se_file = paste0(sim_para,"/SE_YG.out.xls.gz"),
		exposure_pval_file = paste0(sim_para,"/PVAL_XG.exp.xls.gz"),
		outcome_pval_file = paste0(sim_para,"/PVAL_YG.out.xls.gz"),
		ldmat_file = paste0(sim_para,"/genoG_LD.xls.gz")
	)
	# num_mean = c(num_mean, mean(outlier))
	# num_sd[j] = sd(outlier)
	message(paste("Result of", rsqPG, ":", mean(outlier)))
	# save(file=paste0("gsmr/rsqPG_", rsqPG, ".RData"), outlier)
}


# marker = rep("X_Y", times=11)

# testData = data.frame(c(ratioList,ratioList), marker, num_mean, num_sd)
# colnames(testData) = c("ratio","dir","mean","sd")


# ## Plot grouped barplot
# limits <- aes(ymax=testData$mean+testData$sd,
#               ymin=testData$mean-testData$sd)

# p <- ggplot(data=testData, aes(x=factor(testData$ratio), y=testData$mean, fill=factor(testData$dir)))

# p + geom_bar(stat = "identity", position = 'dodge', width = 0.5) +
#   geom_errorbar(limits, position = 'dodge', width = 0.5) +
#   labs(x = "ratio", y = "Numbers") +
#   ggtitle("Number of outlier SNP") +
#   scale_fill_discrete(name = "Direction")