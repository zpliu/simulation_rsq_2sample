library(argparse)
library(flexmix)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doMC)

###########################################
# use GWAS threshold to filter SNPs
###########################################
## Parse argument from command line
parser <- ArgumentParser()
parser$add_argument("--thread", type="integer", help="multi-thread number")
parser$add_argument("--gwas_thresh", default=5e-8, type="double", help="GWAS P threshold to filter the SNPs: 5e-8 (default)")
args <- parser$parse_args()

registerDoMC(args$thread)
gwas_thresh <- args$gwas_thresh



#--------------------------------------------------
# Run flexmix
#--------------------------------------------------
get_outlier_SNP <- function(sim_num, rsqPG) {
	num_pleio_total = c()
	num_causal_total = c()
	res_out <- foreach (j=1:sim_num, .combine=rbind) %dopar% {

	 	## Filter SNPs by GWAS P-val
		remain_index <- which(bzx_pval_dat[,j] < gwas_thresh)	  

	 	bzx = bzx_dat[remain_index,j]
	 	bzx_se = bzx_se_dat[remain_index,j]
	 	bzy = bzy_dat[remain_index,j]
	 	bzy_se = bzy_se_dat[remain_index,j]
	 	
	 	x = bzx/bzx_se
	 	y = bzy/bzy_se
	 	
	 	snp_type <- ifelse(remain_index <= snp_num, "c", "p")
	 	data <- data.frame(x=as.numeric(x), class=as.factor(snp_type))
	 	remain_causal_num <- length(which(remain_index <= snp_num))
	 	remain_pleio_num <- length(which(remain_index > snp_num))
	 	
	 	mo1 <- FLXMRglm(family = "gaussian")
	 	mo2 <- FLXMRglm(family = "gaussian")
	 	flexfit <- flexmix(y ~ x, data = data, k = 2, model = list(mo1, mo2))
	 	# table(clusters(flexfit))
	 	
	 	r1 <- length(snp_type[clusters(flexfit)==1])
	 	r2 <- length(snp_type[clusters(flexfit)==2])

	 	## Prediction after mixture model	 	
	 	predict_total_num <- ifelse(r1<r2, r1, r2)
	 	predict_true_pleio_num <- ifelse(r1<r2, 
				 	           length(grep("p", snp_type[clusters(flexfit)==1])), 
				 	           length(grep("p", snp_type[clusters(flexfit)==2])))
	 	predict_false_pleio_num <- predict_total_num - predict_true_pleio_num  # i.e. actually causal SNPs, but predicted as pleio

	 	c(predict_true_pleio_num, remain_pleio_num, predict_false_pleio_num, remain_causal_num) 
	}
	return(res_out)
}



#-------------------------
# simulation parameters
#-------------------------
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
pleio_snp_num = round(snp_num*as.numeric(pleio_ratio), digits=0)

rsqPG_list = c("0.03","0.04","0.05","0.06","0.07","0.08","0.09","0.1")
# rsqPG_list = c("0.1")

power_ave = c()
power_sd = c()
fp_ave = c()
fp_sd = c()
res_out = list()
for (rsqPG in rsqPG_list){
	message(rsqPG)
	sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"Ref",sample_ref,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)

	bzx_dat = read.table(paste0(sim_para,"/BETA_XG.exp.xls.gz"), head=TRUE, sep="\t")
	bzx_se_dat = read.table(paste0(sim_para,"/SE_XG.exp.xls.gz"), head=TRUE, sep="\t")
	bzx_pval_dat = read.table(paste0(sim_para,"/PVAL_XG.exp.xls.gz"), head=TRUE, sep="\t")

	bzy_dat = read.table(paste0(sim_para,"/BETA_YG.out.xls.gz"), head=TRUE, sep="\t")
	bzy_se_dat = read.table(paste0(sim_para,"/SE_YG.out.xls.gz"), head=TRUE, sep="\t")

	res_out_tmp <- get_outlier_SNP(sim_num, rsqPG) 
	colnames(res_out_tmp) <- c("true_pleio", "remain_pleio", "false_pleio", "remain_causal")

	power_ave <- c(power_ave, mean(res_out_tmp[,1]/res_out_tmp[,2]))
	fp_ave <- c(fp_ave, mean(res_out_tmp[,3]/res_out_tmp[,4]))

	res_out[[rsqPG]] <- res_out_tmp

	message(paste(rsqPG, mean(res_out_tmp[,1]/res_out_tmp[,2]), mean(res_out_tmp[,3]/res_out_tmp[,4])))
}

# out = data.frame(power_ave, power_sd, fp_ave, fp_sd)
out = data.frame(rsqPG_list, power_ave, fp_ave)
save(file="flexmix_res/flexmix_percentage.RData", out)
save(file="flexmix_res/flexmix_num.RData", res_out)


# ---------------------
#   Plot results
# ---------------------
title_para = paste0("beta=",beta, " | maf=",maf, " | sim=",sim_num, 
					"\nsizeX=",sample_X, " | sizeY=",sample_Y, " | sizeRef=",sample_ref, 
					"\ncausalSNP=",snp_num," | pleioSNP=", pleio_snp_num, 
					"\nrsqGX=",rsqGX," | rsqXY=",rsqXY," | rsqC=",rsqC,
					"\nGWAS_thresh=", gwas_thresh)

# limits <- aes(ymax=out$power_ave+out$power_sd,
#               ymin=out$power_ave-out$power_sd)
p1 <- ggplot(data=out, aes(x=factor(rsqPG_list), y=out$power_ave)) +
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.1), limits=c(0,1)) +
  # geom_errorbar(limits, position = 'dodge', width = 0.2) + 
  labs(x = "Variance explained by pleiotropic SNPs", y = "true positive rate") +
  ggtitle(title_para)
  # scale_fill_discrete(name = "Direction")


# limits <- aes(ymax=out$fp_ave+out$fp_sd,
#               ymin=out$fp_ave-out$fp_sd)
p2 <- ggplot(data=out, aes(x=factor(rsqPG_list), y=out$fp_ave)) +
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.1), limits=c(0,1)) +
  # geom_errorbar(limits, position = 'dodge', width = 0.2) + 
  labs(x = "Variance explained by pleiotropic SNPs", y = "false positive rate")
  # ggtitle(title_para)
  # scale_fill_discrete(name = "Direction")


ggsave("flexmix_res/plot_flexmix.pdf", plot = grid.arrange(p1, p2, nrow=2), 
	dpi = 300, units="cm")



