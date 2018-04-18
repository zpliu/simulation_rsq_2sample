args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
# library(argparse)
library(flexmix)

############################################# 
## Rscript script.R rsqPG gwas_thresh
## use a GWAS threshold;
## plot uvw
#############################################
## Parse argument from command line
# parser <- ArgumentParser()
# parser$add_argument("--sim_slt", default=6, type="integer", help="Select one simulation results to plot: 6 (default)")
# parser$add_argument("--rsq_PG", default=0.1, type="double", help="Variance explained by pleiotropic SNPs in both exposure and outcome: 0.1 (default)")
# parser$add_argument("--gwas_p_thresh", default=5e-8, type="double", help="GWAS P threshold to filter the SNPs: 5e-8 (default)")
# args <- parser$parse_args()


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

# rsqPG = args$rsq_PG
# j = args$sim_slt  # selected simulation
# GWAS_thresh = args$gwas_p_thresh

rsqPG = args[1]
j = 3  # selected simulation
GWAS_thresh = as.numeric(args[2])


# ---------------------------------
#  Read files
# ---------------------------------
sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"Ref",sample_ref,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)
bzx_dat = read.table(paste0(sim_para,"/BETA_XG.exp.xls.gz"), head=TRUE, sep="\t")
bzx_se_dat = read.table(paste0(sim_para,"/SE_XG.exp.xls.gz"), head=TRUE, sep="\t")
bzx_pval_dat = read.table(paste0(sim_para,"/PVAL_XG.exp.xls.gz"), head=TRUE, sep="\t")

bzy_dat = read.table(paste0(sim_para,"/BETA_YG.out.xls.gz"), head=TRUE, sep="\t")
bzy_se_dat = read.table(paste0(sim_para,"/SE_YG.out.xls.gz"), head=TRUE, sep="\t")
# bzy_pval_dat = read.table(paste0(sim_para,"/PVAL_YG.out.xls.gz"), head=TRUE, sep="\t")

remain_index <- which(bzx_pval_dat[,j] < GWAS_thresh)

bzx = bzx_dat[remain_index,j] 
bzx_se = bzx_se_dat[remain_index,j]
bzy = bzy_dat[remain_index,j]
bzy_se = bzy_se_dat[remain_index,j]


# ---------------------------------
# Try flexmix using remaining SNPs
# ---------------------------------
x = bzx/bzx_se
y = bzy/bzy_se

# snp_type <- c(rep("causal", snp_num), rep("pleio", pleio_snp_num))
snp_type <- ifelse(remain_index <= snp_num, "causal", "pleio")
data <- data.frame(x=as.numeric(x), class=as.factor(snp_type))
# stop(length(snp_type))
mo1 <- FLXMRglm(family = "gaussian")
mo2 <- FLXMRglm(family = "gaussian")
flexfit <- flexmix(y ~ x, data = data, k = 2, model = list(mo1, mo2))
# table(clusters(flexfit))

r1 <- length(snp_type[clusters(flexfit)==1])
r2 <- length(snp_type[clusters(flexfit)==2])

pleio_predict <- ifelse(r1<r2, 1, 2)
cls_res <- clusters(flexfit)

color_cls <- c()
for (i in 1:length(remain_index)){
	if (snp_type[i]=='causal' & cls_res[i]==pleio_predict){
		color_cls[i] <- 'causal, but wrongly predicted as pleio'  # False positive
	} else if (snp_type[i]=='causal' & cls_res[i]!=pleio_predict) { 
		color_cls[i] <- 'causal, and rightly predicted as causal'  # True negative
	} else if (snp_type[i]=='pleio' & cls_res[i]==pleio_predict) {
		color_cls[i] <- 'pleio, and rightly predicted as pleio'  # True positive
	} else if (snp_type[i]=='pleio' & cls_res[i]!=pleio_predict) {
		color_cls[i] <- 'pleio, but wrongly predicted as causal'  # False negative
	}
}


# ---------------------
# Plot results
# ---------------------
title_para = paste0("beta=",beta, " | maf=",maf, " | sim_selected=",j, 
					"\nsizeX=",sample_X, " | sizeY=",sample_Y, " | sizeRef=",sample_ref, 
					"\ncausalSNP=",snp_num," | pleioSNP=", pleio_snp_num, 
					"\nrsqGX=",rsqGX, " | rsqPG=",rsqPG, " | rsqXY=",rsqXY," | rsqC=",rsqC,
					"\nGWAS_P_thresh=", GWAS_thresh)

df1 <- data.frame(x = x, y = y)
p1 <- ggplot(df1, aes(x, y)) + 
	 geom_point(aes(shape=snp_type, color=color_cls), size=2) +
	 scale_x_continuous(expand = c(0, 0), limits=c(min(x), max(x))) + 
	 scale_y_continuous(expand = c(0, 0), limits=c(min(y), max(y))) +
	 scale_shape_manual(values=c(16, 17)) +
	 scale_colour_manual(values=c("grey65", "red", "limegreen", "mediumorchid3")) +
	 labs(x = "Z_exposure", y = "Z_outcome") +
	 geom_hline(yintercept=0) + geom_vline(xintercept=0) +
	 ggtitle(title_para)

ggsave(paste0("flexmix_res/plot_uvw_",rsqPG,".pdf"), plot = p1, 
	dpi = 300, units="cm")


