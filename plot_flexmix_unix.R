library(flexmix)
library(ggplot2)
library(gridExtra)
library(foreach)
# library(doParallel) # For both windows and unix
library(doMC)
registerDoMC(10)

# rm(list=ls())

get_outlier_SNP <- function(sim_num, rsqPG) {
	num_pleio_total = c()
	num_causal_total = c()
	res_out <- foreach (i=1:sim_num, .combine=rbind) %dopar% {
	  
	  # sim_col = which(single_all_dat$SIM==i)
	  # snp_num = length(sim_col) - 2
	  # bxy_alldat = single_all_dat[sim_col,]$b
	  # bxy = bxy_alldat[1:snp_num]
	  
	  bzx = bzx_dat[,i]
	  # bzx_se = bzx_se_dat[,i]^2
	  bzx_se = bzx_se_dat[,i]
	  bzy = bzy_dat[,i]
	  # bzy_se = bzy_se_dat[,i]^2
	  bzy_se = bzy_se_dat[,i]
	  
	  x = bzx/bzx_se
	  y = bzy/bzy_se
	  
	  data <- data.frame(x=as.numeric(x), class=as.factor(class))
	  
	  mo1 <- FLXMRglm(family = "gaussian")
	  mo2 <- FLXMRglm(family = "gaussian")
	  flexfit <- flexmix(y ~ x, data = data, k = 2, model = list(mo1, mo2))
	  # flexfit <- flexmix(y ~ x, data = data, k = 2, model = mo1)
	  # table(clusters(flexfit))
	  
	  r1 <- length(class[clusters(flexfit)==1])
	  r2 <- length(class[clusters(flexfit)==2])
	  
	  num_total = ifelse(r1<r2, r1, r2)
	  num_pleio = ifelse(r1<r2, 
	             length(grep("p", class[clusters(flexfit)==1])), 
	             length(grep("p", class[clusters(flexfit)==2])))
	  num_causal = num_total - num_pleio
	  
	  # num_pleio_total = c(num_pleio_total, num_pleio)
	  # num_causal_total = c(num_causal_total, num_causal)
	  c(num_pleio, num_causal) 
	}
	# out <- c(mean(num_pleio_total), sd(num_pleio_total), mean(num_causal_total), sd(num_causal_total))
	return(res_out)
}



##############

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
	# single_all_dat = read.table(paste0("statSim_",sim_para,"/X1_1.0_Y1_1.0/MR_G1_X_Y/G1_X_Y.single.txt"), head=TRUE, sep="\t")

	bzx_dat = read.table(paste0(sim_para,"/BETA_XG.exp.xls.gz"), head=TRUE, sep="\t")
	bzx_se_dat = read.table(paste0(sim_para,"/SE_XG.exp.xls.gz"), head=TRUE, sep="\t")

	bzy_dat = read.table(paste0(sim_para,"/BETA_YG.out.xls.gz"), head=TRUE, sep="\t")
	bzy_se_dat = read.table(paste0(sim_para,"/SE_YG.out.xls.gz"), head=TRUE, sep="\t")

	class <- c(rep('c', snp_num), rep('p', pleio_snp_num))

	res_out_tmp <- get_outlier_SNP(sim_num, rsqPG) #
	colnames(res_out_tmp) <- c("true_pleio", "false_pleio")

	power_ave <- c(power_ave, mean(res_out_tmp[,1])/pleio_snp_num)
	fp_ave <- c(fp_ave, mean(res_out_tmp[,2])/snp_num)
	# power_sd <- c(power_sd, res_out_tmp[2]/pleio_snp_num)
	# fp_sd <- c(fp_sd, res_out_tmp[4]/snp_num)
	res_out[[rsqPG]] <- res_out_tmp

	message(paste(rsqPG, mean(res_out_tmp[,1])/pleio_snp_num, mean(res_out_tmp[,2])/snp_num))
	# save(file=paste0("flexmix_res/flexmix_rsq", rsqPG, ".RData"), res_out_tmp)
}

# stop("DONE")
# out = data.frame(power_ave, power_sd, fp_ave, fp_sd)
out = data.frame(power_ave, fp_ave)
save(file="flexmix_res/flexmix_percentage.RData", out)
save(file="flexmix_res/flexmix_num.RData", res_out)




## Plot barplot

title_para = paste0("beta=",beta, " | maf=",maf, " | sim=",sim_num, "\nsizeX=",sample_X, " | sizeY=",sample_Y, " | sizeRef=",sample_ref, "\ncausalSNP=",snp_num," | pleioSNP=", pleio_snp_num, "\nrsqGX=",rsqGX," | rsqXY=",rsqXY," | rsqC=",rsqC)

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


############################################
# v <- .05
# causal_num = 100
# pleio_num = 20

# u0 <- rnorm(causal_num, 0, 1)
# e0 <- rnorm(causal_num, 0, 1)
# y0 <- sum(u0)*v + e0

# u1 <- rnorm(pleio_num, 0, 1)
# w <- rnorm(pleio_num, 0, 1)
# e1 <- rnorm(pleio_num, 0, 1)
# y1 <- sum(u1)*v + sum(w) + e1

# # a <- 0.2
# # y = a*y0 + (1-a)*y1

# plot(u0, y0, xlim=c(-3,3), ylim=c(-10,10))
# plot(u1, y1)
# u <- c(u0, u1)
# y <- c(y0, y1) 
# plot(u, y, xlim=c(-3,3), ylim=c(-10,10))

# x <- c(u0, u1)
# class <- c(rep('c', causal_num), rep('p', pleio_num))
# data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))

# mo1 <- FLXMRglm(family = "gaussian")
# mo2 <- FLXMRglm(family = "gaussian")
# flexfit <- flexmix(y ~ u, data = data, k = 2, model = list(mo1, mo2))
# print(table(clusters(flexfit), data$class))

###############
