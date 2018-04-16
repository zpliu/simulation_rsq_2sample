suppressMessages(library(TwoSampleMR))
library(ggplot2)
library(argparse)
library(foreach)
# library(doParallel)
library(doMC)
library(data.table)

###########################################################
## Run TwoSampleMR for all simluations, output to one file
###########################################################

## Parse argument from command line
parser <- ArgumentParser()
parser$add_argument("--beta_exp", help="Beta values of genotype and exposure")
parser$add_argument("--beta_out", help="Beta values of genotype and outcome")
parser$add_argument("--se_exp", help="se values of genotype and exposure")
parser$add_argument("--se_out", help="se values of genotype and outcome")
parser$add_argument("--pval_exp", help="pval of genotype and exposure")
parser$add_argument("--pval_out", help="pval of genotype and outcome")
parser$add_argument("--sample_X", type="integer", help="sample size for exposure")
parser$add_argument("--sample_Y", type="integer", help="sample size for outcome")
parser$add_argument("--sim_num", type="double", help="simulation number")
parser$add_argument("--outprefix", help="output prefix, including the DIR")
parser$add_argument("--thread", type="integer", help="multi-thread number")
parser$add_argument("--plot", action="store_true", help="whether or not plot the results: FALSE (default)")
args <- parser$parse_args()

registerDoMC(args$thread)  #change the 2 to your number of CPU cores

message("MR: Reading files")
beta_exp = fread(paste("gzip -dc", args$beta_exp), header=FALSE, skip=1, sep="\t")
beta_out = fread(paste("gzip -dc", args$beta_out), header=FALSE, skip=1, sep="\t")
se_exp = fread(paste("gzip -dc", args$se_exp), header=FALSE, skip=1, sep="\t")
se_out = fread(paste("gzip -dc", args$se_out), header=FALSE, skip=1, sep="\t")
pval_exp = fread(paste("gzip -dc", args$pval_exp), header=FALSE, skip=1, sep="\t")
pval_out = fread(paste("gzip -dc", args$pval_out), header=FALSE, skip=1, sep="\t")

## Using an existing data frame
SNP = beta_exp[[1]]
snp_num = length(SNP)
sim_num = args$sim_num

message("MR: Prepare data")
dat <- foreach (j=1:sim_num) %dopar% {
	exp_dat <- data.frame(
		SNP,
		beta = beta_exp[[j+1]],
		se = se_exp[[j+1]],
		effect_allele = rep("G", times=snp_num),
		other_allele = rep("C", times=snp_num),
		pval = pval_exp[[j+1]],
		samplesize = args$sample_X
	)

	out_dat <- data.frame(
		SNP,
		beta = beta_out[[j+1]],
		se = se_out[[j+1]],
		effect_allele = rep("G", times=snp_num),
		other_allele = rep("C", times=snp_num),
		pval = pval_out[[j+1]],
		samplesize = args$sample_Y
	)

	exp_dat <- suppressMessages(format_data(exp_dat, type="exposure"))
	out_dat <- suppressMessages(format_data(out_dat, type="outcome"))

	suppressMessages(harmonise_data(
			exposure_dat = exp_dat, 
			outcome_dat = out_dat,
			action = 1 # 1=Assume all ref alleles are on the positive strand (REV example in COG, ANA and SCZ: rs3131966; rs1498232)
		))

}


### Method list
# methodlist <- c(
# 			"mr_wald_ratio",
# 			"mr_two_sample_ml",
# 			"mr_egger_regression",
# 			"mr_egger_regression_bootstrap",
# 			"mr_simple_median",
# 			"mr_weighted_median",
# 			"mr_penalised_weighted_median"  "mr_ivw",
# 			"mr_simple_mode",
# 			"mr_weighted_mode",
# 			"mr_weighted_mode_nome",
# 			"mr_simple_mode_nome"
# )

message("MR: running")
out_res <- foreach (i=1:sim_num, .combine=c) %dopar% {
	resIVW <- suppressMessages(mr(dat[[i]], method_list="mr_ivw"))
	resEggerCausal <- suppressMessages(mr(dat[[i]], method_list="mr_egger_regression") )
	# invisible(capture.output(resEggerBootstrapCausal <- suppressMessages(mr(dat[[i]], method_list="mr_egger_regression_bootstrap") )))
	resEggerIntercept <- suppressMessages(mr_pleiotropy_test(dat[[i]]) )
	# ## heterogeneity test
	resHET <- suppressMessages(mr_heterogeneity(dat[[i]], method_list="mr_ivw"))
	# ## Directionality test
	Dir_test <- suppressMessages(directionality_test(dat[[i]])$correct_causal_direction)
	Dir_steiger <- suppressMessages(directionality_test(dat[[i]])$steiger_pval)


	paste(i,
		snp_num,
		resIVW$b,
		resIVW$se,
		resIVW$pval,
		resEggerCausal$b,
		resEggerCausal$se,
		resEggerCausal$pval,
		# resEggerBootstrapCausal$b,
		# resEggerBootstrapCausal$se,
		# resEggerBootstrapCausal$pval,
		resEggerIntercept$egger_intercept,
		resEggerIntercept$se,
		resEggerIntercept$pval,
		resHET$Q,
		resHET$Q_pval,
		Dir_test,
		Dir_steiger,
		sep="\t")
}

message("MR: single SNP analysis")
out_single <- foreach (i=1:sim_num, .combine=rbind) %dopar% {
	res <- mr_singlesnp(dat[[i]])
	as.matrix(cbind(rep(i,times=snp_num+2), res[c(6,7,8,9)]))
	# rbind(res_single, res)
}


## ------------------
## Text Output
## ------------------
out_header <- paste("SIM", 
				"N_SNP",
				"IVW_beta",
				"IVW_se",
				"IVW_Pval",
				"EggerCausal_beta",
				"EggerCausal_se",
				"EggerCausal_Pval",
				# "EggerBootstrapCausal_beta",
				# "EggerBootstrapCausal_se",
				# "EggerBootstrapCausal_Pval",
				"EggerIntercept",
				"EggerIntercept_se",
				"EggerIntercept_Pval",
				"HET_IVW_Q",
				"HET_IVW_Q_Pval",
				"Dir_test",
				"Dir_steiger_Pval",
				sep="\t")
out <- c(out_header, out_res)

out_single_header = c("SIM", "SNP", "b", "se", "P")
message("MR: Output results")
write.table(file=paste0(args$outprefix,'.MR.txt'), out, quote=F, row.names=F, col.names=F)
write.table(file=paste0(args$outprefix,'.single.txt'), out_single, quote=F, row.names=F, col.names=out_single_header, sep="\t")


if (args$plot){
	message("HAHa plot")
# 	### Plot results
	# plot_method <- 'mr_ivw'
# 	res <- mr(dat, method_list=plot_method)

# 	### Scatter plot
# 	p1 <- mr_scatter_plot(res, dat)
# 	# length(p1)
# 	# p1[[1]]

# 	### Forest plot
# 	p2 <- mr_forest_plot(res_single)
# 	# p2[[1]]

# 	### Funnel plot
# 	p4 <- mr_funnel_plot(res_single)
# 	# p4[[1]]

# 	### Leave-one-out plot
# 	# res_loo <- mr_leaveoneout(dat)
# 	# p3 <- mr_leaveoneout_plot(res_loo)
# 	# p3[[1]]

# 	### save plot
# 	ggsave(p1[[1]], file=paste0(outprefix, ".scatter.pdf"), width=4, height=4)
# 	ggsave(p2[[1]], file=paste0(outprefix,".forest.pdf"), width=4, height=4)
# 	ggsave(p4[[1]], file=paste0(outprefix,".funnel.pdf"), width=4, height=4)
}