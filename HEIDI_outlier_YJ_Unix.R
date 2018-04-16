# NOTE: this was edited from GSMR original script to
#  just keep the pleiotropic SNPs
library(foreach)
library(doMC)

registerDoMC(10)

# rm(list=ls())

eps = 1e-6;
# ************************************************** #
#       variance-covariane matrix of bXY             #
# ************************************************** #
cov_bXY <- function(bzx, bzx_se, bzy, bzy_se, ldrho) {
  bXY = bzy / bzx
  zscoreZX = bzx / bzx_se
  nsnp = dim(ldrho)[1]
  covbXY = diag(nsnp)
  if(nsnp>1) {
    zszxij = zscoreZX%*%t(zscoreZX)
    bxyij = bXY%*%t(bXY)
    sezyij = bzy_se%*%t(bzy_se)
    bzxij = bzx%*%t(bzx)
    covbXY = ldrho*sezyij/bzxij + ldrho*bxyij/zszxij 
  }
  return(covbXY)
}

# ************************************************** #
#                 HEIDI test                         #
# ************************************************** #
#' @importFrom survey pchisqsum
heidi <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, 
                  heidi_thresh = pchisq(10, 1, lower.tail=F),
                  nSNPs_thresh=10, maxid=integer(0) ) {

  ## add by myself
  # bzx = bzx[c(refid,i)]
  # bzx_se = bzx_se[c(refid,i)]
  # bzx_pval = bzx_pval[c(refid,i)]
  # bzy = bzy[c(refid,i)] 
  # bzy_se = bzy_se[c(refid,i)]
  # ldrho = ldrho[c(refid,i), c(refid,i)]
  # heidi_thresh = 0.01
  # nSNPs_thresh = 2
  # maxid = 1
  ## add end

  remain_index = seq(1, length(bzx))
  m = length(remain_index)
  
  # recaculate the zscore_ZX for vector has been changed
  zscore_ZX = bzx / bzx_se
  bXY = bzy / bzx
  seSMR = sqrt( (bzy_se^2*bzx^2 + bzx_se^2*bzy^2) / bzx^4 )
  # remap the maxid according to filtered data
  maxid = which(remain_index==maxid)
  if(length(maxid) != 1) {
    stop("The top SNP for the HEIDI analysis is missing.")
  }
  # diff = bXY_top - bXY_-i
  dev = bXY[maxid] - bXY[-maxid]
  # v matrix
  covbXY = as.matrix(cov_bXY(bzx, bzx_se, bzy, bzy_se, ldrho))
  tmp1 = diag(covbXY)[maxid]
  tmp2 = covbXY[-maxid, -maxid]
  tmp3 = covbXY[maxid, -maxid]
  vdev = tmp1 + tmp2 - tmp3
  vdev = t(t(vdev) - tmp3)
  diag(vdev) = diag(vdev) + eps
  # variance of diff
  vardev = diag(covbXY)[-maxid] + tmp1 - 2*tmp3
  vardev = vardev + eps
  chisq_dev = dev^2 / vardev
  if(m>2) {
    # correlation matrix of diff
    corr_dev = diag(m-1)
    # more than 2 instruments
    for( i in 1 : (m-2) ) {
      for( j in (i+1) : (m-1) ) {
        corr_dev[i,j] = corr_dev[j,i] =
          vdev[i,j] / sqrt(vdev[i,i]*vdev[j,j])
      }
    }
    
    # estimate the p value
    lambda = eigen(corr_dev, symmetric=TRUE, only.values=TRUE)$values
    t = length(lambda)
    pHet = pchisqsum(sum(chisq_dev)+eps, df=rep(1,t), a=lambda, method="sadd", lower.tail=F)
    
  } else {
    pHet = pchisq(chisq_dev, 1, lower.tail=F)
  }
  
  return(list(pheidi=pHet, nsnps=m, used_index=remain_index))
}

# ************************************************** #
#          Iterations for HEIDI-ouliter              #
# ************************************************** #

heidi_outlier_iter <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh, heidi_thresh, remain_index) {
  # remove pleiotropic instruments
  m = length(remain_index)
  # grab reference instrument
  bxy = bzy/bzx
  bxy_q = quantile(bxy, probs = seq(0, 1, 0.2))
  min_bxy = as.numeric(bxy_q[3]); max_bxy = as.numeric(bxy_q[4]);
  slctindx = which(bxy <= max_bxy & bxy >= min_bxy)
  if(length(slctindx)==0) {
    stop("The top SNP for the HEIDI-outlier analysis is missing. None SNPs are in the third quintile of the distribution of bxy.");
  }
  refid = slctindx[which.min(bzx_pval[slctindx])]
  pheidi = as.numeric()
  for( i in 1 : m ) {
    if( i==refid ) next
    heidi_result = heidi(bzx[c(refid,i)], bzx_se[c(refid,i)], bzx_pval[c(refid,i)],
                         bzy[c(refid,i)], bzy_se[c(refid,i)],
                         ldrho[c(refid,i), c(refid,i)],
                         heidi_thresh, 2, 1)
    pheidi[i] = as.numeric(heidi_result$pheidi)
  }
  remain_index = remain_index[sort(c(refid,which(pheidi>=heidi_thresh)))]
  return(remain_index)
}


# ******************************************************** #
#   HEIDI-ouliter analysis: return true pleio percentage   #
#  This subfunction is edited by myself
# ******************************************************** #
heidi_outlier_modified <- function(bzx_dat, bzx_se_dat, bzx_pval_dat, bzy, bzy_se, ldrho, causal_snp_num, pleio_snp_num) {
	heidi_outlier_thresh = 0.01
 #  true_pleio_num = list()
	# false_pleio_num = list()

	num <- foreach (j=1:sim_num, .combine=rbind) %dopar% {
		# cat(j,"\r")
		bzx = bzx_dat[,j]
		bzx_se = bzx_se_dat[,j]
		bzx_pval = bzx_pval_dat[,j]

		bzy = bzy_dat[,j]
		bzy_se = bzy_se_dat[,j]

		total_snp_num = length(rownames(bzx_dat))

		gwas_thresh = 1.0
		# heidi_outlier_thresh = outlier_pval_thresh
		remain_index = c(1:total_snp_num)

		heidi_causal_snp <- heidi_outlier_iter(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh, heidi_outlier_thresh, remain_index)
		# true_pleio_num <- c(true_pleio_num, (pleio_snp_num - length(which(heidi_causal_snp > causal_snp_num))))
		true_pleio_num <- pleio_snp_num - length(which(heidi_causal_snp > causal_snp_num))
    false_pleio_num <- causal_snp_num - length(which(heidi_causal_snp <= causal_snp_num))

    c(true_pleio_num, false_pleio_num)
	}

	return(num)
}


#######################

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
fp_ave = c()
res_out = list()
for (rsqPG in rsqPG_list) {
  message(rsqPG)
	sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"Ref",sample_ref,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)

	bzx_se_dat = read.table(paste0(sim_para,"/SE_XG.exp.xls.gz"), sep="\t", header=T)
	bzx_dat = read.table(paste0(sim_para,"/BETA_XG.exp.xls.gz"), sep="\t", header=T)
	bzx_pval_dat = read.table(paste0(sim_para,"/PVAL_XG.exp.xls.gz"), sep="\t", header=T)

	bzy_se_dat = read.table(paste0(sim_para,"/SE_YG.out.xls.gz"), sep="\t", header=T)
	bzy_dat = read.table(paste0(sim_para,"/BETA_YG.out.xls.gz"), sep="\t", header=T)

	ldrho = read.table(paste0(sim_para,"/genoG_LD.xls.gz"), sep="\t", header=T)
	
	res_out_tmp <- heidi_outlier_modified(bzx_dat, bzx_se_dat, bzx_pval_dat, bzy_dat, bzy_se_dat, ldrho, snp_num, pleio_snp_num)
  colnames(res_out_tmp) <- c("true_pleio", "false_pleio")

	power_ave <- c(power_ave, mean(as.numeric(res_out_tmp[,1]))/pleio_snp_num)
  fp_ave <- c(fp_ave, mean(as.numeric(res_out_tmp[,2]))/snp_num)
  res_out[[rsqPG]] <- res_out_tmp
  # save(file=paste0("flexmix_res/HEIDI_rsq", rsqPG, ".RData"), res_out)
}

out = data.frame(power_ave, fp_ave)
save(file="flexmix_res/HEIDI_percetage.RData", out)
save(file="flexmix_res/HEIDI_num.RData", res_out)



# 3rd qunitile
# 0.03: 0.28108
# 0.04: 0.2081
# 0.05: 0.25078
# 0.06: 0.3853
# 0.07: 0.36474
# 0.08: 0.34662
# 0.09: 0.33374
# 0.1: 0.4309

# min P value
# 0.03: 0.28136
# 0.04: 0.20826
# 0.05: 0.25102
# 0.06: 0.38538
# 0.07: 0.39532
# 0.08: 0.34714
# 0.09: 0.33434
# 0.1: 0.43078
