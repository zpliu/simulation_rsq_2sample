# options(error=recover)
library(argparse)
library(foreach)
# library(doParallel)  # both Windows and Unix
library(doMC) # Unix only

#########################################################################
### MR Simulation: control Rsq 
###	Ref: YJ Nature Communication GSMR paper 

### Current function in this script:
###   + Simulate in Matrix format.
###   + Parallel computation;
###   + Include pleiotropic SNPs;
###     + Control Rsq(causal), Rsq(pleiotropy), Rsq(C) in exposure;
###     + Control Rsq(XY), Rsq(pleiotropy) in outcome
###   + Two-sampel MR context:
###     +Two different sample sizes for exposure and outcome
###     + Build a LD matrix using reference samples
###   + Standardiztion:
###     + Use standardised genotype to generate exposure, then use non-standardised exposure to generate outcome; 
###     + Use standardised genptype and standardised phenotype to calculate association; 
###   + Output files:
###     + Compress outputs as gz files
###     + Results in output files for genotpye and phenotypes are NOT standardised (phenoype values are limited to 4 decimal)
###   + random seed for genotypes and effects, i.e. as long as the sample size is the same during the simulation, the genotypes will be the same.
###       This is for comparison among different rsqPG.

### X = g + pgX + c + eX, where
###   g = Wu, W is the causal SNP matrix, u is the effect (i.e. gamma in the script);
###   pgX = W(p)u(p), W(p) is the pleiotropic SNP matrix, u(p) is the effect;
###   SNP follows ~ Binomial(2,f), where f ~ Uniform(MAF_min, MAF_max);
###   each SNP effect u~ N(0, 1), i.e. both causal and pleiotropic SNP
###   c is the latent non-genetic confounding variable, c~N(0,var(c)), var(c)=var(g)*Rsq(c)/Rsq(g)
###   eX~N(0, var(eX)), var(eX) = var(g + pgX + c)*[1/(Rsq(g) + Rsq(pgX) + Rsq(c)) - 1]

### Y = bX + pgY + c + eY, where 
###   b is the causal effect (i.e. beta in the script)
###   pgY is the effect of pleiotropic SNP on outcome
###   c is the same latent non-genetic confounding variable for the same individual
###   eY~N(0, var(eY)), var(eY) = var(bX + pgY)*[1/(Rsq(XY) + Rsq(pgY) - 1] - var(c)*(1+2*beta), 
###           here, Rsq(pgY) is NOT equal to Rsq(pgX), Rsq(pgY) = var(pgY)*Rsq(XY)/var(bX)
#########################################################################

## Parse argument from command line
parser <- ArgumentParser()
parser$add_argument("--exposure_symbol", default="X", help="Symbol for exposure trait")
parser$add_argument("--outcome_symbol", default="Y", help="Symbol for outcome trait")
parser$add_argument("--geno_symbol", default="G", help="Symbol for genotype")
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--sample_X", type="integer", help="Sample size in exposure")
parser$add_argument("--sample_Y", type="integer", help="Sample size in outcome")
parser$add_argument("--sample_ref", type="integer", help="Reference sample size for LD matrix")
parser$add_argument("--sim_num", type="integer", help="Number of simulations")
parser$add_argument("--snp_num", type="integer", help="Number of causal SNPs in the simulation")
parser$add_argument("--min_maf", default=0.01, type="double", help="Minimum MAF:0.01 (default)")
parser$add_argument("--max_maf", default=0.5, type="double", help="Maximum MAF:0.5 (default)")
parser$add_argument("--rsq_GX", type="double", help="Variance explained by causal SNPs in exposure")
parser$add_argument("--rsq_PG", default=0, type="double", help="Variance explained by pleiotropic SNPs in both exposure and outcome: 0.0 (default)")
parser$add_argument("--rsq_XY", type="double", help="Variance explained by exposure in outcome")
parser$add_argument("--rsq_C", type="double", help="Variance explained by non-genetic confounding variable C in exposure")
parser$add_argument("--beta", type="double", help="Causal effect of exposure on outcome")
parser$add_argument("--thread", type="integer", help="multi-thread number")
parser$add_argument("--seed", default=123456, type="integer", help="random seed No.")
parser$add_argument("--pleio_ratio", default=0, type="double", help="% Percentage of pleiotropic SNPs: 0.0 (default)")
args <- parser$parse_args()

registerDoMC(args$thread)  #change to your number of CPU cores

### Parsing parameters
myseed = args$seed
exposure_symbol = args$exposure_symbol   
outcome_symbol = args$outcome_symbol  
geno_symbol = args$geno_symbol  

simuDIR = args$outdir
beta = args$beta

sample_X = args$sample_X 
sample_Y = args$sample_Y
sample_ref = args$sample_ref

sim_num = args$sim_num 
snp_num = args$snp_num 

MAF_min = args$min_maf # 0.01 # Allele frequency threshold
MAF_max = args$max_maf # 0.5

rsq_GX = args$rsq_GX
rsq_XY = args$rsq_XY
rsq_C = args$rsq_C
rsq_PG = args$rsq_PG
pleio_ratio = args$pleio_ratio

if (pleio_ratio==0 & rsq_PG!=0){
  stop("Check: pleio_ratio==0 but rsq_PG!=0")
} else if (pleio_ratio!=0 & rsq_PG==0) {
  stop("Check: pleio_ratio!=0 & rsq_PG==0")
}


# -------------------------------------------
# Generate the causal and pleio SNPs effect
# -------------------------------------------
set.seed(myseed); maf <- runif(snp_num, MAF_min, MAF_max)
var_gamma_X <- 1/(2*maf*(1-maf))
set.seed(myseed); gamma_X <- rnorm(snp_num, 0, sqrt(var_gamma_X))
gamma_X <- round(gamma_X, digits=4)

pleio_snp_num = 0
pleio_geno = c()
pleio_maf = c()
pleio_gamma_X = 0   # effect on exposure X
pleio_alpha_Y = 0   # effect on outcome Y
if (rsq_PG != 0) {
  pleio_snp_num <- round(snp_num*pleio_ratio, digits=0)
  set.seed(myseed); pleio_maf <- runif(pleio_snp_num, MAF_min, MAF_max)

  var_pleio_effect <- 1/(2*pleio_maf*(1-pleio_maf))
  set.seed(myseed); pleio_gamma_X <- rnorm(pleio_snp_num, 0, sqrt(var_pleio_effect))
  set.seed(myseed); pleio_alpha_Y <- rnorm(pleio_snp_num, 0, sqrt(var_pleio_effect))

  pleio_gamma_X <- round(pleio_gamma_X, digits=4)
  pleio_alpha_Y <- round(pleio_alpha_Y, digits=4)
} 


# -----------------------------------------------------------------------------------------------------------
### This is important, to remove the effect of seed on the following random functions
###   Otherwise, will be the same as long as the sample size did not change, e.g. rbinom(sample_ref, 2, p=maf[x])
### However, for a comparion among different rsqPG, this is fine, because it will show how the pleiotropic 
###   effect increase as the rsqPG increase (as the genotypes are the same for all rsqPG simulations, only
###   phenotypes changes)

# rm(.Random.seed, envir=.GlobalEnv)
# -----------------------------------------------------------------------------------------------------------


# -------------------------------------------
#  Generate reference samples for LD matrix
# -------------------------------------------
message("Generate LD matrix using reference sample")
geno_ref <- sapply(1:snp_num, function(x) rbinom(sample_ref, 2, p=maf[x]))
pleio_geno_ref = c()
if (rsq_PG != 0) { 
  pleio_geno_ref <- sapply(1:pleio_snp_num, function(x) rbinom(sample_ref, 2, p=pleio_maf[x]))
}
ldmat <- cor(cbind(geno_ref, pleio_geno_ref))


# -------------------------------------------
#  subfunction: Simulate exposure
# -------------------------------------------
sim_exposure <- function(sample_size) {
  res_out = list()
  # Generate causal SNPs
  geno_X <- sapply(1:snp_num, function(x) rbinom(sample_size, 2, p=maf[x]))
  geno_std_X <- sapply(1:snp_num, function(x) (geno_X[,x] - 2*maf[x])/sqrt(2*maf[x]*(1 - maf[x])))
  geno_std_X <- round(geno_std_X, digits=4)

  ### Generate the pleiotropic SNPs
  if (rsq_PG != 0) {
  	pleio_geno <- sapply(1:pleio_snp_num, function(x) rbinom(sample_size, 2, p=pleio_maf[x]))
  	pleio_geno_std <- sapply(1:pleio_snp_num, function(x) (pleio_geno[,x] - 2*pleio_maf[x])/sqrt(2*pleio_maf[x]*(1 - pleio_maf[x])))
    pleio_geno_std <- round(pleio_geno_std, digits=4)
  } 
    
  ### Generate exposure X for all individuals
  g <- geno_std_X%*%gamma_X
  g <- round(g, digits=4)
  var_g <- sd(g)^2     # var(zbzx) in GSMR paper

  pleio_g_X <- rep(0, times=sample_size)
  var_pleio_g_X <- 0
  if (rsq_PG != 0) {
  	pleio_g_X <- pleio_geno_std%*%pleio_gamma_X
    pleio_g_X <- round(pleio_g_X, digits=4)
  	var_pleio_g_X <- sd(pleio_g_X)^2
  }
    
  var_c <- var_g*rsq_C/rsq_GX  ## This is important

  ### Obtain weight for pleiotropic effects
  pleio_weight_X = 0
  if (rsq_PG != 0) { pleio_weight_X = sqrt((rsq_PG*var_g)/(rsq_GX*var_pleio_g_X)) }

  var_epsilon_X <- (var_g + pleio_weight_X^2*var_pleio_g_X + var_c)*(1/(rsq_GX + rsq_PG + rsq_C) - 1) # from GSMR paper
  ## Equivalent to this: 
  ##  (1) var_epsilon_X <- (var_g+var_pleio_g_X)*(1/(rsq_GX+rsq_PG)-1) - var_c
  ##  (2) var_epsilon_X <- (var_g/rsq_GX)*(1-(rsq_GX+rsq_PG+rsq_C)) 
    
  epsilon_X <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_size, 0, sqrt(var_epsilon_X)) }
  cx <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_size, 0, sqrt(var_c)) }
  phenoX <- foreach (j=1:sim_num, .combine=cbind) %dopar% { g + pleio_weight_X*pleio_g_X + cx[,j] + epsilon_X[,j] }
  phenoX <- round(phenoX, digits=4)

  res_out[["var_c"]] <- var_c
  # res_out[["cx"]] <- cx
  res_out[["pleio_weight_X"]] <- pleio_weight_X

  res_out[["phenoX"]] <- phenoX
  res_out[["phenoX_std"]] <- round(scale(phenoX), digits=4); #rm(phenoX);
  res_out[["geno"]] <- geno_X
  res_out[["geno_std"]] <- geno_std_X
  res_out[["pleio_geno"]] <- pleio_geno
  res_out[["pleio_geno_std"]] <- pleio_geno_std

  ## Check Rsq in exposure:
  chk_var_phenoX <- sd(phenoX)^2
  res_out[["chk_rsqGX"]] <- round(var_g/chk_var_phenoX, digits=3)
  res_out[["chk_rsqPX"]] <- round(sd(pleio_weight_X*pleio_g_X)^2/chk_var_phenoX, digits=3)
  res_out[["chk_rsqCX"]] <- round(sd(cx)^2/chk_var_phenoX, digits=3)

  return(res_out)
}


# -------------------------------------------
#  Simulate exposure and outcome
# -------------------------------------------
### Generate exposure X for all individuals
message("Simulate phenotype X")
phenoX_s1 <- sim_exposure(sample_X)

### Generate outcome Y for all individuals
message("Simulate phenotype Y")
phenoX_s2 <- sim_exposure(sample_Y)

pleio_g_Y <- rep(0, times=sample_Y)
var_pleio_g_Y <- 0
if (rsq_PG != 0) {
	pleio_g_Y <- phenoX_s2[["pleio_geno_std"]]%*%pleio_alpha_Y
	var_pleio_g_Y <- sd(pleio_g_Y)^2
}

X_Y <- beta*phenoX_s2[["phenoX"]]
var_X_Y <- sd(X_Y)^2

### Obtain weight for pleiotropic effects
pleio_weight_Y = 0
if (rsq_PG != 0) { pleio_weight_Y = sqrt((rsq_PG*var_X_Y)/(rsq_XY*var_pleio_g_Y)) }

var_c <- phenoX_s2[["var_c"]]
var_epsilon_Y <- (var_X_Y + pleio_weight_Y^2*var_pleio_g_Y)*(1/(rsq_XY + rsq_PG) - 1) - var_c*(1+2*beta)

epsilon_Y <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_epsilon_Y)) }
cy <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_c)) }
phenoY <- foreach (j=1:sim_num, .combine=cbind) %dopar% { X_Y[,j] + pleio_weight_Y*pleio_g_Y + cy[,j] + epsilon_Y[,j] }
phenoY <- round(phenoY, digits=4)

## Check Rsq in outcome:
chk_var_phenoY <- sd(phenoY)^2
chk_rsqPY <- round(sd(pleio_weight_Y*pleio_g_Y)^2/chk_var_phenoY, digits=3)
chk_rsqXY <- round(var_X_Y/chk_var_phenoY, digits=3)

#-----------------------------------------------
### Get X and its variance for each simulation
#-----------------------------------------------
# X_Y <- beta*phenoX_s2[["phenoX"]]
# var_X_Y <- sapply(1:sim_num, function(x) sd(X_Y[,x])^2)  # var_X_Y is a vector with length of sim_num

# ### Obtain weight for pleiotropic effects
# ### For each simulation, will have its own weight and epsilon to simulate phenotype Y
# pleio_weight_Y = c()
# if (rsq_PG != 0) { pleio_weight_Y = sqrt((rsq_PG*var_X_Y)/(rsq_XY*var_pleio_g_Y)) }

# var_c <- phenoX_s2[["var_c"]]
# var_epsilon_Y <- (var_X_Y + pleio_weight_Y^2*var_pleio_g_Y)*(1/(rsq_XY + rsq_PG) - 1) - var_c*(1+2*beta)

# epsilon_Y <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_epsilon_Y[j])) }
# cy <- foreach (j=1:sim_num, .combine=cbind) %dopar% { rnorm(sample_Y, 0, sqrt(var_c)) }
# phenoY <- foreach (j=1:sim_num, .combine=cbind) %dopar% { X_Y[,j] + pleio_weight_Y[j]*pleio_g_Y + cy[,j] + epsilon_Y[,j] }


# -------------------------------------------
#  Assocation Test
# -------------------------------------------
# For raw genotype output
genoG_X_raw <- cbind(phenoX_s1[["geno"]], phenoX_s1[["pleio_geno"]])
genoG_Y_raw <- cbind(phenoX_s2[["geno"]], phenoX_s2[["pleio_geno"]])

# Redefine phenoX and phenoY: both genotypes and pheontoypes are standardised
genoG_X <- cbind(phenoX_s1[["geno_std"]], phenoX_s1[["pleio_geno_std"]])
genoG_Y <- cbind(phenoX_s2[["geno_std"]], phenoX_s2[["pleio_geno_std"]])
pheno_X <- phenoX_s1[["phenoX_std"]]
pheno_Y <- scale(phenoY); # rm(phenoY); # save memory
total_snp_num <- snp_num + pleio_snp_num

# Define parameters for association summary data
Beta_XG = matrix(nrow = total_snp_num, ncol = sim_num)
seBeta_XG = matrix(nrow = total_snp_num, ncol = sim_num)
rsq_XG = matrix(nrow = total_snp_num, ncol = sim_num)
pval_XG = matrix(nrow = total_snp_num, ncol = sim_num)
# residualSE_XG = matrix(nrow = total_snp_num, ncol = sim_num)

Beta_YG = matrix(nrow = total_snp_num, ncol = sim_num)
seBeta_YG = matrix(nrow = total_snp_num, ncol = sim_num)
rsq_YG = matrix(nrow = total_snp_num, ncol = sim_num)
pval_YG = matrix(nrow = total_snp_num, ncol = sim_num)
# residualSE_YG = matrix(nrow = total_snp_num, ncol = sim_num)

### Association test
message("Association test: G-X")
fit_XG <- foreach (j=1:sim_num) %dopar% {
  fitTMP = list()
  for (i in 1:total_snp_num) {
    fit = summary(lm(pheno_X[,j] ~ genoG_X[,i]))
    fitTMP[["Beta"]][[i]] = fit$coefficients[2,1]
    fitTMP[["seBeta"]][[i]] = fit$coefficients[2,2]
    fitTMP[["rsq"]][[i]] = fit$adj.r.squared
    fitTMP[["pval"]][[i]] = fit$coefficients[2,4]
    # fitTMP[["residualSE"]][[i]] = fit$sigma
  }
  fitTMP
}

for (j in 1:sim_num) {
  for (i in 1:total_snp_num) {
    Beta_XG[i,j] = fit_XG[[j]][["Beta"]][[i]]
    seBeta_XG[i,j] = fit_XG[[j]][["seBeta"]][[i]]
    rsq_XG[i,j] = fit_XG[[j]][["rsq"]][[i]]
    pval_XG[i,j] = fit_XG[[j]][["pval"]][[i]]
    # residualSE_XG[i,j] = fit_XG[[j]][["residualSE"]][[i]]
  }
}

message("Association test: G-Y")
fit_YG <- foreach (j=1:sim_num) %dopar% {
  fitTMP = list()
  for (i in 1:total_snp_num) {
    fit = summary(lm(pheno_Y[,j] ~ genoG_Y[,i]))
    fitTMP[["Beta"]][[i]] = fit$coefficients[2,1]
    fitTMP[["seBeta"]][[i]] = fit$coefficients[2,2]
    fitTMP[["rsq"]][[i]] = fit$adj.r.squared
    fitTMP[["pval"]][[i]] = fit$coefficients[2,4]
    # fitTMP[["residualSE"]][[i]] = fit$sigma
  }
  fitTMP
}

for (j in 1:sim_num) {
  for (i in 1:total_snp_num) {
    Beta_YG[i,j] = fit_YG[[j]][["Beta"]][[i]]
    seBeta_YG[i,j] = fit_YG[[j]][["seBeta"]][[i]]
    rsq_YG[i,j] = fit_YG[[j]][["rsq"]][[i]]
    pval_YG[i,j] = fit_YG[[j]][["pval"]][[i]]
    # residualSE_YG[i,j] = fit_YG[[j]][["residualSE"]][[i]]
  }
}





# ------------------------------------------------------
#  Output results: genotype and phenotype simulation
# ------------------------------------------------------
message("Writing results: genotype and phenotype simulation")
# Output log file
loginfo = paste(
                paste0("exposure_symbol=", exposure_symbol),
                paste0("outcome_symbol=", outcome_symbol),
                paste0("geno_symbol=", geno_symbol),
                paste0("sample_X=", sample_X),
                paste0("sample_Y=", sample_Y),
                paste0("sample_ref=", sample_ref, " # For LD matrix"),
                paste0("sim_num=", sim_num),
                paste0("snp_num=", snp_num, "  # direct causal snp"),
                paste0("pleiotropic_snp_num=", pleio_snp_num),
                paste0("MAF_min=", MAF_min),
                paste0("MAF_max=", MAF_max),
                paste0("rsq_GX=", rsq_GX, " # variance explained by causal SNPs in exposure"),
                paste0("     check Rsq: S1 exp=", phenoX_s1[["chk_rsqGX"]], " | ", "S2 exp=", phenoX_s2[["chk_rsqGX"]]),
                paste0("rsq_PG=", rsq_PG, " # variance explained by pleiotropic SNPs in both exposure and outcome"),
                paste0("     check Rsq: S1 exp=", phenoX_s1[["chk_rsqPX"]], " | ", "S2 exp=", phenoX_s2[["chk_rsqPX"]], " | ", "S2 out=", chk_rsqPY),
                paste0("rsq_XY=", rsq_XY, " # variance explained by exposure in outcome"),
                paste0("     check Rsq: S2 out=", chk_rsqXY),
                paste0("rsq_C=", rsq_C, " # Variance explained by non-genetic confounding variable C"),
                paste0("     check Rsq: S1 exp=", phenoX_s1[["chk_rsqCX"]], " | ", "S2 exp=", phenoX_s2[["chk_rsqCX"]]),
                paste0("beta=", beta, " # Causal effect"),
                paste0("Effec size of SNPs were generated from N(0,1/2f[1-f])"),
                sep="\n"
)
write.table(file=paste0(simuDIR, "/", geno_symbol,"_",exposure_symbol,"_",outcome_symbol,".log"), loginfo, col.names=FALSE, row.names=FALSE, quote=FALSE)
                

## Naming columns and rows
causal_snp_symbol = paste0(geno_symbol, "_", c(1:snp_num))
pleio_snp_symbol = c()
if (rsq_PG != 0) { pleio_snp_symbol = paste0(geno_symbol, "_pleio_", c(1:pleio_snp_num)) }

all_snps = c(causal_snp_symbol, pleio_snp_symbol)
all_simulations = paste0("sim", c(1:sim_num))

all_samples_X = c(1:sample_X)
all_samples_Y = c(1:sample_Y)

## Output LD matrix
write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_LD.xls"), ldmat, col.names=all_snps, row.names=all_snps, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/geno", geno_symbol, "_LD.xls")))

## Output genotype
write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_raw.exp.xls"), genoG_X_raw, col.names=all_snps, row.names=all_samples_X, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_raw.out.xls"), genoG_Y_raw, col.names=all_snps, row.names=all_samples_Y, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/geno", geno_symbol, "_raw.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/geno", geno_symbol, "_raw.out.xls")))

# write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_std.exp.xls"), genoG_X, col.names=all_snps, row.names=all_samples_X, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/geno", geno_symbol, "_std.out.xls"), genoG_Y, col.names=all_snps, row.names=all_samples_Y, quote=FALSE, sep="\t")
# system(paste("gzip", paste0(simuDIR, "/geno", geno_symbol, "_std.exp.xls")))
# system(paste("gzip", paste0(simuDIR, "/geno", geno_symbol, "_std.out.xls")))

## Output phenotype
write.table(file=paste0(simuDIR, "/pheno", exposure_symbol, "_raw.exp.xls"), phenoX_s1[["phenoX"]], col.names=all_simulations, row.names=all_samples_X, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/pheno", outcome_symbol, "_raw.out.xls"), phenoY, col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/pheno", outcome_symbol, "_raw.exp2out.xls"), phenoX_s2[["phenoX"]], col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/pheno", exposure_symbol, "_raw.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/pheno", outcome_symbol, "_raw.out.xls")))
system(paste("gzip", paste0(simuDIR, "/pheno", outcome_symbol, "_raw.exp2out.xls")))

# write.table(file=paste0(simuDIR, "/pheno", exposure_symbol, "_std.exp.xls"), pheno_X, col.names=all_simulations, row.names=all_samples_X, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/pheno", outcome_symbol, "_std.out.xls"), pheno_Y, col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/pheno", outcome_symbol, "_std.exp2out.xls"), phenoX_s2[["phenoX_std"]], col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")
# system(paste("gzip", paste0(simuDIR, "/pheno", exposure_symbol, "_std.exp.xls")))
# system(paste("gzip", paste0(simuDIR, "/pheno", outcome_symbol, "_std.out.xls")))
# system(paste("gzip", paste0(simuDIR, "/pheno", outcome_symbol, "_std.exp2out.xls")))


## Output c
# write.table(file=paste0(simuDIR, "/c.exp.xls"), phenoX_s1[["cx"]], col.names=all_simulations, row.names=all_samples_X, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/c.out.xls"), cy, col.names=all_simulations, row.names=all_samples_Y, quote=FALSE, sep="\t")
# system(paste("gzip", paste0(simuDIR, "/c.exp.xls")))
# system(paste("gzip", paste0(simuDIR, "/c.out.xls")))

## Output genotype ref and effects
causal_geno_info = cbind(maf, gamma_X)
pleio_geno_info = cbind(pleio_maf, phenoX_s1[["pleio_weight_X"]]*pleio_gamma_X)
geno_info = rbind(causal_geno_info, pleio_geno_info)
geno_info_col = c("MAF",paste0(exposure_symbol,".exp"))
write.table(file=paste0(simuDIR, "/geno_effect.expSample.xls"), geno_info, col.names=geno_info_col, row.names=all_snps, quote=FALSE, sep="\t")

causal_geno_info = cbind(maf, gamma_X, rep(0, times=snp_num))
pleio_geno_info = cbind(pleio_maf, phenoX_s2[["pleio_weight_X"]]*pleio_gamma_X, pleio_weight_Y*pleio_alpha_Y)
geno_info = rbind(causal_geno_info, pleio_geno_info)
geno_info_col = c("MAF",paste0(exposure_symbol,".exp"), paste0(outcome_symbol,".out"))
write.table(file=paste0(simuDIR, "/geno_effect.outSample.xls"), geno_info, col.names=geno_info_col, row.names=all_snps, quote=FALSE, sep="\t")



# ------------------------------------------------------
#  Output results: association results
# ------------------------------------------------------
message("Writing results: Association test results")

write.table(file=paste0(simuDIR, "/BETA_XG.exp.xls"), Beta_XG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/BETA_YG.out.xls"), Beta_YG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/BETA_XG.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/BETA_YG.out.xls")))

write.table(file=paste0(simuDIR, "/SE_XG.exp.xls"), seBeta_XG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/SE_YG.out.xls"), seBeta_YG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/SE_XG.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/SE_YG.out.xls")))

write.table(file=paste0(simuDIR, "/PVAL_XG.exp.xls"), pval_XG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/PVAL_YG.out.xls"), pval_YG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/PVAL_XG.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/PVAL_YG.out.xls")))

write.table(file=paste0(simuDIR, "/adj.RSQ_XG.exp.xls"), rsq_XG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
write.table(file=paste0(simuDIR, "/adj.RSQ_YG.out.xls"), rsq_YG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
system(paste("gzip", paste0(simuDIR, "/adj.RSQ_XG.exp.xls")))
system(paste("gzip", paste0(simuDIR, "/adj.RSQ_YG.out.xls")))

# write.table(file=paste0(simuDIR, "/RESIDUALSE_XG.exp.xls"), residualSE_XG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
# write.table(file=paste0(simuDIR, "/RESIDUALSE_YG.out.xls"), residualSE_YG, col.names=all_simulations, row.names=all_snps, quote=FALSE, sep="\t")
# system(paste("gzip", paste0(simuDIR, "/RESIDUALSE_XG.exp.xls")))
# system(paste("gzip", paste0(simuDIR, "/RESIDUALSE_YG.out.xls")))





# ------------------------------------------------------
## For test: Check is the Rsq is controlled
# ------------------------------------------------------

testbutton = 'NO'
if (testbutton == 'YES') {
sample_X = 2000
sample_Y = 1500
sample_ref = 500
sim_num = 1000
snp_num = 50
MAF_min = 0.01
MAF_max = 0.5
rsq_GX = 0.4
rsq_XY = 0.05
rsq_C = 0.2
rsq_PG=0.06
beta=1
pleio_ratio=0.1

exposure_symbol="X"
outcome_symbol="Y"
geno_symbol="G"
# simuDIR=""



# Parallel test: 
PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% {  
  PctExpTMP = 0
  for (i in 1:snp_num) {
    fit = summary(lm(phenoX[,j] ~ geno_std_X[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
  PctExpTMP = 0
  for (i in 1:pleio_snp_num) {
    fit = summary(lm(phenoX[,j] ~ pleio_geno_std[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
  PctExpTMP = 0
  for (i in 1:pleio_snp_num) {
    fit = summary(lm(phenoY[,j] ~ pleio_geno_std[,i]))
    PctExpTMP = PctExpTMP+fit$adj.r.squared
  }
  PctExpTMP
}
mean(as.numeric(PctExp))


PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoY[,j] ~ phenoX[,j]))
    fit$adj.r.squared
}
mean(as.numeric(PctExp))

PctExp = c()
PctExp <- foreach (j=1:sim_num) %dopar% { 
    fit = summary(lm(phenoX[,j] ~ cx[,j]))
    fit$adj.r.squared
}
mean(as.numeric(PctExp))

## Test without parallel
# PctExp = c()
# for (j in 1:sim_num) { 
#   PctExpTMP = 0
#   for (i in 1:snp_num) {
#     fit = summary(lm(phenoX1[,j] ~ genoG1[,i]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExp = c(PctExp, PctExpTMP)
# }
# mean(PctExp)


# PctExp = c()
# for (j in 1:sim_num) { 
#   PctExpTMP = 0
#   for (i in 1:pleio_snp_num) {
#     k = i + snp_num
#     fit = summary(lm(phenoX1[,j] ~ genoG1[,k]))
#     PctExpTMP = PctExpTMP+fit$adj.r.squared
#   }
#   PctExp = c(PctExp, PctExpTMP)
# }
# mean(PctExp)
}