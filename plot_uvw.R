args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
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

rsqPG = args[1]
j = 3  # selected simulation

sim_para = paste0("beta",beta,"_maf",maf,"_snp",snp_num,"_sizeX",sample_X,"Y",sample_Y,"Ref",sample_ref,"_sim",sim_num,"_pleio", pleio_ratio,"_rsqGX",rsqGX, "_rsqPG",rsqPG,"_rsqXY",rsqXY,"_rsqC", rsqC)
bzx_dat = read.table(paste0(sim_para,"/BETA_XG.exp.xls.gz"), head=TRUE, sep="\t")
# bzx_var_dat = read.table(paste0(sim_para,"/SE_XG.exp.xls.gz"), head=TRUE, sep="\t")
bzy_dat = read.table(paste0(sim_para,"/BETA_YG.out.xls.gz"), head=TRUE, sep="\t")
# bzy_var_dat = read.table(paste0(sim_para,"/SE_YG.out.xls.gz"), head=TRUE, sep="\t")

causal_zx = c()
for (i in 1:snp_num){ causal_zx = c(causal_zx, bzx_dat[i,j]) }
pleio_zx = c()
for (i in 1:pleio_snp_num){ pleio_zx = c(pleio_zx, bzx_dat[i+snp_num,j]) }

causal_zy = c()
for (i in 1:snp_num){ causal_zy = c(causal_zy, bzy_dat[i,j]) }
pleio_zy = c()
for (i in 1:pleio_snp_num){ pleio_zy = c(pleio_zy, bzy_dat[i+snp_num,j]) }

zx = c(causal_zx, pleio_zx)
zy = c(causal_zy, pleio_zy)
snp_type = c(rep("causal", snp_num), rep("pleio", pleio_snp_num))

# plot(zy, zx, col=as.factor(snp_type))

title_para = paste0("beta=",beta, " | maf=",maf, " | sim_selected=",j, "\nsizeX=",sample_X, " | sizeY=",sample_Y, " | sizeRef=",sample_ref, "\ncausalSNP=",snp_num," | pleioSNP=", pleio_snp_num, "\nrsqGX=",rsqGX, " | rsqPG=",rsqPG, " | rsqXY=",rsqXY," | rsqC=",rsqC)
df <- data.frame(x = zx, y = zy)
p <- ggplot(df, aes(x, y, colour=snp_type)) + geom_point()
p <- p + scale_x_continuous(expand = c(0, 0), limits=c(min(zx), max(zx))) + 
	 scale_y_continuous(expand = c(0, 0), limits=c(min(zy), max(zy))) +
	 scale_colour_manual(values=c("black", "red")) +
	 labs(x = "SNP-exposure coefficient", y = "SNP-outcome coefficient") +
	 geom_hline(yintercept=0) + geom_vline(xintercept=0) +
	 ggtitle(title_para)

ggsave("flexmix_res/plot_uvw.pdf", plot = p, 
	dpi = 300, units="cm")