## Compare the results of HEIDI-outlier and Mixture model
library(ggplot2)
library(gridExtra)

rsqPG_list = c("0.03","0.04","0.05","0.06","0.07","0.08","0.09","0.1")

load("flexmix_res/flexmix_percentage.RData")
flmr_res <- out

load("flexmix_res/HEIDI_percetage.RData")
heidi_res <- out

power_ave <- c(flmr_res[,"power_ave"], heidi_res[,"power_ave"])
fp_ave <- c(flmr_res[,"fp_ave"], heidi_res[,"fp_ave"])

data <- data.frame(
            rsqPG_ratio = c(rsqPG_list, rsqPG_list),
            power = power_ave, 
            fp = fp_ave,
            class = c(rep("mixture-model", length(rsqPG_list)),
                      rep("HEIDI-outlier", length(rsqPG_list)))
        )



 ## Plot true positive results
p1 <- ggplot(data=data, aes(x=factor(data$rsqPG_ratio), y=data$power, fill=factor(data$class))) +
  theme(legend.title=element_blank()) +
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) +
  # geom_errorbar(limits, position = 'dodge', width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.1), limits=c(0,1)) +
  labs(x = "Variance explained by pleiotropic SNPs", y = "True positive rate") +
  ggtitle("mixture-model v.s. HEIDI-outlier")


## Plot false positive results
p2 <- ggplot(data=data, aes(x=factor(data$rsqPG_ratio), y=data$fp, fill=factor(data$class))) +
  theme(legend.title=element_blank()) +
  geom_bar(stat = "identity", position = 'dodge', width = 0.5) +
  # geom_errorbar(limits, position = 'dodge', width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.1), limits=c(0,1)) +
  labs(x = "Variance explained by pleiotropic SNPs", y = "False positive rate")
  # ggtitle("mixture-model v.s. HEIDI-outlier")


ggsave("flexmix_res/plot_HEIDIvMIX.pdf", plot = grid.arrange(p1, p2, nrow=2), 
	dpi = 300, units="cm")