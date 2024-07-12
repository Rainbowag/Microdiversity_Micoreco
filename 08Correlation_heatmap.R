rm(list=ls())

library(microeco)
library(ggplot2)


# data_source  Fungal or Bacteria
data_source <- "Bacteria"

# taxonomy type: Kingdom	Phylum	Class	Order	Family	Genus	Species
tax_type <- "Genus"

# read data
# Bacterial
if(data_source=="Bacteria")
{
  otu_table <- read.csv('Data/Bacteria_OTU.csv', row.names = 1)
  group_table <- read.csv('Data/Bacteria_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Bacteria_taxonomy.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}else if(data_source=="Fungal")
{
  # Fungal
  otu_table <- read.csv('Data/Fungal_OTU_4Genus.csv', row.names = 1)
  group_table <- read.csv('Data/Fungal_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Fungal_taxonomy_4Genus.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}

micro_data <- microtable$new(sample_table = group_table,
                             otu_table = otu_table,
                             tax_table = tax_table)

t1 <- trans_env$new(dataset = micro_data, add_data = env_table[, 1:9])

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1$cal_cor(use_data = tax_type, p_adjust_method = "fdr", p_adjust_type = "Env")

t1$plot_cor()

# filter genera that donot have at least one ***
t1$plot_cor(filter_feature = c("", "*", "**"))


# first create trans_diff object as a demonstration
t2 <- trans_diff$new(dataset = micro_data, method = "rf", group = "group", taxa_level = tax_type ,p_adjust_method = "none")
# then create trans_env object
t1 <- trans_env$new(dataset = micro_data, add_data = env_table[, 1:9])
# use other_taxa to select taxa you need
tax_len <- length(t2[["res_diff"]][["Taxa"]])
t1$cal_cor(use_data = "other", p_adjust_method = "none", other_taxa = t2$res_diff$Taxa[1:tax_len])
p1<-t1$plot_cor()

# clustering heatmap; require pheatmap package
# Let's take another color pallete
p1<-t1$plot_cor(pheatmap = TRUE, color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))

# t1
data_file_name1 <- paste0("Figures/Correlation/",data_source,"_Correlation_heatmap_",tax_type,".csv")
data_abundance1 <- data.frame(Correlation = t1[["res_cor"]])
write.csv(data_abundance1, data_file_name1, quote=F)

data_file_name2 <- paste0("Figures/Correlation/",data_source,"_Correlation_heatmap_",tax_type,".pdf")
ggsave(p1,filename = data_file_name2, width = 11.69, height = 8.27)

# ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Figure 2 diversity heatmap
# t1 <- trans_env$new(dataset = micro_data, add_data = env_table[, 1:9])
# # use add_abund_table parameter to add the extra data table
# 
# t1$cal_cor(add_abund_table = micro_data$alpha_diversity,p_adjust_method = "none")
# # try to use ggplot2 with clustering plot
# # require ggtree and aplot packages to be installed (https://chiliubio.github.io/microeco_tutorial/intro.html#dependence)
# t1$plot_cor(cluster_ggplot = "both")

