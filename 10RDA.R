rm(list=ls())
##----------------------------------------------------------------- prepare data
library(microeco)
library(ggplot2)
theme_set(theme_classic())


data_source <- "Fungal"
# taxonomy type: Kingdom	Phylum	Class	Order	Family	Genus	Species
tax_level <- "Phylum"
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
  otu_table <- read.csv('Data/Fungal_OTU.csv', row.names = 1)
  group_table <- read.csv('Data/Fungal_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Fungal_taxonomy.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}

micro_data <- microtable$new(sample_table = group_table,
                             otu_table = otu_table,
                             tax_table = tax_table)

# As an example, use 10000 sequences in each sample
micro_data$rarefy_samples(sample.size = 10000)

##-------------------------------------------------------------- db-RDA Analysis
t1 <- trans_env$new(dataset = micro_data, add_data = env_table[, 1:9])

# use bray-curtis distance for dbRDA
t1$cal_ordination(method = "dbRDA", use_measure = "bray")
# t1$res_rda is the result list stored in the object
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "group")

# the function cal_ordination_anova is implemented to check the significance of 
# the ordination model instead of the encapsulation in cal_ordination. 
# Furthermore, the function cal_ordination_envfit can be used to get the 
# contribution of each variables to the model.
t1$cal_ordination_anova(method = "dbRDA", use_measure = "bray")
t1$cal_ordination_envfit(method = "dbRDA", use_measure = "bray")

# use Genus

# RDA
t1$cal_ordination(method = "RDA", taxa_level = tax_level)
# # dbRDA
# t1$cal_ordination(method = "dbRDA", use_measure = "bray", taxa_level = tax_level)

# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
p1 <- t1$plot_ordination(plot_color = "group")

p1 

data_file_name2 <- paste0("Figures/RDA/",data_source,"_RDA_",tax_level,".pdf")
ggsave(p1,filename = data_file_name2, width = 11.69, height = 8.27)

# t1
data_file_name1 <- paste0("Figures/RDA/",data_source,"_RDA_data_",tax_level,".csv")
data_abundance1 <- data.frame(RDA = t1[["res_ordination_terms"]])
write.csv(data_abundance1, data_file_name1, quote=F)