rm(list=ls())
##----------------------------------------------------------------- prepare data
library(microeco)
library(ggplot2)
theme_set(theme_classic())

data_source <- "Bacteria"
# read data
# Bacterial
if(data_source=="Bacteria")
  {
  otu_table <- read.csv('Data/Bacteria_OTU_4Family.csv', row.names = 1)
  group_table <- read.csv('Data/Bacteria_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Bacteria_taxonomy_4Family.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
  }else if(data_source=="Fungal")
    {
    # Fungal
    otu_table <- read.csv('Data/Fungal_OTU_4Family.csv', row.names = 1)
    group_table <- read.csv('Data/Fungal_group.csv', row.names = 1)
    tax_table <- read.csv('Data/Fungal_taxonomy_4Family.csv', row.names = 1)
    env_table <- read.csv('Data/Env.csv', row.names = 1)
    }


micro_data <- microtable$new(sample_table = group_table,
                             otu_table = otu_table,
                             tax_table = tax_table)

# # As an example, use 10000 sequences in each sample
# micro_data$rarefy_samples(sample.size = 10000)


##----------------------------------------------------- Figure 2 Abundant Family
t3 <- trans_abund$new(dataset = micro_data, taxrank = "Family", ntaxa = 20)
                      # prefix = "\\|")

t3$plot_bar(bar_type = "full", 
                 use_alluvium = TRUE, 
                 clustering = TRUE, 
                 xtext_angle = 0, 
                 xtext_size = 8, 
                 color_values = RColorBrewer::brewer.pal(8, "Set2"))

t5 <- trans_abund$new(dataset = micro_data, taxrank = "Family", ntaxa = 20, 
                        high_level = "Phylum")

t5$plot_bar(facet = c("group", "SampleID"),ggnested = TRUE,xtext_keep = FALSE,xtext_angle = 0)

# g1 <- t5$plot_bar(facet = c("group", "SampleID"),ggnested = TRUE,xtext_keep = FALSE,xtext_angle = 0)
# 
# g1 <- g1+theme(legend.position="bottom", legend.box = "horizontal")
# 
# g1

# data_file_name5 <- paste0(data_source,"_PyhlmandFamily_Abundance.csv")
# data_abundance5 <- data.frame(Abundance = t5[["data_abund"]][["Abundance"]],
#                               Taxonomy = t5[["data_abund"]][["Taxonomy"]],
#                               SampleID = t5[["data_abund"]][["SampleID"]])
# write.csv(data_abundance5 ,data_file_name5,quote=F)