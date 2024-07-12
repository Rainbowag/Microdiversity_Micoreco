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

# # As an example, use 10000 sequences in each sample
# micro_data$rarefy_samples(sample.size = 10000)


##------------------------------------------------------ Figure 1 Abundant Phyla
# select top 10 abundant Phyla.
t1 <- trans_abund$new(dataset = micro_data, taxrank = "Phylum", ntaxa = 10)

t1$plot_bar(others_color = "grey70", 
            facet = c("group", "SampleID"), 
            xtext_keep = FALSE, 
            legend_text_italic = FALSE, 
            barwidth = 1)
# t1
data_file_name <- paste0(data_source,"_Phylum_Abundance.csv")
data_abundance <- data.frame(Abundance = t1[["data_abund"]][["Abundance"]],
                             Taxonomy = t1[["data_abund"]][["Taxonomy"]],
                             SampleID = t1[["data_abund"]][["SampleID"]])
write.csv(data_abundance ,data_file_name,quote=F)


# The groupmean parameter can be used to obtain the group-mean barplot.
t2 <- trans_abund$new(dataset = micro_data, taxrank = "Phylum", ntaxa = 10, 
                      groupmean = "group")

g1 <- t2$plot_bar(legend_text_italic = TRUE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))

##---------------------------------------------------- Figure 3 Phylum and Family
t4 <- trans_abund$new(dataset = micro_data, taxrank = "Family", ntaxa = 20, 
                         high_level = "Phylum")
t4$plot_bar(ggnested = TRUE,xtext_angle = 0, facet = c("group", "SampleID"))

t5 <- trans_abund$new(dataset = micro_data, taxrank = "Family", ntaxa = 20, 
                      high_level = "Phylum",
                      groupmean = "group", group_morestats = TRUE)

t5$plot_bar(ggnested = TRUE,xtext_angle = 0)

##------------------------------------------------------------- Figure 4 Heatmap
# show 40 taxa at Genus level
t6 <- trans_abund$new(dataset = micro_data, taxrank = "Genus", ntaxa = 40)
t6$plot_heatmap(facet = "group", xtext_keep = FALSE, withmargin = FALSE)

# t6
data_file_name6 <- paste0(data_source,"_Genus_Abundance.csv")
data_abundance6 <- data.frame(Abundance = t6[["data_abund"]][["Abundance"]],
                             Taxonomy = t6[["data_abund"]][["Taxonomy"]],
                             SampleID = t6[["data_abund"]][["SampleID"]])
write.csv(data_abundance6 ,data_file_name6,quote=F)