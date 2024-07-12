rm(list=ls())
##----------------------------------------------------------------- prepare data
library(microeco)
library(ggplot2)
theme_set(theme_classic())


data_source <- "Fungal"

for (landuse_type in c("SU","SF","GL","PF","OD")){
# landuse type SU SF GL PF OD ALL
# landuse_type <- "SF"

# read data
# Bacterial
if(data_source=="Bacteria")
{
  otu_table <- read.csv('Data/Bacteria_OTU_4Species.csv', row.names = 1)
  group_table <- read.csv('Data/Bacteria_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Bacteria_taxonomy_4Species.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}else if(data_source=="Fungal")
{
  # Fungal
  otu_table <- read.csv('Data/Fungal_OTU_4Species.csv', row.names = 1)
  group_table <- read.csv('Data/Fungal_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Fungal_taxonomy_4Species.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}

data.extract <- function(out_table, landuse_type){
  switch(landuse_type,
         SU =  otu_table[1:3],
         SF =  otu_table[4:6],
         GL =  otu_table[7:9],
         PF =  otu_table[10:12],
         OD =  otu_table[13:15],
         ALL = otu_table[1:15])
}

otu_table = data.extract(otu_table,landuse_type)

micro_data <- microtable$new(sample_table = group_table,
                             otu_table = otu_table,
                             tax_table = tax_table)

# As an example, use 10000 sequences in each sample
# micro_data$rarefy_samples(sample.size = 10000)

##-----------------------------------------------------------------
# The parameter cor_method in trans_network is used to select correlation calculation method.
# default pearson or spearman correlation invoke R base cor.test, a little slow
# t1 <- trans_network$new(dataset = micro_data, cor_method = "pearson", filter_thres = 0.0001,taxa_level = "Phylum")
# return t1$res_cor_p list, containing two tables: correlation coefficient table and p value table

# require WGCNA package
# if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
t1 <- trans_network$new(dataset = micro_data, cor_method = "pearson", 
                        use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)


# construct network; require igraph package
# t1$cal_network(COR_p_thres = 0.01, COR_optimization = FALSE,add_taxa_name = c("Phylum"))
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6,add_taxa_name = c("Phylum"))
# return t1$res_network

# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_walktrap")

t1$cal_network_attr()
t1$res_network_attr


# t1$get_node_table(node_roles = TRUE)
# 
# mean(t1$res_node_table$degree)
# mean(t1$res_node_table$betweenness)

# require rgexf package to be installed
t1$save_network(filepath = 
                  paste0("Figures/Network/Network_File/",data_source,"_network_",landuse_type,".gexf"))
}