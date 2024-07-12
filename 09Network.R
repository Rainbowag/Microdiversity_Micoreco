rm(list=ls())
##----------------------------------------------------------------- prepare data
library(microeco)
library(ggplot2)
theme_set(theme_classic())


data_source <- "Fungal"
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

##-----------------------------------------------------------------
# The parameter cor_method in trans_network is used to select correlation calculation method.
# default pearson or spearman correlation invoke R base cor.test, a little slow
t1 <- trans_network$new(dataset = micro_data, cor_method = "spearman", filter_thres = 0.001)
# return t1$res_cor_p list, containing two tables: correlation coefficient table and p value table

# require WGCNA package
# if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
t1 <- trans_network$new(dataset = micro_data, cor_method = "pearson", 
                        use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)


# construct network; require igraph package
t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t1$res_network

# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")

# require rgexf package to be installed
t1$save_network(filepath = paste0("network_",data_source,".gexf"))


# calculate network attributes
t1$cal_network_attr()
t1$res_network_attr


# get node properties
t1$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
t1$plot_taxa_roles(use_type = 1)

# plot node roles with phylum information
t1$plot_taxa_roles(use_type = 2)

# Now, we show the eigengene analysis of modules. The eigengene of a module, 
# i.e. the first principal component of PCA, represents the main variance of the
# abundance in the species of the module.
t1$cal_eigen()
# return t1$res_eigen

# Then we perform correlation heatmap to show the associations between 
# eigengenes and environmental factors.
# create trans_env object
t2 <- trans_env$new(dataset = micro_data, add_data = env_table[, 1:6])
# calculate correlations
t2$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap
t2$plot_cor()





# default parameter represents using igraph plot.igraph function
t2$plot_network()
# use ggraph method; require ggraph package
# If ggraph is not installed; first install it with command: install.packages("ggraph")
t2$plot_network(method = "ggraph", node_color = "Phylum")
# use networkD3 package method for the dynamic network visualization in R
# If networkD3 is not installed; first install it with command: install.packages("networkD3")
t1$plot_network(method = "networkD3", node_color = "module")
t1$plot_network(method = "networkD3", node_color = "Phylum")


# use_col is used to select a column of t1$res_node_table
tmp <- t1$trans_comm(use_col = "module", abundance = FALSE)
tmp
tmp$otu_table[tmp$otu_table > 0] <- 1
tmp$tidy_dataset()
tmp$cal_abund()
tmp2 <- trans_abund$new(tmp, taxrank = "Phylum", ntaxa = 10)
tmp2$data_abund$Sample %<>% factor(., levels = rownames(tmp$sample_table))
tmp2$plot_line(xtext_angle = 30, 
               color_values = RColorBrewer::brewer.pal(12, "Paired")) + ylab("OTUs ratio (%)")


t1$cal_sum_links(taxa_level = "Phylum")
# require chorddiag package; see https://github.com/mattflor/chorddiag
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10,
                  color_values = RColorBrewer::brewer.pal(10, "Paired"),
                  groupnameFontsize = 15,
                  margin=200
                  )

# t1
data_file_name <- paste0("Figures/Chorddiag/",data_source,"_Chorddiag_link_positive.csv")
data_abundance <- data.frame(positive = t1[["res_sum_links_pos"]])
write.csv(data_abundance, data_file_name, quote=F)