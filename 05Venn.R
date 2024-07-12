##----------------------------------------------------------------- prepare data
library(microeco)
library(ggplot2)
theme_set(theme_classic())


# data_source  Fungal or Bacteria
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
# micro_data$rarefy_samples(sample.size = 10000)

##----------------------------------------------------------------- Figure 1
# merge samples as one community for each group
dataset1 <- micro_data$merge_samples(use_group = "group")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = NULL)
p1 <- t1$plot_venn(text_size = 8,
                   text_name_size = 10)
p1+theme()

##----------------------------------------------------------------- Figure 2
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn(text_size = 3)
# The integer is OTU number
# The percentage data is the sequence number/total sequence number

##----------------------------------------------------------------- Figure 3
# use "Type" column in sample_table
dataset1 <- micro_data$merge_samples(use_group = "group")
t1 <- trans_venn$new(dataset1)
t1$plot_venn(text_size = 3,petal_plot = TRUE, petal_color = RColorBrewer::brewer.pal(8, "Dark2"))
t1$plot_venn(text_size = 3,petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, 
             petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")

##----------------------------------------------------------------- Figure 3
# tmp <- micro_data$merge_samples(use_group = "group")
# tmp
# t1 <- trans_venn$new(dataset = tmp)
# # only show some sets with large intersection numbers
# t1$data_summary %<>% .[.[, 1] > 20, ]
# g1 <- t1$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, 
#                   up_bar_fill = "grey50", left_bar_fill = "grey50", 
#                   bottom_point_color = "black")
# g1
# # g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# # as g1 is comprised of several sub-plots, please adjust the details for each sub-plot
# g1[[1]]
# g1[[2]]
# ggsave("diamonds.png")
# 
# dataset1 <- micro_data$merge_samples(use_group = "group")
# t1 <- trans_venn$new(dataset1)
# ## The details of each venn part is stored in object$data_details ...
# ## The venn summary table used for plot is stored in object$data_summary ...

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)

