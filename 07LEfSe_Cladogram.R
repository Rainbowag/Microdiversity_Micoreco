rm(list=ls())

library(microeco)

pacman::p_load(tidyverse,microeco,magrittr)

# data_source  Fungal or Bacteria
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

head(otu_table)[1:6,1:6]; head(group_table)[1:6, ]; head(tax_table)[,1:6]

# 创建microtable对象
d1 <- microtable$new(sample_table = group_table,
                          otu_table = otu_table, 
                          tax_table = tax_table)

# As an example, use 10000 sequences in each sample
# d1$rarefy_samples(sample.size = 10000)


# 执行lefse分析
t1 <- trans_diff$new(dataset = d1, 
                        method = "lefse", 
                        group = "group", 
                        alpha = 0.05, 
                        lefse_subgroup = NULL,
                        p_adjust_method = "none")

# see t1$res_diff for the result
# 查看分析结果
head(t1$res_diff)[1:6]

# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)

# 绘制前30个具有最高LDA（log10）的分类单元的差异特征柱状图
t1$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    group_order = c("SU", "SF", "GL","PF","OD"))
  # ggsci::scale_color_npg() +
  # ggsci::scale_fill_npg()

# show part of the table
t1$res_diff[1:5, c(1, 3, 4, 6)]

dir.create('Figures/LEfSe', recursive = TRUE)
# pdf(paste('Figures/LEfSe/','LEfSe_Cladogram_',data_source,'.pdf', sep = ''), width = 11.69, height = 8.27)
# 展示前200个分类单元和前50个特征
# 需要调用ggtree包
p1 <- t1$plot_diff_cladogram(use_taxa_num = 500, 
                          use_feature_num = 50, 
                          clade_label_level = 7, 
                          group_order = c("SU", "SF", "GL","PF","OD"),
                          annotation_shape_size = 8,
                          clade_label_size = 3)
p1+theme(legend.text = element_text(size=12))
# ggsave(p1,filename = paste('Figures/LEfSe/','LEfSe_Cladogram_',data_source,'.pdf', sep = ''), width = 11.69, height = 8.27)
