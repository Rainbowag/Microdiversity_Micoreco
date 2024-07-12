# 
rm(list=ls())#clear Global Environment
setwd('D:/ChaohaoLibrary/Research/11_Microeco')#设置工作路径
#加载包
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
theme_set(theme_bw())


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
                             otu_table = otu_table)

# As an example, use 10000 sequences in each sample
micro_data$rarefy_samples(sample.size = 10000)

micro_data

# PCoA
micro_data$cal_betadiv(unifrac = FALSE)
t1 <- trans_beta$new(dataset = micro_data, 
                     group = "group", 
                     measure = "bray")

t1$cal_ordination(ordination = "PCoA")

# extract the axis scores
tmp <- t1$res_ordination$scores
# use differential test in trans_env class
t2 <- trans_env$new(dataset = micro_data, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "group", method = "anova")


p1 <- t1$plot_ordination(plot_color = "group", 
                         plot_type = c("point", "centroid","ellipse"),
                         point_size = 4, 
                         point_alpha = 0.5, 
                         ellipse_chull_fill = T,
                         centroid_segment_alpha = 0.9, 
                         centroid_segment_size = 1, 
                         centroid_segment_linetype = 1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text = element_text(size=8,colour = "black"))+
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2)

# groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
p2 <- t2$plot_diff(measure = "PCo1", add_sig = T) + 
  theme_bw() + coord_flip() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

p3 <- t2$plot_diff(measure = "PCo2", add_sig = T) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
g <- p1 %>% 
  insert_top(p2, height = 0.2) %>% 
  insert_right(p3, width = 0.2)

g


# use 1.4-fold of the scores as axis ranges
x_lim <- range(tmp[, 1]) * 1.2
y_lim <- range(tmp[, 2]) * 1.2
# limit x and y axis without any extension
p1 <- p1 + scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0))
# limit x axis of upper figure (it's y axis when flipped)
p2 <- p2 + scale_y_continuous(limits = x_lim, expand = c(0, 0))
# limit y axis of right-hand figure
p3 <- p3 + scale_y_continuous(limits = y_lim, expand = c(0, 0))
g <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)
g


p2 <- p2 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3 <- p3 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
g <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)
g
# # save g to computer
# ggsave("test1.pdf", g, width = 7, height= 6)