rm(list=ls())
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
micro_data$rarefy_samples(sample.size = 10000)

##--------------------------------------------- Calculate Alpha Diversity Index
t1 <- trans_alpha$new(dataset = micro_data, group = "group")
# return t1$data_stat
head(t1$data_stat)
## The transformed diversity data is stored in object$data_alpha ...
## The group statistics are stored in object$data_stat ...

##-------------------------------------------- test the differences among groups
# Then, we test the differences among groups using Kruskal-Wallis Rank Sum Test
# (overall test when groups > 2), Wilcoxon Rank Sum Tests (for paired groups), 
# Dunn¡¯s Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2)
# and anova with multiple comparisons.

t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)


t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
head(t1$res_diff)

# more options
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
head(t1$res_diff)
t1$cal_diff(method = "wilcox")
head(t1$res_diff)
t1$cal_diff(method = "t.test")

##------------------------------------------------ Choose Anova and Plot Figure
t1$cal_diff(method = "anova")
# return t1$res_diff
head(t1$res_diff)

##  Observed
t1$res_diff %<>% base::subset(Significance != "ns")
p1 <- t1$plot_alpha(measure = "Observed", boxplot_add = "dotplot", xtext_size = 15)

##  CHAO1
t1$res_diff %<>% base::subset(Significance != "ns")
p2 <- t1$plot_alpha(measure = "Chao1", boxplot_add = "dotplot", xtext_size = 15)

##  ACE
t1$res_diff %<>% base::subset(Significance != "ns")
p3 <- t1$plot_alpha(measure = "ACE", boxplot_add = "dotplot", xtext_size = 15)

##  Shannon
t1$res_diff %<>% base::subset(Significance != "ns")
p4 <- t1$plot_alpha(measure = "Shannon", boxplot_add = "dotplot", xtext_size = 15)

##  Simpson
t1$res_diff %<>% base::subset(Significance != "ns")
p5 <- t1$plot_alpha(measure = "Simpson", boxplot_add = "dotplot", xtext_size = 15)

# ##  InvSimpson
# t1$res_diff %<>% base::subset(Significance != "ns")
# p6 <- t1$plot_alpha(measure = "InvSimpson", boxplot_add = "dotplot", xtext_size = 15)

##  Fisher
t1$res_diff %<>% base::subset(Significance != "ns")
p7 <- t1$plot_alpha(measure = "Fisher", boxplot_add = "dotplot", xtext_size = 15)

# ##  Coverage
# t1$res_diff %<>% base::subset(Significance != "ns")
# p8 <- t1$plot_alpha(measure = "Coverage", boxplot_add = "dotplot", xtext_size = 15)

grid.arrange(p1,p2,p3,p4,p5,p7)

# t1
dir.create('Figures/Diversity', recursive = TRUE)
data_file_name1 <- paste0('Figures/Diversity/', data_source,"_Alpha_Diversity_stat.csv")
data_abundance1 <- data.frame(Diversity_stat = t1[["data_stat"]])
write.csv(data_abundance1 ,data_file_name1,quote=F)

data_file_name2 <- paste0('Figures/Diversity/', data_source,"_Alpha_Diversity_data.csv")
data_abundance2 <- data.frame(Diversity_stat = t1[["data_alpha"]])
write.csv(data_abundance2 ,data_file_name2,quote=F)
