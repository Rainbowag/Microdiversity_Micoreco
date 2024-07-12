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

##--------------------------------------------- Calculate Beta Diversity Index
micro_data$cal_betadiv()

# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = micro_data, group = "group", measure = "bray")


# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "group", plot_shape = "group", plot_type = c("point", "ellipse"))

##-------------------------------- Then we plot and compare the group distances.
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "anova")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean")

##-------------------------- calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "anova")
t1$plot_group_distance(boxplot_add = "mean")
