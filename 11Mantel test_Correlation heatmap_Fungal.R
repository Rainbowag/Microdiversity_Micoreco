rm(list=ls())
# prepare data
library(microeco)
library(magrittr)

data_source <- "Fungal"
# read data

# Fungal
otu_table <- read.csv('Data/Fungal_OTU.csv', row.names = 1)
group_table <- read.csv('Data/Fungal_group.csv', row.names = 1)
tax_table <- read.csv('Data/Fungal_taxonomy.csv', row.names = 1)
env_table <- read.csv('Data/Env.csv', row.names = 1)


micro_data <- microtable$new(sample_table = group_table,
                             otu_table = otu_table,
                             tax_table = tax_table)

# As an example, use 10000 sequences in each sample
micro_data$rarefy_samples(sample.size = 10000)


micro_data$sample_table <- data.frame(micro_data$sample_table, 
                                      env_table[rownames(micro_data$sample_table), ])
# extract two phyla to show the steps
d1 <- clone(micro_data)
d1$tax_table <- d1$tax_table[d1$tax_table$Phylum == "p__Ascomycota", ]
d1$tidy_dataset()
d1$cal_betadiv()

d2 <- clone(micro_data)
d2$tax_table <- d2$tax_table[d2$tax_table$Phylum == "p__Basidiomycota", ]
d2$tidy_dataset()
d2$cal_betadiv()

d3 <- clone(micro_data)
d3$tax_table <- d3$tax_table[d3$tax_table$Phylum == "p__Mortierellomycota", ]
d3$tidy_dataset()
d3$cal_betadiv()

d4 <- clone(micro_data)
d4$tax_table <- d4$tax_table[d4$tax_table$Phylum == "p__Glomeromycota", ]
d4$tidy_dataset()
d4$cal_betadiv()

# 
# d5 <- clone(micro_data)
# d5$tax_table <- d5$tax_table[d5$tax_table$Phylum == "p__Mucoromycota", ]
# d5$tidy_dataset()
# d5$cal_betadiv()

# first perform mantel test
t1 <- trans_env$new(dataset = d1, env_cols = 3:11)
t1$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

t2 <- trans_env$new(dataset = d2, env_cols = 3:11)
t2$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

t3 <- trans_env$new(dataset = d3, env_cols = 3:11)
t3$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

t4 <- trans_env$new(dataset = d4, env_cols = 3:11)
t4$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

# t5 <- trans_env$new(dataset = d5, env_cols = 3:11)
# t5$cal_mantel(use_measure = "bray", partial_mantel = TRUE)

# extract a part of the results 
x1 <- data.frame(spec = "Ascomycota", t1$res_mantel) %>% .[, c(1, 3, 6, 8)]
x2 <- data.frame(spec = "Basidiomycota", t2$res_mantel) %>% .[, c(1, 3, 6, 8)]
x3 <- data.frame(spec = "Mortierellomycota", t3$res_mantel) %>% .[, c(1, 3, 6, 8)]
x4 <- data.frame(spec = "Glomeromycota", t4$res_mantel) %>% .[, c(1, 3, 6, 8)]
# x5 <- data.frame(spec = "Mucoromycota", t5$res_mantel) %>% .[, c(1, 3, 6, 8)]
# rename columns
colnames(x1) <- colnames(x2) <- colnames(x3) <- colnames(x4) <- c("spec", "env", "r", "p.value")
# generate interval data
x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x3 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x4 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# x5 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
#                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# cobine two tables
plot_table <- rbind(x1, x2,x3,x4)
# install ggcor package from (https://github.com/mj163163/ggcor-1)
# or follow the steps (https://chiliubio.github.io/microeco_tutorial/intro.html#github-packages)
library(ggplot2)
library(ggcor)
set_scale()

g1 <- quickcor(t1$data_env, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = plot_table) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#A2A2A288", "#1B9E77")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

g1

# t1
data_file_name1 <- paste0("Figures/MentalTest_Correlation/",data_source,"_MentalTest_correlation_Phylum.csv")
data_abundance1 <- data.frame(X1 = plot_table)
write.csv(data_abundance1, data_file_name1, quote=F)

data_file_name2 <- paste0("Figures/MentalTest_Correlation/",data_source,"_MentalTest_correlation_Phylum.pdf")
ggsave(g1,filename = data_file_name2, width = 11.69, height = 8.27)