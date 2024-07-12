# 清除变量
rm(list=ls())
library(vegan)
library(ggplot2)
library(ggrepel)

# sname <- 'Funge'
# # read data
# df<-read.delim('真菌PCoA.txt',row.names=1,sep='\t',header = TRUE, check.names = FALSE)
# # read group
# group<-read.delim('真菌分组.txt',sep='\t',stringsAsFactors = FALSE,check.names = FALSE)


data_source <- "Funal"
# read data
# Bacterial
if(data_source=="Bacteria")
{
  df <- read.csv('Data/Bacteria_OTU.csv', row.names = 1)
  group <- read.csv('Data/Bacteria_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Bacteria_taxonomy.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}else if(data_source=="Funal")
{
  # Fungal
  df <- read.csv('Data/Fungal_OTU.csv', row.names = 1)
  group <- read.csv('Data/Fungal_group.csv', row.names = 1)
  tax_table <- read.csv('Data/Fungal_taxonomy.csv', row.names = 1)
  env_table <- read.csv('Data/Env.csv', row.names = 1)
}


# data transform
df1<-t(df)
# distance
distance <- vegdist(df1, method = 'bray')

# anosim
anosim.result <- anosim(distance, group$group,permutations = 999)

summary(anosim.result)

pdf(paste('Figures/anosim_two/',data_source,'_anosim.pdf', sep = ''), width = 7, height = 5)
plot(anosim.result,col=c('gray','red','green','blue','orange','purple'))
dev.off()

# 

##ANOSIM 分析（使用循环处理，进行小分组间比较，如两组间）
#推荐使用 OTU 丰度表作为输入数据，每次筛选分组后重新计算样本距离，避免由于样本数减少可能导致的距离变动而造成误差
group_name <- unique(group$group)

dir.create('Figures/anosim_two', recursive = TRUE)
anosim_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, group %in% c(group_name[i], group_name[j]))
    otu_ij <- df[group_ij$SampleID]
    otu_ij <- t(otu_ij)
    otu_ij_distance <- vegdist(otu_ij, method = 'bray')
    anosim_result_otu_ij <- anosim(otu_ij_distance, group_ij$group, permutations = 999)     #Bray-Curtis 距离测度，基于 999 次置换
    anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif))
    
    #每次循环输出图片
    pdf(paste('Figures/anosim_two/',data_source,'_anosim.', group_name[i], '_', group_name[j], '.pdf', sep = ''), width = 7, height = 5)
    #png(paste('Figures/anosim_two/',data_source,'_anosim.', group_name[i], '_', group_name[j], '.png', sep = ''), width = 600, height = 400)
    plot(anosim_result_otu_ij, col = c('gray', 'red', 'blue'))
    dev.off()
  }
}

#带 R 值和 p 值的表格
anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors = FALSE)
names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value')

#可选添加 p 值校正过程，例如 Benjamini 校正
anosim_result_two$P_value <- as.numeric(anosim_result_two$P_value)
anosim_result_two$P_adj_BH <- p.adjust(anosim_result_two$P_value, method = 'BH')

write.table(anosim_result_two, paste('Figures/anosim_two/',data_source,'_ANOSIM.result_two.txt', sep = ''),row.names = FALSE,sep = '\t',quote = FALSE,na="")
