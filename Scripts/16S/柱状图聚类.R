
#library
library("phyloseq")
library("ggplot2")
library("dplyr")
library(cowplot)
library(dendextend)
library(ggdendro)
require(dplyr)
require(gdata)

theme_set(theme_light())

############################
#0. Data
############################
data("GlobalPatterns")
physeq = GlobalPatterns

Plot_Rank = "Class" #Kingdom Phylum Class Order Family Genus Species

physeq2 <- physeq %>%
  tax_glom(taxrank = Plot_Rank) %>%                   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% # psmelt:melt to long format (for ggploting)                                         
  filter(Abundance > 0.02) %>%   # 去掉每个样本中最低的2%                     
  arrange(!!sym(Plot_Rank)) #按照Class每类别的字母顺序排序


# rank的每一类针对丰度求和，然后取前20的rank 类别
mainAbundance <- physeq2 %>% 
  group_by(!!sym(Plot_Rank)) %>%
  summarise(max_abun = max(Abundance)) %>% 
  top_n(n=20,wt=max_abun) %>%
  arrange(max_abun)


# 仅保留在前20 rank的OTU行
physeq3 <- filter(physeq2, !!sym(Plot_Rank) %in% mainAbundance[[Plot_Rank]])

#clustering聚类
dend <- as.dendrogram(hclust(dist(with(physeq3, tapply(Abundance, Sample, mean)))))

#为physeq3调整样本顺序，与dend一致(以下方法适用于Sample列是有重复的 162行，26unique)
#用于stack bar作图

physeq3$Sample <- reorder.factor(physeq3$Sample, new.order=labels(dend))

#mainAbundance top20排序
physeq3[[Plot_Rank]] <- reorder.factor(physeq3[[Plot_Rank]], new.order=mainAbundance[[Plot_Rank]])

physeq3 = physeq3 %>%
  arrange(desc(!!sym(Plot_Rank))) %>%
  arrange(Sample)


############################
#1. ggplot dendrogram
############################
dend <- dend %>%
  set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) 

ggd1 <- as.ggdend(dend)
p1 <- ggplot(ggd1, horiz = TRUE, label = F) +  #去掉样品label，因为柱状图也会显示同样的label
  theme(panel.border = element_blank()) #去掉最外的灰色边框





############################
#2.1 根据sample 类型分配颜色 （自动按sample类型个数来分配颜色
############################

n_smp_types <- length(unique(physeq3$SampleType))
col_set <- data.frame(SampleType=unique(physeq3$SampleType),color = colorspace::rainbow_hcl(n_smp_types, c = 70, l  = 50))

physeq4 = physeq3 %>% left_join(col_set,by="SampleType") 

sample_list <- physeq4[!duplicated(physeq4$Sample),]

#根据dend的labels将Sample重新排序，然后得到对应color vector
col_branch = sample_list[match(labels(dend), sample_list$Sample),]$color
col_branch <- as.character(col_branch) #cha格式ggplot才会识别为颜色


############################
#3. Plot stack bar
############################
plot_colors <- c(
    "#c8c8c8","#969696","#505050","#010101","#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
  ,"#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248")

p2 <- ggplot(physeq3, aes_string(fill=Plot_Rank, y="Abundance", x="Sample")) + 
  geom_bar(stat="identity", position="fill") + coord_flip()+
  scale_fill_manual(values = plot_colors) + #设置填充颜色
  theme(axis.title.y = element_blank(),axis.text.y = element_text(colour = as.character(col_branch))) #设置y轴样本label颜色

############################
#4. Combine
############################
# scale调整p1和p2相对高度，使每一个sample对应上 
p <- plot_grid(p1, p2, scale = c(1.05,1),align = "h",rel_widths = c(1.3,5))

ggsave(paste(Plot_Rank,"_cluster.pdf",sep=""),p, width=10, height=7, units="in")





