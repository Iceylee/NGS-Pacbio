library("phyloseq")
library("ggplot2")
library("dplyr")
require(gdata)
theme_set(theme_light())

data("GlobalPatterns")
physeq = GlobalPatterns

#Kingdom Phylum Class Order Family Genus Species
Plot_Rank = "Class"

# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

physeq2 <- physeq %>%
  tax_glom(taxrank = Plot_Rank) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                        
# Filter out low abundance taxa
  arrange(!!sym(Plot_Rank))
#sort by rank alphabetically

#get top 20 main abundance Ranks
# rank的每一类针对丰度求和，然后取前20的rank 类别
mainAbundance <- physeq2 %>%
  group_by(!!sym(Plot_Rank)) %>% #将string转成object
  summarise(max_abun = max(Abundance)) %>%
  top_n(n=20,wt=max_abun) %>%
  arrange(max_abun)

#根据上表的前20类别，来filter得到仅含前20类别的表格
physeq3 <- filter(physeq2, !!sym(Plot_Rank) %in% mainAbundance[[Plot_Rank]])

#mainAbundance top20排序
physeq3[[Plot_Rank]] <- reorder.factor(physeq3[[Plot_Rank]], new.order=mainAbundance[[Plot_Rank]])
physeq3 = physeq3 %>%
  arrange(desc(!!sym(Plot_Rank))) 

# Set colors for plotting
plot_colors <- c(
    "#c8c8c8","#969696","#505050","#010101","#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
  ,"#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248")

plot_colors2 <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                  "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

# Plot 
#geom_bar()中的position = “fill”  使得bar补齐到1.如果无此参数，上部会留空。
# aes_string的话， colnames均为string
p <- ggplot(physeq3, aes_string(x = "Sample", y = "Abundance", fill = Plot_Rank)) + 
  geom_bar(stat = "identity",position = "fill") +
  scale_fill_manual(values = plot_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance ( > 2%) \n")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))+
  ggtitle(paste(Plot_Rank," Level Bar Plot",sep="")) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(paste(Plot_Rank,"_Abundance.pdf",sep=""),p, width=10, height=7, units="in")
