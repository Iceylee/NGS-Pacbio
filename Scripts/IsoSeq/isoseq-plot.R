setwd("~/work/2018/0507玉米isoseq/rawData/plot/")
library(ggplot2)
library(dplyr)
curve = read.table(file="saturation_curve.txt", header=T,stringsAsFactors = F)
options(scipen=200) #avoid 1e10
p <- ggplot(data=curve, aes(x=FL_reads_number, y=Isoform_number, group=1)) +
  geom_line()
#   +scale_x_continuous(breaks=seq(0,200000,20000))+
#   scale_y_continuous(breaks=seq(0,120000,10000)) 
pdf(file="reads_isofrom.pdf")
p
dev.off()

p <- ggplot(data=curve, aes(x=FL_reads_number, y=Gene_number, group=1)) +
  geom_line()
#   +scale_x_continuous(breaks=seq(0,200000,20000))+
#   scale_y_continuous(breaks=seq(0,10000,1000)) 
pdf(file="reads_gene.pdf")
p
dev.off()



#splice
#数据只保留五个事件
as = read.table(file="event_count.out", header=F,stringsAsFactors = F)
colnames(as) <- c("AS_Event","number")
as$number =  as.numeric(as$number)

as <- as %>% 
  mutate(per=`number`/sum(`number`)) %>% 
  arrange(desc(AS_Event))
as$label <- scales::percent(as$per)

p <- ggplot(data=as)+
  geom_bar(aes(x="", y=per, fill=AS_Event), stat="identity", width = 1)+
  coord_polar("y", start=0)+
  theme_void()+
  geom_text(aes(x=1, y = cumsum(per) - per/2, label=label))



ggsave("splice.as.pdf",p, width=7, height=5, units="in")


#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置 