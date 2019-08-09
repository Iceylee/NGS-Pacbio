
df1=read.table("Group1.list",stringsAsFactors = FALSE)
Group1=df1[,1]
df2=read.table("Group2.list",stringsAsFactors = FALSE)
Group2=df2[,1]
T<-venn.diagram(list("10w-300"=Group1,"10w-con"=Group2),
                filename=NULL,
                lwd=1,                         #圈线粗度
                lty=1,                         #圈线类型
                #col=c('#0099CC','#FF6666'),    #圈线颜色
                fill=c("red","purple"),   #填充颜色
                #cat.col=c('#0099CC','#FF6666'),#A和B的颜色
                cat.cex = 2,                 #A和B的大小
                #rotation.degree = 45,          #旋转角度
                #main = "A&B",                  #主标题内容
                main.cex = 2,                  #主标题大小
                #sub = "plot : example",        #亚标题内容
                sub.cex = 1,                   #亚标题字大小
                cex=1.5,                       #里面交集字的大小
                alpha = 0.5, 
                cat.fontface = 2,
                reverse=TRUE,
                cat.just=list(c(0,0),c(1,0))) #change label position
pdf("venn_plot.pdf")
grid.draw(T)
dev.off()



####
v <- venn.diagram(list(shams_90d = 1:3, shams_90d_4h = 2:4, sham3__shams_90d = 3:5,
                       sham3_90d__shams = 5:7, sham3_90d__shams_4h = 6:9),
                  fill = c("red", "green", "blue", "yellow", "purple"),
                  alpha = c(0.5, 0.5,0.5, 0.5, 0.5), cex = 1,cat.fontface = 2,
                  lty =1, filename=NULL, cat.cex=0.8, 
                  cat.just=list(c(0.6,1) , c(0,0) , c(0,0) , c(1,1) , c(1,0)))

grid.newpage()
grid.draw(v)


