#slice
slice(mtcars,20)
slice(mtcars,1:20)

#lead  新的vector的每个元素 index+1 ，
b = lead(a)
则 b[i] = a[i+1]

a <- 1:10
b <- lead(a)
a
b


# 将列"a,b,c"取第一个a存入新列
tools <- tools  %>% 
    mutate(work_tools = str_split(WorkToolsSelect,pattern=",")) %>%
    unnest(work_tools)

# 计算每组行数
tool_count <- tool_count  %>% 
    group_by(work_tools)  %>% 
    summarise(n = n())

#作图
# 1.x列为各组名，y列为各组数目即bar高度（这样画bar需要stat="identity"）
# 2.reorder 按大小排序bar
# 3.轴标签转90度
ggplot(tool_count, aes(x=reorder(work_tools,n),y=n)) + 
    geom_bar(stat="identity") +
    theme(axis.text = 
         element_text(angle = 90, vjust = 0.5, hjust = 1))

#case_when 
#grepl：string是否含某个substring，返回boolean
debate_tools <- debate_tools  %>% 
   mutate(language_preference = case_when(
       grepl("R",WorkToolsSelect) & !grepl("Python",WorkToolsSelect) ~ "R",
       !grepl("R",WorkToolsSelect) & grepl("Python",WorkToolsSelect) ~ "Python",
       grepl("R",WorkToolsSelect) & grepl("Python",WorkToolsSelect) ~ "both",
       !grepl("R",WorkToolsSelect) & !grepl("Python",WorkToolsSelect) ~ "neither"))

# row_number和arrange联用，再用filter可取每个分组的n值的top4
recommendations <- recommendations  %>% 
    group_by(language_preference,LanguageRecommendationSelect)  %>% 
    summarise(n=n()) %>%
# Removing empty responses and include the top recommendations
    filter(!is.na(LanguageRecommendationSelect)) %>%
    arrange(-n) %>%
    mutate(row = row_number() )%>%
    filter(row <= 4)

#jupyter R输出图片的大小 设置
library(repr)
options(repr.plot.width=12, repr.plot.height=6)


#ggplot plot data
loess_data <- ggplot_build(p3)$data
dot <- loess_data[[1]]
black_curve <- loess_data[[2]]
green_curve <- loess_data[[3]]

#sliding window
SNP_index3  <- SNP_index3 %>%
  group_by(chrom)%>%
  mutate(avg = roll_mean(SNP_Index, step, na.rm=TRUE, align="left", fill = NA))%>%
  ungroup()

#根据向量筛选data frametop_dest <- flights %>%
top_dest <- flights %>%
  count(dest, sort = TRUE) %>%
  head(10)

flights %>% 
  filter(dest %in% top_dest$dest)

#color & legend 
ggplot(mtcars) + 
  geom_bar(aes(factor(hp), fill=factor(hp))) + 
  scale_fill_manual(values = getPalette(22),
                    guide = guide_legend(nrow=2)) +
  theme(legend.position="bottom")

#根据vector来取特定rows
dt[dt$fct %in% vc,]


#按照tvhours的大小对relig的levels进行排序。这样画图的label就是这样的顺序
relig_summary %>%
  mutate(relig = fct_reorder(relig, tvhours))
ggplot(relig_summary, aes(tvhours, relig)) + geom_point()


#df 提取一列但保留rownames
 A[1,,drop=FALSE]


#string 传为object
Plot_Rank = "Phylum"
arrange(df,!!sym(Plot_Rank))
#等同于
arrange(df,Phylum)

##
df[[Plot_Rank]]

#批量读入文件 
a = list.files("./rawData")                                                       
dir = paste("./rawData/",a,sep="")                                      
n = length(dir)                                                                 
result = read.csv(file = dir[1],header=T,sep="\t",stringsAsFactors=F)   #读入第一个文件内容

for (i in 2:n){
  new.data = read.csv(file = dir[i], header=T, sep="\t",stringsAsFactors=F)
  
  result = left_join(result,new.data,by="geneID")
}

write.table(result,"count.txt" ,row.names = F,quote = F,sep="\t")

#匹配部分文件

path1 = paste(output_path,"3.DiffExprGene",sep="")
sigfiles = list.files(path1,pattern="sig_genes_exprData.txt")

#函数 多参数
sapply(sigfiles,GOKEGG,pSet=0.05)
sapply(sigfiles,GOKEGG,pSet=1)

#bash
system("ls")

#substring
string = c("G1:E001", "G2:E002", "G3:E003")

1) sub

sub(".*:", "", string)
## [1] "E001" "E002" "E003"
2) strsplit

sapply(strsplit(string, ":"), "[", 2)
## [1] "E001" "E002" "E003"
3) read.table

read.table(text = string, sep = ":", as.is = TRUE)$V2
## [1] "E001" "E002" "E003"
4) substring

This assumes second portion always starts at 4th character (which is the case in the example in the question):

substring(string, 4)
## [1] "E001" "E002" "E003"
4a) substring/regex

If the colon were not always in a known position we could modify (4) by searching for it:

substring(string, regexpr(":", string) + 1)
5) strapplyc

strapplyc returns the parenthesized portion:

library(gsubfn)
strapplyc(string, ":(.*)", simplify = TRUE)
## [1] "E001" "E002" "E003"



#得到class中每组的count的前10

go_id_top <- go_id %>%
  group_by(class) %>%
  top_n(n = 10, wt = count) %>%
  arrange(class)


#change F,C,P to MF,CC,BP
levels(go_id$class)[levels(go_id$class)=="C"]<- "CC"
levels(go_id$class)[levels(go_id$class)=="F"]<- "MF"
levels(go_id$class)[levels(go_id$class)=="P"]<- "BP"

ids_du <- ids[!duplicated(ids$UNIGENE),]



conda install r
conda install -y bioconductor-deseq bioconductor-deseq2 bioconductor-edger r-gplots


source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite("DESeq2")
biocLite("edgeR")

#修改levels顺序，按照表格的顺序（逆序），而不是字母表。作图的轴线出现的顺序会受影响	
ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)

p <- p + theme(panel.spacing = unit(0, "lines"))

ggplot
p+theme(axis.text = element_text(size = 8，face="bold"), legend.text = element_text(size = 8)) 





###可以按照表格显示的term顺序作图 ---序列是从小到大的
##如何让表格显示的顺序是BP CC MF呢（level需要修改）
ego_three$term<- factor(ego_three$term, order=TRUE, levels=rev(ego_three$term))
# ego_three$class<- factor(ego_three$class, order=TRUE, levels=c("MF","CC","BP"))



#分开string并计数
p = c('a|b|c|d','b|c|c','x|y')
pp = strsplit(as.character(p),'|',fixed=TRUE)
count = unlist(lapply(pp,length))


##分栏
#1.
p = c('a|b|','b|c','x|y')
pp = strsplit(as.character(p),'|',fixed=TRUE)
ppp = do.call('rbind',pp ) #function name; a list -output a list 函数只能接受两个参数，list中上一个和当前。
foo <- data.frame(ppp)
#2.
within(df, FOO<-data.frame(do.call('rbind', strsplit(as.character(FOO), '|', fixed=TRUE))))
#3.
require(reshape)
df <- data.frame(ID=11:13, FOO=c('a|b','b|c','x|y'))
df = transform(df, FOO = colsplit(FOO, split = "\\|", names = c('a', 'b')))

#
#引号中为一个object，该object是否存在
exists("ego_CC_df")


###研究
#诉求：三个df，检测是否存在。rbind合并存在的几个df
#ego_list <- list(BP=ego_BP_df,CC=ego_BP_CC,MF=ego_BP_MF)
ego_list <- c("ego_BP_df", "ego_MF_df", "ego_CC_df")
TorF <- sapply(ego_list,exists)
real_list <- sapply(ego_list[TorF],as.name)
ppp = do.call('rbind',real_list )


##remove colon:
gene_id <- sapply(gene_id,function (x) {unlist(strsplit(x,split=":"))[2]})


##变量1转列名 变量2转行名 统计变量3
library(reshape2)
rpkm_df <- read.table(file = "AllSamplesRPKMValue.txt",sep = '\t',header =F,stringsAsFactors=F)
dff <- rpkm_df[,1:3] #只能有三列
rpkm_tab <- dcast(dff, V2 ~ V1) #V2为左侧列 V1为最上的行
#转回去是melt和cast

rep(seq(1, 7, by = 2), times = 7)

#正则
# The emails vector has already been defined for you
emails <- c("john.doe@ivyleague.edu", "education@world.gov", "dalai.lama@peace.org",
            "invalid.edu", "quant@bigdatacollege.edu", "cookie.monster@sesame.tv")
# Use grepl() to match for "edu"   [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
grepl("edu",emails)
# Use grep() to match for "edu", save result to hits  [1] 1 2 4 5
grep("edu",emails)   
# Subset emails using hits
emails[grepl("edu",emails)]
emails[grep("edu",emails)]
#匹配@
#匹配任何字符多次 .*
#反义圆点 \\
grepl("@.*\\.edu$",emails)
# Use sub() to convert the email domains to datacamp.edu
sub("@.*\\.edu$","@datacamp.edu",emails)

#报错的情况赋值
qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

# Print the R version details using version
version
# Assign the variable `major` to the major component
major <-  version$major
# Assign the variable `minor` to the minor component
minor <- version$minor

# How long does it take to read movies from CSV?
system.time(read.csv("movies.csv"))
# How long does it take to read movies from RDS?
system.time(readRDS("movies.rds"))
#等号和箭头的区别 在函数内部使用箭头，变量也被赋值。 
median(x <- 1:10); x
#比较函数运行的时间
# Load the microbenchmark package
library(microbenchmark)
# Compare the two functions，10即各运行10次
compare <- microbenchmark(read.csv("movies.csv"), 
                          readRDS("movies.rds"), 
                          times = 10)
# Print compare
print(compare)

##function
#for循环不好 如果df是空 将会报错
df <- data.frame()
1:ncol(df)

for (i in 1:ncol(df)) {
  print(median(df[[i]]))
}
##用seq_along
for (i in seq_along(df)) {
  print(median(df[[i]]))
}
##将df的结果存在output里面
# Create new double vector: output
output = vector("double", length = ncol(df))
#
for (i in seq_along(df)) {
  output[[i]] = median(df[[i]])
}
# Print output
print (output)

##
# Define example vectors x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3,  4)
# Count how many elements are missing in both x and y
# x和y都是NA的，计1.sum求这样的个数
sum(is.na(x) & is.na(y))

rep
seq
