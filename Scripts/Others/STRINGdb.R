STRINGdb$methods()
#查看帮助
STRINGdb$help("map")
#注意版本，看最新支持到哪个版本
#默认可信度是400（属于中等可信度）
string_db <- STRINGdb$new(version = "10", species = 10116, score_threshold = 700, input_directory = "")

#map函数做ID转化
data_mapped <- string_db$map(data, "geneid", removeUnmappedRows = TRUE)

hit <- data_mapped$STRING_id

#还是需要下载string数据库中的互作连接数据，这个比较大，需要等待下的
string_db$plot_network(hit)

#将互作网络的信息导出
info <- string_db$get_interactions(hit)
#将info数据库输出到txt文件中，然后再导入cytoscape即可作图了




