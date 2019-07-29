BEGIN{FS=OFS="\t"}
$9~/structure \"0,1-2\^\"/ {
	i["Exon_Skipping"]++
}
$9~/structure \"1\^,2\^\"/ {
	i["Alternative_Donors"]++
}
$9~/structure \"1-,2-\"/ {
	i["Alternative_Acceptors"]++
}
$9~/structure \"0,1\^2-\"/ {
	i["Intron_Retention"]++
}
$9~/structure \"1-2\^,3-4\^\"/ {
	i["Mutually_Exclusive_Exons"]++
}
$9!~/structure \"0,1-2\^\"/ && $9!~/structure \"1\^,2\^\"/ && $9!~/structure \"1-,2-\"/ && $9!~/structure \"0,1\^2-\"/ && $9!~/structure \"1-2\^,3-4\^\"/ {
	i["Other"]++
}
END{for (event in i){print event,i[event]}}
