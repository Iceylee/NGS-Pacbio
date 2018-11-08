#fastqc是fq.gz和fq两种格式都可以 
for i in `ls *R1_001.fastq`
do
i=${i/R1_001.fastq/}
echo "/opt/FastQC/fastqc -o ../fastqc -t 10 ${i}R1_001.fastq ${i}R2_001.fastq"
done




#脚本
screen -R icey1 #新窗口才会开始跑“step3后开始显示分析进程”
bash /data/script/Step1_FastQC.sh /DataBackUp/2018_02/WangNian/raw_data /DataBackUp/2018_02/WangNian/out

bash /data/script/Step1_FastQC.sh /DataBackUp/2018_02/YanJuan/raw_data /DataBackUp/2018_02/YanJuan/fastqc/raw
Q30_analy.py .




#gff+fasta 转 gbk
python ~/analysis/scripts/gff_to_genbank.py unitig_${i}.gff unitig_${i}.fasta
python ~/Documents/GitHub/python/gff_to_gbk.py 参考序列/Danio_rerio.GRCz10.86.gff3 genome.fa
#加入蛋白信息
awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' Brevibacillus_sp.WF146.pep > Brevibacillus_sp.WF146.1line.pep

cat Brevibacillus_sp.WF146.1line.pep |grep -v ">" > WF146.pep.list

awk '
NR==FNR{i=i+1;a[i]=$1;next}
/CDS/{
	print $0;
	p=p+1;
	print "                     /translation=\"" a[p] "\""
}
!/CDS/{print $0}
'   WF146.pep.list Brevibacillus_sp.WF146.gbk |less




#8个0改成6个0
cat Brevibacillus_sp.WF146.pep.gbk|awk '
/ID=/{
	match($0,/[0-9][0-9][0-9][0-9]\"/,a)
	print "                     /ID=\"Brevibacillus_sp.WF146_00" a[0]  
}
!/ID=/{print $0}
' >Brevibacillus_sp.WF146.pep.gbk.6

awk '
NR==FNR{p[$1]=$2;next}
/ID=/{
	print $0
	match($0,/\".*\"/,a)
	gsub(/"/,"",a[0])
	ii = a[0]
	if(ii in p){
		print "                     /product=\"" p[ii] "\""
	}
	
}
!/ID=/{print $0}
'   gene_nr.list Brevibacillus_sp.WF146.pep.gbk.6 |less

#####

cat Thermoactinomyces_sp.CDF.gbk|awk '
/ID=/{
	match($0,/[0-9][0-9][0-9][0-9]\"/,a)
	print "                     /ID=\"Thermoactinomyces_sp.CDF_00" a[0]  
}
!/ID=/{print $0}
' >Thermoactinomyces_sp.CDF.gbk.6

awk '
NR==FNR{p[$1]=$2;next}
/ID=/{
	print $0
	match($0,/\".*\"/,a)
	gsub(/"/,"",a[0])
	ii = a[0]
	if(ii in p){
		print "                     /product=\"" p[ii] "\""
	}
	
}
!/ID=/{print $0}
'   gene_nr.list Thermoactinomyces_sp.CDF.gbk.6 |less


