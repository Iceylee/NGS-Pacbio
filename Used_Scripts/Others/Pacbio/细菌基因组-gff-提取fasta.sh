cat sRNA.gff|awk '
BEGIN{FS=OFS="\t"}
FNR!=1{
	if ($5=="+") {start=$3;end=$4}
	if ($5=="-") {start=$4;end=$3}
	print $1,"rfam",$2,start,end,"-",$5,$6}
FNR==1{print $1,"Source",$2,$3,$4,"Score",$5,$6}

' > temp.gff

mv temp.gff sRNA.gff 

bedtools getfasta -fi ../assembly.fasta -bed sRNA.gff > sRNA.fasta
