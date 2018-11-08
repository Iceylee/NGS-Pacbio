Rscript /Users/Icey/Documents/GitHub/R/run.R


ln -s CountMatrix4DESeq.csv ID.list


for i in `ls *_GO*.txt`
do
bash ~/Documents/GitHub/R/ID_convert.sh 1_GO_CC_out.txt go_ID.list
done


#if $i_ID_type exist
if [  -f *_ID_type.csv ]; then
    awk '
	NR==FNR{FS=",";a[$1]=$2;next}
	{OFS=",";print a[$1],$2}' *_ID_type.csv ID.list > temp
	mv temp kegg_ID.list
fi

bash ~/Documents/GitHub/R/ID_convert.sh *KEGG*.txt kegg_ID.list


