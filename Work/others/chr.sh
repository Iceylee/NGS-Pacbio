for chr in `seq 1 22` X Y
do

	for smp in A1 A2 A3 A1_input A2_input A3_input
	do
		for i in plus minus
		do

			cat ../${smp}_${i}_depth.txt|awk -v chr="$chr" '
				BEGIN{FS=OFS="\t"}
				$1==chr {print $0}' > temp1

			cat temp1|awk '
				BEGIN{FS=OFS="\t"}
				{loc=int($2/100000);a[loc]=a[loc]+$3}
				END{for (i in a) {print i,a[i]}}' > temp2

			cat temp2 >> ${smp}_${i}_win.txt
		done

	done

done 

#问题：统计应以固定窗口来计算。比如1-1000，1001-2000
#没有输出染色体编号
