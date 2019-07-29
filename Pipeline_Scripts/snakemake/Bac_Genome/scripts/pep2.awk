BEGIN { FS = " # " ; OFS="\t" }
/>/ {

        i = i + 1
        split($1,x,"_")
        contig = "unitig_"x[2]
        start = $2
        end = $3
        if ($4 == 1){
                chain = "+"
        }
        else{
                chain = "-"
        }
        printf (">%s_%06d locus=%s:%d:%d:%s\n",name,i,contig,start,end,chain)}
!/>/ {print $0}
