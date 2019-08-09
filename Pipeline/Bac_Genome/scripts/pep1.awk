BEGIN { FS = " # " ; OFS="\t" }
/>/ {

        i = i + 1
        start = $2
        end = $3
        if ($4 == 1){
                chain = "+"
        }
        else{
                chain = "-"
        }
        printf (">%s_%06d locus=%d:%d:%s\n",name,i,start,end,chain)}
!/>/ {print $0}
