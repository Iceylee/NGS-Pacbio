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
        printf (">PseudomonasPutidaA316_%06d\tlocus=%d:%d:%s\n",i,start,end,chain)}
!/>/ {print $0}