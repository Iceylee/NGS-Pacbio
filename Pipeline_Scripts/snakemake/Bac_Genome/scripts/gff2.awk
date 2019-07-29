BEGIN {FS="\t";OFS="\t" }
!/\#/{
        split($9,x,";")
        i = i + 1

        printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_%08d\n",$1,$2,$3,$4,$5,$6,$7,$8,name,i)
}
/\#/{
        print($0)

}
