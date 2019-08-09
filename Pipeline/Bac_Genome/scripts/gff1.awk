BEGIN {FS="\t";OFS="\t" }
!/\#/{
        split($9,x,";")
        split(x[1],y,"_")

        printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_%08d\n",$1,$2,$3,$4,$5,$6,$7,$8,name,y[2])
}
/\#/{
        print($0)

}
