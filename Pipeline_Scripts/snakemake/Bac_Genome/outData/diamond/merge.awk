BEGIN { FS=OFS="\t" }
{
    curr = $1
    if (curr == prev) {
        rec = rec ";" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }
