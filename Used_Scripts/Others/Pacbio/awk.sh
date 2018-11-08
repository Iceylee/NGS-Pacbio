#3列不为0的输出名字和报酬
awk 'BEGIN{FS=" "}$3 > 0 { print $1, $2 * $3 }' emp.data

