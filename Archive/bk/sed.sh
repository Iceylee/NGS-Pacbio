#打印第一行
cat all.sam|sed -n 1p

#将my替换成your  
sed "s/my/your/g" pets.txt >your_pets.txt
#修改文件内容
sed -i "s/my/your/g" pets.txt
#在每行前加#
sed 's/^/#/g' pets.txt
#在每行最后加
sed 's/$/ --- /g' pets.txt

^ 表示一行的开头。如：/^#/ 以#开头的匹配。
$ 表示一行的结尾。如：/}$/ 以}结尾的匹配。
\< 表示词首。 如：\<abc 表示以 abc 为首的詞。
\> 表示词尾。 如：abc\> 表示以 abc 結尾的詞。
. 表示任何单个字符。
* 表示某个字符出现了0次或多次。
[ ] 字符集合。 如：[abc] 表示匹配a或b或c，还有 [a-zA-Z] 表示匹配所有的26个字符。如果其中有^表示反，如 [^a] 表示非a的字符


#去除html中所有tags
# 如果你这样搞的话，就会有问题
$ sed 's/<.*>//g' html.txt
 Understand?
 
# 要解决上面的那个问题，就得像下面这样。
# 其中的'[^>]' 指定了除了>的字符重复0次或多次。
$ sed 's/<[^>]*>//g' html.txt
This is what I meant. Understand?


#只替换第3行的my
sed "3s/my/your/g" pets.txt
#替换第3到6行的my
sed "3,6s/my/your/g" pets.txt
#替换每行的第一个s
sed 's/s/S/1' my.txt
#替换每行的第二个s
sed 's/s/S/2' my.txt
#替换每行的第3个以后的s
sed 's/s/S/3g' my.txt


#如果我们需要一次替换多个模式，可参看下面的示例：（第一个模式把第一行到第三行的my替换成your，第二个则把第3行以后的This替换成了That）

sed '1,3s/my/your/g; 3,$s/This/That/g' my.txt

sed -e '1,3s/my/your/g' -e '3,$s/This/That/g' my.txt
#我们可以使用&来当做被匹配的变量，然后可以在基本左右加点东西。如下所示：
 sed 's/my/[&]/g' my.txt


#使用圆括号匹配的示例：（圆括号括起来的正则表达式所匹配的字符串会可以当成变量来使用，sed中使用的是\1,\2…）
 sed 's/This is my \([^,&]*\),.*is \(.*\)/\1:\2/g' my.txt


# 其中的1i表明，其要在第1行前插入一行（insert）
$ sed "1 i This is my monkey, my monkey's name is wukong" my.txt


# 其中的1a表明，其要在最后一行后追加一行（append）
$ sed "$ a This is my monkey, my monkey's name is wukong" my.txt

# 注意其中的/fish/a，这意思是匹配到/fish/后就追加一行
$ sed "/fish/a This is my monkey, my monkey's name is wukong" my.txt

#c 命令是替换匹配行
sed "2 c This is my monkey, my monkey's name is wukong" my.txt
sed "/fish/c This is my monkey, my monkey's name is wukong" my.txt

#d 删除匹配行
sed '/fish/d' my.txt
sed '2d' my.txt
sed '2,$d' my.txt

# 使用n参数就好了 p打印
sed -n '/fish/p' my.txt

#从一个模式到另一个模式
sed -n '/dog/,/fish/p' my.txt

#从第一行打印到匹配fish成功的那一行
sed -n '1,/fish/p' my.txt

#删除最后4行
head -n -4 myfile.txt > new.txt
#或者 tac是cat的逆转命令
n=4
tac file.txt | sed "1,$n{d}" | tac





