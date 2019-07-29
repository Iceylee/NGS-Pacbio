git config --global user.name "liumiaocn"
git config --global user.email "liumiaocn@outlook.com"

git status
git add *  #use "git add <file>..." to update what will be committed
git commit -m "message..." #将所有修改内容 提交

git push -u origin master


#查看修改
git log --pretty=oneline test.R
git show HASHID  #具体查看某次commit的修改内容

#冲突
git pull origin master

#分支
#分支尽量会控制在三条以下，其本质则是将持续集成和LeanIT进行扩展，减少中间环节的WIP，尽可能早集成，分支的原则也可以遵循“用后即焚”，删除分支使用–delete即可。
git branch #查看分支
git branch develop #创建分支
git checkout develop #切换到develop分支
git push -u origin develop #
git push origin --delete develop #删除分支