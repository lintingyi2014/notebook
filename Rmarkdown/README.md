Need to install rmarkdown at R  
`install.packages("rmarkdown")`

#YAML   
~~enclose with~~ ---  
title: "Hello"  
author: "Me"  
date: "date"  
rmd_files: ["file.Rmd", "file2.Rmd"]  
output: html_document   
~~Need to enclose the above with~~ ---

# Header 1
## Header 2
### Header 3
End a line with two spaces to start new paragraph.   
*italics*
superscript^2^
~~strikethrough~~
**bold**
[link] (www.) 
> block quote

* unordered list  
  + sub-item  
  
```{r}
code(hi)
```
 #{r eval=TRUE} evaluates code and display results  
 #{r echo=TRUE} display code with results   
 #{r error=TRUE} display error   

#Render 
Run rmarkdown::render("<file path>")

#Configure with Git, usign R console
```
system("git config --global user.name gituser")  
system("git config --global user.email gitemail")  
```
#Link R with git @R console  
`system("git remote add origin https://github.com/[USER]/[repo]")`

#at Rstudio   
`system("git init")`
At the project options popup choose the origin to be ur remote github repo
```
git pull: fetch and merge synchronizing local repo 
git push: make changes to remote repo
