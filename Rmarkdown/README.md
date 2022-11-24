
---
#YAML header: a set of key   
title: "Hello"
author: "Me"
date: date
output: _document (i.e., html_document, pdf_document, word_document)
---

# Header 1
## Header 2
### Header 3
End a line iwth two spaces to start new paragraph.   
*italics*
superscript^2^
~~strikethrough~~
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
