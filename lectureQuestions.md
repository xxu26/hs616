---
title: "LectureQuestions"
author: "Xiangyi Xu (Lydia)"
date: "Saturday, January 31, 2015"
output:
  html_document:
    keep_md: yes
---


## Lecture 1a

Which symbol can be used to slicing and extracting data from a vector?

* `[ ]`
* `[[c( ) ]]`
* `$`
* `[, c( )]`

## Lecture 1b

In the following codes, what is the type of the variable returned?
```
y <- c(5, 6, 7, 8, NA)
is.na(y)
```

* logical
* numeric
* integer
* character

## Lecture 2a

In the following codes, when a and b equal to what values, will the plot show a quarter of a circle?
```r
N <- 10000
x <- runif(N, min=a, max=b)
y <- runif(N, min=a, max=b)
plot(x, y, pch=16, col=ifelse(x^2 + y^2<1, "red", "blue"))
```
![plot of chunk throwing_darts](01b_literate_programming-figure/throwing_darts.png) 

* a=0,    b=0
* a=-1.0, b=1.0
* a=-2.0, b=2.0
* a=-3.0, b=3.0




## Lecture 2b



