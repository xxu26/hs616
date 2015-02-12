---
title: "LectureQuestions"
author: "Xiangyi Xu (Lydia)"
date: "Saturday, January 31, 2015"
output:
  html_document:
    keep_md: yes
---


### Lecture 1a

Which symbol can be used to slicing and extracting data from a vector?

a. `[ ]`
b. `[[c( ) ]]`
c. `$`
d. `[, c( )]`

### Lecture 1b

In the following codes, what is the type of the variable returned?
```
y <- c(5, 6, 7, 8, NA)
is.na(y)
```

a. logical
b. numeric
c. integer
d. character

### Lecture 2a

In the following codes, when m and n equal to what values, will the plot show a quarter of a circle?
```r
N <- 10000
x <- runif(N, min=m, max=n)
y <- runif(N, min=m, max=n)
plot(x, y, pch=16, col=ifelse(x^2 + y^2<1, "red", "blue"))
```

a. m=0,    n=0
b. m=-1.0, n=1.0
c. m=-2.0, n=2.0
d. m=-3.0, n=3.0


### Lecture 2b

What does the following function return?

```r
f <- function(x) {
  f <- function(x) {
    f <- function(x) {
      x ^ 2
    }
    f(x) + 1
  }
  f(x) * 2
}
f(10)
```
a. 202
b. 441
c. 40
d. 200


### Lecture 3a

What is the name of the following distribution? 

```{r uniform_distribution, echo=FALSE}
x <- seq(0,1,length=200)
plot(c(-0.1, 0, x, 1, 1.1), c(0,0,dunif(x),0,0), type='l', ylab="P")
```
a. Uniform distribution
b. Normal distribution
c. Poisson distribution
d. Binominal distribution


### Lecture 3b

persp() can be used to produce a three-dimentional plot. What is the function of its arguments theta and phi?

```r
  persp(x, y, fa, theta=30, phi=20)
```
a. Control the angles at which the plot is viewed. 
b. Control the width and height of the plot
c. Control the size and dimension of the plot
d. Control the length and width of the plot
