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

In the following codes, when m and n equal to what values, will the plot show a quarter of a circle?
```r
N <- 10000
x <- runif(N, min=m, max=n)
y <- runif(N, min=m, max=n)
plot(x, y, pch=16, col=ifelse(x^2 + y^2<1, "red", "blue"))
```

* m=0,    n=0
* m=-1.0, n=1.0
* m=-2.0, n=2.0
* m=-3.0, n=3.0

## Lecture 2b

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

* 202
* 441
* 40
* 200

## Lecture 3a

What is the name of the following distribution? 
```{r uniform_distribution, echo=FALSE}
x <- seq(0,1,length=200)
plot(c(-0.1, 0, x, 1, 1.1), c(0,0,dunif(x),0,0), type='l', ylab="P")
```

* Uniform distribution
* Normal distribution
* Poisson distribution
* Binominal distribution

## Lecture 3b

persp() can be used to produce a three-dimentional plot. What is the function of its arguments theta and phi?
```r
  persp(x, y, fa, theta=30, phi=20)
```

* Control the angles at which the plot is viewed. 
* Control the width and height of the plot
* Control the size and dimension of the plot
* Control the length and width of the plot

## Lecture 4a

Consider the following function, what value will be returned?
```{r function environment}
f <- function(x){
        g <- function(y){
                y+z
        }
        z<- 4
        x + g(x)
}

z<-10
f(3)
```

* 10 			
* 7 			
* 16 	
* 4

## Lecture 4b

The following code will produce a warning in R. Please explain why?
```r
x <- 1:10
if(x > 5) {
        x <- 0
}
```

* 'x' is a vector of length 10 and 'if' can only test a single logical statement. 	 
* use 'x' is a vector and 0 is a scalar. 			
* There are no elements in 'x' that are greater than 5 			
* The expression uses curly braces. 

## Lecture 5a

Which of the following commands is equivalent to with(x, f(z))?

* x$f(x$z)
* f(x$z)
* x$f(z)
* f(z)
* It depends.

## Lecture 5b

