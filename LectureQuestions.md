---
title: "LectureQuestions"
author: "Xiangyi Xu (Lydia)"
date: "Saturday, January 31, 2015"
output:
  html_document:
    keep_md: yes
---

## Lecture 1a

Which symbol can be used to slicing and extracting data from a vector in R?

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

x is a data frame and z is a feature of x. Which of the following commands is equivalent to with(x, f(z))?

* f(x$z)
* x$f(z)
* f(z)
* It depends.

## Lecture 5b

Simulated coin-tossing can be done using different methods. Which of the following will NOT work?

* coin <- sample(c("H", "T"), 10, replace = F)
* rbinom(10, 1, .5)
* ifelse(rbinom(10, 1, .5) == 1, "H", "T")
* c("H", "T") [1 +rbinom(10, 1, .5)]

## Lecture 6a

A vector x <- 1:10, which of the following choice will NOT insert 1.23 between x[7] and x[8]?

* z <- rbind(x, 1.23, after = 7)
* z <- append(x, 1.23, after = 7)
* z <- c(x[1:7], 1.23, x[8:10]) 
* v <- 1.23; k <- 7; i <- seq(along = x); z <- c(x[i <= k], v, x[i > k])
  
## Lecture 6b

What is the name of the following distribution? 
```{r snippet code}
set.seed(1)
rpois(5, 2)
```

* A vector with the numbers 2, 2, 3, 5, 2
* A vector with the numbers 3.3, 2.5, 0.5, 1.1, 1.7
* It is impossible to tell because the result is random
* A vector with the numbers 1, 4, 1, 1, 5

## Lecture 7a

What R function can be used to generate standard Normal random variables?

* rnorm
* pnorm
* dnorm
* qnorm

## Lecture 7b

When simulating data, why is using the set.seed() function important?

* It can be used to specify which random number generating algorithm R should use, ensuring consistency and reproducibility.
* It can be used to generate non-uniform random numbers.
* It ensures that the random numbers generated are within specified boundaries.
* It ensures that the sequence of random numbers is truly random.

## Lecture 8a

The principal() function will perform a principal componets analysis in R, starting with a matrix. The format is as the following. Which of the following decription is NOT correct regarding the parmaters?
```r
  principal(r, nfactors=, rotate=, scores=)
```

* r is a covariance matrix or a raw data matrix 
* nfactors specifies the number of principal components to extract (1 by default)
* rotate indicates the rotation to be applied (varimax by default)
* scores specifies whether or not to calculate principal component scores (false by default)

## Lecture 8b

Which criteria is correct for deciding how many components to retain in a PCA? 

* All of these choices.
* Basing the number of components on prior experience and theory. 
* Selecting the number of components needed to account for some threshold cumulative amount of variance in the variables.
* Selecting the number of components to retain by examing the eigenvalues of the k*k correlation matrix among the variables.

## Lecture 9a

Assume that `library(ggplot2)`has been loaded and mtcars is its built-in database. Which of the following code will NOT achieve the purpose as the other three?  

* plot(wt~mpg, data=mtcars) 
* plot(mtcars$wt, mtcars$mpg)
* qplot(mtcars$wt, mtcars$mpg)
* ggplot(mtcars, aex(x=wt, y=mpg)) + geom_point()

## Lecture 9b

Assume that `library(ggplot2)` has been loaded and database pressure is built-in. Which of the following 2 codes are equivalent?
```{r graphs}
1. qplot(temperature, pressure, geom="lines")
2. ggplot(presure, aes(x=temperature, y=pressure)) + geom_line()
3. qplot(temperature, pressure, data=pressure, geom=c("line", "point"))
4. ggplot(pressure, aex(x=temperature, y=pressure)) + geom_line() + geom_point()
```

* 1 and 2
* 1 and 3
* 2 and 2
* 2 and 4
