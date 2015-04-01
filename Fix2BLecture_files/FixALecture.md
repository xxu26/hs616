# Fix-A-Lecture
Xiangyi Xu  
Thursday, March 26, 2015  

##Objectives


* Describe what a probability distribution function represents, and what it can be used for.
* For the following probability distributions, say whether they are continuous or discrete, and describe the kinds of systems they can be used to model:
    - `Binomial`
    - `Poisson`
    - `Normal`
    - `Uniform`
* Call an R function to generate a vector of random numbers from a normal distribution with a given mean and standard deviation.
* Describe the conditions under which the Normal distribution closely approximates the binomial or Poisson distributions.







##Questions to think About:


* Name two commonly used `Probability Density Functions` 
  (`PDF` describes the relative likelihood for a `CONTINUOUS` random variable to take a given value. )
  
* Name two commonly used `Probability Mass Functions`
  (`PMF` gives the probability that a `DISCRETE` random variable is exactly equal to some value.)
  
  





##Answer:


Distribution | Graph  | Type | Parameter 1 | Parameter 2
------------ | ------ | ---- | ----------- | ------------
Normal | ![normal](https://docs.oracle.com/cd/E12825_01/epm.111/cb_user/images/graphics/dist_normal.gif) | continuous | mean value | Standard Deviation Value
Uniform | ![uniform](https://docs.oracle.com/cd/E12825_01/epm.111/cb_user/images/graphics/dist_uniform.gif) |   continuous | minimum value | maximum value
Binomial | ![binomial](https://docs.oracle.com/cd/E12825_01/epm.111/cb_user/images/graphics/dist_binomial.gif) | discrete | probability (between 0 and 1) | trials 
Poisson | ![poisson](https://docs.oracle.com/cd/E12825_01/epm.111/cb_user/images/graphics/dist_poisson.gif)|   discrete | rate value | NA







##R Random


####Generation of normals:


```r
xn <- rnorm(20) #You will generate 20 standard normals
xn
```

```
##  [1]  0.37219590  2.54127440 -0.47280381  1.85010413 -1.03705456
##  [6] -1.72274577  0.73269936  0.26983711  0.49344343 -0.38928618
## [11] -0.19472442  0.41972985 -0.39308796 -0.24275188  0.36008322
## [16] -0.01413667 -1.61463166  1.34272451 -0.69434328 -0.85772717
```
  
  

####Two types of uniform:


#####**1. Continuous uniform:**


You can have a distribution that has all numbers in some range to be equally likely - a continuous uniform. There are additional arguments to change the range between 0~1.





```r
x_con <- runif(20)
x_con
```

```
##  [1] 0.67944997 0.95976766 0.27795541 0.33371281 0.47875004 0.56323360
##  [7] 0.06623327 0.44593199 0.98745675 0.74940094 0.16168211 0.77478204
## [13] 0.64491569 0.59863940 0.38702334 0.29603127 0.12918492 0.95139099
## [19] 0.85952607 0.09300235
```



#####**2. Discrete uniform:**


You can have a distribution that is equally likely for some finite set of objects, such as a range of integers - a discrete uniform. sample() can generate uniformaly from some set of integers (or other types of objects).





```r
x_dis <- sample(1:20, 5, replace=TRUE)
x_dis
```

```
## [1]  1  5 10  5  9
```







###In R, probability functions take the form 


######`[dpqr]distribution_abbreviation()` 
where the first letter refers to the aspect of the distribution returned:

Letter | Distribution
------ | ------------
d | density
p | cumulative probability function
q | quantile function
r | random generation(random deviates)





Name | Distribution Functions
---- | ----------------------
Binomial | rbinom dbinom pbinom qbinom
Poisson | rpois dpois ppois qpois
Normal | rnorm dnorm pnorm qnorm
Uniform | runif dunif punif qunif






##The Central Limit Theorem


_I know of scarcely anything so apt to impress the imagination as the wonderful form of cosmic order expressed by the "Law of Frequency of Error". The law would have been personified by the Greeks and deified, if they had known of it._                            
--[Sir Francis Galton, 1889](http://en.wikipedia.org/wiki/Central_limit_theorem#History)






##The Galton Box

![Galton Box](http://upload.wikimedia.org/wikipedia/commons/7/7f/Quincunx_%28Galton_Box%29_-_Galton_1889_diagram.png)





##The Normal Distribution

$$ \frac{1}{\sigma \sqrt{2 \pi}} e^{\frac{-(x - \mu)^2}{2 \sigma^2}} $$


```r
normal_pdf <- function(x, mu=0, sigma=1)
    1/(sigma * sqrt(2*pi)) * exp(-(x - mu)^2/(2*sigma^2))

integrate(normal_pdf, -Inf, Inf)
```

```
## 1 with absolute error < 9.4e-05
```





##The Normal Distribution


```r
library(manipulate)
x <- seq(-20, 20, by=0.1)
manipulate(plot(x, dnorm(x, mean=m, sd=s), type="l", col="red"), 
           m = slider(-20.0, 20.0, initial=0),
           s = slider(0.01, 10, initial=1))
```





##Normal Approximation to the Binomial Distribution


```r
n <- 20
k <- 0:n
paths <- choose(n,k)
plot(k, paths/sum(paths))
p <- 0.5
lines(k, dbinom(k, n, p), col="blue", lty=3, lwd=2)
x <- seq(0, n, length=100)
lines(x, dnorm(x, mean=n * p, sd=sqrt(n * p * (1-p))), col="red")
```

![](FixALecture_files/figure-html/binomial_normal_curve-1.png) 





##DeMoivre - Laplace Theorem

[de Moivre, A. _The Doctrine of Chances_, 1738](http://en.wikipedia.org/wiki/De_Moivre%E2%80%93Laplace_theorem)

$${n \choose k}p^kq^{n-k}\simeq \frac{1}{\sqrt{2\pi npq}}e^{-\frac{(k-np)^2}{2npq}}$$

For the biomial distribution:

$$\mu = np$$
$$\sigma = \sqrt{npq}$$ where q = (1-p)





##Bell Curves


```r
x <- seq(-3, 3, length=100)
plot(x, exp(-(x^2)/2))
```

![](FixALecture_files/figure-html/bell_curve-1.png) 





##The Bell Curve as a Probability Distribution Function


```r
integrate(function(x) exp(-(x^2)/2), -Inf, Inf)
```

```
## 2.506628 with absolute error < 0.00023
```

```r
sqrt(2 * pi)
```

```
## [1] 2.506628
```





##The Normal Distribution


```r
x <- rnorm(1e5)
hist(x)
```






```r
dfrm <- data.frame( 
  x=c(x, x * 3 + 10), 
  distro=rep(c("Normal", "m10sd3"), each=length(x))
)
require("ggplot2")
```

```
## Loading required package: ggplot2
```

```r
g <- ggplot(dfrm, aes(x=x, col=distro)) + geom_density()
g
```

![](FixALecture_files/figure-html/adjusting_mean_and_sd-1.png) 





## Poisson Distribtuion

The Poisson distribution is the probability distribution of independent event occurrences in an interval. Its important property is that the sum of independent Poisson random variables is a futher Poisson random varialbe - with mean equal to the sum of the individual means. It is the typical distribution of a count. 

A discrete random variable X  is said to have a Poisson distribution with parameter &lambda;>0, if, for k = 0, 1, 2,..., the probability mass function of X  is given by:               

$$ f(k; \lambda)= \Pr(X{=}k)= \frac{\lambda^k e^{-\lambda}}{k!} $$


For the Poisson distribution with parameter lambda, probabilities and cumulative probablities are given in R by:

* Pr(X=k) = dpois(k, lambda)
* Pr(X<=K) = ppois(k, lambda)


```r
pk = dpois(0:20, 0.5)
pk
```

```
##  [1] 6.065307e-01 3.032653e-01 7.581633e-02 1.263606e-02 1.579507e-03
##  [6] 1.579507e-04 1.316256e-05 9.401827e-07 5.876142e-08 3.264523e-09
## [11] 1.632262e-10 7.419371e-12 3.091405e-13 1.189002e-14 4.246435e-16
## [16] 1.415478e-17 4.423370e-19 1.300991e-20 3.613864e-22 9.510169e-24
## [21] 2.377542e-25
```

```r
barplot(pk, names=0:20, xlab="x", ylab="probability function")
```

![](FixALecture_files/figure-html/unnamed-chunk-5-1.png) 

Note that while the Poisson distribution assigns strictly positive probability to `ANY` positive integer, in the case where the mean &lambda;= 0.5, nearly all the probability is concentrated on the integers 0 ~ 6. 

To verity this by looking at the successive values of the cumulative distribution function with:


```r
ppois(0:20, 0.5)
```

```
##  [1] 0.6065307 0.9097960 0.9856123 0.9982484 0.9998279 0.9999858 0.9999990
##  [8] 0.9999999 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
```


i.e. if X is a Poisson random variable with the mean of 0.5, Pr(X<=6) = 0.9999. We can also verify that the mean of this distribution is indeed 0.5 with 


```r
sum((0:20)*pk)
```

```
## [1] 0.5
```




###How to use poisson in a practical question?


If there are twelve cars crossing a bridge per minute on average, what is the probability of having twenty or more cars crossing the bridge in a particular minute?





###Solution:


The probability of having sixteen or less cars crossing the bridge in a particular minute is given by the function ppois.



```r
ppois(19, lambda=12)   # lower tail 
```

```
## [1] 0.9787202
```



Hence the probability of having twenty or more cars crossing the bridge in a minute is in the upper tail of the probability density function.



```r
ppois(19, lambda=12, lower=FALSE)   # upper tail 
```

```
## [1] 0.02127977
```






##Review Questions

Name this distribution:

![](FixALecture_files/figure-html/uniform_distribution-1.png) 




##Review Questions

Name this distribution:

![](FixALecture_files/figure-html/normal_distribution-1.png) 




##Review Questions

Which distribution might this be?:

![](FixALecture_files/figure-html/poisson_distribution-1.png) 

