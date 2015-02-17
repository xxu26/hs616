# GaussJordon_Practice
Xiangyi Xu (Lydia)  
Monday, February 16, 2015  

#### Matrix Operation--Practice Matrix Indexing, Identity Matrix


```r
random_invertible_matrix <- function(n, seed=NA){
if (!is.na(seed)){
if(is.character(seed)){
require(digest)
seed <- strtoi(paste0("0x",substr(digest(seed),1,4)))
}
set.seed(seed)
}

# Initialize with a singular matrix
A <- matrix( 1:n^2, nrow=n)

# Choose random integer matrixes until we find a non-singular one
while (det(A) == 0)
A <- matrix( sample(1:99, n^2, replace=T), nrow=n)
A
}
```


```r
random_invertible_matrix(3, "jan32.lydia@gmail.com")
```

```
## Loading required package: digest
```

```
##      [,1] [,2] [,3]
## [1,]   13   82   58
## [2,]   70   39   17
## [3,]   98   52   12
```

```r
n <- 3
A <- random_invertible_matrix(n, "jan32.lydia@gmail.com")
A
```

```
##      [,1] [,2] [,3]
## [1,]   13   82   58
## [2,]   70   39   17
## [3,]   98   52   12
```

```r
# Augment the square matrix by binding an identity matrix on the #right.
Aaug <- cbind(A, diag(n))

# Perform the row operations necessary to complete the Gauss-Jordan elimination.

Aaug[1,] <- Aaug[1,] / Aaug[1,1]
Aaug[2,] <- Aaug[2,] - Aaug[1,] * Aaug[2,1]
Aaug[2,] <- Aaug[2,] / Aaug[2,2]

Aaug[3,] <- Aaug[3,] - Aaug[1,] * Aaug[3,1]
Aaug[3,] <- Aaug[3,] - Aaug[2,] * Aaug[3,2]
Aaug[3,] <- Aaug[3,] / Aaug[3,3]

Aaug[2,] <- Aaug[2,] - Aaug[3,] * Aaug[2,3]
Aaug[1,] <- Aaug[1,] - Aaug[2,] * Aaug[1,2]
Aaug[1,] <- Aaug[1,] - Aaug[3,] * Aaug[1,3]

Ainv <- Aaug[,4:6]
Ainv
```

```
##              [,1]        [,2]        [,3]
## [1,] -0.008035852  0.03925205 -0.01676711
## [2,]  0.015955803 -0.10678411  0.07415778
## [3,] -0.003515685  0.14217277 -0.10108561
```

```r
# To check using built-in solver
solve(A)
```

```
##              [,1]        [,2]        [,3]
## [1,] -0.008035852  0.03925205 -0.01676711
## [2,]  0.015955803 -0.10678411  0.07415778
## [3,] -0.003515685  0.14217277 -0.10108561
```

```r
# check by multiplication to get Identity matrix
A %*% Ainv
```

```
##              [,1]          [,2]          [,3]
## [1,] 1.000000e+00  0.000000e+00  0.000000e+00
## [2,] 9.159340e-16  1.000000e+00 -8.881784e-16
## [3,] 1.075529e-15 -1.776357e-15  1.000000e+00
```

```r
round(A %*% Ainv,12) # round to 12 decimal places
```

```
##      [,1] [,2] [,3]
## [1,]    1    0    0
## [2,]    0    1    0
## [3,]    0    0    1
```
