#Week1
v <- c("hello", "world")
v2 <- 1:20

var1 <- seq(1,10,2) #var1 <- seq(1,11,2)
names(var1) <- letters[1:4] #names(var1) <- letters[1:6]
names(var1)
names(var1) <- letters[1:5] #same:names(var1) <- letters[1:6]
names(var1)
length(var1)

yORn <- c(T, F)
length(yORn)
df1 <- data.frame(var1,yORn) # Error: var1 change as commented above
df1
nchar(names(var1[1]))

is.numeric(var1)
var1[1] <- 2.3
var1[2] <- FALSE   # type?


class(df1["yORn"])
df1[1:2]
df1$yORn
df1[c("yORn","var1")]
# Create var2 and mult by var1
df1[,]
df1[1:2,2]
gender <- factor(c("male", "female"), levels=c("male", "female"))
class(gender) 
df3 <- subset(df2,v2>13)
data.frame(df2$v2>13)

v1 <- 1:6
v2 <- 11:16
c1 <-  c("hello","world")
df <- data.frame(v1,v2,c1)
print(df)
df[3,"c1"]<df[4,"c1"]
df[3,"c1"] 
df2 <- data.frame(v1,v2,c1,stringsAsFactors = F)
df2[3,"c1"]<df2[4,"c1"]




#Week2
#read file from file store
babies <-read.csv("babies.txt") #look at Environment !!
babies <-read.csv("babies.txt", sep="")
rm(babies)
babies <-read.table("babies.txt",header=TRUE)
write.csv(babies,"babies.csv")

babies2<- read.xlsx("babies.xlsx", sheetIndex =1)


#download file from web
#install.packages("downloader") 
library(downloader)
         url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
         filename <- basename(url)
         download(url, destfile=filename)
         babies <- read.table("babies.txt", header=TRUE) #alt to read.csv

#read file from web
theLink <- "http://www.jaredlander.com/data/Tomato%20First.csv"
tomatoes <- read.csv(theLink)
tomatoes <- read.table(theLink, header=TRUE, sep=",") # alt

# make sure the directory data exists
save(tomatoes,file="data/tomatoes.rdata")
rm (tomatoes)
load(file="data/tomatoes.rdata")
load("data/tomatoes.rdata")

#write to csv file
write.table(tomatoes,"data/tomatoes.csv", sep=",")

#read data from excel format:  Not recommended: 
#   preferred to save as csv and read from that file
#install.packages("xlsx")
library("xlsx") 
tomatoes2<- read.xlsx("data/tomatoes.xlsx", sheetIndex =1)

#load file from R  
data(mtcars)
hist(mtcars$mpg,20, main="miles per gallon",xlab="mpg")
 

fit <- lm(mpg ~ hp, data = mtcars)
summary(fit)
plot(mpg ~ hp, data =mtcars)
abline(reg = fit) # we see the max point has lots of leverage

# Note that the lm will pass through the avgx,avgy
abline(a = 30.1, b = -0.06823, col = 5)
max(mtcars$hp)
any(is.na(mtcars$hp))
sum(is.na(mtcars$hp))
is.na(mtcars$hp) <- which(mtcars$hp==335)
fit2 <- lm(mpg ~ hp, data = mtcars)
abline(reg = fit2, col = 2) #we have a new avgx,avgy as the pivot of our plot
mean(mtcars$hp,na.rm=TRUE) # use na.rm=TRUE
mean(mtcars$mpg,na.rm=TRUE)
summary(fit2)
plot(fit2)

# library (ggplot2)  




#Data Visualization with ggplot2 cheat sheet: 
#https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf

setwd("/Users/Pat/Documents/R/hs616/")
library(ggplot2)

#install.packages("downloader")
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

#load file from R  
data(mtcars)
hist(mtcars$mpg,10,xlim = c(0,35), main="mtcars", xlab = "mpg")
hist(mtcars$mpg,10,freq = F, main="mtcars", xlab = "mpg")

#Save the dataset to plot, then generate  density plot and histogram
g <- ggplot(mtcars, aes(x=mpg))
g + geom_histogram(aes(color=5))
g + geom_density(aes(color=2))

#boxplots
boxplot(mtcars$mpg)
boxplot(babies$gestation)

plot(mtcars)

#Side-by-side compare small to large cars
small.cars <- subset(mtcars, mtcars$cyl<6)
large.cars <- subset(mtcars, cyl>=6)
hist(mtcars$mpg,20) 
par(mfrow=c(1,2)) #plots: 1 row 2 columns
hist(small.cars$mpg,20,ylim=c(0,3),xlim=c(10,35),col = "red")
hist(large.cars$mpg,20,ylim=c(0,3),xlim=c(10,35),col="blue")
#above viz is not ideal for  comparison, why not?
par(mfrow=c(1,2))
# Do better:
hist(small.cars$mpg,20,xlim=c(10,36),ylim=c(0,4), col="blue",
     xlab="mpg",
     main="Small Cars")  
hist(large.cars$mpg,20,xlim=c(10,36),ylim=c(0,4), col="red",
     xlab="mpg",
     main="Large Cars")  


par(mfrow=c(1,1))
fit.mpg.hp <- lm (mpg~ hp, data = mtcars)
summary(fit.mpg.hp)
plot(fit.mpg.hp,1)

fit.mpg.disp <- lm (mpg~ disp, data = mtcars)
c <-coef(fit.mpg.disp)
summary(fit.mpg.disp)
plot(fit.mpg.disp,1)

fitsmall<- lm (mpg~ disp, data = small.cars)
c.small <-coef(fitsmall)
summary(fitsmall)
plot(fitsmall,1)

fitlarge<- lm (mpg~ disp, data = large.cars)
c.large <-coef(fitlarge)
summary(fitlarge)
plot(fitlarge,1)

g <- ggplot(mtcars,aes(x=disp, y=mpg))

g + geom_point(aes(color=cyl)) + geom_abline(intercept=c[1],slope=c[2])

g + geom_point(aes(color=cyl)) + geom_abline(intercept=c[1],slope=c[2]) +
  geom_abline(intercept=c.small[1],slope=c.small[2], color=5) + 
  geom_abline(intercept=c.large[1],slope=c.large[2], color=3)





library(ggplot2)


setwd("C:\\Users\\lydia\\Documents\\R\\HS616\\Week3")
#install.packages("downloader")
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)


cnames <- names(babies)
# for loop


# The 3 components of a function are its body, arguments, and environment.
# aes doesn't know your function's environment and it only looks within global environment
# So the variable dat declared within the function is not visible to ggplot2's aes function
# unless you pass it explicitly as:
showplot1<-function(indata, inx, iny) {
  p <- ggplot(indata, aes(x=indata[,inx], y=indata[,iny]), environment = environment())
  p <- p + geom_point(size=4, alpha = 0.5)
  print(p) # comment not part of body
}
#http://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r






#If the environment isn't displayed, when you print a function
#  it means that the function was created in the global environment.
showplot1
environment(showplot1)
body(showplot1)
showplot1(babies,"bwt", "gestation")
showplot1(babies2,"bwt", "gestation")
# Caution: the assignment forms of body(), formals(), and environment() can be used to modify functions.








#Alternative with ggplot2 and facet_wrap 
library(ggplot2)
library(reshape2)
b <- melt(babies[,])
head(b)
tail(b)
ggplot(b,aes(x = value)) + 
  # comparing diverse rather than similar data, so scales should range freely to display data best
  facet_wrap(~variable,scales = "free") +  # also try scales = "free"
  geom_histogram(color = 2,fill=5)


g <- ggplot(babies, aes(x=1, fill=group))
g + geom_histogram(position="dodge", binwidth=1)




library(ggplot2)

#load file from R  
data(mtcars)

#add a variable mtcars$engine.size to the data table
mtcars$engine.size <- ifelse(mtcars$cyl>=6, "large", "small")

#boxplots
# A basic box with the fill by engine.size
ggplot(mtcars, aes(x=engine.size, y=mpg, fill=engine.size)) + geom_boxplot()

# The above adds a redundant legend. Here the legend is removed:
ggplot(mtcars, aes(x=engine.size, y=mpg, fill=engine.size)) + geom_boxplot() +
  guides(fill=FALSE)

# With flipped axes
ggplot(mtcars, aes(x=engine.size, y=mpg, fill=engine.size)) + geom_boxplot() + 
  guides(fill=FALSE) + coord_flip()

cnames <- names(mtcars)
# for loop
for (c in cnames){
  print (c);
  hist(mtcars[[c]]) #we can't take the hist of a dataframe, so we need to retrieve column as vector
  #q<-ggplot(mtcars, aes_string(c)) + geom_histogram(fill=5) # ****
  #q<-ggplot(NULL, aes(x=mtcars[[c]])) + geom_histogram(color=2, fill=5)
  #q<- qplot(mtcars[[c]], geom="histogram")  # qplot can plot histogram of vector
  print (q) # draw the plot each iteration
}






#From: Advanced R by Hadley Wickham
# http://adv-r.had.co.nz/Functions.html#lexical-scoping
# The advanced topic linked above is for interest only, not required 

# What does the following code return?
# The last expression evaluated in a function becomes the
#  return value, the result of invoking the function
x <- 10
f1 <- function(x) {
  function() {
    x + 10
  }
}

#These functions were created with their own environment, rather than the global environment.
#So the environment is displayed when you print a function
f1(5) # return value is a function that contains x in the environment
f1(5)() # function call
f1(1) # return value is a function that contains x in the environment

# pure functions have no side effects: they don't affect the state of the world
# in any way apart from the value they return.
# R protects you from one type of side effect:
#  most R objects have copy-on-modify semantics. So 
#  modifying a function argument does not change the original value
# But note that there are two important exceptions to the
#  copy-on-modify rule: environments and reference classes.
#  These can be modified in place, so extra care is needed when working with them.
f1 <- function(x) {
  x$a <- 2      
  x
}
x <- list(a = 1)
x
f1(x) #function returns last expression that was evaluated
x # was global x modified?

#Write a function that adds an element named b with value of 2 to x
f2 <- function(x) {
  x$b <- 2  #adds an element **** remove line ************
  x
}
x <- list(a = 1)
x
f2(x) #function returns last expression that was evaluated
x # was global x modified?

f <- function() {
  x <- 1
  y <- 2
  c(x, y)
}
f()
rm(f)

x <- 2
g <- function() {
  y <- 1
  c(x, y)
}
g()
rm(x, g)

# If a function is defined inside another function: look inside the current function,
#  then where that function was defined, and so on, all the way up to 
#  the global environment, and then on to other loaded packages. 
# Run the following code in your head, then confirm the output by running the R code.
x <- 1
h <- function() {
  y <- 2
  i <- function() {
    z <- 3
    c(x, y, z)
  }
  i()
}
h()
rm(x, h)

# The same rules apply to closures, functions created by other functions. 
# Closures will be described in more detail in functional programming:
#  http://adv-r.had.co.nz/Functional-programming.html#functional-programming
# Here we'll just look at how they interact with scoping. 
# The following function, j(), returns a function. 
# What class will this function return when we call it?
j <- function(x) {
  y <- 2
  function() {
    c(x, y)
  }
}
k <- j(1)  # k contains environment for the returned function
k()    # so x and y are known when the function is called
rm(j, k)
6
# every time a function is called, 
# a new environment is created to host execution.






# From Advanced R by Hadley Wickham:  Memory
# http://adv-r.had.co.nz/memory.html

install.packages("pryr")
library(pryr)

# While determining that copies are being made is not hard, 
# preventing such behaviour is. If you find yourself resorting to 
# exotic tricks to avoid copies, it may be time 
# to rewrite your function in C++, as described in Rcpp.

#For loops in R have a reputation for being slow. 
#Often that slowness is because you're modifying a copy instead of modifying in place. 
#Consider the following code. 
# It subtracts the median from each column of a large data frame:
x <- data.frame(matrix(runif(100 * 1e4), ncol = 100))
medians <- vapply(x, median, numeric(1))

for(i in seq_along(medians)) {
  x[, i] <- x[, i] - medians[i]
}

#You may be surprised to realise that every iteration of the loop 
#copies the data frame. We can see that more clearly by 
# using tracemem or address() and refs() for a small sample of the loop:

untracemem(x)
x <- data.frame(matrix(runif(100 * 1e4), ncol = 100))
tracemem(x)
medians <- vapply(x, median, numeric(1)) # median of each column, 100 columna
for (i in 1:5) {
  x[, i] <- x[, i] - medians[i] # ith column subtract the ith median from each element
  #print(c(address(x), refs(x)))
}

#For each iteration, x is moved to a new location so refs(x) is always 2. 
#This occurs because [<-.data.frame is not a primitive function, 
# so it always increments the refs. We can make the function substantially
# more efficient by using a list instead of a data frame. 
#Modifying a list uses primitive functions, so the refs are 
# not incremented and all modifications occur in place:
y <- as.list(x)
for(i in 1:5) {
  y[[i]] <- y[[i]] - medians[i]
  print(c(address(y), refs(y)))
}

#This behaviour was substantially more problematic prior to R 3.1.0,
# because every copy of the data frame was a deep copy. 
#This made the motivating example take around 5 s, compared to 0.01 s today.







# Adapted with minor modifications from Bob Horton's script beer.R 
# https://www.youtube.com/watch?v=5Dnw46eC-0o

# Beer consumption increases human attractiveness to malaria mosquitoes
# http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0009546
# http://dx.plos.org/10.1371/journal.pone.0009546
# http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/20209056/

# LefÃ¨vre T, Gouagna L-C, DabirÃ© KR, et al. Beer Consumption Increases Human Attractiveness to Malaria Mosquitoes. Tregenza T, ed. PLoS ONE 2010;5(3):e9546. doi:10.1371/journal.pone.0009546.


beer <- c( 27, 20, 21, 26, 27, 31, 24, 21, 20, 19, 23, 24, 28, 19, 24, 29, 18, 20, 17, 31, 20, 25, 28, 21, 27 )

water <- c( 21, 22, 15, 12, 21, 16, 19, 15, 22, 24, 19, 23, 13, 22, 20, 24, 18, 20 )

library("ggplot2")

results <- data.frame(
	mosquitoes = c(beer, water),
	group = c(rep("B", length(beer)), rep("W", length(water)))
)

g <- ggplot(results, aes(x=group, y=mosquitoes, fill=group))
g + geom_boxplot() 

par(mfrow=c(1,2))
boxplot(subset(results$mosquitoes,results$group=="B"), ylim=c(12,32)) 
boxplot(subset(results$mosquitoes,results$group=="W"), ylim=c(12,32)) 

g <- ggplot(results, aes(x=group, y=mosquitoes, fill=group))
g + geom_boxplot() + geom_jitter(position=position_jitter(w=0.1))

# Interleaved histogram
g <- ggplot(results, aes(x=mosquitoes, fill=group))
g + geom_histogram(position="dodge", binwidth=1)

# Overlaid histogram
g <- ggplot(results, aes(x=mosquitoes, fill=group))
g + geom_histogram(position="identity", binwidth=1, alpha=.5) # alpha sets opacity

# Density plot
g <- ggplot(results, aes(x=mosquitoes, col=group))
g + geom_density()

# Density plot with semi-transparent fill
g <- ggplot(results, aes(x=mosquitoes, fill=group))
g + geom_density(alpha=.2)

            



#install.packages("dplyr")
library(dplyr)
library(ggplot2)

### available on github
#The mammals sleep data set contains the 
# sleeptimes and weights for a set of mammals 
# column name	Description
# name	common name
# genus	taxonomic rank
# vore	carnivore, omnivore or herbivore?
# order	taxonomic rank
# conservation	the conservation status of the mammal
# sleep_total	total amount of sleep, in hours
# sleep_rem	rem sleep, in hours
# sleep_cycle	length of sleep cycle, in hours
# awake	amount of time spent awake, in hours
# brainwt	brain weight in kilograms
# bodywt	body weight in kilograms

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- "msleep_ggplot2.csv"
if (!file.exists(filename)) download(url,filename)
msleep <- read.csv("msleep_ggplot2.csv")
head(msleep)

#apply family vectorized functions


# dplyr tutorial adapted from:
# http://genomicsclass.github.io/book/pages/dplyr_tutorial.html
# The dplyr package uses efficient data storage backends, is fast
# dplyr functions take a data frame and return a data frame:
# select (return new df that is subset of cols: like select in SQL)
# filter (return  new df that is subset of rows that meet logical condition: like where clause)
# arrange: (return new df that has reordered rows: like  orderby, can use desc)
# rename:  (return new df that with variable(s) renamed) 
# mutate: (return  new df with transformed or new variables)  

# More infofrom: Introduction to dplyr
# https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html


#select columns: select name, sleep_total from msleep
# (return new df that is subset of cols: like select in SQL)

#select columns 3:7 from msleep


#select all columns *except* genus


#select all columns *except* genus, vore, order
###select(msleep, -(genus, vore, order))  ## wrong

#select all columns *except* columns 3:7 
select(msleep, -vore:sleep_rem) ## wrong

#select all columns that start with (starts_with) the character string "sl"

# distinct:  select   the distinct orders

#filter the rows for mammals that sleep between 12 and 15 hours inclusive
# Here are 2 base R ways to do this, replace with filter
msleep[msleep$sleep_total >= 12 & msleep$sleep_total<=15, ]
subset(msleep, sleep_total >= 12 & msleep$sleep_total<=15)


#filter() works similarly to subset()
# you can give it any number of filtering conditions, 
# which are joined together with & (not &&). 
# You can also use other boolean operators

#filter the rows for mammals that sleep <12 hours or are awake>12 hours


# use pipe to filter as above then select columns sleep_total and awake, display head


#investigate relationship between sleep_total and awake

all(msleep$sleep_total+msleep$awake ==24)
which(msleep$sleep_total+msleep$awake !=24)

#Filter the rows for mammals that are in the Perissodactyla and Primates taxonomic orders


#selectcolumns  name and order, arrange (ie: re-order) rows by taxonomic order, display head


#select columns name, order, sleep_total from msleep, 
# arrange rows by taxonomic order, then by sleep_total, display head


# # Similar to SQL:
# select name, order, sleep_total
#   from msleep
#   order by order,sleep_total
#   limit 6

#Same as above except filter the rows for mammals that 
# sleep for 16+ hours instead of displaying head,
#  arrange rows by taxonomic order, then by descending sleep_total


#Create a new column named rem_proportion which is 
# the ratio of rem sleep to total amount of sleep.

#SQL: alter table <tablename> add column <columnname> type

# NB: mutate returns a  new df with transformed or new variables 
#  mutate does not mutate the original data frame in place


# summarize values


#group_by: takes an existing tbl and converts it into a grouped tbl
# where operations are performed "by group"


#group by order, summarize mean, min, max, std and totals of sleep_total by group


#preliminaries:rnorm; cbind and rbind
# rnorm randomly generates numbers from the normal distribution
rnums <- rnorm(10)
mean(rnums)
sd(rnums)
# What happens when we increase sample size? Rerun with 10000 rows:


rnums1 <- rnorm(10000)
hist(rnums1, ylim=c(0,3000))

hist(rnums2, xlim=c(50,150),ylim=c(0,3000))

hist(rnums3, xlim=c(50,150),ylim=c(0,3000))

rnums_cols<-
        head(rnums_cols)
rnums_rows<-
        head(rnums_rows)

# apply family of functions: 
# apply operates on every row or every column of an array (including matrix, dataframe)
#   will coerce all elements in the row of column to be of the same type
# lapply operates on a list, vector or data frame and returns a list
# sapply operates on a list,vector or data frame and attempts to return a vector
#  sapply is a wapper function for lapply, to simplify the output when possible

# review matrix positional args: data, nrow, ncol
matrix(1:12, 4, 2)# positional args: row before column
m1 <- matrix(1:12, 4) # if only 2 positional args: data, nrows
m1
m2 <- matrix(1:12, ncol=4)# keyword arg
m2

# apply mean to rows 


#what happpens with lapply or sapply?
lapply(m2,mean) #operates on every element of list or vector, returns a list
sapply(m2,mean) #operates on every element of list or vector, returns a vector
# show  matrix can be coerced to list or vector: above is equivalent to:


# matrix with randomly generated numbers:
A <-rnorm(10)
B <-rnorm(10,5)
C <-rnorm(10,5,3)
# create matrix having columns of A,B,C  *

#apply the mean to columns  *

v <- 1:10 # vector: apply functions


a<-array(1:9,dim=c(3,3)) #array


lst1 <- list("hello", c(T, F), 1:10) #list


lst2 <- list(A = 1:7, B = matrix(1:12, 4), D = 12)#list with named elements
#lapply operates on each item of a list and returns a list


# apply with a user-defined function:
f <- function(x){
        x^2 + 6*x - 1
}


# data frame
data (mtcars)
min(mtcars$mpg)
which(mtcars$mpg == min(mtcars$mpg))
# cell 1  to access the first occurrence of the min value:
which(mtcars$mpg == min(mtcars$mpg))[1]

# do above, but apply to each column of the data frame:
apply(mtcars, 2,min) 
apply (mtcars, 2, function(x) which(x == min(x))[1])

# break apart the above function:
f2<-function(x){
        which(x == min(x)) [1] 
} 
apply (mtcars, 2, f2) # applies f2 to each column of mtcars

library(dplyr)
mtcars2<- mtcars # copy mtcars
#explore
mtcars2 %>% filter(cyl>=6, disp<150 )
mtcars2 %>% filter(cyl<6, disp>145 )
mtcars2$engine_size <- ifelse(mtcars$disp<145, "small", "large")
# make engine_size ordered with small < large


#compare apply functions


#preliminaries:rnorm; cbind and rbind
# rnorm randomly generates numbers from the normal distribution
rnums <- rnorm(10)
hist(rnums)
mean(rnums)
sd(rnums)
# What happens when we increase sample size? Rerun with 10000 rows:


rnums1 <- rnorm(10000)
hist(rnums1, ylim=c(0,3000))
mean(rnums1)
sd(rnums1)
rnums2 <- rnorm(10000,100)
hist(rnums2, xlim=c(50,150),ylim=c(0,3000))
rnums3 <- rnorm(10000,100,10)
hist(rnums3, xlim=c(50,150),ylim=c(0,3000))

rnums_cols<-cbind(rnums1,rnums2, rnums3)
head(rnums_cols)
rnums_rows<-rbind(rnums1,rnums2, rnums3)
head(rnums_rows)

# apply family of functions: 
# apply operates on every row or every column of an array (including matrix, dataframe)
#   will coerce all elements in the row of column to be of the same type
# lapply operates on a list, vector or data frame and returns a list
# sapply operates on a list,vector or data frame and attempts to return a vector
#  sapply is a wapper function for lapply, to simplify the output when possible

# review matrix positional args: data, nrow, ncol
matrix(1:12, 4, 2)# positional args: row before column
m1 <- matrix(1:12, 4) # if only 2 positional args: data, nrows
m1
m2 <- matrix(1:12, ncol=4)# keyword arg
m2

# apply mean to columns 
apply(m1,2,mean)# operates row-wise or column-wise on a matrix
apply(m2,2,mean)

#what happpens with lapply or sapply?
lapply(m2,mean) #operates on every element of list or vector, returns a list
sapply(m2,mean) #operates on every element of list or vector, returns a vector
# show  matrix can be coerced to list or vector: above is equivalent to:
sapply(as.vector(m2),mean) 

# matrix with randomly generated numbers:
A <-rnorm(10) 
B <-rnorm(10,5) 
C <-rnorm(10,5,3) 
# create matrix having columns of A,B,C  *
m<-cbind(A,B,C)	
m <- matrix(data=cbind(A,B,C), nrow=10, ncol=3)#fully specified form

#apply the mean to columns  
apply(m,2,mean)

v <- 1:10  # vector:  # not done in class #
apply(v,1,mean) # Error: vectors have neither columns nor rows to apply against
sapply(v,mean) #operates on a list or vector, returns a vector
lapply(v,mean) #operates on a list or vector, returns a list

a<-array(1:9,dim=c(3,3)) #array	  # not done in class #
apply(a,2,mean)
sapply(a,mean) 
lapply(a,mean)

lst1 <- list("hello", c(T, F), 1:10) #list 
apply(lst1,2,mean) # Error: lists have neither columns nor rows to apply against
lapply(lst1,min)
sapply(lst1,min)

lst2 <- list(A = 1:7, B = matrix(1:12, 4), D = 12)#list with named elements
#lapply operates on each item of a list and returns a list
lapply(lst2, min)     # not done in class #
lapply(lst2, max)     # not done in class #
lapply(lst2, mean)    # not done in class #
lapply(lst2, median)  # not done in class #
sapply(lst2, median)  # not done in class #

# apply with a user-defined function:
f <- function(x){
        x^2 + 6*x - 1
}
sapply(1:10,f)
lapply(lst2,f)	   # not done in class #
sapply(lst2,f)	   # not done in class #

# data frame	   # all below not done in class #
data (mtcars)
min(mtcars$mpg)
which(mtcars$mpg == min(mtcars$mpg))
# cell 1  to access the first occurrence of the min value:
which(mtcars$mpg == min(mtcars$mpg))[1]

# do above, but apply to each column of the data frame:
apply(mtcars, 2,min) 
apply (mtcars, 2, function(x) which(x == min(x))[1])

# break apart the above function:
f2<-function(x){
        which(x == min(x)) [1] 
} 
apply (mtcars, 2, f2) # applies f2 to each column of mtcars

library(dplyr)
mtcars2<- mtcars # copy mtcars
#explore
mtcars2 %>% filter(cyl>=6, disp<150 )
mtcars2 %>% filter(cyl<6, disp>145 )
mtcars2$engine_size <- ifelse(mtcars$disp<145, "small", "large")
# make engine_size ordered with small < large
mtcars2$engine_size <- factor(mtcars2$engine_size, ordered=T,levels=c("small", "large"))

#compare apply functions
apply(mtcars2, 2,min)
lapply(mtcars2,min) # returns what type?
sapply(mtcars2, min) # returns what type?





  #install.packages("dplyr")
library(dplyr)
library(ggplot2)

# Variable Descriptions:
# bwt Birth weight in ounces 
# gestation Length of pregnancy in days 
# parity 0 = first born, 1 = otherwise 
# age mother's age in years 
# height mother's height in inches 
# weight Mother's pre-pregnancy weight in pounds 
# smoke Smoking status of mother: 0 = not now, 1 = yes now

#install.packages("downloader")
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

# Substitute NA for the codes used to indicate missing data
babies$gestation <- ifelse(babies$gestation==999, NA,babies$gestation)
babies$age <- ifelse(babies$age == 99, NA,babies$age)
babies$height <- ifelse(babies$height == 99, NA,babies$height)
babies$weight <- ifelse(babies$weight == 999, NA,babies$weight)
babies$smoke <- ifelse(babies$smoke == 9, NA,babies$smoke)
max(babies$age, na.rm = TRUE) # max age is now reasonable

babies2 = babies
# Leave out outliers  x<220 and x>330
babies2$gestation <- ifelse(babies2$gestation<220 | babies2$gestation> 330, NA,babies2$gestation) # alt: now 20 NA

babies2$age.level <- factor( rep("early.20s",nrow(babies2)), ordered=T,levels =c("teens", "early.20s", "late.20s", "early.30s", "late.30s", "forty.plus"))
babies2$age.level[babies2$age<20] <- "teens"
babies2$age.level[babies2$age>=25 & babies2$age<30] <- "late.20s"
babies2$age.level[babies2$age>=30 & babies2$age<35] <- "early.30s"
babies2$age.level[babies2$age>=35 & babies2$age<40] <- "late.30s"
babies2$age.level[babies2$age>=40] <- "forty.plus"
babies2$age.level#optional, shows ordering


#1 save babies2 as babies2.Rdata
#   save babies2 as babies2.csv


#2
#select returns new df with subset of columns: like select in SQL
# 2a select the 5 columns of babies2 from age to age.level inclusive
# 2b select columns parity, smoke, age.level from  babies2 
# 2c select all columns except those 5 you selected in 2a
# 2d select all columns except those 3 you selected in 2b


#3
#filter returns  new df that is subset of rows that meet logical condition: like where in SQL
#filter babies2 to include all rows with age between 28 and 32 inclusive



#4
#arrange: return new df with rows reordered:: like  orderby in SQL, can use desc
# Create babies3 by reordering babies2 by age level, then by descending gestation



#5
#rename: returns new df with variable(s) renamed *****************
# 5a: Rename parity to is.firstborn, rename smoke to is.smoker
# 5b: Use a non dplyr way: rename bwt to birthweight



#6
#mutate: return new df with transformed or new variables 
# Mutate babies3 so that the mean height is subtracted from each value for height
# Look at your babies3 and make sure you have modifed it but haven't rendered the data unusable



#7 new variable:
#Mutate babies2to add a new factor variable named birth.size
# such that bwt>120 is large and  bwt <= 120 is small



#8
# group_by:takes an existing tbl and converts it into a grouped tbl
#  where operations are performed "by group". 
#  Returns the converted grouped tbl
# 8a: Create group babies2 by age.level and store in new variableresultset 
# 8B: Summarize resultset: provide for each group the med.weight and 
#  med.gest, which are output names for median weight and median gestation.
# Look at your output and make sure you do not have NAs.
# 8c: Now use the same command to summarize babies2 
# 8d: write a brief sentence on the difference in summarizing babies2 and resultset



#9 pipeline operator %>% 
#Use a piped expression to selct the following information:
# select birthweight , gestation 
#  from babies3 
#  where weight > median(weight) and is.smoker == T and age.level==early.20s
#   order by gestation



#10 boxplots
#Create 3 boxplot for the variable gestation as follows:
# 10a: create a plot containing side-by-side boxplots for each of the 6 age.levels
# 10b: create a plot containing side-by-side boxplots for the 2 smoke levels 
# 10c: create a plot containing side-by-side boxplots for the 2 parity levels
ggplot(babies2, aes(x=age.level,y=gestation, fill=age.level)) + geom_boxplot()
ggplot(babies2, aes(x=as.factor(smoke),y=gestation, fill=as.factor(smoke))) + geom_boxplot()
ggplot(babies2, aes(x=as.factor(parity),y=gestation, fill=as.factor(parity))) + geom_boxplot()




# debugging: reproducibility set seed for pseudo-random numbers
set.seed=1
# if commented out, seed is taken from the clock, different each time

r1<- rnorm(10000)
hist (r1, 10, ylim=c(0,4000), xlim=c(-10,30))  #####note the spread of the function#####

# pnorm is the cumulative distribution function
# pnorm gives you the probablity that a random variable (RV)
#  x is less than the first argument 
pnorm(q=0)
pnorm(1) #probability that random x < 1
pnorm(2) 
pnorm(3) 
qnorm(.5)  # expect 50% of randomly generated values to be below this value
qnorm(0.8413447)
qnorm(0.9772499)
qnorm(0.9986501)

#How to we get the prob that RV x is 1 standard deviation from the mean?

#How to we get the prob that RV x is 2 standard deviation from the mean?

#How to we get the prob that RV x is 3 standard deviation from the mean?


#We can pass in mean and/or sd rather than use default of 0 and 1 
pnorm(0,0,1)
pnorm(q=0,mean=0,sd=1)
pnorm(0,10,4)
hist(rnorm(10000,10,4), 40, ylim=c(0,4000), xlim=c(-10,30))  #####note the spread of the function#####

# dnorm returns the height of the normal curve at a RV.
dnorm(0)
dnorm(0,10,4)


# Create a data frame of the size you need, fill in as you go
data.frame(matrix(NA, nrow = 2, ncol = 3))

# melt from wide to long format (ex: for faceted plotting)
# id columns 
#default:  all columns that are not id variables will be measured variables
# measured variables are stacked into a single column of data
library (reshape2)
b <- melt(babies3)
head(b)
tail(b)
ggplot(b,aes(x = value)) + 
  # comparing diverse rather than similar data, so scales should range freely to display data best
  facet_wrap(~variable,scales = "free") +  # also try scales = "free"
  geom_histogram(color = 2,fill=5)

b <- melt(babies2) 
b <- melt(babies2, id.vars = "birth.size") 


# Histograms and density plots
g <- ggplot(babies3, aes(x=birthweight, fill=age.level))

# Overlaid histograms with position="identity"
g + geom_histogram(binwidth=10, alpha=.5, position="identity")

# Interleaved histograms with position="dodge"
g + geom_histogram(binwidth=10, position="dodge")

# Density is used compare distributions with very different counts
g + geom_density()

# Density plots with semi-transparent fill
g + geom_density(alpha=.4)


# look over assignment
#if time
runif()

lapply(1:4, runif)
lapply(1:4, runif, min=-1,max=1)







data (mtcars) 

# creates a vector of labels for 1000 elements 
rep("label1",10)

# combines two vectors made up of repetitions into one
c((rep("label1",10)),(rep("label2",5)))

rnorm1 <- rnorm(1000, mean=100, sd=15)
hist(rnorm1, probability=TRUE)
#we want to draw a line to approximate the density function
# generate 100 xxvalues from min(x) to max(x)
xx <- seq(min(rnorm1), max(rnorm1), length=100)
yy <- dnorm(xx, mean=100, sd=15)
# draws a line on the most recent plot (similarly to abline)
lines(xx, yy)# creates a line from the (xx,yy) points


# Plots for the standard normal distribution
# source: ___________________
set.seed(3000)
xseq<-seq(-4,4,.01)
densities<-dnorm(xseq, 0,1)
cumulative<-pnorm(xseq, 0, 1) # always increasing
randomdeviates<-rnorm(1000,0,1)

par(mfrow=c(1,3), mar=c(3,4,4,2))
plot(xseq, densities, col="darkgreen",xlab="", ylab="Density", type="l",lwd=2, cex=2, main="PDF of Standard Normal", cex.axis=.8)
plot(xseq, cumulative, col="darkorange", xlab="", ylab="Cumulative Probability",type="l",lwd=2, cex=2, main="CDF of Standard Normal", cex.axis=.8)
hist(randomdeviates, main="Random draws from Std Normal", cex.axis=.8, xlim=c(-4,4))



binom1 <- rbinom(1000, 100, .2)
# Since this is a discrete function we have a probability function (PF) 
#     (also known as probability mass function)
#  rather than a probablility density function (PDF)
# dbinom looks up P(X = 20) when X is drawn from 
#    binomial(100, 0.2) distribution.
dbinom(20, 100, .2) # probability that number of successes=20
      # over 100 trials, each trial having probability of success = .2
hist(binom1) # does this surprise you?
# Will a different histogram be more informative about dbinom?

dbinom(27, size=100, prob=0.25)




# Shapiro-wilk test for normality 
# non-missing values must be between 3 and 5000
# p-value gives us probablilty that the distribution is normal
# tests the NULL hypothesis that the samples came from a Normal distribution
#  This means that if your p-value <= 0.05, then reject the NULL hypothesis 
#   and conclude your sample is likely to be from a non-normal distribution
r1 <- rnorm(5000,10,4)
shapiro.test(r1)
#Normal Q-Q plot: plot against theoretical quantiles
#  qqline draw straight line, color=red
qqnorm(r1); qqline(r1, col = 2,lwd=2)
par(mfrow=c(1,1))
qqnorm(rnorm(1000,10,4)); qqline(rnorm(1000,10,4), col = 2,lwd=2)
qqnorm(runif(1000)); qqline(runif(1000), col = 2,lwd=2)
qqnorm(rpois(1000,1)); qqline(rpois(1000,1), col = 2,lwd=2) 
qqnorm(rbinom(1000,30, .5)); qqline(rbinom(1000,30, .5), col = 2,lwd=2) 


#A poisson distribution is for count data, 
# lambda is both mean and variance.
pois1<- rnorm(100,1)
pois2<- rnorm(100,10)
par(mfrow=c(1,2))
hist(pois1)
hist(pois2)






#from Pengvideolecture: Lapply (v2)-HD 720p.mov: 
# extract first column of each matrix in lst
lst = list(a=matrix(1:4,2,2), b=matrix(1:6,3,2))
lapply(lst,function(elem){elem[,1]}) #using anonymous function

#generate a factor variable with 3 levels, 10 values each
gl(3, 10)

# split (returns a list) followed by lapply to operate on the lists
split_babies <- split(babies,babies3$age.level)
#colMeans, rowMeans, colSums, rowSums operate on columns only,
#  so can be applied to lists of arrays, matrices or data frames
sapply(split_babies, colMeans, na.rm=T)
sapply(split_babies, function(x) colMeans(x[,c("bwt","gestation", "age")], na.rm=T))

split_babies <- split(babies,list(babies3$age.level,babies3$is.smoker))
sapply(split_babies, function(x) colMeans(x[,c("bwt","gestation", "age")], na.rm=T))

#tapply is typically used on a vector
# results are by group
tapply(babies$bwt,babies3$age.level, mean) 
tapply(babies,babies3$age.level, mean) #Error when used with data frame

#by is a wrapper for tapply for use with dataframe
# results are by group, operates on columns
by(babies[, c(1:2,3:6)], babies3$age.level, colMeans, na.rm=T)
by(babies[, 1:2], babies3$age.level, summary)

# Extract coefficients of linear model by group
#  No need to first subset
tmp <- with(babies3,
            by(babies3, age.level,
               function(x) lm(birthweight ~ gestation, data = x)))
coef.vec <- sapply(tmp, coef)
sapply(tmp, abline(intercept=coef[1],slope=coef[2]))
coef.vec

ggplot(babies3, aes(x=gestation,y=birthweight)) + geom_point(aes(color=age.level)) +
         geom_abline(intercept=coef.vec[1,1],slope=coef.vec[2,1], color=1) +
         geom_abline(intercept=coef.vec[1,2],slope=coef.vec[2,2], color=2) + 
         geom_abline(intercept=coef.vec[1,3],slope=coef.vec[2,3], color=3) +
         geom_abline(intercept=coef.vec[1,4],slope=coef.vec[2,4], color=4) +
         geom_abline(intercept=coef.vec[1,5],slope=coef.vec[2,5], color=5) + 
         geom_abline(intercept=coef.vec[1,6],slope=coef.vec[2,6], color=6)       
      
apply(babies3[,1:7],2, hist)
sapply(babies3[,1:7], hist)






data (mtcars) 

# c takes on each column in mtcars
for (c in (mtcars)){
  print (c)
}

#seq_len(y) takes numeric for y, creates a sequence up to y
# similar to Python: for i in range(len(lst))
# R: generate row numbers (or column numbers)
seq_len(nrow(mtcars))# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
for (i in seq_len(nrow(mtcars))){ 
  print (mtcars[i,]) # i takes on each row index of mtcars
}

# seq_along is generally used for vectors, but
#  if you need column numbers of data frames, this works
# seq_along(x) takes a vector for x, creates a sequence 
#  up to the count of elements in the vector
seq_along(1:10) # [1]  1  2  3  4  5  6  7  8  9 10
seq_along(names(mtcars)) # [1]  1  2  3  4  5  6  7  8  9 10 11
for (i in seq_along(names(mtcars))){
  print (mtcars[,i]) # i takes on each column index of mtcars
}#if you remove names seems to work on data frame, but is intended for vectors


# If you create an empty data frame using matrix, 
# this is fast, and is good for all one type of data
data.frame(matrix(NA, nrow = 2, ncol = 3))

# If your data frame has columns of different types,
#  you can go a slower route and avoid type coercions

#create empty data frame with the same columns as mtcars
my.cars <- mtcars[0,]  

# add rows to a data frame
result<- mtcars %>% filter(105<=hp,hp<=115)
my.cars<- rbind(my.cars, result)   





library(reshape2)
library(ggplot2)

#10
unif1 <- runif(1000, 0, 1)
unif2 <- runif(1000, 20, 30)

#11
#not necessary to wrap up in a function: but makes remaining questions easier
plot_side_by_side <- function(vec1, vec2){ 
  df <- data.frame(vec1,vec2)
  df_long <- melt(df)
  g <- ggplot(df_long, aes(x=value, fill=variable))
  g + geom_histogram(binwidth=1, alpha=.5, position="identity") # Overlaid
}
plot_side_by_side(unif1,unif2) 


#20 Notes
# Shapiro-Wilk test of normality
#  first arg x:	a numeric vector of data values. Missing values are allowed,
#  but the number of non-missing values must be between 3 and 5000
# Tests the null hypothesis that the sample x is from a normal distribution
# We establish an alpha cutoff, often .05
# If p-value<alpha, this suggests that observed data is inconsistent with null hypothesis
#  and we reject null hypothesis that x is normally distributed
# If p-value is high then we affirm the null hypothesis that x is normally distributed
# Hint: use apply to run shapiro.test on each column of the data frame from #19


# Extra plots: not part of assignment:

# Compare overlaid density plot (the default for density)
#  with stacked density plot for norm1, norm2:
# which gives you a better idea of the distributions?
norm1 <- rnorm(1000)
norm2 <- rnorm(1000,5,3) 
norm.dist2 <- data.frame(norm1, norm2)
norm.dist2.melt <- melt (norm.dist2) 

g = ggplot(norm.dist2.melt, aes(x = value, fill=variable))
# Display a smooth density estimate with geom_density
# default is position="identity" (ie: overlaid density)
g + geom_density( alpha=.5)+
  labs(title="Default: Overlaid Density plot: position=identity")

g + geom_density(position="stack")  +
  labs(title="Stacked Density: position=stack")



# Like data frame rand.dist, but 
#  replace pois1 with pois2, norm1 with norm3, unif1 with unif2,
#  keep binom1, take out bern1. *** NB: norm3 is repeated 10 times ***
norm3 <- rnorm(100, 20, 10) 
unif2 <- runif(1000, 20, 30)
pois2 <- rpois(1000, 30)
binom1 <- rbinom(1000, 100, .2)

rand.dist2 <- data.frame(norm3, unif2, pois2, binom1)
rand.dist2.melt <- melt (rand.dist2) 

ggplot(rand.dist2.melt, aes(x = value)) + 
  # comparing diverse rather than similar data, so scales should range freely to display data best
  facet_wrap(~variable,scales = "free") + 
  geom_histogram(color = 2,fill=5)

g = ggplot(rand.dist2.melt, aes(x = value, fill=variable))

# Display a smooth density estimate with geom_density
# default is position="identity" (ie: overlaid density)
g + geom_density( alpha=.5)+
  labs(title="Default: Overlaid Density plot: position=identity")

g + geom_density(position="stack")  +
  labs(title="Stacked Density: position=stack")


#default is stacked histogram with position="stack"
g + geom_histogram(binwidth=1) +
  labs(title="Default: Stacked Histogram: position=stack")

#interleaved histogram with position="dodge"
g + geom_histogram(binwidth=1, position="dodge")+
  labs(title="Interleaved Histogram: position=dodge")

# Overlaid histogram with position="identity"
g + geom_histogram(binwidth=1, alpha=.5, position="identity") +
  labs(title="Overlaid Histogram: position=identity")





#Week 6


norm1 <- rnorm(1000)
norm2 <- rnorm(1000,5,3) 
norm3 <- rnorm(100, mean=20, sd=10) 

df_norm <- data.frame(norm1,norm2, norm3)

# Get summary statistics for df_norm
summary(df_norm)  ######## ######## ######## ######## ########

# Use apply to get the varianceof each distibution
#  NB: that apply must have a function arguement to apply: use variance
apply(df_norm,2,var)   ######## ######## ######## ######## ########

library (reshape2)
babies2 <- read.csv("babies2.csv")
# BTW, summary automatically removes NAs from min, max, mean, etc
summary(babies2)

cor(babies2[2:8]) #default is "everything"
cor(babies2[2:8], use="everything") #NA if either column has any NA
cor(babies2[2:8], use="all.obs") # Error if a single NA in any column
# Next two lines keep only *rows* where every entry is not NA
cor(babies2[2:8], use="na.or.complete")
cor(babies2[2:8], use="complete.obs") # Error if there is no row that is complete

install.packages("GGally")
library (GGally)

# cat vs cat:  bar, box
# cat vs numeric: boxplot , hist
# numeric vs numeric: scatter, cor
babies2$smoke <- factor(babies2$smoke)
ggpairs(babies2[6:9])

# One sample t-test to test null hypothesis that true mean is
# equal to an expected value
# Reject null hypothesis if p-value < alpha
t.test(babies2$age, mu=28) # reject null hypothesis?
t.test(babies2$age, mu=15) # reject null hypothesis?
t.test(babies2$age, mu=27) # reject null hypothesis?



#Wk6 part2
# z-score (standard score)
# Can compare two scores from different *normal* distributions
# We are basically mapping the distribution to the standard normal

# standard normal distribution
norm2 <- rnorm(10000,10,4)
summary(norm2)
m <- mean(norm2)
s <- sd(norm2)
# z-score(1): val is how many SDs above the mean?
val <- 18
z <- (val - m)/s  
z

# normal distribution with mean=30 and sd=3
norm3 <- rnorm(10000, 30,3)
summary(norm3)
m <- mean(norm3)
s <- sd(norm3)
# What value in this distribution is comparable to val of norm1?
# Use formula for z z-score, solve for s
x <- z*s + m
x
# Check:
z <- (x - m)/s  # 33 is how many SDs above the mean?
z


#The p-value is the probability that the observed data could happen,
# under the condition that the null hypothesis is true
#alpha (significance level) the probability of making a type1 error
# type I error is detecting an effect that is not present
# alpha, simply stated, is the probability of making a wrong decision

setwd("/Users/Pat/Documents/R/HS_616/assign")
babies <- read.table("babies.txt", header=TRUE)

babies2 <- read.csv("babies2.csv")
babies2<-babies2[, -1]  
summary(babies2)

s <- sd(babies2$bwt)
stats<-summary(babies2$bwt)
m <- stats[4]  
m

# One sample t-test to test null hypothesis that true mean is
# equal to an expected value
# Reject null hypothesis if p-value < alpha
t.test(babies2$age, mu=28) # reject null hypothesis?
t.test(babies2$age, mu=15) # reject null hypothesis?
t.test(babies2$age, mu=27) # reject null hypothesis?
 
t.test(babies$bwt, babies2$bwt)


#Week7
set.seed(100)

m <- 1
b <- 0
x <- rnorm(20)	
#Generate line with slope m,intercept b, and random noise
y <- m*x + b + rnorm(length(x), mean=0, sd=1)

plot(x,y) # with noise
abline(b, m, lty=2, col="green") # original line

fit<- lm (y~ x)
c <-coef(fit)
summary(fit) # data fits linear model with R^2 0.6696

# Good fit to a linear model
# Not as good a fit to the line data originated from (green)
library (ggplot2)
df <- data.frame(x,y)
g <- ggplot(df,aes(x=x, y=y))
g + geom_point() + 
  geom_abline(intercept=c[1],slope=c[2]) +
  geom_abline(intercept=b, slope=m, lty=2, col="green")


#Week8

library(reshape2)
library(ggplot2)

# Written by Andy Field
logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance 
  nullDev <- LogModel$null.deviance 
  modelN <-  length(LogModel$fitted.values)
  R.l <-  1 -  dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2  ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2        ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2           ", round(R.n, 3),    "\n")
}





#Week9
setwd("/Users/Pat/Documents/R/HS_616/lecture_scripts")

#install.packages("ppcor")
library (ppcor)
library(ggplot2)
library(reshape2)


logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance 
  nullDev <- LogModel$null.deviance 
  modelN <-  length(LogModel$fitted.values)
  R.l <-  1 -  dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2  ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2        ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2           ", round(R.n, 3),    "\n")
}

# PSWQ: degree to which the player worries in general
# Anxious: new measure of the player's anxiety just before the penalty kick
# Previous: percentage of penalties scored by the player over their career
# Scored: outcome: whether the penalty kick scored

penalties <- read.table(file="penalty.txt", header = T)  # tab-delimited file with a header

# Exploratory visualizations
g <- ggplot(data=penalties)
g + geom_histogram(aes(x=PSWQ),binwidth=1, color = 5)  

g + geom_histogram(aes(x=PSWQ, fill=Scored),binwidth=1,position="dodge")  # position="identity" for overlaid
g + geom_histogram(aes(x=Anxious, fill=Scored), binwidth=1,position="dodge")
g + geom_histogram(aes(x=Previous, fill=Scored), binwidth=1,position="dodge")

g + geom_density(aes(x=PSWQ, fill=Scored), alpha=.5)   
g + geom_density(aes(x=Anxious, fill=Scored), alpha=.5)  
g + geom_density(aes(x=Previous, fill=Scored), alpha=.5)   

# Tests for normality of data:
# null-hypothesis is that population is normally distributed
# high p-value -> don't reject null, infer population is normally distributed
shapiro.test(penalties$PSWQ)
shapiro.test(penalties$PSWQ[penalties$Scored=="Scored"])
shapiro.test(penalties$PSWQ[penalties$Scored=="Missed"])
shapiro.test(penalties$Anxious[penalties$Scored=="Scored"]) # normally dist
shapiro.test(penalties$Anxious[penalties$Scored=="Missed"]) # normally dist
shapiro.test(penalties$Previous[penalties$Scored=="Scored"])
shapiro.test(penalties$Previous[penalties$Scored=="Missed"])


# A parameter of the 2-sample t.test is var.equal: indicates
#  whether variances of the 2 samples are equal.
# Could test for equality of variances:
#  var.test if normally distributed, ansari.test if not
# However it is considered safe to use t.test default var.equal=F 
#    even if variances are, in fact, equal.
# So can skip test of variances unless t.test results 
#    turn out to be very near significance level alpha

# T-tests to compare values of scored with missed
t.test(PSWQ ~ Scored, data=penalties) #default var.equal=F
t.test(Anxious ~ Scored, data=penalties)
# Can test variances:
var.test(penalties$Anxious[penalties$Scored=="Missed"],penalties$Anxious[penalties$Scored=="Scored"])
# Since variances are the same, we can redo t.test with var.equal=T:
#   results turn out to be very close to those with var.equal=F
t.test(Anxious ~ Scored, data=penalties, var.equal=T)
t.test(Previous ~ Scored, data=penalties)



# Look for Correlations and partial correlations
cor(penalties[-4])
# accounts for the effect of controlled-for variables on both of the compared vars 
pcor(penalties[-4])     
# can convert Scored to 0,1 before passing to pcor (or cor) 
pcor(as.numeric(penalties)) 

#cor and pcor require that dichotomous variables be coded as 0,1
cor.test(penalties$Scored, penalties$Previous) 

#recode column Scored as 0,1
values <- ifelse(penalties$Scored=="Missed",0,1)
penalties$Scored <-values
cor(penalties)
cor.test(penalties$Scored, penalties$Anxious)
cor.test(penalties$Scored, penalties$PSWQ) # method="spearman" almost same
# r= -0.6676552 so R^2= 0.4457635: PSWQ accounts for 44.6% of the variabliltiy of Scored

# Variable selection: start with all variables, then see how the 
#  model is effected if we leave one out
# AIC (Akaike information criterion): relative estimate of information lost 
#    smaller is good: less informaiton lost
fit_all <- glm(Scored ~ ., family=binomial(link="logit"), data=penalties)
summary(fit_all)  #AIC: 55.416
# log (p(Scored)/p(Missed) ) = -0.25137 PSWQ 
#  for other estimated coef, they are likely to have this fit by chance

# NB: If we decide to use any model we should rerun without the insigificant terms

fit.minus.prev <- glm(Scored ~ . -Previous, family=binomial(link="logit"), data=penalties)
summary(fit.minus.prev) #AIC: 56.074
# log (p(Scored)/p(Missed) ) = -0.2264 PSWQ -0.1190 Anxious , intercept coef is insignificant

fit.minus.anxious <- glm(Scored ~ . -Anxious, family=binomial(link="logit"), data=penalties)
summary(fit.minus.anxious) # AIC: 54.662 # best of these models we ran that have intercept
# log (p(Scored)/p(Missed) ) = -0.23009  PSWQ + 0.06480 Previous, intercept coef is insignificant 

fit.minus.PSWQ <- glm(Scored ~ . -PSWQ, family=binomial(link="logit"), data=penalties)
summary(fit.minus.PSWQ) # AIC: 67.141
# all estimated coef are likely to have this fit by chance

fit.PSWQ <- glm(Scored ~ PSWQ, family=binomial(link="logit"), data=penalties)
summary(fit.PSWQ) # AIC: 64.516
# log (p(Scored)/p(Missed) ) = -0.29397  PSWQ + 4.90010
# Use this model to predict the probablilities of scoring vs PSWQ
pred <- predict(fit.PSWQ, type="response") 
g<- ggplot(penalties, aes(x=PSWQ, y=Scored))
g +  geom_point()  +
  geom_point(aes(x=PSWQ, y= pred), color="red") 

# Look at the distribution of the model residuals
ggplot(fit.PSWQ, aes(x=.resid)) + geom_histogram(binwidth=.2)
# Investigate goodness-of-fit of model using only PSQW
logisticPseudoR2s(fit.PSWQ)

# After fitting the above models  we find fit.minus.anxious is the best
# but the intercept is insignificant so rerun without it:
fit.minus.anxious.inter <- glm(Scored ~ . -Anxious -1, family=binomial(link="logit"), data=penalties)
summary(fit.minus.anxious.inter) # AIC: 53.255  # best of these models
# log (p(Scored)/p(Missed) ) = -0.18223  PSWQ + 0.07639 Previous **********
# Visualize the relationship of the fitted variables on righthand side to Scored 
#  overlaid with the probablilities of scoring as predicted by model
pred_best <- predict(fit.minus.anxious.inter, type="response") 
gbest <- ggplot(penalties, aes(x=I(-0.18223*PSWQ+0.07639*Previous), y=Scored))
gbest +  geom_point()  +
  geom_point(aes(x=I(-0.18223*PSWQ+0.07639*Previous), y= pred_best), color="red") 

# Ex of other models that would take time to investigate
fit.part<- glm(Scored ~ PSWQ + Previous + Anxious:PSWQ + Anxious:Previous, family=binomial(link="logit"), data=penalties)
summary(fit.part) #AIC: 57.797 
# all estimated coef are likely to have this fit by chance

# Use built-in step fucntion (with caution) to select variables 
#  to try differenct models and minimize AIC:
fit.null <- glm(Scored ~ 1, family=binomial(link="logit"), data=penalties)
summary(fit.null) #AIC:

fit.full <- glm(Scored ~ Previous*Anxious*PSWQ, family=binomial(link="logit"), data=penalties)  
summary(fit.full) #AIC: 58.392 
# all estimated coef are likely to have this fit merely by chance

# Caution in using the step function to automate the selection of variables
# Is known to sometimes miss the best model
fit.step <- step(fit.null, 
    scope=list(lower=fit.null, upper=fit.full), direction=  "both")
# best fit: Scored ~ PSWQ + Previous + PSWQ:Previous  AIC=54.25

fit_found <- glm(Scored ~ PSWQ + Previous + PSWQ:Previous,family=binomial(link="logit"), data=penalties)
summary(fit_found)
# log (p(Scored)/p(Missed) ) = -0.584338 PSWQ, other coeficients insignificant

# Pasted from above:
# fit.minus.anxious <- glm(Scored ~ . -Anxious, family=binomial(link="logit"), data=penalties)
# log (p(Scored)/p(Missed) ) = -0.23009  PSWQ + 0.06480 Previous 

# Investigate goodness-of-fit of this model
logisticPseudoR2s(fit.minus.anxious)


# Use this model to estimate the probability of scoring if 
# both PSWQ and Previous are their mean values
newdf = data.frame( PSWQ=mean(penalties$PSWQ), Previous=mean(penalties$Previous), Anxious=mean(penalties$Anxious))
predict(fit.minus.anxious, newdf, type="response") 

# Do same with the insignificant  intercept removed:
predict(fit.minus.anxious.inter, newdf, type="response") 


library (stats)
# test fit.minus.anxious improvement over null   *****************************
modelChi <- fit.minus.anxious$null.deviance - fit.minus.anxious$deviance
chidf <- fit.minus.anxious$df.null - fit.minus.anxious$df.residual

# pchisq is cumulative distribution function for the chi-squared 
#  (chi^2) distribution with chidf degrees of freedom for
#   quantile = modelChi = null deviance - model deviance
# So (1 - pchisq) is the prob of a test stastistic this good or better by chance
chisq.prob <- 1 - pchisq(modelChi, chidf)
modelChi; chidf; chisq.prob # chisq.prob = 1.1533e-12
# indicates improved fit in this model over intercept only (null model) is significant



#Week10
#part1
setwd("/Users/Pat/Documents/R/HS_616/lecture_scripts")

temper <- read.csv("weather.data.csv")
plot(temper$month,temper$upper)
plot(temper$yr,temper$upper)

fit_month0 <- lm(temper$upper~ temper$month)
summary(fit_month0)
plot(fit_month0,1)

temper2 <- temper   ##################
## fill in your code here to correct the problem with this dataset

plot(temper2$month,temper2$upper)
plot(temper2$yr,temper2$upper)

fit_month <- lm(temper2$upper~ temper2$month)
summary(fit_month)
plot(fit_month, 1)

fit_month_1 <- lm(temper2$upper~ temper2$month -1)
summary(fit_month_1)
plot(fit_month_1,1)

fit_month_yr<- lm(temper2$upper~ temper2$month + temper2$yr)
summary(fit_month_yr)
plot(fit_month_yr,1)

fit_month_yr_1<- lm(temper2$upper~ temper2$month + temper2$yr -1)
summary(fit_month_yr_1)
plot(fit_month_yr_1,1)


# fill in your exploratory code here


#part2

# set the working directory
setwd("/Users/Pat/Documents/R/HS_616/lecture_scripts")

# There is a MySQL database for public access at genome-mysql.cse.ucsc.edu.
# This server allows MySQL access to the same set of data currently available on the public UCSC Genome Browser site

#install.packages("RMySQL")  # Run the first time
install.packages("sqldf")
library(RMySQL) # Database Interface and 'MySQL' Driver for R
library(sqldf)  # Perform SQL selects on R data frames

# Adapted from: http://playingwithr.blogspot.com/2011/05/accessing-mysql-through-r.html
#Establish a connection to the UCSC genone browser
con = dbConnect(MySQL(), user='genome', dbname='hg19', host='genome-mysql.cse.ucsc.edu')
# Return a list of the tables in our connection
dbListTables(con)
# Return a list of the fields in a specific table
dbListFields(con, 'knownGene')
dbListFields(con, 'refGene')

#Run a query
# To retrieve results a chunk at a time, use dbSendQuery, dbFetch, then dbClearResult
# Alternatively, if you want all the results (and they'll fit in memory) use dbGetQuery
#  which sends, fetches and clears for you.
resultSet <- dbSendQuery(con, 'SELECT * FROM refGene')

# Fetch records from a previously executed query
#  and save as a data frame object. 
#  n specifies the number of records to retrieve, n=-1 retrieves all pending records
hg19_refgene = dbFetch(resultSet,n=-1, stringsAsFactors=F)
str(hg19_refgene)
head(hg19_refgene)

# hg19_refgeneF = dbFetch(resultSet,n=-1, stringsAsFactors=T)
# str(hg19_refgeneF)
# head(hg19_refgeneF)


# A data frame is used for storing data tables. It is a list of vectors of equal length
# You can think of the vectors as columns in a database or excel spreadsheet
typeof(hg19_refgene)
colnames(hg19_refgene)
dbClearResult(resultSet)

# Disconnect when done with the mysql database
dbDisconnect(con)


# Add a new column to hold transcription start site (tss) relative to the strand
#  on the + strand, this is the left end of the range (txStart), on the - strand it is the right end (txEnd)
hg19_refgene$tss <- ifelse(hg19_refgene$strand == '+', hg19_refgene$txStart, hg19_refgene$txEnd)

write.table(hg19_refgene, file = "hg19_refgene2.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

# Create a subset of the data
hg19 <- hg19_refgene[, c("name", "chrom", "strand", "tss")]
summary(hg19)



#Diabetes_post

#VIF values greater than 10 may indicate that multicollinearity is unduly influencing your regression results. 
#If you see high VIF values, you may want to remove some of the correlated predictors from your model.

#install.packages("foreign")
library(foreign)
brca <- read.arff("Breast_Cancer_0123_PFL.arff")

logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance 
  nullDev <- LogModel$null.deviance 
  modelN <-  length(LogModel$fitted.values)
  R.l <-  1 -  dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2  ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2        ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2           ", round(R.n, 3),    "\n")
}

sum_na <-function(x)  {return (sum(is.na(x)))}
apply (brca, 2, sum_na)

sum(brca$Cyp1B1_48_loc3 ==3) #5
sum(brca$Cyp1B1_119_loc4 ==3) #5
sum(brca$GSTM1_loc8 ==3) #9
names(brca)
is.na(brca$Cyp1B1_48_loc3) <- which(brca$Cyp1B1_48_loc3 ==3)
is.na(brca$Cyp1B1_119_loc4) <- which(brca$Cyp1B1_119_loc4 ==3)
is.na(brca$GSTM1_loc8) <- which(brca$GSTM1_loc8 ==3)

apply (brca, 2, sum_na)


# family=binomial(link="log") indicates fit to poisson rather than logistic
epirr <- read.arff("DataSets_2loc/EPIRR.34-caco10.arff")
fit_epirr <- glm(class ~ location5 + location10 , family=binomial(link="log"), data=epirr)
summary(fit_epirr)
logisticPseudoR2s(fit_epirr)  # Nagelkerke R^2  0.281, with link=log: R^2=0.316, AIC=2242.5

fit2_epirr <- glm(class ~ location5*location10 , family=binomial(link="log"), data=epirr)
summary(fit2_epirr)
logisticPseudoR2s(fit2_epirr)# Nagelkerke R^2   0.346, with link=log: R^2=0.346, AIC=2189.1

fit3_epirr <- glm(class ~ location5:location10 , family=binomial(link="log"), data=epirr)
summary(fit3_epirr) # deviances and AIC are same as fit2_epirr, but coef are now sig
logisticPseudoR2s(fit3_epirr)# Nagelkerke R^2   0.346
# but location52:location102 not defined because of singularities

#when we use -1 to remove the constant term, singularities are gone,
#  we get a coef for location52:location102
fit4_epirr <- glm(class ~ location5:location10 -1 , family=binomial(link="logit"), data=epirr)
summary(fit4_epirr)# should be interactions only
logisticPseudoR2s(fit4_epirr)# Nagelkerke R^2   0.346

plot(fit_epirr,1)
plot(fit2_epirr,1)
plot(fit3_epirr,1)
plot(fit4_epirr,1)


#Exploration
sum(epirr$location5==2 & epirr$location10==2 ) #873  ** striking
sum(epirr$location5==2 & epirr$location10==1 ) #289
sum(epirr$location5==2 & epirr$location10==0 ) #66
sum(epirr$location5==1 & epirr$location10==2 ) #263
sum(epirr$location5==0 & epirr$location10==2 ) #61

sum(epirr$location5==1 & epirr$location10==1 ) #288
sum(epirr$location5==0 & epirr$location10==0 ) #18


#geneInterPost

#VIF values greater than 10 may indicate that multicollinearity is unduly influencing your regression results. 
#If you see high VIF values, you may want to remove some of the correlated predictors from your model.


setwd("/Users/Pat/Documents/R/HS_617/geneInters")

#install.packages("foreign")
library(foreign)
brca <- read.arff("Breast_Cancer_0123_PFL.arff")


logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance 
  nullDev <- LogModel$null.deviance 
  modelN <-  length(LogModel$fitted.values)
  R.l <-  1 -  dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2  ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2        ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2           ", round(R.n, 3),    "\n")
}

sum_na <-function(x)  {return (sum(is.na(x)))}
apply (brca, 2, sum_na)

sum(brca$Cyp1B1_48_loc3 ==3) #5
sum(brca$Cyp1B1_119_loc4 ==3) #5
sum(brca$GSTM1_loc8 ==3) #9
names(brca)
is.na(brca$Cyp1B1_48_loc3) <- which(brca$Cyp1B1_48_loc3 ==3)
is.na(brca$Cyp1B1_119_loc4) <- which(brca$Cyp1B1_119_loc4 ==3)
is.na(brca$GSTM1_loc8) <- which(brca$GSTM1_loc8 ==3)

apply (brca, 2, sum_na)


# family=binomial(link="log") indicates fit to poisson rather than logistic
epirr <- read.arff("DataSets_2loc/EPIRR.34-caco10.arff")
fit_epirr <- glm(class ~ location5 + location10 , family=binomial(link="log"), data=epirr)
summary(fit_epirr)
logisticPseudoR2s(fit_epirr)  # Nagelkerke R^2  0.281, with link=log: R^2=0.316, AIC=2242.5

fit2_epirr <- glm(class ~ location5*location10 , family=binomial(link="log"), data=epirr)
summary(fit2_epirr)
logisticPseudoR2s(fit2_epirr)# Nagelkerke R^2   0.346, with link=log: R^2=0.346, AIC=2189.1

fit3_epirr <- glm(class ~ location5:location10 , family=binomial(link="log"), data=epirr)
summary(fit3_epirr) # deviances and AIC are same as fit2_epirr, but coef are now sig
logisticPseudoR2s(fit3_epirr)# Nagelkerke R^2   0.346
# but location52:location102 not defined because of singularities

#when we use -1 to remove the constant term, singularities are gone,
#  we get a coef for location52:location102
fit4_epirr <- glm(class ~ location5:location10 -1 , family=binomial(link="logit"), data=epirr)
summary(fit4_epirr)# should be interactions only
logisticPseudoR2s(fit4_epirr)# Nagelkerke R^2   0.346

plot(fit_epirr,1)
plot(fit2_epirr,1)
plot(fit3_epirr,1)
plot(fit4_epirr,1)



#Exploration
sum(epirr$location5==2 & epirr$location10==2 ) #873  ** striking
sum(epirr$location5==2 & epirr$location10==1 ) #289
sum(epirr$location5==2 & epirr$location10==0 ) #66
sum(epirr$location5==1 & epirr$location10==2 ) #263
sum(epirr$location5==0 & epirr$location10==2 ) #61

sum(epirr$location5==1 & epirr$location10==1 ) #288
sum(epirr$location5==0 & epirr$location10==0 ) #18


#Week11 SVM

#Practical session: Introduction to SVM in R
#Jean-Philippe Vert
#https://escience.rpi.edu/data/DA/svmbasic_notes.pdf

#install.packages("kernlab")
#install.packages("ROCR")
library(kernlab)
library(kernlab)
n <- 150 # number of data points
p <- 2 # dimension

sigma <- 1 # variance of the distribution
meanpos <- 0 # centre of the distribution of positive examples
meanneg <- 3 # centre of the distribution of negative examples
npos <- round(n/2) # number of positive examples
nneg <- n-npos # number of negative examples

# Generate the positive and negative examples
xpos <- matrix(rnorm(npos*p,mean=meanpos,sd=sigma),npos,p)
xneg <- matrix(rnorm(nneg*p,mean=meanneg,sd=sigma),npos,p)
x <- rbind(xpos,xneg)

# Generate the labels
y <- matrix(c(rep(1,npos),rep(-1,nneg)))

# Visualize the data
plot(x,col=ifelse(y>0,1,2))
legend("topleft",c('Positive','Negative'),col=seq(2),pch=1,text.col=seq(2))

#Now we split the data into a training set (80%) and a test set (20%):
## Prepare a training and a test set ##
ntrain <- round(n*0.8) # number of training examples
tindex <- sample(n,ntrain) # indices of training samples
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]
ytest <- y[-tindex]
istrain=rep(0,n)
istrain[tindex]=1

# Visualize
plot(x,col=ifelse(y>0,1,2),pch=ifelse(istrain==1,1,2))
legend("topleft",c('Positive Train','Positive Test','Negative Train','Negative Test'),
       col=c(1,1,2,2),pch=c(1,2,1,2),text.col=c(1,1,2,2))

# load the kernlab package
library(kernlab)

# train the SVM
svp <- ksvm(xtrain,ytrain,type="C-svc",kernel='vanilladot',C=100,scaled=c())

# General summary
svp

# Attributes that you can access
attributes(svp)
# For example, the support vectors
alpha(svp) # alpha vector (possibly scaled)
# indices of support vectors in data matrix
#  after the possible effect of na.omit and subset
alphaindex(svp)
coef(svp) #The corresponding coefficients times the training labels
# Use the built-in function to pretty-plot the classifier
plot(svp,data=xtrain)
#              plot(scale(x), col=y+2, pch=y+2, xlab="", ylab="")
w <- colSums(coef(svp)[[1]] * x[unlist(alphaindex(svp)),])
b <- b(svp)

# Predict labels on test
ypred = predict(svp,xtest)
table(ytest,ypred)
# Compute accuracy
sum(ypred==ytest)/length(ytest)
# Compute at the prediction scores
ypredscore = predict(svp,xtest,type="decision")










#######Schedule

week 0 (before first class meeting get up and runningon R and R studio): 
ch 1-3 from Lander text 
 Getting R
 The R Environment
 R packages

week 1 ch 3-4 from Lander text 
 Basics of R 
 Advanced Data Structures (include data tables)

week 2  ch 6-7 from Lander text 
 Introduction to R Markdown
 Reading Data into R
 Statistical Graphics

week 3 ch 8-9, gloss over 10 from Lander text   
 Writing R Functions  
 Control Statements
 Loops: The Un-R Way to Iterate (gloss over, not on quizzes)

week 4 ch 11-12 from Lander text 
 Group Manipulation
 Data Reshaping

week 5  ch 13-14 from Lander text 
 Manipulating Strings 
 Probability Distributions

week 6 ch 15 from Lander text 
 Basic Statistics
 Brief Intro to ANOVA 

week 7  ch 16 from Lander text  
 Linear Models (in depth)  

week 8 ch 17-18 from Lander text 
 Genearalized Linear Models  (17.1) 
 Survivial Analysis  (17.4)
 (gloss over other Genearalized Linear Models)
 Model Diagnostics - basic linear algebra, ANOVA for model comparison
 Residuals
 Bootstrap

week 9
 Using SQL queries with R
 Data Simulation

week 9
 Machine Learning   SVM
 Learning curves
 ROC curves: sensitivity/specificity tradeoff

week 10
 RF
 Regularization and Shrinkage (ch 16 from Lander text)

week 11
 Dimension reduction
 PCA
 Clustering (ch 22 from Lander text)

week 12
 Nonlinear Models (not on quiz) (ch 20 from Lander text)

week 13
 Using SQL databases with R

week 14
 Visualization

week 15
Quick topics:
 Time Series and Auto-Correlation (ch 21 from Lander text)
 Eigenvectors + solving systems of equations + Gauss-Jordan elimination
 Web scraping
 Interactive analysis









