# Data Simulation
Xiangyi Xu  
Thursday, April 09, 2015  

####Pastor-Valero, M. (2013). Fruit and vegetable intake and vitamins C and E are associated with a reduced prevalence of cataract in a Spanish Mediterranean population. BMC Ophthalmology, 13:52. 



####Background
Cataract is among the major causes of vision impairment and blindness 
worldwide. Epidemiological studies support the role of antioxidants in the etioloty of cataract, but the evidence of one specific antioxidant over another is inconsistent. Few studies have examined the association of cataract with fruit and vegetable intake with inconclusive results. In this study, the relationship between cataract and fruit and vegetable intake and dietary and blood levels of carotenoids, vitamins C and E were examined in a Spanish Mediterranean population. 



####Data Simulation
The dataset is simulated from the above study. Among 433 elderly with cataract or cataract extraction, 54% are women and 46% are men, their average age iss around 72. The study shows that increasing quartiles  `(13~83, 83~107, 107~143, 143~408 mg/d)` of dietary intakes from `107mg/d` of vitamin C indicating a significant decreasing association with prevalence of cataract or cataract extraction. In the simulated dataset, with enough data, it will show that both age and intake of vitamin C matter much in having cataract or not. Noise such as zodiac signs are added in the study too. 
Hint: log of VitaminC for the analysis of data



####Features
Feature1: age    
Feature2: obesity- Y, N     
Feature3: sex- M, F     
Feature4: marital_status- Y, N     
Feature5: diabetes- Y, N    
Feature6: smoking- Never, Past, Current    
Feature7: VitaminC     
Feature8: zodiac signs     
Outcome:  cataract Y, N      



####Manipulate function is used to generate coefficients.

library(manipulate)
manipulate(with(df,{  
score <- 0.73 + 2^a*(age - 72) - 2^b*log(VitaminC)
prob <- logistic(score)
hist(prob, breaks=50)
}), a=slider(-9, 9, step=0.1, initial = 0), b=slider(-9, 9, step=0.1, initial = 0))



####generate dataset



```r
generate_dataset <- function(N=100){
        
        age <- runif(N, min=72-10, max=72+10)
        
        obesity <- sample(c("Y", "N"), N, replace=TRUE, prob=c(.36, .64))
        
        sex <- sample(c("M", "F"), N, replace=TRUE, prob=c(.46, .54))
        
        marital_status <- sample(c("Y", "N"), N, replace=TRUE, 
                                  prob=c(.70, .30))
                                        
        diabetes <- sample(c("Y", "N"), N, replace=TRUE, prob=c(.26, .74))
        
        smoking <- sample(c("Never", "Past", "Current"), N, replace=TRUE, 
                             prob=c(.52, .31, .17))
                             
        VitaminC <- rnorm(N, mean=107, sd=15)

        

        zodiac <- c("Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
                      "Libra","Scorpio", "Sagittarius", "Capricorn","Aquarius", 
                        "Pisces") 
        
        sign <- sample(zodiac, N, replace=TRUE)
        
        simulate_cataract <- function(age, VitaminC){
                logistic <- function(t) 1 / (1 + exp(-t))
                score <- 0.73 +  2^(-1.2)*(age - 72) - 2^(-1.2)*log(VitaminC)  
        
                prob <- logistic(score)
        
                result <- ifelse(runif(length(prob)) < prob, "Y", "N")
        } 
                
        cataract <- simulate_cataract(age, VitaminC)

        data.frame(age, sex, obesity, marital_status, smoking, VitaminC,
                         sign, cataract)
        

}

df <- generate_dataset(2e5)

head(df)
```

```
##        age sex obesity marital_status smoking  VitaminC     sign cataract
## 1 77.35373   F       N              Y    Past 120.34500   Gemini        Y
## 2 62.34716   F       N              N Current 124.10608   Pisces        N
## 3 73.75635   M       N              N Current  97.33514   Taurus        N
## 4 64.24056   M       N              Y Current 115.89233    Virgo        N
## 5 78.03271   F       N              Y    Past 118.38592    Aries        Y
## 6 74.21270   F       Y              Y   Never  88.19312 Aquarius        Y
```

```r
fit <- glm(cataract ~ I(age-72) + log(VitaminC), data=df, family="binomial")
summary(fit)
```

```
## 
## Call:
## glm(formula = cataract ~ I(age - 72) + log(VitaminC), family = "binomial", 
##     data = df)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.5613  -0.4805  -0.1574   0.4896   3.3777  
## 
## Coefficients:
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept)    0.403480   0.219439   1.839    0.066 .  
## I(age - 72)    0.434358   0.001888 230.047  < 2e-16 ***
## log(VitaminC) -0.367157   0.047062  -7.802 6.12e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 260015  on 199999  degrees of freedom
## Residual deviance: 139226  on 199997  degrees of freedom
## AIC: 139232
## 
## Number of Fisher Scoring iterations: 6
```

```r
plot(fit)
```

![](DataSimulation_files/figure-html/unnamed-chunk-1-1.png) ![](DataSimulation_files/figure-html/unnamed-chunk-1-2.png) ![](DataSimulation_files/figure-html/unnamed-chunk-1-3.png) ![](DataSimulation_files/figure-html/unnamed-chunk-1-4.png) 

```r
fit2 <- glm(cataract ~ age + log(VitaminC) + obesity + marital_status, data=df, family="binomial")
summary(fit2)
```

```
## 
## Call:
## glm(formula = cataract ~ age + log(VitaminC) + obesity + marital_status, 
##     family = "binomial", data = df)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.5643  -0.4805  -0.1574   0.4893   3.3754  
## 
## Coefficients:
##                   Estimate Std. Error  z value Pr(>|z|)    
## (Intercept)     -30.871666   0.258064 -119.628  < 2e-16 ***
## age               0.434364   0.001888  230.047  < 2e-16 ***
## log(VitaminC)    -0.367240   0.047063   -7.803 6.03e-15 ***
## obesityY         -0.014089   0.014086   -1.000    0.317    
## marital_statusY   0.009141   0.014743    0.620    0.535    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 260015  on 199999  degrees of freedom
## Residual deviance: 139225  on 199995  degrees of freedom
## AIC: 139235
## 
## Number of Fisher Scoring iterations: 6
```















