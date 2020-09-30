#Exercise 4.2 generalized linear models

timedeath = c(65,156,100,134,16,108,121,4,39,143,56,26,22,1,1,5,65)
whitecell = c(3.36 ,2.88 ,3.63 ,3.41 , 3.78 , 4.02 , 4.00 , 4.23 , 3.73 , 3.85 , 3.97 , 4.51 , 4.54 , 5.00 , 5.00 , 4.72 , 5.00)

leukemia = data.frame(timedeath, whitecell)

library(ggplot2)
#4.2 a

#Plot
ggplot(leukemia, aes(timedeath, whitecell))+ 
  geom_point() +
  ggtitle("Survival Time vs Initial White Blood Cell Count") +
  xlab("Survival Time (Weeks)") + ylab("Initial White Blood Cell Count")+
  geom_smooth(method = "lm")

# as y decreses, x increases (y appears to decrease exponentially)  

#4.2 b <  log, because log is a function that tranforms the equation to be linear
# The link function is the function of y that gives you Xb on the right hand side.

# 4.2 c (see pdf with latex)

# 4.2 d
# log on both sides E[log(Y)] = b1 + b_2 X

summary(lm(log(timedeath) ~ (whitecell), data=leukemia))

# b1 = 11.0738 b2= -1.8829

#usar glm para familia exponencial
#exp_model <- glm(timedeath~whitecell, family=Gamma(link="log"), data=leukemia)


# 4.2 e

#Predicted y-hats
leukemia["Predicted"] = exp(11.0738 - 1.8829 *(leukemia[2]))

# Linear models time of death with predicted
leukpred = lm((timedeath) ~ Predicted, data=leukemia) 
#Standarized Residuals
leukpred.stand = rstandard(leukpred)
#Plot
plot(leukemia$timedeath, leukpred.stand, ylab="Standardized Residuals", xlab="Waiting Time", main="Leukemia") 
abline(0, 0)                  # the horizon

#Model fits the data well.


#5.4 a (not sure)

exp_model <- glm(timedeath~whitecell, family=Gamma(link="log"), data=leukemia)

summary(exp_model)

beta_1 <- exp_model$coefficients[1]
beta_2 <- exp_model$coefficients[2]

#extracting parameter variances from variance-covariance matrix
betaVar <- c(vcov(exp_model)[1,1], vcov(exp_model)[2,2])

#using N(0,1) to derive confidence interval of 95%
CI_1 <- c(beta_1 - 1.96*std_error_1, beta_1 + 1.96*std_error_1)
CI_2 <- c(beta_2 - 1.96*std_error_2, beta_2 + 1.96*std_error_2)

CI_1
CI_2


"
firstmodel = lm(log(timedeath) ~ (whitecell), data=leukemia)
confint(firstmodel, level=0.90)
WaldTest = function(L,thetahat,Vn,h=0) # H0: L theta = h
  # Note Vn is the asymptotic covariance matrix, so it's the
  # Consistent estimator divided by n. For true Wald tests
  # based on numerical MLEs, just use the inverse of the Hessian.
{
  WaldTest = numeric(3)
  names(WaldTest) = c("W","df","p-value")
  r = dim(L)[1]
  W = t(L%*%thetahat-h) %*% solve(L%*%Vn%*%t(L)) %*%
    (L%*%thetahat-h)
  W = as.numeric(W)
  pval = 1-pchisq(W,r)
  WaldTest[1] = W; WaldTest[2] = r; WaldTest[3] = pval
  WaldTest
} # End function WaldTest
WaldTest()
"

#5.4 b

dev_diff <- exp_model$null.deviance - exp_model$deviance
dev_diff


#OTHER CODE
cor(leukemia$timedeath, leukemia$whitecell)
        
