library(data.table)
library(MASS)
library(gtools)
library(mvtnorm)

c9 <- as.data.table(read.table(file = "creatinine.csv", sep = ",", header = T))
# Coercing the dates to the correct format
c9[, result_date := as.Date(result_date, origin = "1900-01-01")]

# Choosing patients having more than 5 test results
serious.patients <- c9[, .N >= 5, by = Person_ID][V1 == T, Person_ID]
# Removing un-serious patients from the table
c9 <- c9[Person_ID %in% serious.patients]

monthsDiff <- function(d1, d2){
    # Calculate months from an origin data
    monthsFromOrigin <- function(d){
        posixData <- as.POSIXlt(d)
        posixData$year * 12 + posixData$mon
    }
    
    monthsFromOrigin(d1) - monthsFromOrigin(d2)
}

# result_date = Date - Date[1] <in months>, to make all records start from month 0
c9 <- c9[, result_date_months := as.numeric(monthsDiff(result_date, result_date[1])),
             by = Person_ID]
print("Cleaned Data Set Loaded")
###############################################################################################
source('~/740 Project/makeCachePolyFits.R')
source('~/740 Project/calculatePolyFits.R')
source('~/740 Project/givePBetas.R')

degree <- 4
# To determine the mean and the covariance matrix for the normal distribution of the 
# Betas, the dataset is sampled from again and again, and on the resultant samples the lm is 
# fit to get a variety of Beta values.
poly.fits.cacheFunctions <- makeCachePolyFits()

if(!file.exists("~/740 Project/polyFits.Rda")){
    poly.fits <- calculatePolyFits(c9, degree)
    poly.fits.cacheFunctions$setCache(poly.fits)
}
# Calculated tonnes of Beta values from the data
poly.fits <- poly.fits.cacheFunctions$getCache()
# Computing the covariance matrix for the poly.fits
betas.covariance.matrix <- cov(poly.fits)
# Calculating the means of the Betas
betas.mean <- apply(poly.fits, 2, mean)
# Using a Multivariate normal distribution to get the intial values of the betas
betas <- mvrnorm(n = 12, betas.mean, betas.covariance.matrix)

multinomial.pi <- as.vector(rdirichlet(n = 1, alpha = rep(2, 12)))

# Outer Loop from here
multinomial.probs <- as.vector(rmultinom(n = 1, size = length(unique(c9$Person_ID)), prob = multinomial.pi))

person_id <- unique(c9$Person_ID)
values1 <- list()
values2 <- list()
expectation.probabilities <- matrix(nrow = length(person_id), ncol = 12)

# Expectation Step    
for(i in 1:length(unique(c9$Person_ID))){
     Xs <- c9[Person_ID == person_id[i], result_date_months]
     phiMat <- data.table(X1 = Xs)
     phiMat <- phiMat[, `:=`(X2 = X1^2, X3 = X1^3, X4 = X1^4)]
     expectation.mean <- as.matrix(phiMat) %*% t(betas)
     s <- ((0.1)^2) * diag(length(Xs))
     
     normal.prob <- numeric(length = 12)
     for(j in 1:12){
        normal.prob[j] <- pmvnorm(lower = -Inf, upper = Xs, mean = expectation.mean[, j], sigma = s)[1]
     }
     expectation.probabilities[i, ] <- normal.prob * multinomial.probs
     normalizing.constant <- sum(expectation.probabilities[i, ])
     expectation.probabilities[i, ] <- expectation.probabilities[i, ] / normalizing.constant
     
     # Values for the Maximization step
     values1[[i]] <- expectation.probabilities[i, ] %*% t(as.matrix(phiMat)) %*% solve(s) %*% phiMat
     values2[[i]] <- expectation.probabilities[i, ] %*% t(as.matrix(phiMat)) %*% solve(s) %*% as.matrix(Xs)
}

alpha <- 2
new.pis <- vector(mode = "numeric", length = 12L)
new.betas <- matrix(nrow = 12, ncol = degree)
sum.value1 <- matrix(data = 0, nrow = degree, ncol = degree)
sum.value2 <- matrix(data = 0, nnrow = degree, ncol = 1)

#Maximization Step
add <- function(x) Reduce("+", x)
sum.value1 <- add(values1)
sum.value2 <- add(values2)
    
for(i in 1:12){
    new.pis[i] <- sum(expectation.probabilities[, i]) + alpha - 1
}
normalizing.constant <- sum(new.pis)
new.pis <- new.pis / normalizing.constant

betas.covariance.inverse <- solve(betas.covariance.matrix)
# For summing the values in Bracket 1
values1 <- list()
values1[[1]] <- betas.covariance.inverse
values1[[2]] <- sum.value1

# For summing the values in Bracket 2
values2 <- list()
values2[[1]] <- betas.covariance.inverse %*% betas.mean
values2[[2]] <- sum.value2

new.betas <- solve(add(values1)) %*% (add(values2))