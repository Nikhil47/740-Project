library(data.table)
library(MASS)
library(gtools)
library(mvtnorm)
library(reshape2)

c9 <- as.data.table(read.table(file = "creatinine.csv", sep = ",", header = T))
# Coercing the dates to the correct format
c9[, result_date := as.Date(result_date, origin = "1900-01-01")]
c9[, group := -1]

# Choosing patients having more than 5 test results
serious.patients <- c9[, .N >= 5, by = Person_ID][V1 == T, Person_ID]
# Removing un-serious patients from the table
c9 <- c9[Person_ID %in% serious.patients]

biWeekDiff <- function(d1, d2){
    # Calculate months from an origin data
#     monthsFromOrigin <- function(d){
#         posixData <- as.POSIXlt(d)
#         posixData$year * 12 + posixData$mon
#     }
#     
#     (monthsFromOrigin(d1) - monthsFromOrigin(d2))
    difftime(d1, d2, units = "weeks")/52
}

# result_date = Date - Date[1] <in months>, to make all records start from month 0
c9 <- c9[, result_date_months := as.numeric(biWeekDiff(result_date, result_date[1])),
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
print("PolyFits Loaded")

g <- 12
alpha <- 5

# Calculated tonnes of Beta values from the data
poly.fits <- poly.fits.cacheFunctions$getCache()
# Computing the covariance matrix for the poly.fits
# betas.covariance.matrix <- cov(poly.fits)
betas.covariance.matrix <- diag(x = 10^2, nrow = degree + 1, ncol = degree + 1)
# Calculating the means of the Betas
# betas.mean <- apply(poly.fits, 2, mean)
beats.mean <- rep(0, degree+1)
# Using a Multivariate normal distribution to get the intial values of the betas
betas <- mvrnorm(n = g, betas.mean, betas.covariance.matrix)

multinomial.pi <- as.vector(rdirichlet(n = 1, alpha = rep(alpha, g)))

multinomial.probs <- as.vector(rmultinom(n = 1, size = length(unique(c9$Person_ID)), prob = multinomial.pi))

person_id <- unique(c9$Person_ID)
values1 <- list()
values2 <- list()
expectation.probabilities <- matrix(nrow = length(person_id), ncol = g)

while(convergenceCondition > 0.001){
    # Expectation Step    
    for(i in 1:length(unique(c9$Person_ID))){
         Xs <- c9[Person_ID == person_id[i], result_date_months]
#          phiMat <- data.table(X1 = Xs)
#          phiMat <- phiMat[, `:=`(X2 = Xs^2, X3 = Xs^3, X4 = Xs^4)]
         phiMat <- data.table(X0 = rep(1, length(Xs)))
         phiMat <- phiMat[, `:=`(X1 = Xs, X2 = Xs^2, X3 = Xs^3, X4 = Xs^4)]
         expectation.mean <- as.matrix(phiMat) %*% t(betas)
         s <- 1 * diag(length(Xs))
#        normal.prob <- rep(1, times = g)
         normal.prob <- numeric(length = g)
         egfrs <- c9[Person_ID == person_id[i], eGFR]
         for(j in 1:g){
            normal.prob[j] <- dmvnorm(x = egfrs,
                                      mean = as.vector(expectation.mean[, j]), sigma = s, log = T)
         }
         
         normalizing.constant <- sum(normal.prob)
         normal.prob <- normal.prob / normalizing.constant
         
         expectation.probabilities[i, ] <- normal.prob * multinomial.probs
         normalizing.constant <- sum(expectation.probabilities[i, ])
         expectation.probabilities[i, ] <- expectation.probabilities[i, ] / normalizing.constant
         
         # Values for the Maximization step
         values1[[i]] <- t(as.matrix(phiMat)) %*% solve(s) %*% as.matrix(phiMat)
         values2[[i]] <- t(as.matrix(phiMat)) %*% solve(s) %*% as.matrix(egfrs)
    }
    print("Expectation Step Completed")
    
    new.pis <- vector(mode = "numeric", length = g)
    new.betas <- matrix(nrow = g, ncol = degree+1)
#     sum.value1 <- matrix(data = 0, nrow = degree, ncol = degree)
#     sum.value2 <- matrix(data = 0, nrow = degree, ncol = 1)
    
    #Maximization Step
    add <- function(x) Reduce("+", x)
    
    betas.covariance.inverse <- solve(betas.covariance.matrix)
    # For summing the values in Bracket 1
    Bracket1 <- list()
    Bracket1[[1]] <- betas.covariance.inverse
    
    # For summing the values in Bracket 2
    Bracket2 <- list()
    Bracket2[[1]] <- betas.covariance.inverse %*% betas.mean
    
    for(j in 1:g){
        temp.list1 <- list()
        temp.list2 <- list()
        for(i in 1:length(unique(c9$Person_ID))){
            temp.list1[[i]] <- expectation.probabilities[i, j] * values1[[i]]
            temp.list2[[i]] <- expectation.probabilities[i, j] * values2[[i]]
        }
        
        Bracket1[[2]] <- add(temp.list1)
        Bracket2[[2]] <- add(temp.list2)
        diagMat <- 0.1 * diag(nrow = degree + 1, ncol = degree + 1)
        #new.betas[j, ] <- solve(add(Bracket1) + diagMat) %*% (add(Bracket2))
        new.betas[j, ] <- solve(add(Bracket1)) %*% (add(Bracket2))
    }
    
    # Calculating Pis
    for(i in 1:g){
        new.pis[i] <- sum(expectation.probabilities[, i]) + alpha - 1
    }
    normalizing.constant <- sum(new.pis)
    new.pis <- new.pis / normalizing.constant
    
    print("Old Betas")
    print(betas)
    print("New Betas")
    print(new.betas)
    
    # New values for the next iteration
    betas <- new.betas
    betas.mean <- apply(betas, 2, mean)
    betas.covariance.matrix <- cov(betas)
    convergenceCondition <- sum(abs(new.betas - betas)/betas)
    multinomial.probs <- new.pis
    
    # Generate Plot for Beta values
    j <- 1
    for(i in unique(c9$Person_ID)){
        maximum <- max(expectation.probabilities[j, ])
        c9[person_id == i, group := which(maximum == expectation.probabilities[j, ])]
        j <- j + 1
    }
    # 60 : months for 5 years
    gcolors <- c("red", "blue", "black", "green", "purple", "cyan", "yellow", "brown", "orange",
                 "orchid", "navy")
    xvalues <- data.frame(x = seq(from = 1, to = 7))
    betas.plot <- ggplot(data = xvalues, aes(x = x))
    plines <- function(x, row){ 
        eq <- 0
        for(i in 1:degree+1)
            eq <- eq + betas[row, i] * (x^i)
        return(eq)
    }
    
    for(i in 1:g){
        betas.plot <- betas.plot + 
            stat_function(fun = plines, args = list(row = i),
                          aes(colour = paste("line", i, sep = " ")))
    }
    betas.plot <- betas.plot + 
                    scale_color_manual("Groups", values = gcolors)
    print(betas.plot)
    # the print statement can be modified to save the plot to disk
}

dest <- paste0("~/", gsub(pattern = "(\\s|:)", replacement = "-", x = Sys.time()))
dir.create(path = dest)
print("Folder Created")
for(i in 1:g){
    
    groupset <- c9[, .SD[group == i], by = group]
    group.plot <- ggplot()
    
    for(j in head(unique(groupset$Person_ID))){
        group.plot <- group.plot + 
            stat_function(fun = plines, data = as.data.frame(groupset[j, eGFR]), args = list(row = i))
    }
    ggsave(filename = paste("Group", i, sep = " "), path = dest)
}