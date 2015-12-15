library(data.table)
library(MASS)
library(gtools)
library(mvtnorm)
library(reshape2)

dest <- paste0("~/", gsub(pattern = "(\\s|:)", replacement = "-", x = Sys.time()))
dir.create(path = dest)
setwd(dest)

c9 <- as.data.table(read.table(file = "creatinine.csv", sep = ",", header = T))
# Coercing the dates to the correct format
c9[, result_date := as.Date(result_date, origin = "1900-01-01")]
c9[, group := -1]

# Choosing patients having more than 5 test results
serious.patients <- c9[, .N >= 10, by = Person_ID][V1 == T, Person_ID]
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
betas.mean <- rep(0, degree+1)
# Using a Multivariate normal distribution to get the intial values of the betas
betas <- mvrnorm(n = g, betas.mean, betas.covariance.matrix)

multinomial.pi <- as.vector(rdirichlet(n = 1, alpha = rep(alpha, g)))

person_id <- unique(c9$Person_ID)
values1 <- list()
values2 <- list()
expectation.probabilities <- matrix(nrow = length(person_id), ncol = g)

convergenceCondition <- 1
iterations <- 0
logsumexp <- function(x)
{
    m <- max(x)
    m + log(sum(exp(x - m)))
}

while(convergenceCondition > 0.001 && iterations < 20){
    # Expectation Step
    multinomial.probs <- as.vector(rmultinom(n = 1, size = length(unique(c9$Person_ID)), prob = multinomial.pi))
    
    for(i in 1:length(unique(c9$Person_ID))){
         Xs <- c9[Person_ID == person_id[i], result_date_months]
         phiMat <- data.table(X0 = rep(1, length(Xs)))
         phiMat <- phiMat[, `:=`(X1 = Xs, X2 = Xs^2, X3 = Xs^3, X4 = Xs^4)]
         expectation.mean <- as.matrix(phiMat) %*% t(betas)
         s <- 1 * diag(length(Xs))

         normal.prob <- numeric(length = g)
         egfrs <- c9[Person_ID == person_id[i], eGFR]
         for(j in 1:g){
            normal.prob[j] <- dmvnorm(x = egfrs,
                                      mean = as.vector(expectation.mean[, j]), sigma = s, log = T)
         }
         
         expectation.probabilities[i, ] <- normal.prob + log(multinomial.probs) -
                                        logsumexp(normal.prob + log(multinomial.probs))
         expectation.probabilities[i, ] <- exp(expectation.probabilities[i, ])
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
        #diagMat <- 0.1 * diag(nrow = degree + 1, ncol = degree + 1)
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
    convergenceCondition <- sum(abs((new.betas - betas)/betas))
    betas <- new.betas
    #betas.mean <- apply(betas, 2, mean)
    #betas.covariance.matrix <- cov(betas)
    multinomial.pi <- new.pis
    iterations <- iterations + 1
    
    # Generate Plot for Beta values
    j <- 1
    for(i in unique(c9$Person_ID)){
            maximum <- max(expectation.probabilities[j, ])
            c9[Person_ID == i, group := which(maximum == expectation.probabilities[j, ])]
    }
    # 60 : months for 5 years
    gcolors <- c("red", "blue", "black", "green", "purple", "cyan", "yellow", "brown", "orange",
                 "orchid", "navy")
    
    plines <- function(x, row){ 
        eq <- 0
        for(i in 1:(degree+1))
            eq <- eq + betas[row, i] * (x^(i-1))
        return(eq)
    }
    
    xvalues <- data.table(x = seq(from = 0, to = 3, by = 0.1))
    
    for(i in 1:g){
        expr <- parse(text = paste("Group", i, " := sapply(xvalues$x, plines, i)", sep = ""))
        xvalues[, eval(expr)]
    }
    xvalues <- melt(xvalues, measure.vars = names(xvalues)[2:length(names(xvalues))], variable.name = "Groups")
    ggplot(data = xvalues, aes(x = x, y = value, colour = variable)) + geom_line() + 
        ylab("eGFRs") + xlab("years") + ggtitle("Beta plot")
    ggsave("Beta Curve.png")
    # the print statement can be modified to save the plot to disk
}

print("Folder Created")
for(i in 1:g){
    
    groupset <- c9[, .SD[group == i], by = group]
    xvalues <- data.frame(Person_ID = factor(),
                          result_date_months = numeric(),
                          eGFR = numeric(),
                          eGFR_hat = numeric())
    
    for(j in head(unique(groupset$Person_ID), n = 6)){
        tempset <- groupset[Person_ID == j, .(Person_ID, result_date_months, eGFR)]
        tempset <- tempset[, eGFR_hat := sapply(result_date_months, plines, i)]
        
        xvalues <- rbindlist(list(xvalues, tempset))
    }
    group.plot <- ggplot() + 
        geom_point(data = xvalues, aes(x = result_date_months, y = eGFR, colour = Person_ID)) +
        geom_line(data = xvalues, aes(x = result_date_months, y = eGFR_hat, colour = Person_ID)) +
        ggtitle(paste("Group ", i, sep = "")) + xlab("Years") +
        xlim(0, 3) + ylim(0, 200)
    filename <- paste("Group", i, sep = " ")
    filename <- paste(filename, ".png")
    #print(group.plot)
    ggsave(filename = filename, plot = group.plot)
    
    group.plot <- group.plot + geom_density()
    filename <- paste("GroupWithIndividualLines", i, sep = " ")
    filename <- paste(filename, ".png")
    ggsave(filename = filename, plot = group.plot)
}