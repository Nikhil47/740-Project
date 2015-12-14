# Divide data in 12 parts and fit polynomial curve of degree 3 to it
givePBetas <- function(c9, parts, poly.fits, iteration, degree){
    for(i in 1:12){
        fitting.data <- c9[Person_ID == parts[[i]]]
        fit <- lm(data = fitting.data, 
                  formula = fitting.data$eGFR ~ poly(fitting.data$result_date_months, 
                                                     degree = degree, raw = T))
        poly.fits[12*(iteration - 1) + i, ] <- as.vector(summary(fit)$coef[, 1])
    }
    poly.fits
}