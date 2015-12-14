calculatePolyFits <- function(c9, degree = 3, g){
    poly.fits = matrix(nrow = 12000, ncol = degree)
    for(i in 1:1000){
        # Divide unique ids in 12 parts
        ids <- sample(x = unique(c9$Person_ID), replace = F)
        parts <- split(ids, ceiling(seq_along(ids)/2977))
        
        poly.fits <- givePBetas(c9, parts, poly.fits, i, degree)
    }
    poly.fits
}