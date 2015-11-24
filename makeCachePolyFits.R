makeCachePolyFits <- function(){
    
    #resetCache <- function() saveRDS(NULL, file = "polyFits.Rda")
    setCache <- function(poly.fits) saveRDS(poly.fits, file = "~/740 Project/polyFits.Rda")
    getCache <- function() readRDS(file = "~/740 Project/polyFits.Rda")
    
    list(setCache = setCache, getCache = getCache)
}