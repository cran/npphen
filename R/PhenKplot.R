#' @title PhenKplot
#' @description Plot the most probable vegetation greenness values.
#' @encoding UTF-8
#' @param x	Numeric vector with greenness values
#' @param dates	Vector with dates at which the greenness values were recorded
#' @param h	Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h=1 if the vegetation is in the Northern Hemisphere (season starting at January 1st), h=2 if it is in the Southern Hemisphere (season starting at July 1st)
#' @param nGS Numeric. Number of greenness values within a single growing season. For example, nGS=23 for MODIS	Vegetation Index 16-days composites
#' @param xlab Character vector (or expression) giving plot title in x axis label
#' @param ylab Caracter vector (or expression) giving plot title in y axis label
#' @param rge  A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @details It is the graphical version of the Phen function. It calculates and plot a likelihood map of the vegetation-greenness – time space using a numeric vector of vegetation canopy greenness values (e.g. Leaf Area Index (LAI) or greenness proxies values such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). Also a vector with dates for the greenness values is required. This function calculates the confidence areas on a per year basis. Functions for confidence intervals on per day basis are under development. This function is partially based on the ci2d function on package gplots.
#' @seealso \code{\link{Phen}}
#' @examples
#' \donttest{
#' library(lubridate)
#'
#' ## Testing North Hemisphere data. Raster data from Slovenia (EVI index), h=1 ##
#'
#' # Load data
#' sl.path<-system.file("extdata/HN_slovenia",package="npphen")
#' sl_rasters<-list.files(path=sl.path, pattern=glob2rx("slovenia*.tif"), full.names=TRUE)
#' Slovenia_rasters<-stack(sl_rasters)
#' sl_dates<-read.csv(system.file("extdata/date_tables/Slovenia_dates.csv", package="npphen"))
#' Slovenia_dates <- as.Date(sl_dates$date, format='%d/%m/%Y')
#'
#' # Generate a Raster time series using a raster stack and a date database from Slovenia
#' # Obtain data from a particular pixel generating a time series
#' sl_pixel<-cellFromXY(Slovenia_rasters,c(474368,5096979))
#' sl_pixelts<-as.numeric(Slovenia_rasters[sl_pixel])
#' plot(Slovenia_dates,sl_pixelts, type='l')
#'
#' # Phenology for the given pixel
#' PhenKplot(x=sl_pixelts,dates=Slovenia_dates,h=1,nGS=23, xlab="DOY",
#' ylab="EVI", rge=c(0,10000))
#'
#'
#' ## Testing South Hemisphere data. Raster data from Chile (EVI index), h=2 ##
#'
#' # Load data
#' ay.path<-system.file("extdata/HS_aysen",package="npphen")
#' ayrasters<-list.files(path=ay.path, pattern=glob2rx("aysen*.tif"), full.names=TRUE)
#' Aysen_rasters<-stack(ayrasters)
#' ay_dates<-read.csv(system.file("extdata/date_tables/Aysen_dates.csv", package="npphen"))
#' Aysen_dates <- as.Date(ay_dates$date, format='%d/%m/%Y')
#'
#' # Generate a Raster time series using a raster stack and a date database from Aysen
#' # Obtain data from a particular pixel generating a time series
#' ay_pixel<-cellFromXY(Aysen_rasters,c(228373,4806975))
#' ay_pixelts<-as.numeric(Aysen_rasters[ay_pixel])
#' plot(Aysen_dates,ay_pixelts, type = 'l')
#'
#' # Phenology for the given pixel
#' PhenKplot(x=ay_pixelts,dates=Aysen_dates,h=2,nGS=23, xlab="DOY",
#' ylab="EVI", rge=c(0,10000))
#'}
#' @export

PhenKplot <-
  function(x,dates,h,nGS,xlab,ylab,rge) {

    # a.Preparing dataset
    if(length(rge)!=2){stop("rge must be a vector of length 2")}
    if(rge[1]>rge[2]){stop("rge vector order must be minimum/maximum")}
    if(length(dates)!=length(x)){stop("N of dates and files do not match")}
    if (all(is.na(x))) {
      return(rep(NA,nGS))
    }
    DOY <- yday(dates)
    D1<-cbind(DOY,x)
    if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*length(D1))) {
      return(rep(NA,nGS))
    }

    # b. Kernel calculation using the reference period (D1)
    if(h!=1 && h!=2){stop("Invalid h")}
    DOGS<-cbind(seq(1,365),c(seq(185,365),seq(1,184)))
    if(h==2){
      for(i in 1:nrow(D1)){
        D1[i,1]<-DOGS[which(DOGS[,1]==D1[i,1],arr.ind=TRUE),2]}}

    Hmat<-Hpi(na.omit(D1))
    Hmat[1,2]<-Hmat[2,1]
    K1<-kde(na.omit(D1),H=Hmat,xmin=c(1,rge[1]),xmax=c(365,rge[2]),gridsize=c(365,500))
    K1Con<-K1$estimate
    for(j in 1:365){
      K1Con[j,]<-K1$estimate[j,]/sum(K1$estimate[j,])}
    MAXY<-apply(K1Con,1,max)
    for(i in 1:365){
      MAXY[i]<-median(K1$eval.points[[2]][which(K1Con[i,]==MAXY[i],arr.ind=TRUE)])}

    # Plot
    h2d <- list()
    h2d$x <- seq(1,365)
    h2d$y <- seq(rge[1],rge[2],len=500)
    h2d$density <- K1Con/sum(K1Con)
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumProbs <- cumsum(uniqueVals)
    names(cumProbs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]

    image(h2d$x, h2d$y, h2d$cumDensity, xlab=xlab, ylab=ylab, font.lab=2, breaks = c(0,0.5, 0.75, 0.9,0.95),col=heat.colors(n=4,alpha=0.6))
    contour(h2d$x, h2d$y, h2d$cumDensity, levels = c(0,0.5, 0.75, 0.9,0.95),add=T,col=grey(0.25),labcex=1)
    lines(seq(1,365),MAXY,lwd=3,col='dark red')
}
