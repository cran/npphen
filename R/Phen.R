#' @title Phen
#' @description Estimates the annual phenological cycle from a time series of vegetation greenness.
#' @param x	Numeric vector with greenness values
#' @param dates	Vector with dates at which the greenness values were recorded
#' @param h	Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h=1 if the vegetation is in the Northern Hemisphere (season starting at January 1st), h=2 if it is in the Southern Hemisphere (season starting at July 1st)
#' @param nGS Numeric. Number of greenness values within a single growing season. For example, nGS=23 for MODIS	Vegetation Index 16-days composites
#' @param rge A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @details Derives the annual phenological cycle for a standard growing season using a numeric vector of vegetation canopy greenness values (e.g. Leaf Area Index, LAI) or satellite based greenness proxies such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). A vector with dates for the greenness values is also required.
#' @return A numeric vector of length = nGS, where each value represents the expected greeness at that date
#' @seealso \code{\link{PhenMap}}
#' @examples
#' \dontshow{
#' ## Testing function with time series of Slovenian data (EVI)
#' # Load data
#' phents<-read.table(system.file("extdata/date_tables/datats",package="npphen"),
#' dec='.',sep='\t',header=TRUE)
#' # Phenology for the given data
#' Phen(x=as.vector(phents$x),dates=phents$dates,h=1,nGS=23,rge=c(0,10000))
#' }
#' \donttest{
#' library(rts)
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
#' sl_ts<-rts(Slovenia_rasters,Slovenia_dates)
#'
#' # Obtain data from a particular pixel generating a time series
#' sl_pixel<-cellFromXY(sl_ts,c(474368,5096979))
#' sl_pixelts<-extract(sl_ts,sl_pixel)
#' plot(sl_pixelts)
#'
#' # Phenology for the given pixel
#' Phen(x=as.vector(sl_pixelts),dates=Slovenia_dates,h=1,nGS=23,rge=c(0,10000))
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
#' ay_ts<-rts(Aysen_rasters,Aysen_dates)
#'
#' # Obtain data from a particular pixel generating a time series
#' ay_pixel<-cellFromXY(ay_ts,c(228373,4806975))
#' ay_pixelts<-extract(ay_ts,ay_pixel)
#' plot(ay_pixelts)
#'
#' # Phenology for the given pixel
#' Phen(x=as.vector(ay_pixelts),dates=Aysen_dates,h=2,nGS=23,rge=c(0,10000))
#' }
#' @import raster
#' @import ks
#' @import grDevices
#' @import methods
#' @import rgdal
#' @import snow
#' @importFrom lubridate yday
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom graphics abline
#' @export

Phen <-
  function(x,dates,h,nGS, rge) {

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

    # c. Getting the most probable values for a standard year (nGS time steps)
    first.DOGS <- which.min(D1[,1])
    first.DOGS <- first.DOGS[1]
    last.DOGS <- first.DOGS+nGS-1
    Ref<-rep(NA,nGS)
    for(i in first.DOGS:last.DOGS){
      Ref[i]<-MAXY[D1[i,1]]}
    plot(seq(1,365,round(365/nGS)),Ref[first.DOGS:last.DOGS],xlab='DGS',ylab='VI',font.lab=2,type='l') # VI=vegetation index, DGS=day of the growing season
    Ref[first.DOGS:last.DOGS]
}

