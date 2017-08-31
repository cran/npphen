#' @title PhenMap
#' @description Estimates annual Land Surface Phenology (LSP) using time series of  a vegetation greenness raster stack.
#' @param s	Raster stack with greenness (e,g. NDVI or EV) values
#' @param dates	Vector with dates at which the greenness values were recorded
#' @param h	Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h=1 if the vegetation is in the Northern Hemisphere (season starting at January 1st), h=2 if it is in the Southern Hemisphere (season starting at July 1st)
#' @param nGS Numeric. Number of greenness values within a single growing season. For example, nGS=23 for MODIS	Vegetation Index 16-days composites
#' @param nCluster	Numeric. Number of CPU cores to be used for computational calculations
#' @param outname	Character vector with the output path and filename with extension or only the filename and extension if work directory was set. For example outname="output_phen.tif". See \code{\link{writeRaster}}
#' @param format	Character. Output file type. See \code{\link{writeFormats}}
#' @param datatype	Character. Output data type. See \code{\link{dataType}}
#' @param rge  A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @details Derives the annual Land Surface Phenological (LSP) cycle for a standard growing season using a raster stack of satellite based greenness values such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). The LSP cycle is calculated for all pixels of the input raster stack in the same way as for the Phen function. The output is a multiband raster where every band is the expected greeness value at a given time step of the standard growing season. For example, for MODIS Vegetation Index 16-days composites the number of time steps of the growing season (nGS) is 23 , and therefore, the output raster will have 23 bands. A vector with dates for the greenness values is also required.
#' @return RasterStack
#' @seealso \code{\link{Phen}}
#' @examples
#' \donttest{
#' ##DEPENDING ON HARDWARE, THIS PROCESS CAN BE HIGHLY TIME CONSUMING##
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
#' # Making the LSP raster, n bands = 23
#'
#' library(snow)
#'
#' # Define the number of cores to be use. In this example we use 1
#' nc1<-1
#'
#' PhenMap(s=Slovenia_rasters,dates=Slovenia_dates,h=1,nGS=23, nCluster=nc1,
#' outname="phen_slov.tif", format="GTiff", datatype="FLT4S",rge=c(0,10000))
#' map1<-raster("phen_slov.tif")
#' plot(map1)
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
#' # Making the LSP raster, n bands = 23
#' # Define the number of cores to be use. In this example we use 1
#' nc1<-1
#'
#' PhenMap(s= Aysen_rasters,dates=Aysen_dates,h=2,nGS=23, nCluster=nc1,
#' outname="phen_aysen.tif", format="GTiff", datatype="FLT4S",rge=c(0,10000))
#' map2<-raster("phen_aysen.tif")
#' plot(map2)
#' }
#' @export

PhenMap <-
  function(s,dates,h,nGS,nCluster,outname,format,datatype,rge) {
    ff <- function(x) {

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
      #plot(Ref[first.DOGS:last.DOGS],xlab='Time',ylab='Anomaly',font.lab=2,type='l')
      Ref[first.DOGS:last.DOGS]
    }

    #----------------------------------------------------------------------------------------
    # Calculating the annual phenological curve using n clusters
    beginCluster(n=nCluster) # write 'beginCluster(n=3)' for using e.g. 3 cores, default uses all available cores)
    dates<<-dates
    clusterR(x=s,calc, args=list(ff),export=c('dates'),filename=outname,format=format,datatype=datatype,overwrite=T)
    endCluster()
}
