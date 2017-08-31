#' @title PhenAnoMap
#' @description Calculates anomalies with respect to the regular phenological cycle using time series (raster) of vegetation greenness.
#' @param s	Raster stack with greenness (e,g. NDVI or EV) values
#' @param dates	Vector with dates at which the greenness values were recorded
#' @param h	Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h=1 if the vegetation is in the Northern Hemisphere (season starting at January 1st), h=2 if it is in the Southern Hemisphere (season starting at July 1st)
#' @param refp Numeric vector with the correlative number of dates to be used as reference period. For example, refp=c(1:393) for MODIS Vegetation Index 16-days composites (18/02/2000 – 06/06/2017)
#' @param anop Numeric vector with the correlative number of dates for the period in which the anomalies will be calculated. For example refp=c(21:43) for the first complete year for MODIS Vegetation Index 16-days composites (01/01/2001 – 19/12/2001). anop y refp can be overlapped
#' @param nCluster	Numeric. Number of CPU cores to be used for computational calculations
#' @param outname	Character vector with the output path and filename with extension or only the filename and extension if work directory was set. For example outname="output_phen.tif". See \code{\link{writeRaster}}
#' @param format	Character. Output file type. See \code{\link{writeFormats}}
#' @param datatype	Character. Output data type. See \code{\link{dataType}}
#' @param rge  A vector containing minimum and maximum values of the response variable used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge =c(0,10000)
#' @details Similar to PhenAnoma, it calculates phenological anomalies but using a raster stack  instead of a numeric vector of vegetation canopy greenness values (e.g. Leaf Area Index, LAI) or satellite based greenness proxies such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). For this purpose, it divides the time series (raster stack) of vegetation greeness into 2: the reference period, from which the annual phenological cycle is calculated (same as the Phen function), and the observation period, for which we want to calculate anomalies with respect to the annual phenological cycle. Negative anomalies correspond to observed values lower than the reference and positive anomalies to values higher than the reference. It delivers a raster stack with anomalies per date.
#' @return RasterStack
#' @seealso \code{\link{PhenAnoma}}
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
#' # Calculating the anomalies for the last growing season [343:365], refp [1:342]
#' # In this case refp and anop do not overlap
#'
#' library(snow)
#'
#' # Define the number of cores to be use. In this example we use 1
#' nc1<-1
#'
#' PhenAnoMap(s=Slovenia_rasters,dates=Slovenia_dates,h=1,refp=c(1:342), anop=c(343:365),
#' nCluster=nc1,outname="ano_slov.tif", format="GTiff", datatype="FLT4S", rge=c(0,10000))
#' map_an1<-raster("ano_slov.tif")
#' plot(map_an1)
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
#' # Making the LSP raster, n bands = 23
#' # Define the number of cores to be use. In this example we use 1
#' nc1<-1
#'
#' PhenAnoMap(s=Aysen_rasters,dates=Aysen_dates,h=2,refp=c(1:354), anop=c(309:331),
#' nCluster=nc1,outname="ano_aysen.tif", format="GTiff", datatype="FLT4S",rge=c(0,10000))
#' map_an2<-raster("ano_aysen.tif")
#' plot(map_an2)
#'
#'}
#' @export

PhenAnoMap <-
  function(s,dates,h,refp,anop,nCluster,outname,format,datatype,rge) {
    ff <- function(x) {

      # a.Preparing dataset

      if(length(rge)!=2){stop("rge must be a vector of length 2")}
      if(rge[1]>rge[2]){stop("rge vector order must be minimum/maximum")}
      if(length(dates)!=length(x)){stop("N of dates and files do not match")}

      ref.min <- min(refp)
      ref.max <- max(refp)
      ano.min <- min(anop)
      ano.max <- max(anop)
      ano.len <- ano.max-ano.min+1

      if(ref.min>=ref.max | ano.min>=ano.max){stop("for refp or anop, lower value > upper value")}

      if (all(is.na(x))) {
        return(rep(NA,ano.len))
      }

      DOY <- yday(dates)
      D1<-cbind(DOY[ref.min:ref.max],x[ref.min:ref.max])
      D2<-cbind(DOY[ano.min:ano.max],x[ano.min:ano.max])

      if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*length(D1))) {
        return(rep(NA,ano.len))
      }

      if (all(is.na(D2[,2]))) {
        return(rep(NA,ano.len))
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

      # c. Anomaly calculation (D2)

      if(h==2){
        for(i in 1:nrow(D2)){
          D2[i,1]<-DOGS[which(DOGS[,1]==D2[i,1],arr.ind=TRUE),2]}}

      Anoma<-rep(NA,ano.len)
      for(i in 1:nrow(D2)){
        Anoma[i]<-D2[i,2]-MAXY[D2[i,1]]}
      plot(Anoma[1:ano.len],xlab='Time',ylab='Anomaly',font.lab=2,type='l')
      abline(h=0,col="red")
      Anoma[1:ano.len]
    }

    #----------------------------------------------------------------------------------------
    # Calculating the annual phenological curve using n clusters
    beginCluster(n=nCluster) # write 'beginCluster(n=3)' for using e.g. 3 cores, default uses all available cores)
    dates<<-dates
    clusterR(x=s,calc, args=list(ff),export=c('dates'),filename=outname,format=format,datatype=datatype,overwrite=T)
    endCluster()
}

