#' Binning beacon or beagle data
#'
#'Binning beacon/ beagle data (or any csv in format: tagX, tagY, time) and applying by default a default 50% minimum signal per bin on the data to count the bin as signal.
#' @param dataCSV csv in format: tagX, tagY, time
#' @param binsize can be anything >1
#' @return returns the data as X, Y, bin indicator (e.g., bin indicator = 1 for all connections made during bin 1)
#' @export
#'
binningBeaconData<- function(dataCSV,binsize) {
  dat<-utils::read.csv(dataCSV)

  #only keep unique rows (there are repeated rows with same ID1, ID2, time information)
  dat <- dat[!duplicated( dat [ ,c(1,2,4)]),] #unique (ignoring the power variable)
  colnames(dat)<- c("tagX", "tagY", "power", "time")

  #add timestamp info:
  dat<-cbind(dat,rep(1:length(tapply (dat$tagX,dat$time,length)), tapply (dat$tagX,dat$time,length)))
  colnames(dat)<- c("tagX", "tagY", "power", "time", "timeStamp")


  groupsize<-length(unique(c(dat[,1],dat[,2])))

  #get unique interactions per timestamp (1:timestamps)
  reX<-plyr::mapvalues(dat$tagX, from = sort(unique(c(dat$tagY,dat$tagX))),
                 to = c(1:groupsize))
  reY<-plyr::mapvalues(dat$tagY, from = sort(unique(c(dat$tagY,dat$tagX))),
                 to = c(1:groupsize))

  dat<- unique(cbind(dat,reX,reY))

  timestamps<- dat$timeStamp [length(dat$timeStamp)]

  #bin data
  indicatorHelp<- rep(c(1:ceiling(timestamps/binsize)), each=binsize)

  binIndicator<-plyr::mapvalues(dat$timeStamp, from = c(1:timestamps), to = indicatorHelp[1:timestamps])


  datBinned <-cbind(dat,binIndicator)

  #clean data to get unique interactions
  uniquecounts<-plyr::count(datBinned,c("reX","reY","binIndicator"))

  #here: keep only those that have the following combination:
  #smoothing at binsize/2 i.e. at least binsize/2 times a signal per bin:
  keep<- uniquecounts[which(uniquecounts$freq > floor(binsize/2)), c(1,2,3)]
  data_p <- do.call("paste", datBinned[,c(6,7,8)])
  keep_p <- do.call("paste", keep)
  datBinned<- cbind(datBinned,data_p)
  out<-datBinned[data_p %in% keep_p, ]

  out<- out[!duplicated(out[ ,9]),]

  return(out)
}










