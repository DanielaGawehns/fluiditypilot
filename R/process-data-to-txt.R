#' Binning beacon or beagle data
#'
#'Binning beacon/ beagle data (or any csv in format: tagX, tagY, time) and applying by default a default 50% minimum signal per bin on the data to count the bin as signal.
#' this function needs as input the imported csv file from the beacons.
#' it will return data that are very similar to when you are using read.responder and tempSmooth on the Extractor files.

#'Variable names in output
#' tagX       tagY power       time timeStamp reX reY binIndicator   data_p
#' power here is RSSI in the extractor file
#' time/ timeStamp correspond (in theory) with time Local on extractor file
#' binIndicator is calculated in the same way as for extractor file
#' data_p is a helper variable and can be ignored

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

##############################################################################


#' Extracting connections per-child
#'
#'beacon Data is given in overall overview, here, the information is broken into individual chunks, representing individual children and their
# connections. The output is similar to the data from the responders/ extractors.
#' @param beaconDat Takes input from binningBeaconData function
#' @return a list of connections per child
#' @export

breakBeaconData<- function(beaconDat) {

  beaconData <- beaconDat[order(beaconDat[,1]),]

  #check which Y is not in X
  blockData <- list()
  for (i in 1:length(unique(c(beaconData[,1], beaconData[,2])))) { # check all unique tags

    tagOwn<-unique(c(beaconData[,1], beaconData[,2])) [i]

    ownData<- beaconData [which(beaconData[,1] == tagOwn), ]
    tagMe <- ownData[,1]
    tagOther<- ownData[,2]
    ownDataTwo<- beaconData [which(beaconData[,2] == tagOwn), ]

    #now change order of tagX and tagY in ownDataTwo block for formatting
    tagMe <- c(ownData[,1], ownDataTwo[,2])
    tagOther<- c(ownData[,2], ownDataTwo[,1])
    ownData<- rbind(ownData,ownDataTwo)

    ownData [,1] <- tagMe
    ownData [,2] <- tagOther

    blockData [[i]] <- ownData
    #print(i)
  }
  return(blockData)
}


####################################################################
#output wanted is:
#ID I , ID II, Time Start, Time End
#all in a .txt file with tabs between values and no header/first row information

#Information can be extracted via breakBeaconData () and ContactlengthsBlocks ()

### here, we can specify a parameter to only use contacts that are above a certain min threshold:
#default is 1
#minlength<-5
#filterLengths<- lapply(lengthsStarts, function (x) {which ( x[[1]] > minlength)})

### here, we can specify that between last time end and next time start a minimum lag needs to be:
# default is 1

#############################################################################
#Helperfunction

#' Find all newstart moments
#' If tapplied, it accounts for newstart moments depending on tagY (or whichever tapply dependence)
#'
#' @param x Depends on how the function is tapplied, which ever subset of data is fed. Data needs to be formated as the output from the binningBeaconData function
#' @param minlength indicates the minimum length of an interaction
#' @param minlag indicates how much lag between interactions is needed/wanted
#'
#' @return filtered data, where minlength and minlag are applied on all connections
#' @export

newstart<- function (x, minlength,minlag){
  differences<- diff (x)

  lengths<- which(differences >1)

  #include first start:
  startsIndices<- c(1,lengths+1)

  #include last interaction length:
  lastlength<- rle(differences)$lengths [length(rle(differences)$lengths)] +1
  lengths<-c(lengths, lastlength)

  #return indices that are above minlength
  minlengthIndex<- lengths >= minlength

  #return indices that are at least minlag apart
  #minlagIndex <- get difference between starts with diff and compare to min lag
  minlagIndex<- diff(startsIndices) > minlag
  #watch out --> first entry is always one (first starting point needs to be included)
  minlagIndex<- c(TRUE, minlagIndex)

  #combine both Indices into one by adding them up and asking if combined is == 2 - if so keep entry
  minLengthLagIndex <- minlengthIndex + minlagIndex == 2

  fullReturn<- list(lengths, x [startsIndices], startsIndices)
  filteredReturn<- lapply(fullReturn, function (x) x[minLengthLagIndex])
  return(filteredReturn)
}

###############################################################################
#Helperfunction
## turns list elements that newstart returns when tapplied into useful output:

#' Extract and format start and end info
#'
#' @param listelement as from the output of newstart helper-function
#'
#' @return a cbind startinfo, endinfo combination
#' @export

turnListintoOutput <- function (listelement) {
  startInfo<- listelement [[2]]
  endInfo <- startInfo + listelement [[1]]
  ifelse(length(startInfo)>1, output<-cbind(startInfo,endInfo), output<-c(startInfo,endInfo))

  return(output)
}

###############################################################################

#turn this into the format: ID1, ID2, start, end
#go through broken file and keep track of which "tagX" / "me" are already covered

#' Turn raw data into format for txt as needed in jupyter notebook
#'
#' @param brokenfile as from breakBeaconData function
#' @param minlength indicates the minimum length of an interaction
#' @param minlag indicates how much lag between interactions is needed/wanted
#'
#' @return matrix with tagX, tagY, start, end
#' @export

getOutputStartEnds<- function (brokenfile, minlength,minlag) {
  datOutput<- matrix(ncol = 4)
  accounting<-numeric()

  for (i in 1:length(brokenfile)) {
    #remove those interactions that are already accounted for
    filterDone<- which(brokenfile[[i]]$tagY %in% accounting)
    ifelse(length(filterDone) >= 1 , dat<- brokenfile [[i]] [- c( filterDone),], dat <- brokenfile [[i]] )

    if(length(nrow(dat)) == 0) {next}

    ifelse (nrow(dat) > 1,
            #use newstart to get all info, then filter on minlength,then get start/end info and combine
            #newstart needs to be tapplied to dat$binIndicator, otherwise it just takes all/any tagY and mixes them up
            {lengthStarts<- tapply(dat$binIndicator,dat$tagY,newstart, minlength=minlength, minlag=minlag);
            #deal with list element that tapply returns
            #only keep lists with at least one entry
            keepIndex<- sapply(lengthStarts, function (x) length(x[[1]]) != 0 );
            lengthStarts <- lengthStarts[keepIndex];
            #turn every list entry into proper Output object
            output<- lapply(lengthStarts, turnListintoOutput)

            #output can be empty:
            if(all(sapply(sapply(output,nrow),length) ==0 )) {next}

            #add Xtag and Ytag info to output:
            Ytags<- as.numeric(names(output))
            newOutput<-numeric()
            for (j in 1:length(output)) {ifelse(length(nrow(output[[j]])) != 0,
                                                output [[j]] <- cbind(rep(Ytags[j], nrow(output[[j]])), output [[j]] ),
                                                output [[j]] <- c( Ytags[j], output [[j]]))
              newOutput<-rbind(newOutput,output[[j]])

            }

            ifelse(length(nrow(newOutput)) != 0,
                   newOutput <- cbind(rep(dat$tagX[1], nrow(newOutput)), newOutput),
                   newOutput <- c (dat$tagX[1],newOutput)  )

            #get info and compile info
            #startInfo<- lengthStarts [[2]]; endInfo <- startInfo + lengthStarts [[1]];
            #newOutput<-cbind(dat$tagX[lengthStarts[[3]]], dat$tagY[lengthStarts[[3]]], startInfo, endInfo)
            },

            {startInfo<- dat$binIndicator; endInfo<- startInfo +1; newOutput<- c(dat$tagX, dat$tagY, startInfo, endInfo)  })

    accounting[i]<-brokenfile [[i]]$tagX[1]
    datOutput <- rbind (datOutput, newOutput)
    print(i)
  }

  return (datOutput[-1,])
}

##############################################################################################

#### Write txt file
### Add option to change filename to the name of the school/school and assessment

#' Write txt for python script
#'
#' Exports txt file needed for the python scripts/ jupyter notebook to run
#' Keeps track of binsize, minlag, minlength parameter when naming the output files
#'
#' @param output the object that you want to store as txt
#' @param binsize size of the window or bin to be binned over
#' @param minlag indicates how much lag between interactions is needed/wanted
#' @param minlength indicates the minimum length of an interaction
#' @param filename filename that indicates which data we are working with (e.g., School5_morning or similar)
#' @param group depending on how txt are supposed to be named, this is an option to group data by location, provenance or other info
#'
#' @return .txt file in the folder ~/fluiditypilot/analysis/data/derived_data/textdata_","/ add group if you would like.
#' Here, all is written into the derived_data folder
#' @export

write.output <- function (output,binsize,minlag,minlength, filename, group) {
  # write.csv(output, file = paste0("~/fluiditypilot/analysis/data/derived_data/textdata_",group,"/",filename,".txt"), sep= "\t")
  utils::write.table(output, file =  paste0("~/fluiditypilot/analysis/data/derived_data/",
                                     filename,"_", binsize,"_",minlength,"_",minlag,"_",".txt"), sep = "\t",
              row.names = FALSE, col.names = FALSE)
}


#######################################################################

#' Apply csv to txt script on one group
#'
#' Function to collect info on interaction lengths.
#' @param binsize size of the window or bin to be binned over
#' @param minlength indicates the minimum length of an interaction
#' @param minlag indicates how much lag between interactions is needed/wanted
#' @param group depending on how txt are supposed to be named, this is an option
#' to group data by location, provenance or other info
#'
#' @return .txt files in the folder ~/fluiditypilot/analysis/data/derived_data/textdata_",group,"/ ; also returns a list with two elements
#' 1) list of interaction lengths; 2) vector with info on how long the recordings are
#'
#' @export

get.txt.files <- function (binsize, minlength, minlag,group) {

  files <- list.files(pattern="*.csv", full.names=FALSE, recursive=FALSE)
  filenames<- sapply(strsplit(files, "[.]"), function (x) x[1])

  lengthsInfo<- list()
  lengthrecord<-numeric()

  for (i in 1:length(filenames)) {
    filename <- filenames [i]
    print(filename)
    binnedData<- binningBeaconData(paste0(filename,".csv"),binsize)
    brokenData<- breakBeaconData(binnedData)
   # print(brokenData)
    outputData<-getOutputStartEnds(brokenData,minlength,minlag)
    #collect info on lengths of interactions
    lengthsInfo[[i]] <- outputData[,4] - outputData[,3]
    start<-min(c(outputData[,4], outputData[,3]))
    end<- max(c(outputData[,4], outputData[,3]))

    #length of the recording in terms of bin Indicator -- i.e. take the value *binsize to get actual length in seconds
    lengthrecord[i]<- c(end-start)

    #exports/ writes .txt files
    write.output(outputData,binsize,minlag,minlength,filename,group)
  }
  return(list(lengthsInfo,lengthrecord))
}

#################################################################################




