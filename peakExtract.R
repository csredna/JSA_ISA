library("pdist")

setwd("C:/Users/Minoxsta/Dropbox/ISA2016/data")

loadMA <- function(name) {
  ma <- read.table(paste("manual/", name, sep = "", collapse = NULL), sep="\t", header = TRUE)
  return(ma)
}

getPeakAnno <- function(ma, x) {
  peakAnno <- ma[x, c("measurement_name", "peak_name", "signal", "index_t", "index_r")]
  return(peakAnno)
}

loadRaw <- function(name) {
  raw <- read.csv(paste("raw/", name, sep = "", collapse = NULL))
  return(raw)
}

# Given a peak (information from annotation file), a name of a raw data file and a radius
# returns a dataframe or matrix of the square of "pixels" around the peak,
# the "Peak Window".
getPW <- function(peak, raw, rad) {
  pw <- NULL
  row <- peak[1,"index_t"]
  col <- peak[1,"index_r"]
  
  # Usually a window is copied directly from the raw file.
  if ((row-rad) >= 133 & (row+rad) <= nrow(raw) & (col-rad) >= 3 & (col+rad) <= ncol(raw)) {
    pw <- raw[(row-rad):(row+rad),(col-rad):(col+rad)]
  } else {
    # Near edges "pixels" are put in the resulting window one at a time
    pw <- matrix(0, nrow = (2*rad+1), ncol = (2*rad+1)) # Create a matrix the right size, but with all zeros
    m <- 1
    n <- 1
    for (i in (row-rad):(row+rad)) {
      for (j in (col-rad):(col+rad)) {
        if (i >= 133 & i <= nrow(raw) & j >= 3 & j <= ncol(raw)) {
          pw[m,n] <- raw[i,j]
        } 
        n <- n+1
      }
      n <- 1
      m <- m+1
    }
  }
  return(pw)
}

getPeaksAnno <- function(ma) {
  peaksAnno <- NULL
  for (i in 1:nrow(ma)) {
    peakAnno <- getPeakAnno(ma, i)
    peaksAnno <- rbind(peaksAnno, peakAnno)
  }
  return(peaksAnno)
}

# Given an annotation-file name and a peak window radius,
# returns a dataframe with a row for each peak, 
# and coloumns containing important information from annotation file, 
# and a value for each "pixel" in the peak window.
getPeaks <- function(ma, rd, peaksAnno, pwRad) {
  peaksFul <- NULL
  
  # For each peak get the important data from the annotation file,
  # and the "pixel"-values from the raw file,
  # then put it in a single row/vector,
  # and append/bind it to the full list of peaks.
  for (i in 1:nrow(ma)) {
    peakAnno <- peaksAnno[i,]
    peakDat <- getPW(peakAnno, rd, pwRad)
    peakFul <- peakAnno
    
    # Add fields for the "pixel"-values
    for (i in 1:((2*pwRad+1)^2)) {
      peakFul[,paste("P", i, sep = "")] <- 0
    }
    
    # Insert "pixel"-values into the fields
    k <- 1
    for (i in 1:nrow(peakDat)) {
      for (j in 1:ncol(peakDat)) {
        peakFul[,paste("P", k, sep = "")] <- peakDat[i,j]
        k <- k+1
      }
    }
    peaksFul <- rbind(peaksFul, peakFul)
  }
  
  return(peaksFul)
}

getFakePeaks <- function(ma, raw, rad) {
  dist <- 3*2*rad
  fakesFul <- NULL
  fakesAnno <- NULL
  peaksAnno <- getPeaksAnno(ma)
  amount <- nrow(ma)
  
  # print("part 1 worked")
  
  for (i in 1:amount) {
    fp <- peaksAnno[i,]
    space <- FALSE
    while (!space) {
      fp[, "index_t"] <- sample(133:nrow(raw), 1)
      fp[, "index_r"] <- sample(3:ncol(raw), 1)
      x <- pdist(fp[,c("index_t", "index_r")], peaksAnno[,c("index_t", "index_r")])
      y <- as.matrix(x)
      if (min(y[-ncol(y)] > dist)) {
        space <- TRUE
      }
    }
    
    fakesAnno <- rbind(fakesAnno, fp)
  }
  fakesFul <- getPeaks(ma, raw, fakesAnno, rad)
  return(fakesFul)
}

# Create file with full peaks
createTrainingSet <- function() {
  for (i in dir("manual", recursive = FALSE)) {
    print(i)
  }
}