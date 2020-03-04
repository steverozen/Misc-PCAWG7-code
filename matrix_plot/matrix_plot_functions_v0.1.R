###########################################################################
###########################################################################
# Matrix plot functions
# 
# v0.1
# 
# An alpha version
# 
# Copyright (C) 2018 Steven G. Rozen and Mi Ni Huang
# 
# The code is released under GPL-3
# https://www.gnu.org/licenses/gpl-3.0.en.html
# you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# Contact: steverozen@gmail.com
###########################################################################
###########################################################################

# Dependencies
library(plotrix) ## for color legend
library(RColorBrewer)
library(shape)

### define the color scheme
colPal <- colorRampPalette(c('#FFC651', "#FF6B6B", "#6767F7", "#003D7C"))
col.ref <- colPal(10)

### match col to mutation rates
rate.to.col <- function(mut.rate){
  ### get the index of the color 
  idx <- max(which(mut.rate.ref < mut.rate))
  color <- col.ref[idx]
  color
}

### define size and percentage breaks
percent.ref <- seq(0, 0.9, 0.1)
size.ref <- percent.ref*1.5 + 0.5

### match percentage to sizes
percent.to.size <- function(percent){
  idx <- max(which(percent.ref < percent))
  size <- size.ref[idx]
  size
}

## read contribution matrix function
## exp.data: number of muts (samples|cohort * sigs)
## returns -- 
## mat1: proportion of samples w/ sigs  
## mat2: average (proportion of) mutations due to sigs 
## mat3: average mutations due to sigs 
## mat4: number of samples w/ sigs
## mat5: median number of mutations due to sigs 
exp.data.to.mats <- function(expr.data){
  Cancer.Type <- gsub(':.*', '', expr.data$Sample.Names)
  cancer.names <- unique(Cancer.Type)
  sig.names <- names(expr.data[,3:ncol(expr.data)])
  
  ## initiate five matrixes
  mat1 <- matrix(0, length(sig.names), length(cancer.names))
  rownames(mat1) <- sig.names
  colnames(mat1) <- cancer.names
  mat2 <- mat1
  mat3 <- mat1
  mat4 <- mat1
  mat5 <- mat1
  ## vec1: number of samples in each cancer type
  vec1 <- table(Cancer.Type)
  
  for (cancer in cancer.names){
    ## mat4: number of samples w/ sigs
    n.pos.samples <- apply(expr.data[Cancer.Type==cancer, 3:ncol(expr.data)], 
                           2, function(x){sum(x!=0)})
    mat4[, cancer] <- n.pos.samples[match(rownames(mat1), names(n.pos.samples))]
    
    ## mat1: proportion of samples w/ sigs
    n.samples <- sum(Cancer.Type==cancer)
    percent.samples <- n.pos.samples/n.samples
    mat1[, cancer] <- percent.samples[match(rownames(mat1), names(percent.samples))]

    ## mat3: average mutations due to sigs
    total.mut <- colSums(expr.data[Cancer.Type==cancer, 3:ncol(expr.data)])
    per.sample.mut <- total.mut/n.pos.samples
    per.sample.mut[is.na(per.sample.mut)] <- 0
    mat3[, cancer] <- per.sample.mut[match(rownames(mat3), names(per.sample.mut))]
    
    ## mat2: average proportion of mutations due to sigs
    ## mat5: median number of mutations due to sigs (except 0s)
    for (sig in sig.names){
      percent.mut <- ifelse(total.mut[sig] == 0, 0, 
                            total.mut[sig]/sum(expr.data[Cancer.Type==cancer&expr.data[,sig]!=0, 3:ncol(expr.data)]))
      mat2[sig, cancer] <- percent.mut
      median.mut <- ifelse(total.mut[sig] == 0, 0, 
                           median(expr.data[Cancer.Type==cancer&expr.data[,sig]!=0, sig]))
      mat5[sig, cancer] <- median.mut
    }
  }
  
  output = list(mat1, mat2, mat3, mat4, mat5, vec1)
  names(output) <- c('sample.proportion.mat', 'mut.proportion.mat', 
                     'mut.load.mat', 'sample.n.mat', 'mut.median.load.mat', 'n.samples')
  output
}

## cluster cancer types
## mat1: average (median) mutations due to sigs ／ proportion of samples have the sig
## binary: present or not (else use the numbers)
## n: number of clusters
cluster.by.types<-function(mat1, binary=FALSE, n=10){
  if (binary){
    mat1[mat1!=0] <- 1
  }
  #distance <- as.dist(1-cor(mat1, method="pearson"))
  distance <- dist(t(mat1), method = "euclidean")
  hc <- hclust(distance, "complete")
  #plot(hc, hang = -1, xlab='', sub='', yaxt='n', ylab='')
  #rect.hclust(hc, k=n, border="red")
  groups <- cutree(hc, k=n)
  groups <- groups[hc$order]
  groups <- table(groups)[unique(groups)]
  output <- list(order=hc$order, groups=groups)
  output
}

### keep > 0.90 accuracy samples with certain number of mutations cut-off
filter.samples <- function(indata, cutoff=100){
  indata <- indata[indata$Accuracy>0.9, ]
  indata <- indata[rowSums(indata[,4:ncol(indata)])>=cutoff, ]
  sample.ids <- paste(indata$Cancer.Type, indata$Sample.Names, sep=':')
  outdata <- indata[,2:ncol(indata)]
  outdata$Sample.Names <- sample.ids
  outdata
}

## Green matrix plot function
## mat1: proportion of samples w/ sigs --- size 
## mat2: average mutations due to sigs --- color
## mat1 and mat2 should have the same 0 elements
## vec1: number of samples in each cancer type
## vec2: etiology of each sig
matrix.plot<-function(mat1, mat2, vec1, vec2=NULL, mut.rate.ref, clustering=FALSE, genome.length=3200, 
                      cancer.names=NULL, no.names=FALSE, size.legend=TRUE, color.legend='bottom'){
  ### sort two mats by clustering method
  if (clustering) {
    clusters <- cluster.by.types(mat1, n=9)
    mat1 <- mat1[,clusters$order]
    mat2 <- mat2[,clusters$order]
  }
  
  n.sigs <- nrow(mat1)
  n.cancers <- ncol(mat1)
  sig.names <- rownames(mat1)
  cancer.types <- colnames(mat1)
  
  if (!is.null(cancer.names)){
    cancer.types <- cancer.names
    n.cancers <- length(cancer.names)
  }
  
  ## sort cancer type names
  cancer.types <- sort(cancer.types)
  
  ## match number of samples
  n.samples <- vec1[match(cancer.types, names(vec1))]
  n.samples[is.na(n.samples)] <- 0
  ## sort to mats by cancer type names
  mat1 <- mat1[,match(cancer.types, colnames(mat1))]
  mat2 <- mat2[,match(cancer.types, colnames(mat2))]
  mat1[is.na(mat1)] <- 0
  mat2[is.na(mat2)] <- 0
  
  x <- rep(1:n.cancers, each=n.sigs)
  y <- rep(n.sigs:1, n.cancers)
 
  percents <- as.vector(mat1)
  mut.load <- as.vector(mat2)
  
  ## plot background color --- add two lines on top for numbers
  x.plot.region <- rep(1:n.cancers, each=n.sigs+2)
  y.plot.region <- rep((n.sigs+2):1, n.cancers)  
  plot(x.plot.region, y.plot.region, 
       type = "n", axes = FALSE, xlab='', ylab='') ## no axes
  ## bg scheme: vertical stripes
  for (i in 1:length(x)){
    if (i %% n.sigs == 0) {
      bg.col <- ifelse((i/n.sigs) %% 2 == 0, 'grey90', 'grey95')
    } else {
      bg.col <- ifelse(floor(i/n.sigs) %% 2 == 0, 'grey95', 'grey90')
    }
    rect(x[i]-0.5, y[i]-0.5, x[i]+0.5, y[i]+0.5, border = "white", col = bg.col)
  }
  
  ## drop 0s in the positions
  x <- x[percents != 0]
  y <- y[percents != 0]
  percents <- percents[percents != 0]
  mut.load <- mut.load[mut.load != 0]
  
  ## use intervals for dot size: 0-0.1, 0.1-0.2 ...
  # sizes <- sapply(percents, percent.to.size)
  ## use continous numbers for dot size 
  sizes <- percents*1.5 + 0.5
  ## mutation load to mutation per MB
  mut.rate <- mut.load/genome.length
  colors <-sapply(mut.rate, rate.to.col)

  ## plot points 
  points(x, y, pch=19, cex=sizes, col=colors)
  
  ## add cancer names
  if (no.names==FALSE){
  text(x=1:n.cancers - 0.1, y=n.sigs+2, labels=cancer.types, 
       srt=45, cex=0.75, xpd=T, adj=0)
  }
  ## add total number of tumors
  for (i in 1:n.cancers){
    roundrect(c(i, n.sigs+1), radx=0.2, rady=0.45, rx=0.25, dr = 0.01, 
              col = "#CCE5FF", lcol = NA, lwd = 1, angle = 0, xpd=T)
  }
  text(x=1:n.cancers, y=n.sigs+1, labels=n.samples, srt=0, cex=0.55, xpd=T, adj=0.5)
  
  ## adjust sig name size
  # sig.cex <- ifelse(n.sigs > 30, 0.8, 1)
  sig.cex <- 0.75
  ## add sig names 
  print.sig.names <- gsub('Signature\\.', '', sig.names)
  print.sig.names <- gsub('Subs', 'SBS', print.sig.names)
  print.sig.names <- gsub('Indels', 'ID', print.sig.names)
  print.sig.names <- gsub('Dinucs', 'DBS', print.sig.names)
  print.sig.names <- gsub('\\.0', '', print.sig.names)
  print.sig.names <- gsub('\\.', '', print.sig.names)
  text(x=-0.01, y=n.sigs:1, labels=print.sig.names, srt=0, cex=sig.cex, xpd=T, adj=1)

  ## add color legend
  n.col <- length(col.ref)
  if (color.legend=='bottom'){
    x.col <- seq(n.cancers+1.5,n.cancers+n.col*2+0.5,2)
    y.col <- rep(-1, 10)
    color.legend(x.col[1]-1, -0.5, x.col[n.col]+1, -1.5, c('',''), col.ref, 
                 align="rb", gradient="x")
    text(x.col-1, y.col-1, labels=names(mut.rate.ref), cex=0.85, xpd=T)
    text(x=(x.col[1]+x.col[n.col])/2, y=-3.5, 
         labels='Median mutations/Mb due to signature\n(among tumors with the signature)', 
         cex=0.8, xpd=T)
  } else if (color.legend=='right'){
    x.col <- rep(n.cancers+21, 10)
    y.col <- seq(2, 11, 1)
    color.legend(x.col[1]+0.5, y.col[1]-0.5, x.col[n.col]+1.5, y.col[n.col]+0.5, c('',''), col.ref, 
                 align="rb", gradient="y")
    text(x.col+1.6, y.col-0.5, labels=names(mut.rate.ref), cex=0.75, xpd=T, adj=0)
    text(x=x.col[1]+4, y=(y.col[1]+y.col[n.col])/2, 
         labels='Median mutations/Mb due to signature\n(among tumors with the signature)', 
         cex=0.8, xpd=T, srt=90)
  }
  
  if (size.legend==TRUE){
    ## add size legend
    x.size <- 1:10
    y.size <- rep(-1, 10)
    ## size legend background
    for (i in 1:length(x.size)){
      bg.col <- ifelse(i %% 2 == 0, 'grey95', 'grey90')
      rect(x.size[i]-0.5, y.size[i]-0.5, x.size[i]+0.5, y.size[i]+0.5, 
           border = "white", col = bg.col, xpd=T)
    }
    points(x.size, y.size, pch=19, cex=size.ref, col='grey60', xpd=T)
    n.size <- 2
    x.pos <- seq(x.size[1]-0.5, x.size[length(x.size)]+0.5, 10/n.size)
    text(x.pos, rep(y.size[1]-1, n.size+1), 
         labels=seq(0,1,1/n.size), cex=0.85, xpd=T)
    text(x=mean(c(x.size[1], x.size[length(x.size)])), y=-3.5, 
         labels='Proportion of tumors\nwith the signature', cex=0.8, xpd=T)
    
    ## add total number of tumors legend
    roundrect(c(15+2, -1), radx=0.2, rady=0.45, rx=0.25, dr = 0.01, 
              col = "#CCE5FF", lcol = NA, lwd = 1, angle = 0, xpd=T)
    text(x=15+2, y=-1, labels='n', cex=0.75, font=4, xpd=T)
    text(x=15+3, y=-1, labels='Number of tumors', cex=0.8, xpd=T, srt=0, adj=0)
  }
  
  ## add etiology
  if (!is.null(vec2)){
    ## match etiology of sigs
    etiology <- vec2[match(sig.names, names(vec2))]
    ## add background color
    for (i in 1:n.sigs){
      # bg.col <- ifelse(i %% 2 == 0, 'grey95', 'grey90')
      bg.col <- 'grey95'
      len.etiology <- 20
      rect(n.cancers+0.5, n.sigs+1-i-0.5, n.cancers+0.5+len.etiology, 
           n.sigs+1-i+0.5, border = "white", col = bg.col, xpd=T)
    }
    ## add text
    text(x=n.cancers+0.5+len.etiology/2, y=(n.sigs+1):1, 
         labels=c('Proposed etiology', etiology), cex=0.70, xpd=T)
  }
  
  ## add clustering lines
  if (clustering) {
    x.line = 0.5
    for (i in clusters$groups){
      x.line <- x.line + i
      ## subs
      #segments(x.line-0.1, n.sigs+2, x.line+2.25-0.1, n.sigs+2+3, lwd = 1.5, xpd=T)
      ## indels
      #segments(x.line-0.1, n.sigs+2, x.line+2.55-0.1, n.sigs+2+2.8, lwd = 1.5, xpd=T)
      ## dinucs
      segments(x.line-0.1, n.sigs+2, x.line+2.5-0.1, n.sigs+2+2.8, lwd = 1.5, xpd=T)
    }
  }
}

