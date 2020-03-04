###########################################################################
###########################################################################
# Gaddy gram plot functions
# 
# v0.1
# 
# An alpha version; code used to plot the "snake plots" a.k.a.
# "Gaddy grams" as seen for example at 
# https://cancer.sanger.ac.uk/cosmic/signatures/SBS/SBS1.tt
# 
# Copyright (C) 2018-2020 Steven G. Rozen and Mi Ni Huang
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


merge.cancers <- function(expr.data){
  sample.names <- expr.data$Sample.Names
  sample.names <- gsub('Bone-Epith', 'Bone-Other', sample.names)
  sample.names <- gsub('Bone-Benign', 'Bone-Other', sample.names)
  sample.names <- gsub('Breast.*:', 'Breast:', sample.names)
  sample.names <- gsub('Cervix.*:', 'Cervix:', sample.names)
  sample.names <- gsub('Myeloid-MDS', 'Myeloid-MDS/MPN', sample.names)
  sample.names <- gsub('Myeloid-MPN', 'Myeloid-MDS/MPN', sample.names)
  sample.names <- gsub('AdenoCa', 'AdenoCA', sample.names)
  sample.names <- gsub('^AML:', 'Myeloid-AML:', sample.names)
  sample.names <- gsub('Lung-Small', 'Lung-SCC', sample.names)
  sample.names <- gsub('Transitional-cell-carcinoma', 'Bladder-TCC', sample.names)
  expr.data$Sample.Names <- sample.names
  expr.data
}


gaddy.gram <- function(sig.data, cancer.type, sig.name, n=NULL){
  ## get total counts
  counts <- table(cancer.type)
  # ## remove n <= 10
  ## for indels and dinus: do not use the 10 cut-off
  # counts <- counts[counts>10]
  
  ## remove 0s
  sig.Cancer.Type <- cancer.type[sig.data!=0]
  sig.data <- sig.data[sig.data!=0]
  
  cancer.names <- unique(sig.Cancer.Type)
  cancer.names <- cancer.names[cancer.names%in%names(counts)]
  # ## remove n <= 10
  # cancer.type.to.remove <- names(table(sig.Cancer.Type))[table(sig.Cancer.Type) <= 10]
  # cancer.names <- cancer.names[!cancer.names%in%cancer.type.to.remove]
  n.cancers <- length(cancer.names)
  
  ### show all types ...
  if (is.null(n)){
    n.bg <- length(counts)
  } else {
    n.bg <- n
  }
  ymax <- ceiling(max(sig.data)/100)*100
  ymin <- floor(min(sig.data)*1000000)/1000000
  ## plot background color
  x <- 1:n.bg
  y <- c(log10(ymin/1.5), log10(ymax))
  y <- c(floor(y[1]), ceiling(y[2]))
  x.plot.region <- rep(x, each=2)
  y.plot.region <- rep(y, n.bg)  
  plot(x.plot.region, y.plot.region, 
       type = "n", axes = FALSE, xlab='', ylab='') ## no axes
  ## bg scheme: vertical stripes
  for (i in 1:n.cancers){
    bg.col <- ifelse(i %% 2 == 0, 'grey90', 'grey95')
    rect(x[i]-0.5, y[1], x[i]+0.5, y[2], border = "white", col = bg.col, xpd=T)
  }
  ## y axis label
  axis.pos <- seq(y[1], y[2], by = 1)
  axis.label <- 10^axis.pos
  axis(side = 2, at = axis.pos, labels = axis.label, las=2, pos = 0.5, cex.axis=0.8)
  # mtext("# mutations / Mb", side=2, line = 1.5)
  y.lab.pos <- 2.1
  text(-y.lab.pos, (y[1]+y[2])/2, labels="# mutations / Mb", srt=90, cex=1, xpd=T)
  
  ## sort cancers by medians 
  med <- rep(0, n.cancers)
  names(med) <- cancer.names
  for (cancer in cancer.names){
    med[cancer] <- median(sig.data[sig.Cancer.Type == cancer])
  }
  med[is.na(med)] <- 0
  med <- sort(med)
  counts <- counts[match(names(med), names(counts))]
  
  present.counts <- counts
  ## per cancer plot
  for (i in 1:n.cancers){
    cancer <- names(med)[i]
    present.counts[cancer] <- sum(sig.Cancer.Type==cancer)
    if (sum(sig.Cancer.Type==cancer)==0) next
    sig.cancer.data <- sort(sig.data[sig.Cancer.Type==cancer])
    sig.cancer.data <- sig.cancer.data[sig.cancer.data!=0]
    seg.len <- 0.2
    if (length(sig.cancer.data) == 1) x.plot <- x[i]
    else x.plot <- seq(x[i]-seg.len,x[i]+seg.len,seg.len*2/(length(sig.cancer.data)-1))
    y.plot <- log10(sig.cancer.data)
    points(x.plot, y.plot, pch=16, cex=0.5, xpd=T)
    med.plot <- log10(med[i])
    segments(x[i]-seg.len, med.plot, x[i]+seg.len, med.plot, lwd=2, col='Red')
  }
  
  ## add cancer names
  text(x=1:n.cancers, y=y[2], labels=names(med), 
       srt=90, cex=0.8, xpd=T, adj=0)
  
  ## add sample numbers
  y.len <- (y[2] - y[1])/25
  text(x=1:n.cancers, y=y[1] - y.len, labels=present.counts, cex=0.62, xpd=T)
  text(x=1:n.cancers, y=y[1] - y.len*4, labels=counts, cex=0.62, xpd=T, adj=0.5)
  ## for labels
  for (i in 1:n.cancers){
    seg.len <- 0.4
    segments(x[i] - seg.len, y[1] - y.len*2.5, x[i] + seg.len, y[1] - y.len*2.5, lwd=1, col='blue', xpd=T)
  }
}

get.width <- function(sig.data, cancer.type){
    ## get total counts
    counts <- table(cancer.type)
    # ## remove n <= 10
    ## for indels and dinus: do not use the 10 cut-off
    # counts <- counts[counts>10]
    
    ## remove 0s
    sig.Cancer.Type <- cancer.type[sig.data!=0]
    sig.data <- sig.data[sig.data!=0]
    
    cancer.names <- unique(sig.Cancer.Type)
    cancer.names <- cancer.names[cancer.names%in%names(counts)]
    n.cancers <- length(cancer.names)
    a <- 0.82 ## left margin
    b <- 0.3 ## right margin
    w <- (11.2677 - a - b)/42 ## per cancer length in inch 
    if (n.cancers > 1) {
      width.page <- n.cancers*w + a + b 
      width.left <- a + w*(0.5-0.04*(n.cancers-1))
      width.right <- b + w*(0.5-0.04*(n.cancers-1))
    }
    else {
      width.page <- 1*0.864 + a + b
      width.left <- a + 1*0.432
      width.right <- b + 1*0.432
    } 
    return(c(width.page, width.left, width.right))
}



## use assignment from combined file?? already # mutations/Mb
exp.file = '~/Desktop/PCAWG7-calls/2018-03-13-ludmil-assignments/combined_green_matrix.csv'

### reading data
data <- read.csv(exp.file, stringsAsFactors = F)
data <- merge.cancers(data)
Cancer.Type <- gsub(':.*', '', data$Sample.Names)
## filter rare cancer types
filtered.cancers <- names(table(Cancer.Type))[table(Cancer.Type) < 10]
data <- data[!Cancer.Type%in%filtered.cancers, ]
Cancer.Type <- gsub(':.*', '', data$Sample.Names)
sig.names <- names(data[,3:ncol(data)])

### per sig plot
for (sig in sig.names){
  print.sig.name <- gsub('Signature\\.', '', sig)
  print.sig.name <- gsub('Subs', 'SBS', print.sig.name)
  print.sig.name <- gsub('\\.0', '', print.sig.name)
  print.sig.name <- gsub('\\.', '', print.sig.name)
  output.name <- paste0('~/Desktop/PCAWG7-calls/2018-03-13-ludmil-assignments/gaddy-gram/SBS/',print.sig.name,'.png')
  sig.data <- data[, sig]
  if (sum(sig.data!=0) == 0) next
  # width.pdf <- get.width(sig.data, Cancer.Type)
  # cairo_pdf(output.name, width=width.pdf[1], height=3.1)
  png(output.name, width=11.2677, height=3.1, units = 'in', res = 300)
  par(mai=c(0.32,0.82,1.16,0.3), omi=c(0,0,0,0))
  gaddy.gram(sig.data, Cancer.Type, sig)
  dev.off()
}



## INDELs
exp.file1 = '~/Desktop/PCAWG7-calls/2018-03-13-ludmil-assignments/PCAWG_sigProfiler_ID_signatures_in_samples_filtered.csv'
exp.file2 = '~/Desktop/PCAWG7-calls/2018-03-13-ludmil-assignments/TCGA_WES_sigProfiler_ID_signatures_in_samples_filtered.csv'

### reading data
data <- read.csv(exp.file1, stringsAsFactors = F)
data <- merge.cancers(data)
data.E <- read.csv(exp.file2, stringsAsFactors = F)
data.E <- merge.cancers(data.E)
## calculate mut per Mb
data <- cbind(data[,1:2], data[,3:ncol(data)]/3200)
data.E <- cbind(data.E[,1:2], data.E[,3:ncol(data)]/50)

data <- rbind(data, data.E)
Cancer.Type <- gsub(':.*', '', data$Sample.Names)
sig.names <- names(data[,3:ncol(data)])

### per sig plot
for (sig in sig.names){
  print.sig.name <- gsub('Signature\\.', '', sig)
  print.sig.name <- gsub('Indels', 'ID', print.sig.name)
  print.sig.name <- gsub('\\.0', '', print.sig.name)
  print.sig.name <- gsub('\\.', '', print.sig.name)
  output.name <- paste0('~/Desktop/PCAWG7-calls/2018-03-13-ludmil-assignments/gaddy-gram/ID/',print.sig.name,'.png')
  sig.data <- data[, sig]
  if (sum(sig.data!=0) == 0) next
  width.pdf <- get.width(sig.data, Cancer.Type)
  # cairo_pdf(output.name, width=width.pdf[1], height=3.1)
  png(output.name, width=11.2677, height=3.1, units = 'in', res = 300)
  par(mai=c(0.32,0.82,1.16,0.3), omi=c(0,0,0,0))
  gaddy.gram(sig.data, Cancer.Type, sig)
  dev.off()
}






#############################
## DBS ### no change yet....
exp.file = '~/Desktop/PCAWG7-calls/2018-01-11-ludmil-assignments/All_PCAWG_WGS_DINUCS_signatures_in_samples_filtered.csv'
### reading data
data <- read.csv(exp.file, stringsAsFactors = F)
data <- merge.cancers(data)
Cancer.Type <- gsub(':.*', '', data$Sample.Names)
sig.names <- names(data[,3:ncol(data)])

### per sig plot
for (sig in sig.names){
  print.sig.name <- gsub('Signature\\.', '', sig)
  print.sig.name <- gsub('Dinucs', 'DBS', print.sig.name)
  print.sig.name <- gsub('\\.0', '', print.sig.name)
  print.sig.name <- gsub('\\.', '', print.sig.name)
  output.name <- paste0('~/Desktop/PCAWG7-calls/2018-01-11-ludmil-assignments/gaddy-gram/DBS/',print.sig.name,'.png')
  sig.data <- data[, sig]
  sig.data <- sig.data/3200 ## for indels and dinucs....
  if (sum(sig.data!=0) == 0) next
  width.pdf <- get.width(sig.data, Cancer.Type)
  # cairo_pdf(output.name, width=width.pdf[1], height=3.1)
  png(output.name, width.pdf[1], height=3.1, units = 'in', res = 300)
  par(mai=c(0.32,width.pdf[2],1.16,width.pdf[3]), omi=c(0,0,0,0))
  gaddy.gram(sig.data, Cancer.Type, sig)
  dev.off()
}


