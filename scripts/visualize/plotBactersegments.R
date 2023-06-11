library(stringr)

plot_detected_segments <- function(paths,ylab,main) {
  plot(x=0,xlim=c(0,144551),ylim=c(1,length(paths)),type="n",xaxt="n",yaxt="n",
       xlab="Site no.",ylab="Run",mgp=c(2,1,0))
  title(main=main,line=1)
  axis(side=1,at=seq(from=10000,to=140000,by=10000),
       labels=format(seq(from=10000,to=140000,by=10000),scientific=T),
       cex.axis=0.8)
  axis(side=2,at=1:length(paths),labels=ylabs,cex.axis=0.8,las=1)
  for (i in 1:length(paths)) {
    # 35 is a magic number! It works only for the NEXUS files for set1c.
    nexusinfo_str <- readLines(paths[i])[35]
    start_end_pts <- data.frame(str_match_all(nexusinfo_str,"region=\\{(?<start>\\d+),(?<end>\\d+)\\}")[[1]])
    segments(x0=as.numeric(start_end_pts$start),y0=i,
             x1=as.numeric(start_end_pts$end),y1=i,lwd=2)
  }
}

# posterior support threshold 50%
main <- expression("Posterior support" >= "50%")
ylabs <- c("1","2","3","4","5","4a","5a","5b")
path <- "results/bhv/bacter/experimentMC3_rep%s/combined_mc3_test.set1c.summary copy.tree"
paths <- sapply(X=ylabs,FUN=function(ylab) sprintf(path,ylab))
paths[1] <- "results/bhv/bacter/experimentMC3/combined_mc3_test.set1c.summary copy.tree"

pdf(file='figures/bacter/detectedsegments_support50.pdf',width=8,height=8)
plot_detected_segments(paths,ylab,main)
dev.off()

# posterior support threshold 20%
main <- expression("Posterior support" >= "20%")
ylabs <- c("1","2","3","4","5","4a","5a","5b")
path <- "results/bhv/bacter/experimentMC3_rep%s/combined_mc3_test.set1c.summary20 copy.tree"
paths <- sapply(X=ylabs,FUN=function(ylab) sprintf(path,ylab))
paths[1] <- "results/bhv/bacter/experimentMC3/combined_mc3_test.set1c.summary20 copy.tree"

pdf(file='figures/bacter/detectedsegments_support20.pdf',width=8,height=8)
plot_detected_segments(paths,ylab,main)
dev.off()