# Extract from BIC scores and log-likelihoods for best k-reticulation network
# (k = 0,1,2, ...) from NetRAX logfile.
# path: path to logfile; run: run number (distinct start trees correspond to
# distinct runs)
retVSbic <- function(path,run) {
  # unaesthetic and hack-ish, but functional
  s <- system(paste("sed -n '/n_reticulations/,$p'",path,"| sed '$d'"),intern=T)
  df <- read.table(text=paste(s[1:length(s)]),header=T)[,1:3]
  colnames(df) <- c("n_reticulations","logl","bic")
  df$n_reticulations <- as.integer(substr(df$n_reticulations,1,
                                          nchar(df$n_reticulations)-1))
  df$logl <- as.numeric(substr(df$logl,1,nchar(df$logl)-1))
  df$bic <- as.numeric(substr(df$bic,1,nchar(df$bic)-1))
  df$run <- run
  return(df)
}

retVSbic_runs <- vector("list",15) # list of dataframes, one for each start tree
for (i in 0:14) {
  path <- sprintf("results/bhv/netrax/multistart-titanium_subruns/run_%d_logfile.txt",i)
  retVSbic_runs[[i+1]] <- retVSbic(path,i)
}

###############################################################################

pdf(file="figures/BIC_15genes_blksize10000.pdf",width=12,height=6)
# Plot BIC score against no. of reticulations for each run
par(oma=c(2,0,0,0),mar=c(3,5,3,0)+0.1)
layout(matrix(c(1,2),1,2,byrow=T))
# min(vapply(retVSbic_runs,FUN=function(x) min(x$bic),FUN.VALUE=1))
# max(vapply(retVSbic_runs,FUN=function(x) max(x$bic),FUN.VALUE=1))
with(retVSbic_runs[[1]],
     plot(n_reticulations,bic,type="l",ylab="BIC score",xlab="",col="deeppink",
          lwd=2,yaxt="n",ylim=c(517827.1,522207.5),cex.axis=1,cex.lab=1.5))
axis(2,at=seq(from=518000,to=522000,by=2000),
     labels=format(seq(from=518000,to=522000,by=2000),scientific=T))
for (i in c(2:12,14,15)) {
  with(retVSbic_runs[[i]],lines(n_reticulations,bic,col=rgb(0,0,0,0.5)))
}
with(retVSbic_runs[[13]],lines(n_reticulations,bic,col="darkorange",lwd=2))
legend(x="topright",legend=c("lowest BIC","most reticulations"),
       col=c("darkorange","deeppink"),
       lty=c("solid","solid"),bty="n",cex=1.5,lwd=3)

# Plot log-likelihood against no. of reticulations for each run
# min(vapply(retVSbic_runs,FUN=function(x) min(x$logl),FUN.VALUE=1))
# max(vapply(retVSbic_runs,FUN=function(x) max(x$logl),FUN.VALUE=1))
with(retVSbic_runs[[1]],
     plot(n_reticulations,logl,type="l",ylab="Log-likelihood",xlab="",
          col="deeppink",lwd=2,yaxt="n",ylim=c(-260064.1,-257706.4),cex.axis=1,
          cex.lab=1.5))
axis(2,at=seq(from=-260000,to=-258000,by=1000),
     labels=format(seq(from=-260000,to=-258000,by=1000),scientific=T))
for (i in c(2:12,14,15)) {
  with(retVSbic_runs[[i]],lines(n_reticulations,logl,col=rgb(0,0,0,0.5)))
}
with(retVSbic_runs[[13]],lines(n_reticulations,logl,col="darkorange",lwd=2))
mtext(text="No. of reticulations",side=1,line=0,outer=T,cex=1.5)
dev.off()
