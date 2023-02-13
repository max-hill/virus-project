# RF-net embedding cost profiles for "15 genes, 10Kbp blocks" partition

# genes defined left to right, RF-net in regular mode
l2r <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted",
                  header=T)
# ... , RF-net in fast mode
l2r_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_f",
                    header=T)
# ..., RF-net in tree-child mode
l2r_tc <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc",
                     header=T)
# ..., RF-net in tree-child mode and fast mode
l2r_tc_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_set1c_ogrooted_tc_f",
                       header=T)

# genes defined right to left, RF-net in regular mode
r2l <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted",
                  header=T)
# ..., RF-net in fast mode
r2l_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_f",
                    header=T)
# ..., RF-net in tree-child mode
r2l_tc <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_tc",
                     header=T)
# ..., RF-net in tree-child mode and fast mode
r2l_tc_f <- read.table(file="results/bhv/rfnet/embedding cost/15genes_blksize10000_rev_set1c_ogrooted_tc_f",
                       header=T)

# compare embedding cost profiles under different RF-net modes
par(mfrow=c(2,2))
plot(x=l2r$r,y=l2r$dist,xlab="No. of reticulations",ylab="embedding cost",
     type="l",col="red")
title(main="l2r vs l2r_f",line=.5)
lines(x=l2r_f$r,y=l2r_f$dist,col="blue")

plot(x=l2r$r,y=l2r$dist,xlab="",ylab="",
     type="l",col="red")
title(main="l2r vs l2r_tc",line=.5)
lines(x=l2r_tc$r,y=l2r_tc$dist,col="blue")

plot(x=l2r$r,y=l2r$dist,xlab="",ylab="",
     type="l",col="red")
title(main="l2r vs l2r_tc_f",line=.5)
lines(x=l2r_tc_f$r,y=l2r_tc_f$dist,col="blue")

par(mfrow=c(2,2))
plot(x=r2l$r,y=r2l$dist,xlab="No. of reticulations",ylab="embedding cost",
     type="l",col="red")
title(main="r2l vs r2l_f",line=.5)
lines(x=r2l_f$r,y=r2l_f$dist,col="blue")

plot(x=r2l$r,y=r2l$dist,xlab="",ylab="",
     type="l",col="red")
title(main="r2l vs r2l_tc",line=.5)
lines(x=r2l_tc$r,y=r2l_tc$dist,col="blue")

plot(x=r2l$r,y=r2l$dist,xlab="",ylab="",
     type="l",col="red")
title(main="r2l vs r2l_tc_f",line=.5)
lines(x=r2l_tc_f$r,y=r2l_tc_f$dist,col="blue")
