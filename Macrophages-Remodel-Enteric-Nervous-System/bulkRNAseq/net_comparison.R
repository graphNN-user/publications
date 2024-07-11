# this code from Lauren for network comparison and calculate p values
#
x = read.table("/Users/user/Documents/Renrichment.csv",sep=",",header=T)

fc = (x[,"overlap"]/x[,"geneset"])/(x[,"subnetwork"]/x[,"network"])
pv = 1-phyper(x[,"overlap"]-1, x[,"subnetwork"], x[,"network"]-x[,"subnetwork"], x[,"geneset"])
x = cbind(x, fc, pv)
write.table(x, "/Users/user/Documents/Renrichmentresults.csv",sep=",",quote=F)
