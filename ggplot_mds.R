
# adding coloumns in passport. FID1 and Lgabs
mds<- read.table("plink04_04_16.mds", header=TRUE, quote="\"")
passport <- read.csv("passport.csv")
mds$FID1<-gsub ('\\-[0-9]{1}','',mds[,1])
mds.passport <- merge(passport, mds, by="FID1")
colnames(mds.passport)[7]<-'Country'
library(ggplot2)

### plot
pmdscountry <- ggplot(mds.passport, aes(C1, C2)) + geom_point(aes(colour = Country))
pmdscountry  + ggtitle("mds by country")
lg <- ggplot(mds.passport, aes(C1, C2)) + geom_point(aes(colour = LG)) + scale_colour_gradient2(mid=("green"), low="black", high="purple")
lg  + ggtitle("mds by LG")
lgabs <- ggplot(mds.passport, aes(C1, C2)) + geom_point(aes(colour = Lgabs)) + scale_colour_gradient2(low="black",mid="green",high="purple")
lgabs  + ggtitle("mds by LGabs")

