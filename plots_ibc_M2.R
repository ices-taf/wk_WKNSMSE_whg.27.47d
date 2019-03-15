setwd(paste("N:/1_Tanja/1_MARLAB/WKNSMSE", sep=""))

data<-read.csv("IBC_M2_corr.csv", sep=",")
data1<-data[data$year%in%c(1996:2007,2009:2017),]

tiff("IBC_catch_1996.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,1))

plot(data1$catch,data1$ibc, bty="l", ylab="IBC (t)", xlab="Catch (t)", xlim=c(0,max(data1$catch)*1.05), ylim=c(0,max(data1$ibc)*1.05),pch=20)
dev.off()


data1<-data[data$year%in%c(1996:2007,2009:2017),]

tiff("IBC_prop_time1996.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,1))

plot(data1$year,100*data1$ibc/data1$catch, bty="l", ylab="IBC (% of total catch)", xlab="Year", xlim=c(1996,max(data1$year)), ylim=c(0,17),pch=20)
dev.off()




#no relationship with SSB
tiff("IBCcatch_SSB_1996.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,5,4,1))

plot(data1$ssb/1000,100*data1$ibc/data1$catch, bty="l", ylab="IBC (% of total catch)", xlab="SSB (1000 t)", ylim=c(0,1.05*max(100*data1$ibc/data1$catch)), xlim=c(50,1.1*max(data1$ssb/1000)),pch=20)
dev.off()


tiff("IBCcatch_rec_1996.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,5,4,1))

plot(data1$rec/1000000,100*data1$ibc/data1$catch, bty="l", ylab="IBC (% of total catch)", xlab="Recruitment (billion)", ylim=c(0,1.05*max(100*data1$ibc/data1$catch)), xlim=c(0,1.1*max(data1$rec/1000000)),pch=20)
dev.off()

tiff("IBCcatch_age1_1996.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,5,4,1))

plot(data1$Nage1/1000000,100*data1$ibc/data1$catch, bty="l", ylab="IBC (% of total catch)", xlab="N age1 (billion)", ylim=c(0,1.05*max(100*data1$ibc/data1$catch)), xlim=c(0,1.1*max(data1$Nage1/1000000)),pch=20)
dev.off()


tiff("SSB_M2.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,5,4,1))

plot(data$ssb/1000,data$M2, bty="l", ylab="Predation mortality (M2, age 1)", xlab="SSB (1000 t)", ylim=c(0.4,1.05*max(data1$M2)), xlim=c(100,1.1*max(data1$ssb/1000)),pch=20)
dev.off()


