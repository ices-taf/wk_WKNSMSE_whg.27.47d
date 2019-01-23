### ------------------------------------------------------------------------ ###
### create FLStock for whiting ####
### ------------------------------------------------------------------------ ###
### load files from package mse for easier debugging

#devtools::install_github("fishfollower/SAM/stockassessment")
#
#devtools::install_github(repo = "flr/mse")
#devtools::install_github("shfischer/FLfse/FLfse")
#install.packages(pkgs = c("FLCore"), repos = "http://flr-project.org/R")
#install.packages(pkgs = c("FLAssess"), repos = "http://flr-project.org/R")
#install.packages(pkgs = c("FLash"), repos = "http://flr-project.org/R")
#install.packages(pkgs = c("FLCore"), repos = "http://flr-project.org/R")
#install.packages(pkgs = c("ggplotFL"), repos = "http://flr-project.org/R")
#devtools::install("D:\\Tanja\\WKNSMSE\\mse-master\\mse-master")

### base on SAM assessment

### check versions of required R packages
if (packageVersion("FLCore") < "2.6.11.9001") 
  stop("please update FLCore")
if (packageVersion("FLfse") < "0.0.0.9003") 
  stop("please update FLfse")
if (packageVersion("stockassessment") < "0.8.1") 
  stop("please update stockassessment")
if (packageVersion("mse") < "0.9.2") 
  stop("please update mse")


### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
### load files from package mse for easier debugging
# devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

setwd(paste0("/home/miethet/MSE_whiting/whiting_new"))

source("a4a_mse_WKNSMSE_funs.R")


### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations/replicates
n <- 1000
### number of years one extra (MSE n_years-1)
n_years <- 20
### last data year
yr_data <- 2018

### ------------------------------------------------------------------------ ###
### fit SAM ####
### ------------------------------------------------------------------------ ###
### use input data provided in FLfse
### recreates the WGNSSK2018 whiting assessment
fit <- FLR_SAM(stk = whg4_stk, idx = whg4_idx, conf = whg4_conf_sam)

is(fit)
fit
plot(fit)

### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit)
sam_initial$logScale <- numeric(0)


### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration
### cod4_stk2 is used as template, i.e. the input values (catch) include
### the catch multiplier, 
### the results (stock numbers & harvest) are used from the real WGNSSK fit
stk <- SAM2FLStock(object = fit, stk = whg4_stk)
summary(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                              "", "", "f", "", ""))
plot(stk)

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### add uncertainty ####
### ------------------------------------------------------------------------ ###
### first approach: use variance-covariance

### add iteration dimension
stk <- FLCore::propagate(stk, n)
dim(stk)

### add uncertainty estimated by SAM as iterations in historic time period
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n, print_screen = FALSE)
### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

### catch noise added later

plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### maximum observed F
max(fbar(stk))
max(harvest(stk))


### ------------------------------------------------------------------------ ###
### extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### special case for NS cod, whiting
### maturity data available for 2018 (based on the IBTS Q1)
### although SAM estimates F in 2018, this is not reported or taken forward into
### forcasts by the WG
# stk_stf2017 <- stf(window(stk, end = 2017), n_years + 1)
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### use all available data
stk_stf <- stk_stf2018

### ------------------------------------------------------------------------ ###
### biological data for OM ####
### ------------------------------------------------------------------------ ###
### Resample maturity and natural mortality from the last 10 years 
### (2008-2017), last 3 years for weights (2015-2017)
### set up an array with one resampled year for each projection year 
### (including intermediate year) and replicate
### use the same resampled year for all biological parameters
### this is the approach used in eqsim for North Sea cod
set.seed(2)

### use last 10 data years to sample biological parameters
sample_yrs <- 2008:2017
### get year position of sample years
sample_yrs_pos <- which(dimnames(stk_stf)$year %in% sample_yrs)

### use last three data years to sample selectivity parameters
sample_yrs1 <- 2015:2017
### get year position of sample years
sample_yrs_pos1 <- which(dimnames(stk_stf)$year %in% sample_yrs1)


### create samples for biological data (weights, etc.)
### the historical biological parameters are identical for all iterations
### and consequently do not need to be treated individually
### (but keep age structure)
### create vector with resampled years
bio_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years + 1) * n, replace = TRUE)
### do the same for selectivity
sel_samples <- sample(x = sample_yrs_pos1, 
                      size = (n_years + 1) * n, replace = TRUE)
### years to be populated
bio_yrs <- which(dimnames(stk_stf)$year %in% 2018:dims(stk_stf)$maxyear)


### insert values
catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1])
mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1])
### maturity data for 2018 exists, re-insert real data
mat(stk_stf)[, ac(2018)] <- mat(stk_orig)[, ac(2018)]
### use different samples for selectivity
harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])

plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### use only data from 1983 and later
sr <- as.FLSR(window(stk_stf, start = 1983), model = "segreg")
sr1 <- as.FLSR(window(stk_stf, start = 2002), model = "segreg")

### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr)))
suppressWarnings(. <- capture.output(sr1 <- fmle(sr1)))

### run in parallel
# library(doParallel)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# sr <- fmle_parallel(sr, cl)


tiff("output/sr1983.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
plot(sr[,,,,,1], cex=0.5)
dev.off()

tiff("output/sr1983_alliters.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
plot(sr, cex=0.5)
dev.off()

tiff("output/sr2002.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
plot(sr1[,,,,,1], cex=0.5)
dev.off()



### check breakpoints
summary(params(sr)["b"]) # median 120141
summary(params(sr1)["b"])# median 116381

### plot model and data
as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb)) %>%
  mutate(age = NULL,
         year = ifelse(qname == "SSB", year + 1, year)) %>%
  tidyr::spread(key = qname, value = data) %>%
  ggplot() +
  geom_point(aes(x = SSB, y = rec, group = iter), 
             alpha = 0.5, colour = "grey", shape = 1) +
  geom_line(aes(x = SSB, y = fitted, group = iter)) +
  theme_bw() + xlim(0, NA) + ylim(0, NA)

### Check extent of autocorrelation
# Some autocorrelation

ACFrec <- acf(window(stock.n(stk_orig)[1], start = 1983))
acfRecLag11 <- round(ACFrec$acf[,,][2],2)

ACFrec1 <-  acf(window(stock.n(stk_orig)[1], start = 2002))
acfRecLag12<- round(ACFrec1$acf[,,][2],2)


# autocorrelation plots
tiff("output/acf1983.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,1))
acf(window(stock.n(stk_orig)[1], start = 1983),lag.max=7, plot=T, main=paste("Autocor. in Rec, Lag1 =",acfRecLag11,sep=" "), cex=0.5)
dev.off()

tiff("output/acf2002.tiff", bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,1))
acf(window(stock.n(stk_orig)[1], start = 2002),lag.max=7, plot=T, main=paste("Autocor. in Rec, Lag1 =",acfRecLag12,sep=" "), cex=0.5)
dev.off()



### Check method proposed for generating recruitment compares with past recruitment estimates
test <- as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb))
test <- mutate(test, age = NULL, year = ifelse(qname == "SSB", year + 1, year))
test <- tidyr::spread(test, key = qname, value = data)
test <- test[complete.cases(test),]
test$res <- test$rres <- rep(NA, nrow(test))

# Generate residuals for future recruitments
foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
        .errorhandling = "pass") %do% {
          
          set.seed(iter_i^2)
          
          ### get residuals for current iteration
          res_i <- c(FLCore::iter(residuals(sr), iter_i))
          res_i <- res_i[!is.na(res_i)]
          
          ### calculate kernel density of residuals
          density <- density(x = res_i)
          ### sample residuals
          mu <- sample(x = res_i, size = length(res_i), replace = TRUE)
          ### "smooth", i.e. sample from density distribution
          rres <- rnorm(n = length(res_i), mean = mu, sd = density$bw)
          
          ### calculate serial correlation
          #res_y <- res_i[-length(res_i)]
          #res_yp1 <- res_i[-1]
          #rho <- sum(res_y * res_yp1) / sqrt(sum(res_y^2) * sum(res_yp1^2))
          ACFrec <- acf(res_i, lag.max=1)
          rho <- ACFrec$acf[,,][2] #lag 1
         
          ### generate autocorrelated residuals
          res <- res_i
          res[2:length(res)] <- 0
          for (r in 2:length(res)){
            res[r] <- rho * res[r-1] + sqrt(1 - rho^2) * rres[r]
          }
          
          test$rres[test$iter==iter_i] <- rres
          test$res[test$iter==iter_i] <- res
          
        }

# Generate future recruitments from past SSBs and generated residuals
test$future <- test$fitted * exp(test$res)

# 10 randomly selected iters for plotting
# should probably increase the number later
i_samp <- sample(seq(dim(sr)[6]), 20, replace=FALSE)



### generate residuals for MSE, three options since 1983 (with AR(1), 2002 (with/without AR(1))
### years with missing residuals
# NW: dimnames produces NULL for me
# yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution

#1983 with AR(1)
res_new_1983 <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %dopar% {
                     
  set.seed(iter_i)
  
  ### get residuals for current iteration
  res_i <- c(FLCore::iter(residuals(sr), iter_i))
  res_i <- res_i[!is.na(res_i)]
  
  ### calculate kernel density of residuals
  density <- density(x = res_i)
  ### sample residuals
  mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
  ### "smooth", i.e. sample from density distribution
  rres <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
  
  ### calculate serial correlation
  #res_y <- res_i[-length(res_i)]
  #res_yp1 <- res_i[-1]
  #rho <- sum(res_y * res_yp1) / sqrt(sum(res_y^2) * sum(res_yp1^2))
  ACFrec <- acf(res_i, lag.max=1)
  rho <- ACFrec$acf[,,][2] #lag 1
         
  ### generate autocorrelated residuals
  res_new <- rep(0, length(yrs_res))
  res_new[1] <- rho * res_i[length(res_i)] + sqrt(1 - rho^2) * rres[1]
  for (r in 2:length(res_new)){
    res_new[r] <- rho * res_new[r-1] + sqrt(1 - rho^2) * rres[r]
  }

  return(res_new)
  
}


# 2002 without autocorrelation
res_new_2002 <- foreach(iter_i = seq(dim(sr1)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %dopar% {
                     
  set.seed(iter_i)
  
  ### get residuals for current iteration
  res_i <- c(FLCore::iter(residuals(sr1), iter_i))
  res_i <- res_i[!is.na(res_i)]
  
  ### calculate kernel density of residuals
  density <- density(x = res_i)
  ### sample residuals
  mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
  ### "smooth", i.e. sample from density distribution
  rres <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
  
  ### calculate serial correlation
  #res_y <- res_i[-length(res_i)]
  #res_yp1 <- res_i[-1]
  #rho <- sum(res_y * res_yp1) / sqrt(sum(res_y^2) * sum(res_yp1^2))
  #ACFrec <- acf(res_i)
  rho <- 0 #ACFrec$acf[,,][2] #lag 1
         
  ### generate autocorrelated residuals
  res_new <- rep(0, length(yrs_res))
  res_new[1] <- rho * res_i[length(res_i)] + sqrt(1 - rho^2) * rres[1]
  for (r in 2:length(res_new)){
    res_new[r] <- rho * res_new[r-1] + sqrt(1 - rho^2) * rres[r]
  }

  return(res_new)
  
}

#since 1983
sr_<-sr
summary(exp(unlist(res_new_1983)))

### insert into model
residuals(sr_)[, yrs_res] <- unlist(res_new_1983)
### exponeniate residuals to get factor
residuals(sr_) <- exp(residuals(sr_))
sr_res <- residuals(sr_)
plot(sr_res)

write.csv(as.array(summary(exp(unlist(res_new_1983)))),file=paste0("output/",n,"summary_res_new1983.csv"))


tiff(paste0("output/",n,"res1983_lag1.tiff"), bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,2))
plot(residuals(sr_))
dev.off()



# since 2002
sr1_0<-sr1
summary(exp(unlist(res_new_2002)))
### insert into model
residuals(sr1_0)[, yrs_res] <- unlist(res_new_2002)
### exponeniate residuals to get factor
residuals(sr1_0) <- exp(residuals(sr1_0))
sr_res10 <- residuals(sr1_0)
plot(sr_res10)

write.csv(as.array(summary(exp(unlist(res_new_2002)))),file=paste0("output/",n,"summary_res_new2002.csv"))


tiff(paste0("output/",n,"res2002_no.tiff"), bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,2))
plot(residuals(sr1_0))
dev.off()



#######################################
# 75% factor, lower residuals, multiplied
#sample periods 1-4 for each iteration

#2002
sr1_0_dip<-sr1_0
year1<-2018
nits<-n
RS <- 20  # to make sure we make enough shifts for total projection period
duration_byIter <- vector()
res2<-residuals(sr1_0) # since 2002
rule<-c(1,0.75)

for (its in 1:nits) {
#print(its)

x<-as.integer(runif(RS)*4+1)   # max period length 4 years
yrx <- x
yrx[1] <- yrx[1] + (year1-1)

for(i in 2:RS){
    yrx[i] <- yrx[i-1] + x[i]
}
yrx <- c((year1-1),yrx)
sxx <- sum(x)
duration_byIter <- c(duration_byIter,sxx)
#print(x)

#forward projection for the differing periods
model<-x #store for recruitment model indices, 2 different scenarios
for(i in 1:RS){
ind<-as.integer(runif(1)*2+1)  # hrules decides how many recruitment scenarios (mdedium, low)
model[i]<-ind
}
#print(model)

for(i in 1:RS){
         for (ii in 1:x[i]){
         
if(as.character(yrx[i]+ii)<=2039) res2[,as.character(yrx[i]+ii),,,,its]<-res2[,as.character(yrx[i]+ii),,,,its]*rule[model[i]]

}}}

residuals(sr1_0_dip)<-res2

tiff(paste0("output/",n,"res2002_dip.tiff"), bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,2))
plot(res2)
dev.off()



#1983
year1<-2018
nits<-n
RS <- 20  
sr_dip<-sr_
rule<-c(1,0.75)
duration_byIter <- vector()
res3<-residuals(sr_)  # since 1983

for (its in 1:nits) {
#print(its)

x<-as.integer(runif(RS)*4+1)   # max period length 4 years
yrx <- x
yrx[1] <- yrx[1] + (year1-1)

for(i in 2:RS){
    yrx[i] <- yrx[i-1] + x[i]
}
yrx <- c((year1-1),yrx)
sxx <- sum(x)
duration_byIter <- c(duration_byIter,sxx)
#print(x)

#forward projection for the differing periods
model<-x #store for recruitment model indices, 2 different scenarios
for(i in 1:RS){
ind<-as.integer(runif(1)*2+1)  # hrules decides how many recruitment scenarios (low/medium, low)
model[i]<-ind
}
#print(model)

for(i in 1:RS){
         for (ii in 1:x[i]){
         
if(as.character(yrx[i]+ii)<=2039) res3[,as.character(yrx[i]+ii),,,,its]<-res3[,as.character(yrx[i]+ii),,,,its]*rule[model[i]]

}}}

residuals(sr_dip)<-res3

tiff(paste0("output/",n,"res1983_dip.tiff"), bg="white",  res=200, width = 900, height = 900)
par(mar=c(4,4,4,2))
plot(res3)
dev.off()



tiff(paste0("output/",n,"res1983_no_hist.tiff"), bg="white",  res=200, width = 900, height = 1100)
par(mar=c(4,4,4,2))
hist(residuals(sr_)[,36:57,,,,],col=rgb(1,0,0,0.2), xlab="",ylim=c(0,6000), main="histogram residuals 1983 series")
legend("topright",legend=c("without level change", "with level change"),bty="n", col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2) ), pch=15)
hist(residuals(sr_dip)[,36:57,,,,],col=rgb(0,0,1,0.2), add=T)
dev.off()

tiff(paste0("output/",n,"res2002_no_hist.tiff"), bg="white",  res=200, width = 900, height =1100)
par(mar=c(4,4,4,2))
hist(residuals(sr1_0)[,17:38,,,,],col=rgb(1,0,0,0.2), xlab="",ylim=c(0,6000), main="histogram residuals 2002 series")
legend("topright",legend=c("without level change", "with level change"),bty="n", col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2) ), pch=15)
hist(residuals(sr1_0_dip)[,17:38,,,,],col=rgb(0,0,1,0.2), add=T)
dev.off()


# Generate future recruitments from past SSBs and generated residuals

rec1 <- as.data.frame( residuals(sr_))
rec2 <- as.data.frame( residuals(sr_dip))

test$future1<-test$fitted*rec1[rec1$year%in%c(2005:2039),"data"]
test$future2<-test$fitted*rec2[rec2$year%in%c(2005:2039),"data"]


# Plot past and future stock recruit pairs for selected iters
tiff("output/sr1983_srpairs_iters.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
ggplot(test[is.element(test$iter, i_samp),]) +
  geom_point(aes(x = SSB, y = rec), 
             alpha = 0.5, colour = "red", shape = 19) +
  geom_point(aes(x = SSB, y = future1), 
             alpha = 0.5, colour = "black", shape = 19) +
  geom_line(aes(x = SSB, y = fitted)) +
  facet_wrap(~iter) +
  theme_bw() + xlim(0, NA) + ylim(0, NA)
dev.off()

# Empirical cumulative distributions for the same iters
tiff("output/sr1983_cumulative_iters.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
ggplot(test[is.element(test$iter, i_samp),]) +
  stat_ecdf(aes(rec), geom = "step", colour = "red") +
  stat_ecdf(aes(future1), geom = "step", colour = "black") +
  facet_wrap(~iter) +
  theme_bw() + xlim(0, NA) + ylim(0, NA)
dev.off()

# Combine previous two plots over iters
# Stock recruit pairs
tiff("output/sr1983_srpairs.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))

ggplot(test) +
  geom_point(aes(x = SSB, y = rec), 
             alpha = 0.5, colour = "red", shape = 19) +
  geom_point(aes(x = SSB, y = future1), 
             alpha = 0.5, colour = "black", shape = 19) +
   ylab("Rec")+           
  theme_bw() + xlim(0, NA) + ylim(0, NA)
dev.off()

# Empirical cumulative distribution

tiff("output/sr1983_cumulative.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))
ggplot(test) +
  stat_ecdf(aes(rec), geom = "step", colour = "red") +
  stat_ecdf(aes(future1), geom = "step", colour = "black") +
  theme_bw() + xlim(0, NA) + ylim(0, NA)
dev.off()


# Combine previous two plots over iters
# Stock recruit pairs
tiff("output/sr1983_compare_srpairs.tiff", bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))

ggplot(test) +
 geom_point(aes(x = SSB, y = future1), 
            alpha = 0.5, colour = "red", shape = 19) +
  geom_point(aes(x = SSB, y = future2), 
             alpha = 0.3, colour = "black", shape = 19) +
      ylab("Rec")+
  theme_bw() + xlim(0, NA) + ylim(0, NA)
  
dev.off()

rm(test, i_samp)

### ------------------------------------------------------------------------ ###
### process noise ####
### ------------------------------------------------------------------------ ###
### create FLQuant with process noise
### this will be added to the values obtained from fwd() in the MSE

### create noise for process error
set.seed(3)
proc_res <- stock.n(stk_stf) %=% 0 ### template FLQuant
proc_res[] <- stats::rnorm(n = length(proc_res), mean = 0, 
                           sd = uncertainty$proc_error)
### the proc_res values are on a normale scale,
### exponentiate to get log-normal 
proc_res <- exp(proc_res)
### proc_res is a factor by which the numbers at age are multiplied

### for historical period, numbers already include process error from SAM
### -> remove deviation
proc_res[, dimnames(proc_res)$year <= 2017] <- 1

### remove deviation for first age class (recruits)
proc_res[1, ] <- 1

### try saving in stock recruitment model ... 
### this gets passed on to the projection module

fitted(sr_dip)<-proc_res
fitted(sr_)<-proc_res


tiff(paste0("output/",n,"process_noise.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(proc_res)
dev.off()



### ------------------------------------------------------------------------ ###
### stf for 2018: assume catch advice is taken ####
### ------------------------------------------------------------------------ ###
c2018 <- 26191
ctrl <- fwdControl(data.frame(year = 2018, quantity = "catch", 
                              val = c2018))

### project forward for intermediate year (2018)
stk_int <- stk_stf
stk_int[] <- fwd(stk_stf, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = TRUE, maxF = 5)[]
### add process noise
stock.n(stk_int) <- stock.n(stk_int) * proc_res
stock(stk_int)[] <- computeStock(stk_int)

### create stock for MSE simulation
stk_fwd <- stk_stf
### insert values for 2018
stk_fwd[, ac(2018)] <- stk_int[, ac(2018)]
### insert stock number for 2019 in order to calculate SSB at beginning of 
### 2019
stock.n(stk_fwd)[, ac(2019)] <- stock.n(stk_int)[, ac(2019)]
stock(stk_fwd)[, ac(2019)] <- computeStock(stk_fwd[, ac(2019)])

#all.equal(window(stk_fwd, end = 2018), window(stk_stf, end = 2018))

### ------------------------------------------------------------------------ ###
### biological data for OEM observation error model ####
### ------------------------------------------------------------------------ ###

### base on OM stock
stk_oem <- stk_fwd

### projection years
proj_yrs <- 2018:range(stk_oem)[["maxyear"]]

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(catch.wt(stk_oem)[, ac(sample_yrs)])
landings.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(landings.wt(stk_oem)[, ac(sample_yrs)])
discards.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(discards.wt(stk_oem)[, ac(sample_yrs)])
stock.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(stock.wt(stk_oem)[, ac(sample_yrs)])
m(stk_oem)[, ac(proj_yrs)] <- yearMeans(m(stk_oem)[, ac(sample_yrs)])
### maturity starts one year later because there is data for 2018
mat(stk_oem)[, ac(proj_yrs[-1])] <- yearMeans(mat(stk_oem)[, ac(sample_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA


tiff(paste0("output/",n,"OM_noise_mat.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(mat(stk_fwd))
dev.off()

tiff(paste0("output/",n,"OM_noise_mort.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(m(stk_fwd))
dev.off()

tiff(paste0("output/",n,"OM_noise_catchwghts.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(catch.wt(stk_fwd))
dev.off()

tiff(paste0("output/",n,"OM_noise_stockwghts.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(stock.wt(stk_fwd))
dev.off()

tiff(paste0("output/",n,"OM_noise_harvest.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(harvest(stk_fwd))
dev.off()



### ------------------------------------------------------------------------ ###
### indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- whg4_idx
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n)

### insert catchability
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
  
}
### create copy of index with original values
idx_raw <- lapply(idx ,index)
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx)

### create deviances for indices
### first, get template
idx_dev <- lapply(idx, index)
### create random noise based on sd
set.seed(4)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                     mean = 0, sd = idx_dev[[idx_i]])
  ### exponentiate to get from normal to log-normal scale
  idx_dev[[idx_i]] <- exp(idx_dev[[idx_i]])
}

### modify residuals for historical period so that index values passed to 
### stock assessment are the ones observed in reality
### IBTS Q1, values up to 2018
idx_dev$"IBTS-Q1"[, dimnames(idx_dev$"IBTS-Q1")$year <= 2018] <- 
  idx_raw$"IBTS-Q1"[, dimnames(idx_raw$"IBTS-Q1")$year <= 2018] /
  index(idx$"IBTS-Q1")[, dimnames(idx$"IBTS-Q1"@index)$year <= 2018]
### IBTS Q3, values up to 2017
idx_dev$"IBTS-Q3"[, dimnames(idx_dev$"IBTS-Q3")$year <= 2017] <- 
  idx_raw$"IBTS-Q3"[, dimnames(idx_raw$"IBTS-Q3")$year <= 2017] /
  index(idx$"IBTS-Q3")[, dimnames(idx$"IBTS-Q3"@index)$year <= 2017]

### compare simulated to original survey(s)

tiff(paste0("output/",n,"Survey.tiff"), bg="white",  res=200, width = 1100, height = 1200)
par(mar=c(4,4,4,2))

as.data.frame(FLQuants(whg4_q1 = index(whg4_idx$"IBTS-Q1"), 
                       whg4_q3 = index(whg4_idx$"IBTS-Q3"),
                       sim_q1 = (index(idx$"IBTS-Q1")),
                       sim_q3 = (index(idx$"IBTS-Q3"))
)) %>%
  mutate(survey = ifelse(grepl(x = qname, pattern = "*_q1$"), "Q1", "Q3"),
         source = ifelse(grepl(x = qname, pattern = "^sim*"), "sim", "data")) %>%
  filter(year <= 2019) %>%
  ggplot(aes(x = year, y = data, colour = source)) +
  facet_grid(paste("age", age) ~ paste("IBTS", survey), scales = "free_y") +
  stat_summary(fun.y = quantile, fun.args = 0.25, geom = "line",
               alpha = 0.5) +
  stat_summary(fun.y = quantile, fun.args = 0.75, geom = "line",
               alpha = 0.5) +
  stat_summary(fun.y = median, geom = "line") +
  theme_bw()
dev.off()

### check survey
# idx0 <- calc_survey(stk = stk_fwd, idx = idx)
# idx0 <- lapply(seq_along(idx0), function(idx_i) {
#   idx_tmp <- idx0[[idx_i]]
#   index(idx_tmp) <- index(idx_tmp) * idx_dev[[idx_i]]
#   return(idx_tmp)
# })
# plot(index(idx0[[2]]))

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

### create noise for catch
set.seed(5)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
### the catch_res values are on a normale scale,
### exponentiate to get log-normal 
catch_res <- exp(catch_res)
### catch_res is a factor by which the numbers at age are multiplied

### for historical period, pass on real observed catch
### -> remove deviation
catch_res[, dimnames(catch_res)$year <= 2017] <- 1

tiff(paste0("output/",n,"catch_noise.tiff"), bg="white",  res=200, width = 900, height = 1200)
par(mar=c(4,4,4,2))
plot(catch_res)
dev.off()

#create extra implementation error noise for total catch

set.seed(6)
catch_res_l <- catch(stk_fwd) %=% 0 ### template FLQuant
catch_res_l[] <- stats::rnorm(n = length(catch_res_l), mean = 0, 
                            sd = fit$sdrep$sd[names(fit$sdrep$value)=="logCatch"])
catch_res_l <- exp(catch_res_l)
catch_res_l[, dimnames(catch_res_l)$year <= 2017] <- 1

tiff(paste0("output/",n,"impl_error.tiff"), bg="white",  res=200, width = 1000, height = 1000)
par(mar=c(4,4,4,2))
plot(catch_res_l, xlab="Year")
dev.off()



### ------------------------------------------------------------------------ ###
### check SAM ####
### ------------------------------------------------------------------------ ###
# 
# stk_tmp <- window(stk_stf, end = 2018)
# catch.wt(stk_tmp)[,ac(2018)] <- landings.wt(stk_tmp)[,ac(2018)] <-
#   discards.wt(stk_tmp)[,ac(2018)] <- NA
# stk_tmp <- stk_tmp[,,,,, 1:10]
# idx_tmp <- window(idx, end = 2018)
# idx_tmp[[2]] <- window(idx_tmp[[2]], end = 2017)
# idx_tmp <- lapply(idx_tmp, FLCore::iter, 1:1)
# 

# fit3 <- FLR_SAM(stk = stk_tmp, idx = idx_tmp, conf = whg4_conf_sam)

# stk3 <- SAM2FLStock(fit3)
# plot(iter(FLStocks(whg4 = stk, sim = stk3), 1))
# 
# all.equal(idx_tmp, whg4_idx)

### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###

### path
input_path <- paste0("input/whg4/", n, "_", n_years, "/")
dir.create(input_path)

### stock
saveRDS(stk_fwd, file = paste0(input_path, "stk.rds"))
### stock recruitment
saveRDS(sr, file = paste0(input_path, "sr.rds"))
saveRDS(sr1, file = paste0(input_path, "sr1.rds"))
saveRDS(sr_, file = paste0(input_path, "sr1983_AR.rds"))
saveRDS(sr1_0, file = paste0(input_path, "sr2002.rds"))
saveRDS(sr_dip, file = paste0(input_path, "sr1983_AR_dip.rds"))
saveRDS(sr1_0_dip, file = paste0(input_path, "sr2002_dip.rds"))

### recruitment residuals
saveRDS(residuals(sr_), file = paste0(input_path, "sr1983_AR1_res.rds"))
saveRDS(residuals(sr1_0), file = paste0(input_path, "sr2002_res.rds"))
saveRDS(residuals(sr_dip), file = paste0(input_path, "sr1983_AR1_dip_res.rds"))
saveRDS(residuals(sr1_0_dip), file = paste0(input_path, "sr2002_dip_res.rds"))

### surveys
saveRDS(idx, file = paste0(input_path, "idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, "idx_dev.rds"))
### catch noise
saveRDS(catch_res, file = paste0(input_path, "catch_res.rds"))
### process error
saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))
### observed stock
saveRDS(stk_oem, file = paste0(input_path, "stk_oem.rds"))
### sam initial parameters
saveRDS(sam_initial, file = paste0(input_path, "sam_initial.rds"))
### sam configuration
saveRDS(whg4_conf_sam, file = paste0(input_path, "whg4_conf_sam"))

### imp error catch noise
saveRDS(catch_res_l, file = paste0(input_path,"catch_res_l.rds"))


rm(sr,sr1)

sr<-sr_dip  # SSR since 1983 (AR1), level shift and implementation error
save.image(file = paste0("input/whg4/base_image/image_OM3_",n,".RData"))


rm(sr,catch_res_l)
sr<-sr_dip  # SSR since 1983 (AR1), level shift and no implementation error
save.image(file = paste0("input/whg4/base_image/image_OM2_",n,".RData"))

rm(sr,catch_res_l)  # 
sr<-sr_  # SSR since 1983 (AR1), no implementation error
save.image(file = paste0("input/whg4/base_image/image_OM1_",n,".RData"))

