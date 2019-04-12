### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(ggplotFL)
library(tidyr)
library(cowplot)
library(dplyr)

library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)

### load additional functions
setwd(paste("/home/miethet/MSE_whiting/", sep=""))
source("a4a_mse_WKNSMSE_funs.R")

iters<-1000
years<-20
om<-1
Blim<-119970
### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- paste0("output/runs/whg4/",iters,"_",years,"/")
files_res <- data.frame(file = list.files(path_res, pattern = "*E.rds"), 
                        stringsAsFactors = FALSE)
files_res <- files_res[files_res$file != "stats.rds",, drop = FALSE]


files_res$OM <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 1), gsub,
                                     pattern = "OM", replacement = ""))
files_res$Ftrgt <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 3), gsub,
                                     pattern = "Ftrgt-", replacement = ""))
                                     
files_res$Btrigger <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 4), gsub,
                                     pattern = "Btrigger-", replacement = ""))
files_res$HCR <- sapply(lapply(strsplit(files_res$file, 
                                        split = "_", fixed = TRUE),
                               "[[", 2), gsub,
                        pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(sapply(lapply(strsplit(files_res$file, 
                                                      split = "_", fixed = TRUE),
                                               "[[", 5), gsub,
                                      pattern = "TACconstr-", replacement = ""))
files_res$BB <- as.logical(
  sapply(lapply(lapply(strsplit(files_res$file, split = "_", fixed = TRUE), 
                       "[[", 6), gsub, pattern = "BB-", replacement = ""), 
         gsub, pattern = ".rds", replacement = ""))

#stats <- readRDS(paste0(path_res, "stats.rds"))
stats <- read.csv(paste0(path_res, "stats.csv"))
stats_new <- merge(stats, files_res, all = TRUE)


stats_new <- stats_new[!stats_new$file %in% stats$file, ]


library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)

res_list <- vector(mode = "list", length = nrow(stats_new))

res_list <- foreach(i = seq(nrow(stats_new))) %dopar% {
  readRDS(paste0(path_res, stats_new$file[i]))
}

### ------------------------------------------------------------------------ ###
### calculate summary statistics ####
### ------------------------------------------------------------------------ ###
### calculate for short- (year 1-5), medium- (year 6-10) and 
### long-term (year 11-20)
### risk 1: proportion of stock below Blim, average over iterations and period
### risk 3: maximum of annual proportions
### catch: mean catch in period (mean over years and iterations)
### iav: inter-annual variation of catch, average over years and iterations


### catch mean
stats_new$catch_long <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2029))
}
stats_new$catch_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2019, end = 2023))
}
stats_new$catch_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2024, end = 2028))
}
### catch median
stats_new$catch_median_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2029))
}
stats_new$catch_median_short <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2019, end = 2023))
}
stats_new$catch_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2024, end = 2028))
}
### risk 1
stats_new$risk1_full <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019) < Blim)
}
stats_new$risk1_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2029) < Blim)
}
stats_new$risk1_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019, end = 2023) < Blim)
}
stats_new$risk1_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2024, end = 2028) < Blim)
}
#risk 3
stats_new$risk3_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2029) < Blim))
}
stats_new$risk3_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < Blim))
}
stats_new$risk3_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) <Blim))
}

### inter-annual variation of catch
stats_new$iav_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2028)), summary_all = median)
}
stats_new$iav_short <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2018, end = 2023)), 
      summary_all = median)
}
stats_new$iav_medium <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2023, end = 2028)), 
      summary_all = median)
}


### inter-annual variation of catch
stats_new$iavTAC_long <- foreach(x = res_list, .packages = "FLCore",
                              .combine = "c") %dopar% {
  iav(object = window(x@tracking["metric.is"], start = 2028, end = 2037), 
      summary_all = median)
}
stats_new$iavTAC_short <- foreach(x = res_list, .packages = "FLCore",
                               .combine = "c") %dopar% {
  iav(object = window(x@tracking["metric.is"], start = 2018, end = 2022), 
      summary_all = median)
}
stats_new$iavTAC_medium <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = window(x@tracking["metric.is"], start = 2023, end = 2027), 
      summary_all = mean)
}



### SSB median
stats_new$ssb_median_long <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2029))
}
stats_new$ssb_median_short <- foreach(x = res_list, .packages = "FLCore",
                                        .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2019, end = 2023))
}
stats_new$ssb_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                         .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2024, end = 2028))
}

### time to recovery
MSYBtrigger <- 166708
stats_new$recovery_proportion <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
  mean(apply(window(ssb(x@stock), start = 2019) >= MSYBtrigger, 6, max))
}
stats_new$recovery_time <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
  median(apply(window(ssb(x@stock), start = 2019)@.Data >= MSYBtrigger, 6, 
               function(x) {
    if (any(x)) {
      which(x)[1]
    } else {
      Inf
    }
  }))
}
### SAM convergence
stats_new$conv_failed <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  sum(x@tracking["conv.est", ac(2018:2037)] != 0)
}
all(stats_new$conv_failed == 0)
### check F maxed (2)
stats_new$F_maxed <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
  sum(window(fbar(stock(x)), start = 2019) >= 2)
}
# realized F
stats_new$F_median_long <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
  median(window(fbar(x@stock), start = 2029))
}

stats_new$F_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
  median(window(fbar(x@stock),  start = 2024, end = 2028))
}
stats_new$F_median_short<- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
  median(window(fbar(x@stock),start = 2019, end = 2023) )
}

#pos <- which(stats_new$F_maxed != 0)
#stats_new[pos, ]
#plot(stock(res_list[[pos]]))
#pos_iter <- which(apply(fbar(stock(res_list[[pos]])), c(1, 6), max) >= 2)
#plot(stock(res_list[[pos]])[,,,,, pos_iter])
### F maxed in THREE scenarios, 
### in each of them in ONE iteration and ONE time only
### HCR B, Btrigger = 110000, Ftrgt = 0.5, iteration 29
### HCR AD, Btrigger = 170000, Ftrgt = 0.5, iteration 442
### HCR AD, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR BE, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR CE, Btrigger = 190000, Ftrgt = 0.5, iteration 442

stats <- rbind(stats, stats_new)
stats$Blim<-119970
stats <- stats[order(stats$file), ]
saveRDS(object = stats, file = paste0(path_res, "stats.rds"))
write.csv(x = stats, file = paste0(path_res, "stats.csv"), row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE, highlight_max=FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_median_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$ssb <- dat[, paste0("ssb_median_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  dat$Btrigger<-dat$Btrigger/1000
# find yield maximum
dat_max <- dat %>% 
    filter(risk <= 0.05) %>%
    filter(catch == max(catch)) %>%
    select(Ftrgt, Btrigger)
 

 p1 <- ggplot() +
    geom_raster(data = dat %>% 
                  filter(risk <= 0.05) %>%
                  filter(catch >= 0.95 * max(catch)),
                aes(x = Btrigger, y = Ftrgt, fill = catch)) +
    scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",

                        high = "green") +
    geom_text(data = dat, 
              aes(x = Btrigger, y = Ftrgt, 
                  label = round(catch), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, "-term catch [t]")) +
    scale_x_continuous(breaks = c(seq(from = 100, to = 280, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
    ### risk
  p2 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 3"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()+
 facet_wrap(~ paste0(time, "-term risk 3")) +
    scale_x_continuous(breaks = c(seq(from = 100, to = 280, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### iav
  p3 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = iav)) +
    geom_text(aes(label = round(iav, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\ninter-annual\ncatch variability"),
                        low = "green", high = "red") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term inter-annual\ catch variability")) +
    scale_x_continuous(breaks = c(seq(from = 100, to = 280, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### SSB
  p4 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = ssb), alpha = 0.75) +
    geom_text(aes(label = round(ssb), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\nSSB [t]"),
                        low = "red", high = "green") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term SSB [t]")) +
    scale_x_continuous(breaks = c(seq(from = 100, to = 280, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk1
  p5 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk1), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk1 <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 1")) +
    scale_x_continuous(breaks = c(seq(from = 100, to = 280, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))


### highlight maximum
  if (isTRUE(highlight_max)) {
    p_add <- geom_tile(data = dat_max, aes(x = Btrigger, y = Ftrgt),
                width = 10, height = 0.01, linetype = "solid",
                alpha = 0, colour = "black", size = 0.3)
    p1 <- p1 + p_add
    p2 <- p2 + p_add
    p3 <- p3 + p_add
    p4 <- p4 + p_add
    p5 <- p5 + p_add
  }
  
  if (isTRUE(add_risk1)) {
    plot_grid(p1, p2, p3, p5, nrow = 2, ncol = 2, align = "hv")
  } else {
    plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")


  }
  
  
}

highlight_max<-TRUE

### A
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == FALSE & TACconstr == FALSE & Ftrgt>0,
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_A_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == FALSE & TACconstr == FALSE & Ftrgt>0,
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_A_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == FALSE & TACconstr == FALSE & Ftrgt>0,
      Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "short")
ggsave(filename =paste0("output/runs/whg4/",iters,"_",years,"/grid_A_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "long", add_risk1 = FALSE,highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_B_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_B_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "short")
ggsave(filename =paste0("output/runs/whg4/",iters,"_",years,"/grid_B_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "long", add_risk1 = FALSE,highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_C_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_C_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == FALSE & TACconstr == FALSE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "short")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_C_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### AD
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "long", add_risk1 = FALSE,highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_AD_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == TRUE & TACconstr == TRUE &
      Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_AD_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "A" & BB == TRUE & TACconstr == TRUE &
      Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "short")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_AD_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### BE
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "long", add_risk1 = FALSE,highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_BE_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_BE_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "B" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "short")
ggsave(filename =paste0("output/runs/whg4/",iters,"_",years,"/grid_BE_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### CE
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "long", add_risk1 = FALSE,highlight_max<-TRUE)
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_CE_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "medium")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_CE_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(OM==1,HCR == "C" & BB == TRUE & TACconstr == TRUE &
       Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "short")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/grid_CE_short.png"), 
width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

  
  
