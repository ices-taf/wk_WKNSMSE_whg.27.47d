### load packages
rm()

#devtools::install_github("flr/ggplotFL")

library(FLfse)
library(stockassessment)
library(FLash)
library(FLAssess)
library(mse)
library(FLCore)
library(ggplot2)
library(ggplotFL)
library(tidyr)
library(cowplot)
library(dplyr)

### load files from package mse for easier debugging
# devtools::load_all("../mse/")


setwd(paste0("/home/miethet/MSE_whiting/"))

source("a4a_mse_WKNSMSE_funs.R")


  library(doParallel)
  cl <- makeCluster(1)
  registerDoParallel(cl)
  cl_length <- length(cl)


. <- foreach(i = seq(cl_length)) %dopar% {
  #devtools::load_all("../mse/")
  setwd(paste("/home/miethet/MSE_whiting", sep=""))

  source("a4a_mse_WKNSMSE_funs.R")
}


iters<-1000
years<-20
Blim<-119970
MSYbtrigger<-166708
Flim<-0.458
Fpa<-0.33
Fmsy<-0.172

datastats<-read.csv(paste0("output/runs/whg4/",iters,"_",years,"/stats.csv"))


data_long<-datastats[,c("OM","Ftrgt","Btrigger","HCR","TACconstr","BB","catch_median_long","ssb_median_long","iav_long","risk3_long","risk1_long","F_median_long","conv_failed","F_maxed")]
data_medium<-datastats[,c("OM","Ftrgt","Btrigger","HCR","TACconstr","BB","catch_median_medium","ssb_median_medium","iav_medium","risk3_medium","risk1_medium","F_median_medium","conv_failed","F_maxed")]
data_short<-datastats[,c("OM","Ftrgt","Btrigger","HCR","TACconstr","BB","catch_median_short","ssb_median_short","iav_short","risk3_short","risk1_short","F_median_short","conv_failed","F_maxed")]

data_medium[,c("catch_median_medium","ssb_median_medium")]<-round(data_medium[,c("catch_median_medium","ssb_median_medium")])
data_medium[,c("iav_medium","risk3_medium","risk1_medium","F_median_medium")]<-round(data_medium[,c("iav_medium","risk3_medium","risk1_medium","F_median_medium")],3)

data_short[,c("catch_median_short","ssb_median_short")]<-round(data_short[,c("catch_median_short","ssb_median_short")])
data_short[,c("iav_short","risk3_short","risk1_short","F_median_short")]<-round(data_short[,c("iav_short","risk3_short","risk1_short","F_median_short")],3)

data_long[,c("catch_median_long","ssb_median_long")]<-round(data_long[,c("catch_median_long","ssb_median_long")])
data_long[,c("iav_long","risk3_long","risk1_long","F_median_long")]<-round(data_long[,c("iav_long","risk3_long","risk1_long","F_median_long")],3)


write.csv(x = data_long, file = paste0("output/runs/whg4/",iters,"_",years,"/stats_long.csv"), row.names = FALSE)
write.csv(x = data_medium, file = paste0("output/runs/whg4/",iters,"_",years,"/stats_medium.csv"), row.names = FALSE)
write.csv(x = data_short, file = paste0("output/runs/whg4/",iters,"_",years,"/stats_short.csv"), row.names = FALSE)



stats<-data_long
table1<-NULL

for(j in 1:3){
combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"), 
                    OM=j,  
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)  


combs <- merge(combs, stats, all.x = TRUE)

table1<-rbind(table1,combs)

}
write.csv(x = data_long, file = paste0("output/runs/whg4/",iters,"_",years,"/HCR_comb_stats_long.csv"), row.names = FALSE)


stats<-data_short
table1<-NULL

for(j in 1:3){
combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"), 
                    OM=j,  
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)  


combs <- merge(combs, stats, all.x = TRUE)

table1<-rbind(table1,combs)

}
write.csv(x = data_short, file = paste0("output/runs/whg4/",iters,"_",years,"/HCR_comb_stats_short.csv"), row.names = FALSE)





om_opt<-"OM3" # vary OM

#combs <- data.frame(name = c("A", "B", "C", "AD", "BE", "CE"),
#                    HCR = c("A", "B", "C", "A", "B", "C"),
#                    BB = c(rep(FALSE, 3), TRUE, TRUE, TRUE),
#                    TACconstr = c(rep(FALSE, 3), TRUE, TRUE, TRUE),
#                    Btrigger = c(190000, 200000, 190000, 180000, 190000,230000),
#                   Ftrgt = c(0.12, 0.16, 0.12, 0.12, 0.15, 0.15),
#                    scenario = 0)
                    
# includes F=0                    
combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"),   
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)                  
 
  
                    
for(i in 1:8){

path1<-paste0(om_opt,
              "_HCR-", combs$HCR[i],
              "_Ftrgt-", combs$Ftrgt[i],
              "_Btrigger-", combs$Btrigger[i],
              "_TACconstr-", combs$TACconstr[i],
              "_BB-", combs$BB[i])

output <- readRDS(paste0("output/runs/whg4/",iters,"_",years,"/", path1,".rds"))

#plot(output)


original <- readRDS(paste0("input/whg4/",iters,"_",years,"/stk.rds"))

new1<-original
new1[, ac(2018:2038)]<-output@stock[,ac(2018:2038)]


#
### save
 saveRDS(object = new1, file = paste0("output/runs/whg4/", iters, "_", years,
                                     "/",path1,"_base_full_stk.rds"))
#### plot
 plot(new1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + 
  xlab("year") + 
  geom_hline(data = data.frame(qname = "SSB", data = Blim),
              aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB", data = combs$Btrigger[i]),
              aes(yintercept = data), linetype = "solid", color="darkgrey") +
  geom_hline(data = data.frame(qname = "SSB", data = MSYbtrigger),
              aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data =Flim),
              aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data =Fmsy),
              aes(yintercept = data), linetype = "solid") +
   geom_vline(xintercept = 2018.5) +    
   theme_bw()+ ylim(0, NA) + theme(legend.position = 0) + 
      facet_wrap(~ qname, ncol = 1, strip.position = "right", scales = "free_y",
                 labeller = as_labeller(c(
                   "Rec" = "Rec [1000]",
                   "SSB" = "SSB [t]",
                   "Catch" = "Catch [t]",
                   "F" = paste0("F (ages ", paste(range(new1)["minfbar"], 
                                      range(new1)["maxfbar"], sep = "-"),
                                ")"))))
ggsave(filename = paste0("output/runs/whg4/", iters, "_", years, "/plots/",path1,"_base_full_stk.png"), 
        width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")



    ### plot iterations
    stk_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(new1),
                                      `SSB [t]` = ssb(new1),
                                      `Catch [t]` = catch(new1),
                                      `F` = fbar(new1)))
    p <- ggplot(data = stk_df, 
           aes(x = year, y = data, group = iter)) +
      geom_line(alpha = 0.025) +
      facet_wrap(~ qname, ncol = 1, strip.position = "right",
                 scale = "free_y") +
      theme_bw() +
      ylim(c(0, NA)) + labs(y = "") + 
      geom_vline(xintercept = 2018.5) +
      geom_hline(data = data.frame(qname = "SSB [t]", data = Blim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "SSB [t]", data = MSYbtrigger),
                 aes(yintercept = data), linetype = "solid") +
      geom_hline(data = data.frame(qname = "F", data = Flim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "F", data = Fmsy),
                 aes(yintercept = data), linetype = "solid")
    
  
  print(p)
  ### save plot
      filename <- paste0("output/runs/whg4/", iters, "_", years, "/plots/stock_plots/",path1,".png")
    ggsave(filename = filename,
           width = 30, height = 20, units = "cm", dpi = 300, 
           type = "cairo")
    
  }



### base OM
### current HCR

combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"),   
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)                  
 combs <- rbind(cbind(combs, OM = 1),
               cbind(combs, OM = 2),
               cbind(combs, OM = 3)
               )

combs <- combs[-c(10, 18), ]  
                    
for(i in 1:length(combs[,1])){

path1<-paste0("OM",combs$OM[i],
              "_HCR-", combs$HCR[i],
              "_Ftrgt-", combs$Ftrgt[i],
              "_Btrigger-", combs$Btrigger[i],
              "_TACconstr-", combs$TACconstr[i],
              "_BB-", combs$BB[i])

output <- readRDS(paste0("output/runs/whg4/",iters,"_",years,"/", path1,".rds"))

#plot(output)


original <- readRDS(paste0("input/whg4/",iters,"_",years,"/stk.rds"))

new1<-original
new1[, ac(2018:2038)]<-output@stock[,ac(2018:2038)]

iters_plot = 1:5
yr_start<-2018.5
plot_iter<-FALSE

 if (!isTRUE(plot_iter)) {
    ### plot percentiles with iterations
    p <- plot(new1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = iters_plot) +
      xlab("year") + geom_vline(xintercept = yr_start) +
        geom_hline(data = data.frame(qname = "SSB", data = Blim),
                   aes(yintercept = data), linetype = "dashed") +
        geom_hline(data = data.frame(qname = "SSB", data = MSYbtrigger),
                   aes(yintercept = data), linetype = "solid") +
        geom_hline(data = data.frame(qname = "F", data = Flim),
                   aes(yintercept = data), linetype = "dashed") +
        geom_hline(data = data.frame(qname = "F", data = Fmsy),
                   aes(yintercept = data), linetype = "solid") +
        theme_bw() + ylim(0, NA) + theme(legend.position = 0) + 
      facet_wrap(~ qname, ncol = 1, strip.position = "right", scales = "free_y",
                 labeller = as_labeller(c(
                   "Rec" = "Rec [1000]",
                   "SSB" = "SSB [t]",
                   "Catch" = "Catch [t]",
                   "F" = paste0("F (ages ", paste(range(new1)["minfbar"], 
                                      range(new1)["maxfbar"], sep = "-"),
                                ")"))))
  } else {
    ### plot all iterations
    stk_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(new1),
                                      `SSB [t]` = ssb(new1),
                                      `Catch [t]` = catch(new1),
                                      `F` = fbar(new1)))
    p <- ggplot(data = stk_df, 
           aes(x = year, y = data, group = iter)) +
      geom_line(alpha = 0.025) +
      facet_wrap(~ qname, ncol = 1, strip.position = "right",
                 scale = "free_y") +
      theme_bw() +
      ylim(c(0, NA)) + labs(y = "") + 
      geom_vline(xintercept = yr_start) +
      geom_hline(data = data.frame(qname = "SSB [t]", data = Blim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "SSB [t]", data = MSYbtrigger),
                 aes(yintercept = data), linetype = "solid") +
      geom_hline(data = data.frame(qname = "F", data = Flim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "F", data = Fmsy),
                 aes(yintercept = data), linetype = "solid")
    
  }

print(p)

    
    filename <-paste0("output/runs/whg4/", iters, "_", years, "/plots/stock_plots/",path1, ifelse(!isTRUE(plot_iter == 0), "", "_iters"),
                       ".png")
    ggsave(filename = filename,
           width = 30, height = 20, units = "cm", dpi = 300, 
           type = "cairo")

}
               








###################################################################################################################################
#F0 
# worm plots
om_opt<-"OM1"

path1<-paste0(om_opt,
              "_HCR-", combs$HCR[1],
              "_Ftrgt-", combs$Ftrgt[1],
              "_Btrigger-", combs$Btrigger[1],
              "_TACconstr-", combs$TACconstr[1],
              "_BB-", combs$BB[1]) 
stkF0 <- readRDS(paste0("output/runs/whg4/",iters,"_",years,"/", path1,"_base_full_stk.rds"))

### iterations
stkF0_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stkF0),
                                  `SSB [t]` = ssb(stkF0),
                                  `Catch [t]` = catch(stkF0),
                                  `F` = fbar(stkF0)))

 
 ggplot(data = stkF0_df[stkF0_df$iter %in% 1:1000, ], 
       aes(x = year, y = data, group = iter)) +
  geom_line(alpha = 0.025) +
  facet_wrap(~ qname, ncol = 1, strip.position = "right",
             scale = "free_y") +
  theme_bw() +
  ylim(c(0, NA)) + labs(y = "") + 
  geom_vline(xintercept = 2018.5) +
  geom_hline(data = data.frame(qname = "SSB [t]", data = Blim),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB [t]", data = MSYbtrigger),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = Flim),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = Fmsy),
             aes(yintercept = data), linetype = "solid")
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/", 
                         "stk_F0_iters_",iters,".png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

 ### ------------------------------------------------------------------------ ###
### summary box plots: compare HCR options for each OM ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations

stats<-datastats
om_<- 1

for(j in 1:3){
combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"), 
                    OM=j,  
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)  


combs <- merge(combs, stats, all.x = TRUE)
combs2 <- gather(data = combs, key = "key", value = "value",
                catch_median_medium, risk3_medium, iav_medium,
                ssb_median_medium, F_median_medium)
combs2$name <- factor(combs2$name, levels = c("F0","A*","A", "B", "C", "AD", "BE", "CE"))



ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/",j, 
                         "compare_HCR_medium.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


combs2 <- gather(data = combs, key = "key", value = "value",
                catch_median_short, risk3_short, iav_short,
                ssb_median_short, F_median_short)
combs2$name <- factor(combs2$name, levels = c("F0","A*","A", "B", "C", "AD", "BE", "CE"))


ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/",j, 
                         "compare_HCR_short.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


combs2 <- gather(data = combs, key = "key", value = "value",
                catch_median_long, risk3_long, iav_long,
                ssb_median_long, F_median_long, recovery_proportion, recovery_time)
combs2$name <- factor(combs2$name, levels = c("F0","A*","A", "B", "C", "AD", "BE", "CE"))


ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/",j, 
                         "compare_HCR_long.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

}
       
### load entire distribution for stats
stats_full <- function(data) {
  #combs_full <- foreach(i = split(data, seq(nrow(data))), 
  #                    .packages = "FLCore", .combine = rbind) %dopar% {
  combs_full<-NULL
  
  for(j in 1:length(data$file)){                      
    
    stk_i <- readRDS(paste0("output/runs/whg4/1000_20/", data$file[j]))
    i<-data[j,]
    MSYBtrigger <- 166708
    Blim<-119970
    res <- rbind(
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_long",
               value = c(window(catch(stk_i@stock), start = 2029))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_medium",
               value = c(window(catch(stk_i@stock), start = 2019, end = 2023))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_short",
               value = c(window(catch(stk_i@stock), start = 2024, end = 2028))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_long",
               value = mean(window(ssb(stk_i@stock), start = 2029) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_medium",
               value = mean(window(ssb(stk_i@stock), 
                                   start = 2024, end = 2028) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_short",
               value = mean(window(ssb(stk_i@stock), 
                                   start = 2019, end = 2023) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_long",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2029) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_medium",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2024, end = 2028) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_short",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2019, end = 2023) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_long",
               value = c(iav(object = catch(window(stock(stk_i), start = 2028))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_medium",
               value = c(iav(object = catch(window(stock(stk_i), 
                                                   start = 2023, end = 2028))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_short",
               value = c(iav(object = catch(window(stock(stk_i), 
                                                   start = 2018, end = 2023))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_long",
               value = c(window(ssb(stk_i@stock), start = 2029))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_medium",
               value = c(window(ssb(stk_i@stock), start = 2024, end = 2028))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_short",
               value = c(window(ssb(stk_i@stock), start = 2019, end = 2023))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "f_long",
               value = c(window(fbar(stk_i@stock), start = 2029))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "f_medium",
               value = c(window(fbar(stk_i@stock), start = 2024, end = 2028))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "f_short",
               value = c(window(fbar(stk_i@stock), start = 2019, end = 2023))),          
               
    data.frame(name = i$name, scenario = i$scenario,
               key = "recovery_proportion",
               value = mean(apply(window(ssb(stk_i@stock), 
                                         start = 2019) >= MSYBtrigger, 6, max))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "recovery_time",
               value = c(apply(window(ssb(stk_i@stock), 
                                      start = 2019)@.Data >= MSYBtrigger, 6,
                               function(x) {
                                 if (any(x)) {which(x)[1]} else {Inf}})))
  )
  if (i$name == "F0") {
    res$value[res$key %in% c("catch_long", "catch_medium", "catch_short",
                             "iav_long", "iav_medium", "iav_short","f_short","f_medium","f_long")] <- 0
  }
  res <- merge(res, i[, c("name", "OM", "HCR", "BB", "TACconstr", "Btrigger",
                          "Ftrgt")])
                          
 combs_full<-rbind(combs_full,res)  
  }
  combs_full$name <- factor(combs_full$name, 
                            levels = c("F0","A*", "A", "B", "C", "AD", "BE", "CE"))
  return(combs_full)
}
############################################################################
### base OM
#iav can be large for catches clos to zero remove
# box plots for OM1


combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"), 
                    OM=1,  
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0,0.172,0.14, 0.16, 0.14, 0.16, 0.16, 0.15),
                    scenario = 0)  

combs <- merge(combs, stats, all.x = TRUE)
combs_base <- stats_full(data = combs)


ggplot(data = combs_base, 
       mapping = aes(x = name, y = value, group = name)) +
       
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)
  
#ggsave(filename = paste0("output/runs/whg4/",iters,"_", years,"/plots/", 
#                         "compare_HCR_boxplots.png"),
#       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

################################################################################################################
# to run for original Ftarget/Btrigger pair
### get median for option A*

combs_base <- left_join(combs_base, 
                         combs_base %>%
  group_by(key, OM, name) %>%
  summarise(value_median = median(value)) %>%
  filter(name == "A*") %>%
    select(-name))
p_catch_long <- ggplot(data = combs_base[combs_base$key == "catch_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "long-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_medium <- ggplot(data = combs_base[combs_base$key == "catch_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "medium-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_short <- ggplot(data = combs_base[combs_base$key == "catch_short", ], 
                         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "short-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_risk1_long <- ggplot(data = combs_base[combs_base$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_blank(data = combs_base[combs_base$key == "risk3_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_base[combs_base$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_base[combs_base$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_base[combs_base$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_base[combs_base$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_base[combs_base$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_base[combs_base$key == "iav_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 3)) + 
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_base[combs_base$key == "iav_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 3)) + 
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_base[combs_base$key == "iav_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 3)) + 
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_base[combs_base$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.75e+06)) +
  labs(x = "", y = "long-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2))
p_ssb_medium <- ggplot(data = combs_base[combs_base$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.7e+06)) +
  labs(x = "", y = "medium-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_ssb_short <- ggplot(data = combs_base[combs_base$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.6e+06)) +
  labs(x = "", y = "short-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))

p_f_long <- ggplot(data = combs_base[combs_base$key == "f_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "", y = "long-term F") 
  
p_f_medium <- ggplot(data = combs_base[combs_base$key == "f_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "", y = "medium-term F") 
  
p_f_short <- ggplot(data = combs_base[combs_base$key == "f_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "", y = "short-term F")  

p_recovery_proportion <- 
  ggplot(data = combs_base[combs_base$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_base[combs_base$key == "recovery_time", ], 
         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "recovery time [years]")


plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long, p_ssb_long,p_f_long,
          align = "hv")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/plots/baseOM_stats/", 
                         "summary_baseOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium, 
          p_ssb_medium,p_f_medium,
          align = "hv")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/plots/baseOM_stats/", 
                         "summary_baseOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short, p_ssb_short,p_f_short,
          align = "hv")
ggsave(filename = paste0("output/runs/whg4/",iters,"_",years,"/plots/baseOM_stats/", 
                         "summary_baseOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_recovery_proportion, p_recovery_time,
          align = "hv")
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### summary plots: base OM, additional scenarios around maximum yield ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
### add: 0.9 & 1.1 * Ftrgt
###      Fmsylower, Fmsyupper
stats<-datastats
combs <- data.frame(name = rep(c("A", "B", "C", "AD", "BE", "CE"), each = 5),
                    HCR = rep(c("A", "B", "C", "A", "B", "C"), each = 5),
                    BB = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    TACconstr = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    Btrigger = rep(c(220000, 200000, 220000, 250000, 210000, 230000), each = 5),
                    Ftrgt = c(0.14 * c(0.9, 1, 1.1), 0.158, 0.172,
                              0.16 * c(0.9, 1, 1.1), 0.158, 0.172,
                              0.14 * c(0.9, 1, 1.1), 0.158, 0.172,
                              0.16 * c(0.9, 1, 1.1), 0.158, 0.172,
                              0.16 * c(0.9, 1, 1.1), 0.158, 0.172,
                              0.15 * c(0.9, 1, 1.1), 0.158, 0.172),
                    scenario = c("0.9*Ftrgt", "Ftrgt", "1.1*Ftrgt",
                                 "Fmsylower", "Fmsyupper"),
                    OM = 1)
combs <- merge(combs, stats)
combs_dat <- stats_full(data = combs)
combs_dat$scenario <- factor(combs_dat$scenario, 
                             levels = c("Fmsylower", "0.9*Ftrgt", "Ftrgt", 
                                        "1.1*Ftrgt", "Fmsyupper"))
ggplot(data = combs_dat, 
       mapping = aes(x = name, y = value, group = interaction(scenario, name), 
                     colour = scenario)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_dat[combs_dat$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]")
p_catch_medium <- ggplot(data = combs_dat[combs_dat$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]")
p_catch_short <- ggplot(data = combs_dat[combs_dat$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]")
p_risk1_long <- ggplot(data = combs_dat[combs_dat$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_dat[combs_dat$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_dat[combs_dat$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_dat[combs_dat$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_dat[combs_dat$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_dat[combs_dat$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_dat[combs_dat$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_dat[combs_dat$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_dat[combs_dat$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_dat[combs_dat$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 5.5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_medium <- ggplot(data = combs_dat[combs_dat$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  coord_cartesian(ylim = c(0, 5.5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_short <- ggplot(data = combs_dat[combs_dat$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 5.5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
                     
 p_f_long <- ggplot(data = combs_dat[combs_dat$key == "f_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term F") +
  theme(legend.direction = "horizontal")
  p_f_medium <- ggplot(data = combs_dat[combs_dat$key == "f_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term F") +
  theme(legend.direction = "horizontal")
  p_f_short <- ggplot(data = combs_dat[combs_dat$key == "f_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term F") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))                    
                     
                     
p_recovery_proportion <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = scenario, colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long, p_f_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_f_long), nrow = 2, rel_heights = c(1, 0.1))
          
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
       
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium, p_f_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_f_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short, p_f_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_f_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")



### ------------------------------------------------------------------------ ###
### summary plots: compare alternative OMs ####
### ------------------------------------------------------------------------ ###

### alternative OMs
### select maximum yield combinations
stats<-datastats

combs <- data.frame(name = c("F0","A*","A", "B", "C", "AD", "BE", "CE"),   
                    HCR = c("A","A","A", "B", "C", "A", "B", "C"),

                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(MSYbtrigger,MSYbtrigger,220000, 200000, 220000, 250000, 210000,230000),
                    Ftrgt = c(0, 0.172, 0.14, 0.16, 0.14, 0.16, 0.16, 0.15),

                    scenario = 0)  

                  
                 
                    
combs_alt <- rbind(cbind(combs, OM = "1"),
                   cbind(combs, OM = "2"),
                   cbind(combs, OM = "3")
                  )
combs_alt <- merge(combs_alt, stats)
combs_alt <- stats_full(data = combs_alt)
ggplot(data = combs_alt, 
       mapping = aes(x = name, y = value, group = interaction(OM, name), 
                     colour = OM)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

#ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
#                         "summary_altOM.png"), 
#       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
       


p_catch_long <- ggplot(data = combs_alt[combs_alt$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]") +
  coord_cartesian(ylim = c(0, 1e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_catch_medium <- ggplot(data = combs_alt[combs_alt$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]") +
  coord_cartesian(ylim = c(0, 1e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_catch_short <- ggplot(data = combs_alt[combs_alt$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]") +
  coord_cartesian(ylim = c(0, 1e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_risk1_long <- ggplot(data = combs_alt[combs_alt$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = OM, 
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_medium <- ggplot(data = combs_alt[combs_alt$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_short <- ggplot(data = combs_alt[combs_alt$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_alt[combs_alt$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = OM,
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_medium <- ggplot(data = combs_alt[combs_alt$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_short <- ggplot(data = combs_alt[combs_alt$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_alt[combs_alt$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_alt[combs_alt$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_alt[combs_alt$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_alt[combs_alt$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  
  coord_cartesian(ylim = c(0, 1e+06))
p_ssb_medium <- ggplot(data = combs_alt[combs_alt$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  coord_cartesian(ylim = c(0, 1e+06))
p_ssb_short <- ggplot(data = combs_alt[combs_alt$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA)) +
  coord_cartesian(ylim = c(0, 6e+05))
  
p_f_long <- ggplot(data = combs_alt[combs_alt$key == "f_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term F") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 0.5))
  
p_f_medium <- ggplot(data = combs_alt[combs_alt$key == "f_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term F") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 0.5))
  
p_f_short <- ggplot(data = combs_alt[combs_alt$key == "f_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term F") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA)) +
  coord_cartesian(ylim = c(0, 0.5))
  
p_recovery_proportion <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = OM, colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long, p_f_long + theme(legend.position = "none"),
                    align = "hv")  ,
          get_legend(p_f_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
       
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium, p_f_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_f_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short, p_f_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_f_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
          align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### compare OM SSB/F with MP SSB/F ####
### ------------------------------------------------------------------------ ###

df <- foreach(OM = c(1, 2, 3),
              .combine = rbind) %do% {
  res <- readRDS(paste0("output/runs/whg4/1000_20/OM", OM, "_HCR-A_Ftrgt-0.14",
                        "_Btrigger-220000_TACconstr-FALSE_BB-FALSE.rds"))
  res_input <- readRDS(paste0("input/whg4/base_run/OM",OM,"_base_run1000.rds"))$om@stock
  tmp <- FLQuant(NA, dimnames = list(
    metric = c("F.om", "SSB.om", "F.est", "SSB.est", "F", "SSB"), 
    year = dimnames(res_input)$year,
    iter = dimnames(res_input)$iter))
  tmp["F.est", ac(an(dimnames(res@tracking)$year) - 1)] <- 
    res@tracking["F.est"]
  tmp["SSB.est", ac(an(dimnames(res@tracking)$year) - 1)] <- 
    res@tracking["B.est"]
  tmp["SSB.om", dimnames(res_input)$year] <- ssb(res_input)
  tmp["F.om", dimnames(res_input)$year] <- fbar(res_input)
  tmp["SSB.om", dimnames(res@stock)$year] <- ssb(res@stock)
  tmp["F.om", dimnames(res@stock)$year] <- fbar(res@stock)
  tmp["SSB"] <- tmp["SSB.est"] / tmp["SSB.om"]
  tmp["F"] <- tmp["F.est"] / tmp["F.om"]
  tmp_out <- cbind(as.data.frame(tmp), OM = OM)
  tmp_out <- tmp_out[!is.na(tmp_out$data) & is.finite(tmp_out$data), ]
  tmp_out
  
}
df <- df %>% group_by(metric, year, OM) %>%
  summarise(X0.05 = quantile(data, probs = 0.05),
            X0.25 = quantile(data, probs = 0.25),
            X0.50 = quantile(data, probs = 0.50),
            X0.75 = quantile(data, probs = 0.75),
            X0.95 = quantile(data, probs = 0.95))
ggplot(data = df %>% filter(metric %in% c("F", "SSB")),
       aes(x = year, y = X0.50)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3, fill = "#F8766D",
              show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, fill = "#F8766D",
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_grid(OM ~ metric) +
  theme_bw() +
  geom_hline(yintercept = 1, alpha = 0.5) +
  labs(y = "MP value / OM value")
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/altOMs_stats/", 
                         "MP_vs_OM.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### plot risk over time ####
### ------------------------------------------------------------------------ ###

### get optimized A
stkA_file <- datastats %>% filter(OM == "1" & Ftrgt == 0.14 & Btrigger == 220000 &
                   TACconstr == FALSE & BB == FALSE & HCR == "A")
### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/runs/whg4/1000_20/", 
                               stkA_file$file))@stock)
### get historical SSB
ssbA <- ssb(readRDS("input/whg4/base_run/OM1_base_run1000.rds")$om@stock)
### combine
ssbA[, dimnames(ssbA_new)$year] <- ssbA_new
### calculate annual risk
riskA <- apply((ssbA < stkA_file$Blim), 2, mean)
### plot
ggplot(data = as.data.frame(window(riskA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "p(SSB<Blim)") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")
ggsave(filename = paste0("output/runs/whg4/1000_20/plots/stock_plots/", 
                         "risk_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")




# script to plot up HCR rules with different Ftrgt and Btrigger values

###---------------------------------------------------###
### set HCR rules
###---------------------------------------------------###

HCR_rule<-function(ssb, data,Blim){
  
  # data should be data frame of sceanrio name, HCR rule to use, Btrigger and Ftrgt
  Btrigger<-data$Btrigger
  Ftrgt<-data$Ftrgt
  rule<-data$HCR
  
  if(rule %in% "A"){
    # A
    HCR<-function(ssb,Btrigger,Ftrgt,Blim){
      #fit striaght line through origin
      slope<-coef(lm(y~x,data=data.frame(x=c(0,Btrigger),y=c(0,Ftrgt))))["x"]
      ifelse(ssb < Btrigger,slope*ssb , Ftrgt)
    }
  }
  if(rule %in% "B"){
    #B
    HCR<-function(ssb,Btrigger,Ftrgt,Blim){
      #fit striaght line through origin
      slope<-coef(lm(y~x,data=data.frame(x=c(0,Btrigger),y=c(0,Ftrgt))))["x"]
      ifelse(ssb < Btrigger,ifelse(ssb<Blim,0.25*Ftrgt,slope*ssb),Ftrgt)
    }
  }
  
  #C
  if(rule %in% "C"){
    HCR<-function(ssb,Btrigger,Ftrgt,Blim){
      #fit striaght line through origin
      slope<-coef(lm(y~x,data=data.frame(x=c(0,Btrigger),y=c(0,Ftrgt))))["x"]
      A<-ifelse(ssb < Btrigger,slope*ssb , Ftrgt)
      B<-ifelse(ssb < Btrigger,ifelse(ssb < Blim,0.25*Ftrgt,slope*ssb),Ftrgt)
      ifelse(A>B,A,B)
    }
  }
  return(HCR(ssb,Btrigger,Ftrgt,Blim))
}
###--------------------------------------------###
### Set F and B pairs
###--------------------------------------------###

# ref points
refpts <- list(msyBtrigger = 166708,
               Blim=119970,
               Flim=0.458,
               Fmsy = 0.172,
               Fpa = 0.33,
               Bpa = 166708)

# optimum pairs
combs <- data.frame(name = c("A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c("A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(refpts$msyBtrigger, 220000, 200000,220000, 250000, 210000, 230000),
                    Ftrgt = c(refpts$Fmsy, 0.14, 0.16, 0.14, 0.16, 0.16, 0.15))

# ssb range for stock
ssb_range<-sort(c(refpts$msyBtrigger,refpts$Blim,refpts$Blim-1,seq(0,300000,10000)))


###------------------------------------------###
### plot rules 
###------------------------------------------###

A_star<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="A*",],Blim<-refpts$Blim)
A<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="A",],Blim<-refpts$Blim)
B<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="B",],Blim<-refpts$Blim)
C<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="C",],Blim<-refpts$Blim)
AD<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="AD",],Blim<-refpts$Blim)
BE<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="BE",],Blim<-refpts$Blim)
CE<-HCR_rule(ssb=ssb_range,data=combs[combs$name=="CE",],Blim<-refpts$Blim)

ylims<-c(0,max(combs$Ftrgt))

#windows()


tiff(paste0("output/final_optHCRs.tiff"), bg="white",  res=200, width = 1200, height = 1100)
par(mar=c(4,4,4,1))


plot(ssb_range,A_star,type="l",col="black", lwd=2,ylim=ylims,xlab="SSB (tonnes)",ylab="Ftrgt", main="HCR scenarios")
lines(ssb_range,A,col="blue",lwd=2)
lines(ssb_range,B,col="red",lwd=2)
lines(ssb_range,C,col="green",lwd=2)
#lines(ssb_range,AD,col="blue",lty="dashed",lwd=2)
#lines(ssb_range,BE,col="red",lty="dashed",lwd=2)
#lines(ssb_range,CE,col="green",lty="dashed",lwd=2)

legend("bottomright",inset=0.02,legend=c("A*","A","B","C"),col=c("black",rep(c("blue","red","green"),1)),
       lty=c(rep(1,4),rep(2,3)),cex=0.8,lwd=2)

dev.off()
#



