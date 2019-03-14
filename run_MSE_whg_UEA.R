### ------------------------------------------------------------------------ ###
### R script to run WKNSMSE cod MSE on HPC ####
### ------------------------------------------------------------------------ ###
### This is designed to be called by a job submission script
### run_mse.qsub for systems using PBS and the qsub commands
### run_mse.bsub for system using LSF and the bsub commands


### ------------------------------------------------------------------------ ###
### load arguments from job script ####
### ------------------------------------------------------------------------ ###

### load arguments
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  
  ### parallelisation environment
  if (!exists("par_env")) par_env <- 1
  if (!exists("n_workers")) n_workers <- 1
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
### load files from package mse for easier debugging
#devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

### load additional functions
#setwd(paste("/home/miethet/MSE_whiting", sep=""))
source("a4a_mse_WKNSMSE_funs_whg.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###
### par_env=1 -> MPI (Rmpi, DoMPI)
### par_env=2 -> DoParallel


if (par_env == 1) {
  
  library(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
  cl_length <- cl$workerCount
  
} else if (par_env == 2) {
  
  library(doParallel)
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  cl_length <- length(cl)
  
}



### load packages and functions into workers
. <- foreach(i = seq(cl_length)) %dopar% {
  #devtools::load_all("../mse/")
  library(mse)
  library(FLash)
  library(FLfse)
  library(stockassessment)
  library(foreach)
  library(doRNG)
  #setwd(paste("/home/miethet/MSE_whiting", sep=""))
  
  source("a4a_mse_WKNSMSE_funs_whg.R")
}

### set random seed for reproducibility
library(doRNG)
registerDoRNG(123)

### ------------------------------------------------------------------------ ###
### load data for MSE ####
### ------------------------------------------------------------------------ ###

### data path
path_data <- paste0("input/whg4/")



### load input objects
input <- readRDS(paste0(path_data,"base_run/OM",om_opt,"_base_run",iters,".rds"))

### modify input for running in parallel
input$genArgs$nblocks <- nblocks

### ------------------------------------------------------------------------ ###
### set up HCR & options ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### set HCR option: A, B, C


if (exists("HCRoption")) {
  
  input$ctrl.mp$ctrl.hcr@args$option <- switch(HCRoption, 
                                               "1" = "A", 
                                               "2" = "B", 
                                               "3" = "C",
                                               "4" = "A",
                                               "5" = "B",
                                               "6" = "C")
  
  cat(paste0("\nSetting custom HCR option: HCRoption = ", HCRoption, 
             " => HCR ", input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
  
} else {
  
  cat(paste0("\nUsing default HCR option: HCR ", 
             input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
  HCRoption <- 0
  
}


### ------------------------------------------------------------------------ ###
### set HCR parameters 
if (HCRoption %in% 1:6) {

  
 ### additional combinations after finding yield maximum
 
  comb_max <- switch(HCRoption, 
                     "1" = c(220000, 0.14), 
                     "2" = c(200000, 0.16), 
                     "3" = c(220000, 0.14),
                     "4" = c(220000, 0.14),
                     "5" = c(210000, 0.16),
                     "6" = c(230000, 0.15))
  hcr_vals <- expand.grid(Ftrgt = c(comb_max[2]*0.9, comb_max[2]*1.1, 0.158, 0.172), Btrigger = comb_max[1])


# for alternative OM 2 and OM3
if(om_opt>1) hcr_vals <- expand.grid(Ftrgt = comb_max[2], Btrigger = comb_max[1])
#

  
}

#input_bckp <- input
#for (HCR_comb in 1:length(hcr_vals$Ftrgt)) {
#  input <- input_bckp
  
  ### implement
  if (exists("HCR_comb")) {
    
    ### set Btrigger
    Btrigger <- hcr_vals[HCR_comb, "Btrigger"]
    input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
    input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
    
    ### set Ftrgt
    Ftrgt <- hcr_vals[HCR_comb, "Ftrgt"]
    input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
    
    cat(paste0("\nSetting custom Btrigger/Ftrgt values.\n",
               "Using HCR_comb = ", HCR_comb, "\n",
               "Ftrgt = ", Ftrgt, "\n",
               "Btrigger = ", Btrigger, "\n\n"))
    
  } else {
    input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
    input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
    
    cat(paste0("\nUsing default Btrigger/Ftrgt values.\n",
               "Ftrgt = ", input$ctrl.mp$ctrl.phcr@args$Ftrgt, "\n",
               "Btrigger = ", input$ctrl.mp$ctrl.phcr@args$Btrigger, "\n\n"))
    
  }
  
  ### ------------------------------------------------------------------------ ###
  ### TAC constraint
  input$ctrl.mp$ctrl.is@args$TAC_constraint <- FALSE
  ### check conditions
  ### either manually requested or as part of HCR options 4-6 
  if (exists("TAC_constraint")) {
    
    if (isTRUE(as.logical(TAC_constraint))) {
      input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
    }
  }
  if (HCRoption %in% 4:6) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
  }
  ### implement
  if (isTRUE(input$ctrl.mp$ctrl.is@args$TAC_constraint)) {
    
    input$ctrl.mp$ctrl.is@args$lower <- 80
    input$ctrl.mp$ctrl.is@args$upper <- 125
    input$ctrl.mp$ctrl.is@args$Btrigger_cond <- TRUE
    
    cat(paste0("\nImplementing TAC constraint.\n\n"))
    
  } else {
    
    cat(paste0("\nTAC constraint NOT implemented.\n\n"))
    
  }
  
  ### ------------------------------------------------------------------------ ###
  ### banking & borrowing
  if(om_opt==3){ input$iem@args$BB<- FALSE } 
  if(!om_opt==3)input$iem <- NULL
  input$ctrl.mp$ctrl.is@args$BB <- FALSE
  
  ### check conditions
  ### either manually requested by BB=TRUE or as part of HCR options 4-6
  if (exists("BB")) {
    
    if (isTRUE(as.logical(BB))) {
      
      
      if(!om_opt==3)input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE, imp.error=FALSE, imp.res=NULL))
      if(om_opt==3){ input$iem@args$BB<- TRUE }  # with IE
      input$ctrl.mp$ctrl.is@args$BB <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
      
    }
    
    
    if (HCRoption %in% 4:6) {
      
      if(!om_opt==3)input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE, imp.error=FALSE, imp.res=NULL))
      if(om_opt==3){ input$iem@args$BB = TRUE } #with IE
      input$ctrl.mp$ctrl.is@args$BB <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_check_hcr <- FALSE
      input$ctrl.mp$ctrl.is@args$BB_check_fc <- FALSE
      input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
      
      
      if (HCRoption %in% 4) {
        
        input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
        
      } else if (HCRoption %in% 5:6) {
        
        input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
        
      }
      
      
    }
    
    if (input$ctrl.mp$ctrl.is@args$BB) { cat(paste0("\nImplementing banking and borrowing.\n\n")) 
      
      
    } else {
      
      cat(paste0("\nBanking and borrowing NOT implemented.\n\n"))
    }
    
    
  }
  
  if(om_opt==3) cat(paste0("\n Implementing Implementation error (IBC).\n\n"))
  
  
  
  
  ### ------------------------------------------------------------------------ ###
  ### run MSE ####
  ### ------------------------------------------------------------------------ ###
  
  #debugonce(mse:::goFish)
  #debugonce(mp)
  
  print(Sys.time())
  
  ### run MSE
  res1 <- mp(om = input$om,
             oem = input$oem,
             iem = input$iem,
             ctrl.mp = input$ctrl.mp,
             genArgs = input$genArgs,
             tracking = input$tracking)
  
  print(Sys.time())
  ### save results
  
  
  path_out <- paste0("output/runs/whg4/", iters, "_", years)
  dir.create(path = path_out, recursive = TRUE)
  file_out <- paste0("OM",om_opt,
                     "_HCR-", input$ctrl.mp$ctrl.hcr@args$option,
                     "_Ftrgt-", input$ctrl.mp$ctrl.phcr@args$Ftrgt,
                     "_Btrigger-", input$ctrl.mp$ctrl.phcr@args$Btrigger,
                     "_TACconstr-", input$ctrl.mp$ctrl.is@args$TAC_constraint,
                     "_BB-", input$ctrl.mp$ctrl.is@args$BB
  )
  
  saveRDS(object = res1, paste0(path_out, "/", file_out, ".rds"))
### }

### ------------------------------------------------------------------------ ###
### combine and plot ####
### ------------------------------------------------------------------------ ###
if(FALSE){output <- readRDS(paste0(path_out, "/", file_out, ".rds"))

original <- readRDS(paste0("input/whg4/",iters,"_",years,"/stk.rds"))

new1<-original
new1[, ac(2018:2038)]<-res1@stock[,ac(2018:2038)]


#tiff(paste0("output/runs/plots_",om_opt,"_",iters,"/mse_stock.tiff"), bg="white",  res=200, width = 900, height = 1200)
#par(mar=c(4,4,4,2))
#plot(new1)
#dev.off()

#
### save
saveRDS(object = new1, file = paste0(path_out, "/", file_out,"_base_full_stk.rds"))
#### plot
plot(new1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + 
  xlab("year") + geom_vline(xintercept = 2018) +
  geom_hline(data = data.frame(qname = "SSB", data = 119970),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB", data = 166708),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = 0.458),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = 0.33),
             aes(yintercept = data), linetype = "solid") +
  theme_bw()
ggsave(filename = paste0(path_out, "/", file_out, ".png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


}
### ------------------------------------------------------------------------ ###
### terminate ####
### ------------------------------------------------------------------------ ###

### close R
# mpi.finalize()
### mpi.finalize() or mpi.quit() hang...
### -> kill R, the MPI processes stop afterwards
quit(save = "no")

