### ------------------------------------------------------------------------ ###
### prepare objects for new a4a standard mse package ####
### ------------------------------------------------------------------------ ###
### https://github.com/flr/mse
rm(list = ls())

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

setwd(paste("/home/miethet/MSE_whiting", sep=""))

source("a4a_mse_WKNSMSE_funs.R")

om_opt<-"OM1"
n<-1000

### save workspace to start from here
# save.image(file = "input/whg4/image_10.RData")
load(file = paste0("input/whg4/base_image/image_",om_opt,"_",n,".RData"))



### reference points
refpts_mse <- list(Btrigger = 166708,
                   Ftrgt = 0.172,
				   Fpa = 0.33,
                   Bpa = 166708,
                   Blim = 119970
        				  
				           )  # vary for different reference points
				   
### some specifications for short term forecast with SAM
whg4_stf_def <- list(fwd_yrs_average = -3:-1,
                     fwd_yrs_rec_start = 2002,
                     fwd_yrs_sel = -3:-1,
                     fwd_yrs_lf_remove = -2:-1, # overwrite recent years landing fraction, does not apply fwd_splitLD=FALSE
                     fwd_splitLD = FALSE)

### some arguments (passed to mp())
genArgs <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
                y0 = dims(stk_fwd)$minyear, ### first data year
                iy = yr_data, ### first simulation (intermediate) year
                nsqy = 3, ### not used, but has to provided
                nblocks = 1, ### block for parallel processing
                seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           projection = mseCtrl(method = fwd_WKNSMSE, 
                                args = list(maxF = 2,
                                            ### process noise on stock.n
                                            proc_res = "fitted"
                                ))
)

### observation (error) model
oem <- FLoem(method = oem_WKNSMSE,
             observations = list(stk = stk_oem, idx = idx), 
             deviances = list(stk = FLQuants(catch.dev = catch_res),   #CATCH RES
                              idx = idx_dev),
             args = list(idx_timing = c(0, -1),
                         catch_timing = -1,
                         use_catch_residuals = TRUE, 
                         use_idx_residuals = TRUE,
                         use_stk_oem = TRUE))
                         
### implementation error model (banking and borrowing)
# iem <- FLiem(method = iem_WKNSMSE, 
#              args = list(BB = FALSE, imp.error=FALSE, imp.res=catch_res_l))  # extra implementation error

### default management
ctrl_obj <- mpCtrl(list(
  ctrl.est = mseCtrl(method = SAM_wrapper,
                     args = c(### short term forecast specifications
                       forecast = TRUE, 
                       fwd_trgt = "fsq", fwd_yrs = 1, 
                       whg4_stf_def,
                       ### speeding SAM up
                       newtonsteps = 0, rel.tol = 0.001,
                       par_ini = list(sam_initial),
                       track_ini = TRUE, ### store ini for next year
                       ### SAM model specifications
                       conf = list(whg4_conf_sam),
                       parallel = FALSE ### TESTING ONLY
                     )),
  ctrl.phcr = mseCtrl(method = phcr_WKNSMSE,
                      args = refpts_mse),
  ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  ctrl.is = mseCtrl(method = is_WKNSMSE, 
                    args = c(hcrpars = list(refpts_mse),
                             ### for short term forecast
                             fwd_trgt = c("fsq", "hcr"), fwd_yrs = 2,
                             whg4_stf_def                            
                             ### TAC constraint
                             #TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,
                             #Btrigger_cond = FALSE,
                             ### banking and borrowing 
                             #BB = TRUE,
                             #BB_check_hcr = FALSE,
                             #BB_check_fc = TRUE,
                             #BB_rho = list(c(-0.1, 0.1))
                    ))#,
  #ctrl.tm = NULL
))
### additional tracking metrics
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")


input <- list(om = om, oem = oem, ctrl.mp = ctrl_obj,
              genArgs = genArgs, tracking = tracking_add)
saveRDS(object = input,file = paste0("input/whg4/base_run/",om_opt,"_base_run",n,".rds"))

