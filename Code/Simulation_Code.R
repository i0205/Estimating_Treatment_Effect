#******************************************************************************
#.  LOAD Functions: <UPDATE PATH> ####
#******************************************************************************
source("<path>/ATE_Functions.R")

#******************************************************************************
#                         SIMULATION RUN                                   ####
#******************************************************************************
  
  obs_month=c(0:24) # Total time points of observation : UPDATE ACCORDINGLY
  M=1000 # No. of simulation run : UPDATE ACCORDINGLY
  set.seed(0707) # Random Seed : UPDATE ACCORDINGLY

  #Storage Matrix
  #*********************************
  #Monthly
    est_mon_cost_iptw_trt_0 <- est_mon_cost_no_iptw_trt_0 <- true_mon_cost_trt_0 <-
    est_mon_cost_iptw_trt_1 <- est_mon_cost_no_iptw_trt_1 <- true_mon_cost_trt_1 <-
    est_mon_cost_iptw_trt_2 <- est_mon_cost_no_iptw_trt_2 <- true_mon_cost_trt_2 <-	
    std_mon_cost_iptw_trt_0 <- std_mon_cost_iptw_trt_1 <- std_mon_cost_iptw_trt_2 <- 
  #Accumulative 
    est_acc_cost_iptw_trt_0 <- est_acc_cost_no_iptw_trt_0 <- true_acc_cost_trt_0 <-
    est_acc_cost_iptw_trt_1 <- est_acc_cost_no_iptw_trt_1 <- true_acc_cost_trt_1 <-
    est_acc_cost_iptw_trt_2 <- est_acc_cost_no_iptw_trt_2 <- true_acc_cost_trt_2 <-	
    std_acc_cost_iptw_trt_0 <- std_acc_cost_iptw_trt_1 <- std_acc_cost_iptw_trt_2 <- 
      
    matrix(NA,nrow=length(obs_month),ncol=M)
  
  #The Simulation loop begins
  #*********************************
  for (m in 1:M){ 
    
    #******************************************************************************
    #.  Running Simulation Function for Data ####
    #******************************************************************************
    data_per_month=sim_fun(n=2000, #sample size in a single simulation
                 months=24, # max study duration in months
                 delta_e=c(0,1,-1,0), #for survival time
                 delta_1=c(2,1,2,1), #for Trt 1 selection
                 delta_2=c(1,2,1,2), #for Trt 2 selection
                 delta_b0= c(100,100,100,100), #for patient level cost variation
                 sigma=25, #for error
                 changpoint = c(3,8), #population level changepoints
                 Delta_t_1=function(t){1.5}, #difference between Trt 1 and Control
                 Delta_t_2=function(t){0.5} #difference between Trt 1 and Control
                 )
# save(data_sim_results, file = "<path>/<file.name>.RData")
    
# Observed Time points
  #**********************
  month=unique(data_per_month$MONTH)
  
  #******************************************************************************
  #.  Monthly Estimates ####
  #******************************************************************************
  
  # RUN 1: Running Monthly Cost Function for estimates with IPTW
  #***************************
  mon_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #data_per_month;
                         iptw=1, #1 or 0;
                         acc_mon='MON' #ACC or MON;
                         )
  est_mon_cost_iptw_trt_0[,m] = mon_iptw$est_trt_0 
  est_mon_cost_iptw_trt_1[,m] = mon_iptw$est_trt_1 
  est_mon_cost_iptw_trt_2[,m] = mon_iptw$est_trt_2 	
  std_mon_cost_iptw_trt_0[,m] = mon_iptw$std_trt_0
  std_mon_cost_iptw_trt_1[,m] = mon_iptw$std_trt_1
  std_mon_cost_iptw_trt_2[,m] = mon_iptw$std_trt_2
    
  # RUN 2: Running Monthly Cost Function for estimates without IPTW
  #***************************
  mon_no_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #data_per_month;
                            iptw=0, #1 or 0;
                            acc_mon='MON' #ACC or MON;
                            )
  est_mon_cost_no_iptw_trt_0[,m] = mon_no_iptw$est_trt_0 
  est_mon_cost_no_iptw_trt_1[,m] = mon_no_iptw$est_trt_1 
  est_mon_cost_no_iptw_trt_2[,m] = mon_no_iptw$est_trt_2
  
  #Observed Monthly Cost
  #**********************
  true_mon_cost_trt_0[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_0, 
                               by=list(data_per_month$MONTH), FUN=mean)$x
  true_mon_cost_trt_1[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_1, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  true_mon_cost_trt_2[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_2, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  
  
  #******************************************************************************
  #.  Accumulative Estimates ####
  #******************************************************************************
  
  # RUN 1:  Running Accumulative Cost Function for estimates with IPTW
  #***************************
  acc_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #data_per_month;
                         iptw=1, #1 or 0;
                         acc_mon='ACC' #ACC or MON;
                         )
  est_acc_cost_iptw_trt_0[,m] = acc_iptw$est_trt_0 
  est_acc_cost_iptw_trt_1[,m] = acc_iptw$est_trt_1 
  est_acc_cost_iptw_trt_2[,m] = acc_iptw$est_trt_2 	
  std_acc_cost_iptw_trt_0[,m] = acc_iptw$std_trt_0
  std_acc_cost_iptw_trt_1[,m] = acc_iptw$std_trt_1
  std_acc_cost_iptw_trt_2[,m] = acc_iptw$std_trt_2
  
  # RUN 2: Running Accumulative Cost Function for estimates without IPTW
  #***************************
  acc_no_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #data_per_month;
                            iptw=0, #1 or 0;
                            acc_mon='ACC' #ACC or MON;
                            )
  est_acc_cost_no_iptw_trt_0[,m] = acc_no_iptw$est_trt_0 
  est_acc_cost_no_iptw_trt_1[,m] = acc_no_iptw$est_trt_1 
  est_acc_cost_no_iptw_trt_2[,m] = acc_no_iptw$est_trt_2
  
  #Observed Accumulative Cost
  #**********************
  true_acc_cost_trt_0[,m] = aggregate(data_per_month$TRUE_ACC_TOT_CHRG_0, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  true_acc_cost_trt_1[,m] = aggregate(data_per_month$TRUE_ACC_TOT_CHRG_1, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  true_acc_cost_trt_2[,m] = aggregate(data_per_month$TRUE_ACC_TOT_CHRG_2, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  #The Simulation loop ends
  #*********************************
  }
  

#SAVING The data for future Use
#**********************
data_sim_results=list( est_mon_cost_iptw_trt_0 ,est_mon_cost_no_iptw_trt_0 ,true_mon_cost_trt_0 ,
      est_mon_cost_iptw_trt_1 ,est_mon_cost_no_iptw_trt_1 ,true_mon_cost_trt_1 ,
      est_mon_cost_iptw_trt_2 ,est_mon_cost_no_iptw_trt_2 ,true_mon_cost_trt_2 ,	
      std_mon_cost_iptw_trt_0 ,std_mon_cost_iptw_trt_1 ,std_mon_cost_iptw_trt_2 ,
      #Accumulative 
      est_acc_cost_iptw_trt_0 ,est_acc_cost_no_iptw_trt_0 ,true_acc_cost_trt_0 ,
      est_acc_cost_iptw_trt_1 ,est_acc_cost_no_iptw_trt_1 ,true_acc_cost_trt_1 ,
      est_acc_cost_iptw_trt_2 ,est_acc_cost_no_iptw_trt_2 ,true_acc_cost_trt_2 ,	
      std_acc_cost_iptw_trt_0 ,std_acc_cost_iptw_trt_1 ,std_acc_cost_iptw_trt_2 
   )
# save(data_sim_results, file = "<path>/<file.name1>.RData")
# load ("<path>/<file.name1>.RData")
   
#******************************************************************************
#                         COST PLOTS                               ####
#******************************************************************************

setwd("/Users/indranil/Desktop/UoL/University Stuff/PhD work/Proj 1/Data/R/OUTPUT")
pdf("Simu_Study_Plot.pdf",height = 8,width = 10)

layout(matrix(c(1:6), ncol=2, byrow = F), heights=c(1,1,1))
# A:
#*****
plot_func(obs_data=apply(true_mon_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_mon_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_mon_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("A: Monthly cost for Control"), 
                ylim=c(0,2000), y_lab="Cost ($)",
                std_iptw=apply(std_mon_cost_iptw_trt_0,1,mean)[-1],
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)
# Add a legend
legend("topright", 
       legend=c("Simulated",
                "Estimate with IPTW",
                "Estimate without IPTW"),
       col=c("black", "red", "blue"),
       lty=c(1,2,3),lwd=c(4,5,5),
       pt.cex = 2,cex=1.25)

# B1:
#*****
plot_func(obs_data=apply(true_mon_cost_trt_1,1,mean)[-1]-
                  apply(true_mon_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_mon_cost_no_iptw_trt_1,1,mean)[-1] -
                  apply(est_mon_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_mon_cost_iptw_trt_1,1,mean)[-1] -
                  apply(est_mon_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("B1: Monthly cost diff: Trt 1 vs Control"), 
                ylim=c(0,2000), y_lab="Cost Diff ($)",
                std_iptw=sqrt((apply(std_mon_cost_iptw_trt_1,1,mean)[-1])^2+(apply(std_mon_cost_iptw_trt_0,1,mean)[-1])^2),
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)
# B2:
#*****
plot_func(obs_data=apply(true_mon_cost_trt_2,1,mean)[-1]-
                  apply(true_mon_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_mon_cost_no_iptw_trt_2,1,mean)[-1] -
                  apply(est_mon_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_mon_cost_iptw_trt_2,1,mean)[-1] -
                  apply(est_mon_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("B1: Monthly cost diff: Trt 2 vs Control"), 
                ylim=c(-1000,0), y_lab="Cost Diff ($)",
                std_iptw=sqrt((apply(std_mon_cost_iptw_trt_2,1,mean)[-1])^2+(apply(std_mon_cost_iptw_trt_0,1,mean)[-1])^2),
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)

# C:
#*****
plot_func(obs_data=apply(true_acc_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_acc_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_acc_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("C: Accumulative cost for Control"), 
                ylim=c(0,12000), y_lab="Cost ($)",
                std_iptw=apply(std_acc_cost_iptw_trt_0,1,mean)[-1],
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)

# D1:
#*****
plot_func(obs_data=apply(true_acc_cost_trt_1,1,mean)[-1]-
                  apply(true_acc_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_acc_cost_no_iptw_trt_1,1,mean)[-1] -
                  apply(est_acc_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_acc_cost_iptw_trt_1,1,mean)[-1] -
                  apply(est_acc_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("D1: Accu cost diff: Trt 1 vs Control"), 
                ylim=c(0,15000), y_lab="Cost Diff ($)",
                std_iptw=sqrt((apply(std_acc_cost_iptw_trt_1,1,mean)[-1])^2+(apply(std_acc_cost_iptw_trt_0,1,mean)[-1])^2),
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)
# D2:
#*****
plot_func(obs_data=apply(true_acc_cost_trt_2,1,mean)[-1]-
                  apply(true_acc_cost_trt_0,1,mean)[-1],
                est_no_iptw=apply(est_acc_cost_no_iptw_trt_2,1,mean)[-1] -
                  apply(est_acc_cost_no_iptw_trt_0,1,mean)[-1],
                est_iptw=apply(est_acc_cost_iptw_trt_2,1,mean)[-1] -
                  apply(est_acc_cost_iptw_trt_0,1,mean)[-1],
                tot_mon=obs_month[-1], title=c("D1: Accu cost diff: Trt 2 vs Control"), 
                ylim=c(-7000,0), y_lab="Cost Diff ($)",
                std_iptw=sqrt((apply(std_acc_cost_iptw_trt_2,1,mean)[-1])^2+(apply(std_acc_cost_iptw_trt_0,1,mean)[-1])^2),
                CI_IPTW=0, # 1 for Yes;
                h_abl=0 # 1 for abline-h;
)

dev.off()


# #######################################################################
#        Performance Matrix Calculation      ####
# #######################################################################

Error_mat=matrix((
  round(rbind(
    #Monthly
    #IPTW
    cbind(mean(apply(abs((est_mon_cost_iptw_trt_1-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((est_mon_cost_iptw_trt_1-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((est_mon_cost_iptw_trt_1-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((est_mon_cost_iptw_trt_2-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((est_mon_cost_iptw_trt_2-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((est_mon_cost_iptw_trt_2-est_mon_cost_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #No IPTW
    cbind(mean(apply(abs((est_mon_cost_no_iptw_trt_1-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((est_mon_cost_no_iptw_trt_1-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((est_mon_cost_no_iptw_trt_1-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_1-true_mon_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((est_mon_cost_no_iptw_trt_2-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((est_mon_cost_no_iptw_trt_2-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((est_mon_cost_no_iptw_trt_2-est_mon_cost_no_iptw_trt_0)-apply(true_mon_cost_trt_2-true_mon_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #Accumulative
    #IPTW
    cbind(mean(apply(abs((est_acc_cost_iptw_trt_1-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((est_acc_cost_iptw_trt_1-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((est_acc_cost_iptw_trt_1-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((est_acc_cost_iptw_trt_2-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((est_acc_cost_iptw_trt_2-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((est_acc_cost_iptw_trt_2-est_acc_cost_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #No IPTW
    cbind(mean(apply(abs((est_acc_cost_no_iptw_trt_1-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((est_acc_cost_no_iptw_trt_1-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((est_acc_cost_no_iptw_trt_1-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_1-true_acc_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((est_acc_cost_no_iptw_trt_2-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((est_acc_cost_no_iptw_trt_2-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((est_acc_cost_no_iptw_trt_2-est_acc_cost_no_iptw_trt_0)-apply(true_acc_cost_trt_2-true_acc_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    )
  ),2)),
  ncol=6,
  dimnames=list(c(paste0(c("No_chg_"),"sig10_",c("Mon_", "Acc_"), rep(c("IPTW","No_IPTW"),2))),
  c(paste0(rep(c("Del1_", "Del2_"),each=3),c("MAE","MBE","RMSE")))))
