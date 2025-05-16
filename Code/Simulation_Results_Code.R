#******************************************************************************
#.  LOAD Functions: <UPDATE PATH> ####
#******************************************************************************
source("<Path>/ATE_Functions.R")

#******************************************************************************
#                         SIMULATION RUN                                   ####
#******************************************************************************
  
  obs_month=c(1:24) # Total time points of observation : UPDATE ACCORDINGLY
  M=1000 # No. of simulation run : UPDATE ACCORDINGLY (Atleast M=2 for plots)
  set.seed(0707) # Random Seed : UPDATE ACCORDINGLY

  #Storage Matrix
  #*********************************
  #Monthly
    #Cost
    est_mon_cost_iptw_trt_0 <- est_mon_cost_no_iptw_trt_0 <- true_mon_cost_trt_0 <-
    est_mon_cost_iptw_trt_1 <- est_mon_cost_no_iptw_trt_1 <- true_mon_cost_trt_1 <-
    est_mon_cost_iptw_trt_2 <- est_mon_cost_no_iptw_trt_2 <- true_mon_cost_trt_2 <-	
    #Diff
    est_mon_diff_iptw_trt_1_vs_0 <- est_mon_diff_iptw_trt_2_vs_1 <- est_mon_diff_iptw_trt_2_vs_0 <- 
    est_mon_diff_no_iptw_trt_1_vs_0 <- est_mon_diff_no_iptw_trt_2_vs_1 <- est_mon_diff_no_iptw_trt_2_vs_0 <- 
  #Accumulative 
    #Cost
    est_acc_cost_iptw_trt_0 <- est_acc_cost_no_iptw_trt_0 <- true_acc_cost_trt_0 <-
    est_acc_cost_iptw_trt_1 <- est_acc_cost_no_iptw_trt_1 <- true_acc_cost_trt_1 <-
    est_acc_cost_iptw_trt_2 <- est_acc_cost_no_iptw_trt_2 <- true_acc_cost_trt_2 <-	
    #Diff
    est_acc_diff_iptw_trt_1_vs_0 <- est_acc_diff_iptw_trt_2_vs_1 <- est_acc_diff_iptw_trt_2_vs_0 <- 
    est_acc_diff_no_iptw_trt_1_vs_0 <- est_acc_diff_no_iptw_trt_2_vs_1 <- est_acc_diff_no_iptw_trt_2_vs_0 <- 

    matrix(NA,nrow=length(obs_month),ncol=M)
    
  #The Simulation loop begins
  #*********************************
  for (m in 1:M){ 
    
  #***********************************************************************
  #.  Running Simulation Function for Data ####
  #***********************************************************************
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
  # save(data_per_month, file = "<Path>/data_per_month.RData")

  # Observed Time points
  #**********************
  month=unique(data_per_month$MONTH)
  
#***********************************************************************
#.  Monthly Estimates ####
#***********************************************************************
  
  # RUN 1: Running Monthly Cost Function for estimates with IPTW
  #***************************
  mon_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #input dataset;
                         iptw=1, #1 or 0;
                         acc_mon='MON', #ACC or MON;
                         knot_gap=2, #smoothing value for knot seq
                         lambda='best' # provide value for use instead
                         )
  est_mon_cost_iptw_trt_0[,m] = mon_iptw$est_trt_0 
  est_mon_cost_iptw_trt_1[,m] = mon_iptw$est_trt_1 
  est_mon_cost_iptw_trt_2[,m] = mon_iptw$est_trt_2 	

  est_mon_diff_iptw_trt_1_vs_0[,m] = mon_iptw$est_diff_1_vs_0 
  est_mon_diff_iptw_trt_2_vs_1[,m] = mon_iptw$est_diff_2_vs_1 
  est_mon_diff_iptw_trt_2_vs_0[,m] = mon_iptw$est_diff_2_vs_0 	

  # RUN 2: Running Monthly Cost Function for estimates without IPTW
  #***************************
  mon_no_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #input dataset;
                            iptw=0, #1 or 0;
                            acc_mon='MON', #ACC or MON;
                            knot_gap=2 , #smoothing value for knot seq
                            lambda='best' # provide value for use instead
                            )
  est_mon_cost_no_iptw_trt_0[,m] = mon_no_iptw$est_trt_0 
  est_mon_cost_no_iptw_trt_1[,m] = mon_no_iptw$est_trt_1 
  est_mon_cost_no_iptw_trt_2[,m] = mon_no_iptw$est_trt_2

  est_mon_diff_no_iptw_trt_1_vs_0[,m] = mon_no_iptw$est_diff_1_vs_0 
  est_mon_diff_no_iptw_trt_2_vs_1[,m] = mon_no_iptw$est_diff_2_vs_1 
  est_mon_diff_no_iptw_trt_2_vs_0[,m] = mon_no_iptw$est_diff_2_vs_0 	

  #Observed Monthly Cost
  #**********************
  #*
  true_mon_cost_trt_0[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_0, 
                               by=list(data_per_month$MONTH), FUN=mean)$x
  true_mon_cost_trt_1[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_1, 
                                by=list(data_per_month$MONTH), FUN=mean)$x
  true_mon_cost_trt_2[,m] = aggregate(data_per_month$TRUE_MON_TOT_CHRG_2, 
                                by=list(data_per_month$MONTH), FUN=mean)$x

 #***********************************************************************
 #.  Accumulative Estimates ####
 #***********************************************************************
  
  # RUN 1:  Running Accumulative Cost Function for estimates with IPTW
  #***************************
  acc_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #input dataset;
                         iptw=1, #1 or 0;
                         acc_mon='ACC', #ACC or MON;
                         knot_gap=2, #smoothing value for knot seq
                         lambda='best' # provide value for use instead
                         )
  est_acc_cost_iptw_trt_0[,m] = acc_iptw$est_trt_0 
  est_acc_cost_iptw_trt_1[,m] = acc_iptw$est_trt_1 
  est_acc_cost_iptw_trt_2[,m] = acc_iptw$est_trt_2 	

  est_acc_diff_iptw_trt_1_vs_0[,m] = acc_iptw$est_diff_1_vs_0 
  est_acc_diff_iptw_trt_2_vs_1[,m] = acc_iptw$est_diff_2_vs_1 
  est_acc_diff_iptw_trt_2_vs_0[,m] = acc_iptw$est_diff_2_vs_0 	

  # RUN 2: Running Accumulative Cost Function for estimates without IPTW
  #***************************
  acc_no_iptw=ATE_IPTW_FUNC(in_data=data_per_month, #input dataset;
                            iptw=0, #1 or 0;
                            acc_mon='ACC', #ACC or MON;
                            knot_gap=2, #smoothing value for knot seq
                            lambda='best' # provide value for use instead
                            )
  est_acc_cost_no_iptw_trt_0[,m] = acc_no_iptw$est_trt_0 
  est_acc_cost_no_iptw_trt_1[,m] = acc_no_iptw$est_trt_1 
  est_acc_cost_no_iptw_trt_2[,m] = acc_no_iptw$est_trt_2
  
  est_acc_diff_no_iptw_trt_1_vs_0[,m] = acc_no_iptw$est_diff_1_vs_0 
  est_acc_diff_no_iptw_trt_2_vs_1[,m] = acc_no_iptw$est_diff_2_vs_1 
  est_acc_diff_no_iptw_trt_2_vs_0[,m] = acc_no_iptw$est_diff_2_vs_0 	
  
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
  
  #***********************************************************************
    #SAVING simulation data for future Use
 #***********************************************************************
  data_sim_results=list(   
  #Monthly
  #Cost
  est_mon_cost_iptw_trt_0= est_mon_cost_iptw_trt_0, 
  est_mon_cost_no_iptw_trt_0= est_mon_cost_no_iptw_trt_0, 
  true_mon_cost_trt_0= true_mon_cost_trt_0,
  est_mon_cost_iptw_trt_1= est_mon_cost_iptw_trt_1, 
  est_mon_cost_no_iptw_trt_1= est_mon_cost_no_iptw_trt_1, 
  true_mon_cost_trt_1= true_mon_cost_trt_1,
  est_mon_cost_iptw_trt_2=  est_mon_cost_iptw_trt_2, 
  est_mon_cost_no_iptw_trt_2= est_mon_cost_no_iptw_trt_2, 
  true_mon_cost_trt_2= true_mon_cost_trt_2,	
  #Diff
  est_mon_diff_iptw_trt_1_vs_0= est_mon_diff_iptw_trt_1_vs_0, 
  est_mon_diff_iptw_trt_2_vs_1= est_mon_diff_iptw_trt_2_vs_1, 
  est_mon_diff_iptw_trt_2_vs_0= est_mon_diff_iptw_trt_2_vs_0, 
  est_mon_diff_no_iptw_trt_1_vs_0= est_mon_diff_no_iptw_trt_1_vs_0, 
  est_mon_diff_no_iptw_trt_2_vs_1= est_mon_diff_no_iptw_trt_2_vs_1, 
  est_mon_diff_no_iptw_trt_2_vs_0= est_mon_diff_no_iptw_trt_2_vs_0, 
  #Accumulative 
  #Cost
  est_acc_cost_iptw_trt_0= est_acc_cost_iptw_trt_0, 
  est_acc_cost_no_iptw_trt_0= est_acc_cost_no_iptw_trt_0, 
  true_acc_cost_trt_0= true_acc_cost_trt_0,
  est_acc_cost_iptw_trt_1= est_acc_cost_iptw_trt_1, 
  est_acc_cost_no_iptw_trt_1= est_acc_cost_no_iptw_trt_1, 
  true_acc_cost_trt_1= true_acc_cost_trt_1,
  est_acc_cost_iptw_trt_2=  est_acc_cost_iptw_trt_2, 
  est_acc_cost_no_iptw_trt_2= est_acc_cost_no_iptw_trt_2, 
  true_acc_cost_trt_2= true_acc_cost_trt_2,	
  #Diff
  est_acc_diff_iptw_trt_1_vs_0= est_acc_diff_iptw_trt_1_vs_0, 
  est_acc_diff_iptw_trt_2_vs_1= est_acc_diff_iptw_trt_2_vs_1, 
  est_acc_diff_iptw_trt_2_vs_0= est_acc_diff_iptw_trt_2_vs_0, 
  est_acc_diff_no_iptw_trt_1_vs_0= est_acc_diff_no_iptw_trt_1_vs_0, 
  est_acc_diff_no_iptw_trt_2_vs_1= est_acc_diff_no_iptw_trt_2_vs_1, 
  est_acc_diff_no_iptw_trt_2_vs_0= est_acc_diff_no_iptw_trt_2_vs_0 
  )
# save(data_sim_results, file = "<Path>/data_sim_results.RData")

#***********************************************************************
#                         COST PLOTS                               ####
#***********************************************************************
# load ("<Path>/data_per_month.RData")
# load ("<Path>/data_sim_results.RData")

setwd("<Path>")
pdf("sim_plot.pdf",height = 8,width = 12)
layout(matrix(c(1:6), ncol=3, byrow = T), heights=c(1,1,1))

# A1:
#*****
plot_func(obs_data=apply(data_sim_results$true_acc_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_acc_cost_no_iptw_trt_0,1,mean),
          est_iptw=apply(data_sim_results$est_acc_cost_iptw_trt_0,1,mean),
          tot_mon=obs_month, title=c("A1: Accu cost for Control"), 
          ylim=c(0,12000), y_lab="Cost ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_acc_cost_iptw_trt_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)

# A2:
#*****
plot_func(obs_data=apply(data_sim_results$true_acc_cost_trt_1,1,mean)-
            apply(data_sim_results$true_acc_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_acc_diff_no_iptw_trt_1_vs_0,1,mean),
          est_iptw=apply(data_sim_results$est_acc_diff_iptw_trt_1_vs_0,1,mean),
          tot_mon=obs_month, title=c("A2: Accu diff: Trt 1 vs Control"), 
          ylim=c(0,15000), y_lab="Cost Diff ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_acc_diff_iptw_trt_1_vs_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)
# A3:
#*****
plot_func(obs_data=apply(data_sim_results$true_acc_cost_trt_2,1,mean)-
            apply(data_sim_results$true_acc_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_acc_diff_no_iptw_trt_2_vs_0,1,mean),
          est_iptw=apply(data_sim_results$est_acc_diff_iptw_trt_2_vs_0,1,mean),
          tot_mon=obs_month, title=c("A3: Accu diff: Trt 2 vs Control"), 
          ylim=c(-6000,0), y_lab="Cost Diff ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_acc_diff_iptw_trt_2_vs_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)

# B1:
#*****
plot_func(obs_data=apply(data_sim_results$true_mon_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_mon_cost_no_iptw_trt_0,1,mean),
          est_iptw=apply(data_sim_results$est_mon_cost_iptw_trt_0,1,mean),
          tot_mon=obs_month, title=c("B1: Monthly cost for Control"), 
          ylim=c(0,2000), y_lab="Cost ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_mon_cost_iptw_trt_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)

# Add a legend
legend("topright", 
       legend=c("True values",
                "Estimate with IPTW",
                "Estimate without IPTW"),
       col=c("black", "red", "blue"),
       lty=c(1,2,3),lwd=c(4,5,5),
       pt.cex = 2,cex=1.5)
# B2:
#*****
plot_func(obs_data=apply(data_sim_results$true_mon_cost_trt_1,1,mean)-
            apply(data_sim_results$true_mon_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_mon_diff_no_iptw_trt_1_vs_0,1,mean),
          est_iptw=apply(data_sim_results$est_mon_diff_iptw_trt_1_vs_0,1,mean),
          tot_mon=obs_month, title=c("B2: Monthly diff: Trt 1 vs Control"), 
          ylim=c(0,2000), y_lab="Cost Diff ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_mon_diff_iptw_trt_1_vs_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)
# B3:
#*****
plot_func(obs_data=apply(data_sim_results$true_mon_cost_trt_2,1,mean)-
            apply(data_sim_results$true_mon_cost_trt_0,1,mean),
          est_no_iptw=apply(data_sim_results$est_mon_diff_no_iptw_trt_2_vs_0,1,mean),
          est_iptw=apply(data_sim_results$est_mon_diff_iptw_trt_2_vs_0,1,mean),
          tot_mon=obs_month, title=c("B3: Monthly diff: Trt 2 vs Control"), 
          ylim=c(-1000,0), y_lab="Cost Diff ($)",
          CI_IPTW=0, # 1 for Yes;
          std_iptw=apply(data_sim_results$est_mon_diff_iptw_trt_2_vs_0,1,sd),
          no_IPTW=1, # 1 for Yes;
          h_abl=0 # 1 for abline-h;
)
dev.off()


#***********************************************************************
#        Performance Matrix Calculation                             ####
#***********************************************************************

Error_mat=matrix((
  round(rbind(
    #Monthly
    #IPTW
    cbind(mean(apply(abs((data_sim_results$est_mon_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((data_sim_results$est_mon_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((data_sim_results$est_mon_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((data_sim_results$est_mon_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((data_sim_results$est_mon_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((data_sim_results$est_mon_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #No IPTW
    cbind(mean(apply(abs((data_sim_results$est_mon_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((data_sim_results$est_mon_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((data_sim_results$est_mon_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_mon_cost_trt_1-data_sim_results$true_mon_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((data_sim_results$est_mon_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((data_sim_results$est_mon_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((data_sim_results$est_mon_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_mon_cost_trt_2-data_sim_results$true_mon_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #Accumulative
    #IPTW
    cbind(mean(apply(abs((data_sim_results$est_acc_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((data_sim_results$est_acc_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((data_sim_results$est_acc_diff_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((data_sim_results$est_acc_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((data_sim_results$est_acc_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((data_sim_results$est_acc_diff_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    ),
    #No IPTW
    cbind(mean(apply(abs((data_sim_results$est_acc_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 1
          abs(mean(apply(((data_sim_results$est_acc_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 1
          mean(sqrt(apply(((data_sim_results$est_acc_diff_no_iptw_trt_1_vs_0)-apply(data_sim_results$true_acc_cost_trt_1-data_sim_results$true_acc_cost_trt_0,1,mean))^2,2,mean))), #RMSE #TRT 1
          
          mean(apply(abs((data_sim_results$est_acc_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean)),2,max)), #MAE #TRT 2
          abs(mean(apply(((data_sim_results$est_acc_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean)),2,mean))), #MBE  #TRT 2
          mean(sqrt(apply(((data_sim_results$est_acc_diff_no_iptw_trt_2_vs_0)-apply(data_sim_results$true_acc_cost_trt_2-data_sim_results$true_acc_cost_trt_0,1,mean))^2,2,mean))) #RMSE #TRT 2
    )
  ),2)),
  ncol=6,
  dimnames=list(c(paste0(rep(c("Mon_", "Acc_"),each=2), rep(c("IPTW", "No_IPTW"),2))),
  c(paste0(rep(c("Del1_", "Del2_"),each=3),c("MAE","MBE","RMSE"))))
)
#Display Error
Error_mat
