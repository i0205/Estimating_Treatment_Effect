library(splines2)
library(glmnet)

#******************************************************************************
#******************************************************************************
#.                        FUNCTION 1: SIMULATION DATASET FUNCTION ####
#******************************************************************************
#******************************************************************************
# OUTPUTS: Dataset for Estimation Function (Function 2)
#******************************************************************************

#Defining the Simulation Function
sim_fun=function(n=2000, #sample size in a single simulation
                 months=24, #study duration in months
                 delta_e=c(0,1,-1,0), # for survival time
                 delta_1=c(2,1,2,1), #for treatment selection
                 delta_2=c(1,2,1,2), #for treatment selection
                 delta_b0= c(100,100,100,100), #for patient level cost variation
                 sigma=10, #for error
                 changpoint = c(3,8), #population level changepoints
                 Delta_t_1=function(t){1.5}, #difference between Trt 1 and Control
                 Delta_t_2=function(t){0.5} #difference between Trt 1 and Control
){
  #******************************************************************************
  #   STEP 1: SIMULATING THE DATA ####
  #******************************************************************************
  #COVARIATES
  #*********************************
  sim_X=data.frame(cbind(
    PATIENT_ID=c(1:n),
    X1=rbinom(n,1,.5)-0.5, X2=rbinom(n,1,.5)-0.5,
    X3=rnorm(n,0,1), X4=rnorm(n,0,1)
  ))
  X_mat=as.matrix(sim_X[,c("X1", "X2", "X3", "X4")])
  
  #Treatment
  #*********************************
  #logit(p)=X*Delta
  del_1_X=exp(X_mat%*%delta_1)
  del_2_X=exp(X_mat%*%delta_2)
  
  p1=del_1_X/(1+del_1_X+del_2_X)
  p2=del_2_X/(1+del_1_X+del_2_X)
  p0=1-p1-p2
  P=cbind(p0,p1,p2)
  
  for (i in 1:n){
    sim_X$TREATMENT[i]=sample(0:2, size = 1, prob = P[i,])
  }
  
  #Months Alive
  #*********************************
  sim_X$MONTH_TOT_con_trt2=round(rexp(n,1/(50+X_mat%*%delta_e )))
  sim_X$MONTH_TOT_trt1=sim_X$MONTH_TOT_con_trt2+ round(rexp(n,1/(5)))
  #Assumption: Trt 1 gives higher life expectancy (+exp(5) months).
  sim_X$MONTH_TOT = ifelse(sim_X$TREATMENT==1, sim_X$MONTH_TOT_trt1, sim_X$MONTH_TOT_con_trt2 ) 
  sim_X$D=ifelse(sim_X$MONTH_TOT <= months + changpoint[2], 1, 0 ) 
  
  #PS WEIGHTS
  #*********************************
  sim_X$no_weight=1
  glm_fit_sim=nnet::multinom(as.factor(TREATMENT)~as.factor(X1)+as.factor(X2)+X3+X4,
                             data = sim_X, trace=F)
  ps=predict(glm_fit_sim,type = "probs", newdata = sim_X)
  sim_X$ps_weights <- rep(0, nrow(sim_X)) #initialize weights
  for (i in levels(as.factor(sim_X$TREATMENT))) {
    sim_X$ps_weights[as.factor(sim_X$TREATMENT) == i] <-
      1/ps[as.factor(sim_X$TREATMENT) == i, i]
  }
  
  #Simulated Monthly Cost
  #*********************************
  #Cost calculation variables to include patient level variation
  beta_0= abs(100 + (X_mat%*%delta_b0))
  #dataset shell
  ACC_TOT_DATA=data.frame(matrix(rep(NA, 9*n*months), nrow=n*months, 
                                 ncol=9, dimnames=list(c(NULL),
                                                       c("PATIENT_ID", "MONTH","MON_TOT_CHRG", 
                                                         "TRUE_MON_TOT_CHRG_0", "TRUE_MON_TOT_CHRG_1",
                                                         "TRUE_MON_TOT_CHRG_2", "MON_TOT_ALIVE", "D", "TREATMENT"
                                                       ))))
  k=0
  for (i in 1:n){
    for (t in 1:months){
      k=k+1
      
      ACC_TOT_DATA[k,"PATIENT_ID"]=i
      ACC_TOT_DATA[k,"MONTH"]=t
      ACC_TOT_DATA[k,"D"]=sim_X$D[i]
      ACC_TOT_DATA[k,"TREATMENT"]=sim_X$TREATMENT[i]
      ACC_TOT_DATA[k,"MON_TOT_ALIVE"]=sim_X$MONTH_TOT[i]
      s_error=rnorm(1,0,sigma)
      
      #Cost Function
      #*********************************
      # Changepoints considered are: changepoint[1] months from trt and changepoint[2] months before death
      Ft= function(x, tot_time) { 
        ind=ifelse(tot_time <= months + changpoint[2], 1, 0 ) #death indicator function
        ifelse(ind==0, 
               #alive
               ifelse(x<=changpoint[1], 20*(beta_0[i])*(0.5)^x, 20*(beta_0[i])*(0.5)^changpoint[1]), 
               #death
               ifelse(x<=changpoint[1], 20*(beta_0[i])*(0.5)^x, 20*(beta_0[i])*(0.5)^changpoint[1]) + 
                 ifelse(x<(tot_time-changpoint[2]), 0, 
                        ifelse(x<=tot_time, (beta_0[i])*(x-(tot_time-changpoint[2])),0))
        )
      }   
      #adjusting for survival
      time=sim_X$MONTH_TOT[i]
      time_cont_trt2= sim_X$MONTH_TOT_con_trt2[i]
      time_trt1= sim_X$MONTH_TOT_trt1[i]
      
      #calculating monthly cost
      ACC_TOT_DATA[k,"MON_TOT_CHRG"]=ifelse(t<=time,
                                            ifelse(sim_X$TREATMENT[i]==0, Ft(t,time) + abs(s_error), #trt=0
                                                   ifelse(sim_X$TREATMENT[i]==1, Delta_t_1(t)*( Ft(t, time) + abs(s_error)), #trt=1 
                                                          Delta_t_2(t)*(Ft(t, time) + abs(s_error)))), #trt=2 
                                            0)
      #calculating true monthly treatment under treatment
      ACC_TOT_DATA[k,"TRUE_MON_TOT_CHRG_0"]=ifelse(t<=time_cont_trt2,Ft(t,time_cont_trt2) + abs(s_error),0) #trt=0
      ACC_TOT_DATA[k,"TRUE_MON_TOT_CHRG_1"]=ifelse(t<=time_trt1,Delta_t_1(t)*Ft(t,time_trt1)+ abs(s_error),0) #trt=1: +6 months of survival
      ACC_TOT_DATA[k,"TRUE_MON_TOT_CHRG_2"]=ifelse(t<=time_cont_trt2,Delta_t_2(t)*Ft(t,time_cont_trt2)+ abs(s_error),0) #trt=2
    }
  }
  
  #Creating all months data
  #*********************************
  dummy_data=sim_X[rep(seq_len(nrow(sim_X)), each=months),]
  dummy_data$MONTH=rep(c(1:months),length(unique(dummy_data$PATIENT_ID)))
  
  sim_monthly_X_=merge(dummy_data,subset(ACC_TOT_DATA, select = -c(TREATMENT)), by=c("PATIENT_ID","MONTH"), 
                      all.x = TRUE)
  sim_monthly_X=sim_monthly_X_[order(sim_monthly_X_$PATIENT_ID, sim_monthly_X_$MONTH),]
  sim_monthly_X[is.na(sim_monthly_X)] <- 0
  
  #Simulated Accumulative Cost
  #*********************************
    for (i in 1:n) {
    for (t in 1:months) {
      l=(i-1)*(months) + t #index for calculation
      
      sim_monthly_X$ACC_TOT_CHRG[l]=ifelse(sim_monthly_X$MONTH[l]==min(sim_monthly_X$MONTH),
                                           sim_monthly_X$MON_TOT_CHRG[l],
                                           sim_monthly_X$MON_TOT_CHRG[l] + sim_monthly_X$ACC_TOT_CHRG[l-1])
      
      sim_monthly_X$TRUE_ACC_TOT_CHRG_0[l]=ifelse(sim_monthly_X$MONTH[l]==min(sim_monthly_X$MONTH),
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_0[l],
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_0[l] + sim_monthly_X$TRUE_ACC_TOT_CHRG_0[l-1])
      
      sim_monthly_X$TRUE_ACC_TOT_CHRG_1[l]=ifelse(sim_monthly_X$MONTH[l]==min(sim_monthly_X$MONTH),
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_1[l],
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_1[l] + sim_monthly_X$TRUE_ACC_TOT_CHRG_1[l-1])
      
      sim_monthly_X$TRUE_ACC_TOT_CHRG_2[l]=ifelse(sim_monthly_X$MONTH[l]==min(sim_monthly_X$MONTH),
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_2[l],
                                                  sim_monthly_X$TRUE_MON_TOT_CHRG_2[l] + sim_monthly_X$TRUE_ACC_TOT_CHRG_2[l-1])
    }
  }
  
  return(sim_mon_cost_data=sim_monthly_X)
  
} #end of simulation function (sim_fun)



#******************************************************************************
#******************************************************************************
#.                        FUNCTION 2: ATE ESTIMATION MODEL ####
#******************************************************************************
#******************************************************************************
# Note: Update return list and treatment specific variables 
#       for change in number of treatments.
#******************************************************************************
# INPUT: The in-dataset must have the following variables:
# PATIENT_ID, MONTH, TREATMENT, ps_weights, no_weight(assign 1 for all values), 
# MON_TOT_CHRG, ACC_TOT_CHRG
#******************************************************************************
# OUTPUTS: Estimates and CI for 3 treatment groups. Can be used in Plot Function (3) 
#******************************************************************************

ATE_IPTW_FUNC=function(in_data=data_per_month, #data_per_month;
                        iptw=1, #1 for Yes ;
                        acc_mon='MON', #ACC or MON;
                       knot_gap=2, #smoothing value for knot seq
                       lambda='best' # provide single value for use instead
                      ){

#******************************************************************************
#.  STEP 1: I/M-SPLINE Model ####
#******************************************************************************
#Knots
knot_start=min(in_data$MONTH)+1
knot_end=max(in_data$MONTH)-1

#Beta I-pline
spline_beta=iSpline(x=in_data$MONTH, df = NULL, 
                    knots = seq(knot_start, knot_end,knot_gap), 
                    degree = 3L,intercept = F)  
spline_beta=data.frame(1, spline_beta)
colnames(spline_beta) = paste0("Beta",(1:ncol(spline_beta)))
spline_est=spline_beta
#Adding Gamma_i I-splines for additional treatment; Control=0
trt_val=sort(unique(in_data$TREATMENT))
  for (i in trt_val[2:length(trt_val)]) {
    spline_val=ifelse(in_data$TREATMENT==i, 1, 0)*spline_beta
    colnames(spline_val) = paste0("Gamma_",i,"_",1:ncol(spline_val))
    spline_est=cbind(spline_est, spline_val)
  }

#******************************************************************************
#. STEP 2: GLMNET Model####
#******************************************************************************

#IPTW/ No Weights
if (iptw==1){WEIGHTS=in_data$ps_weights} else{WEIGHTS=in_data$no_weight}

  cost_var=in_data$ACC_TOT_CHRG
  glm_X=as.matrix(spline_est)

  # Cross-validation to select lambda/provided lambda
  if ((lambda)=='best' ){
    cv_fit <- cv.glmnet(y=cost_var, x=glm_X, alpha = 0, weights = WEIGHTS)
    best_lambda <- cv_fit$lambda.min
  } else { 
    best_lambda<- lambda
    }
  
  # Fit the weighted ridge regression model
  GLMnet_est=glmnet(y=cost_var, x=glm_X, alpha = 0, weights = WEIGHTS, 
                 lambda = best_lambda)

#******************************************************************************
#. STEP 3: Prediction ####
#******************************************************************************

#Total Trt
trt_val=sort(unique(in_data$TREATMENT))
  
#Choosing I-Spline for accumulative prediction
#********************************************

spln=iSpline(x=c(1:max(in_data$MONTH)), df = NULL, 
               knots = seq(knot_start, knot_end,knot_gap), degree = 3L,intercept = F) 
spln=cbind(1,spln)
#Prediction for Control
#**************
# Pred dataset
contrl_est=data.frame(cbind(spln, matrix(data=0,nrow=nrow(spln), 
                      ncol=(length(trt_val)-1)*ncol(spln))))
colnames(contrl_est)=colnames(spline_est)
  
  pred_data=data.frame(predict(GLMnet_est,newx=as.matrix(contrl_est), 
                                          s = best_lambda))
  colnames(pred_data)=c("est_trt_0")
#Prediction for Treatments
#**************
for (i in trt_val[2:length(trt_val)]) {
  # Pred
  trt_est=data.frame(contrl_est)
 trt_est[,(i*ncol(spln)+1):((i+1)*ncol(spln))] =data.frame(spln) #Spline
    
  estimate= data.frame(predict(GLMnet_est,newx=as.matrix(trt_est), 
                               s = best_lambda))
   colnames(estimate)=c(paste0("est_trt_",i))
 
  pred_data=cbind(pred_data,estimate)
}  
#Prediction for Differences between Trts
#****************************
for (i in trt_val[2:length(trt_val)]) {
  for (j in trt_val[1:i]) {
  # Pred dataset
  estimate_diff=data.frame(pred_data[,i+1]-pred_data[,j+1])
  colnames(estimate_diff)=c(paste0("est_diff_",i,"_vs_",j))
  pred_data=cbind(pred_data,estimate_diff)
  }
}
  
#Monthly Cost prediction
#********************************************
  if (acc_mon=='MON'){
#Acc difference is monthly cost
  pred_data_mon=pred_data
  pred_data_mon=pred_data[1,]
  for (i in 2:nrow(pred_data)) {
    pred_data_mon[i,]=pred_data[i,]-pred_data[i-1,]
    }
  pred_data=pred_data_mon
  }
  
  #Outing best lambda
  pred_data=cbind(pred_data, best_lambda=best_lambda)
#end of loop 
return( pred_data=pred_data)
}


#******************************************************************************
#******************************************************************************
#.                        FUNCTION 3: ESTIMATION PLOT FUNCTION ####
#******************************************************************************
#******************************************************************************

#Cost Plot Function
#************************************
plot_func=function(obs_data=obs_cost,
                         est_no_iptw=est_cost_no_IPTW,
                         est_iptw=est_cost_IPTW,
                         tot_mon=c(1:24),title=c("Plot"), 
                         ylim=c(0,2000), y_lab="Cost Diff ($)",
                         CI_IPTW=1, # 1 for Yes;
                         std_iptw=std_cost_IPTW,
                         no_IPTW=1, # 1 for Yes;
                         h_abl=0 # 1 for abline-h;
                         ){
  par(mar = c(5,5,4,2))

  #Observed cost
  plot( tot_mon, obs_data ,type="l",lwd=4,lty=1, col=1,
        main = title,
        xlab="Months",
        ylab=y_lab,
        ylim=ylim,
        yaxt="none",xaxt="none", cex.lab=2, cex.main=2)
  #If cost difference;
    axis(2, at = c(ylim[1],0,ylim[2]), 
       labels =  c(paste0(ylim[1]/1000,"K"), "0K",paste0(ylim[2]/1000,"K")), 
       cex.axis=2)
    axis(1, at = c(1,6,12,18,24), labels =c(1,6,12,18,24), cex.axis=2)
    if (h_abl==1){
      abline(h=0, col="gray",lty=2,lwd=1)
      }
  
  #Estimation with IPTW
  lines(tot_mon,est_iptw, lwd=5,col="red",lty=2) 
  
  # CI for IPTW
  if (CI_IPTW==1){
    ci_L_iptw= est_iptw - 2*std_iptw
    ci_U_iptw= est_iptw + 2*std_iptw
    
    lines(tot_mon,ci_L_iptw, lwd=1,col="black",lty=3)
    lines(tot_mon,ci_U_iptw, lwd=1,col="black",lty=3)
  } 
  if (no_IPTW==1){  #Estimation without IPTW
  lines(tot_mon,est_no_iptw, lwd=5,col="blue",lty=3)
  } 
}

