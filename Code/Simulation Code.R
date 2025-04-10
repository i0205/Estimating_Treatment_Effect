library(MASS)
library(splines2)
library(mgcv)
library(dplyr)
library(stats)
library(weights)
library(reldist)

#******************************************************************************
#                         SIMULATION STUDIES  - SETUP                      ####
#******************************************************************************

#Defining the Simulation Function
sim_fun=function(n=2000, #sample size in a single simulation
                 months=24, #study duration in months
                 delta_e=c(0,1,-1,0), # for survival time
                 delta_1=c(2,1,2,1), #for treatment selection
                 delta_2=c(1,2,1,2), #for treatment selection
                 delta_b0= c(100,100,100,100), #for patient level cost variation
                 sigma=25, #for error
                 changpoint = c(3,8), #population level changepoints
                 Delta_t_1=function(t){1.5}, #difference between Trt 1 and Control
                 Delta_t_2=function(t){0.5}, #difference between Trt 1 and Control
                 M=1 # No. of simulation run
                ){

  #Storage Matrix
  #*********************************
  #Accumulative 
  est_cum_cost_ATE_trt_0 <- est_cum_cost_unweighted_trt_0 <- 	act_cum_cost_trt_0 <- true_cum_cost_trt_0 <-
    est_cum_cost_ATE_trt_1 <- est_cum_cost_unweighted_trt_1 <- 	act_cum_cost_trt_1 <- true_cum_cost_trt_1 <-
    est_cum_cost_ATE_trt_2 <- est_cum_cost_unweighted_trt_2 <- 	act_cum_cost_trt_2 <- true_cum_cost_trt_2 <-	
    est_cum_cost_ATE_diff_1_0 <- est_cum_cost_unweighted_diff_1_0 <- 	act_cum_cost_diff_1_0 <- true_cum_cost_diff_1_0 <-
    est_cum_cost_ATE_diff_2_1 <- est_cum_cost_unweighted_diff_2_1 <- 	act_cum_cost_diff_2_1 <- 	true_cum_cost_diff_2_1 <-
    est_cum_cost_ATE_diff_2_0 <- est_cum_cost_unweighted_diff_2_0 <- 	act_cum_cost_diff_2_0 <- true_cum_cost_diff_2_0 <-
   #Monthly
      est_mon_cost_ATE_trt_0 <- est_mon_cost_unweighted_trt_0 <- 	act_mon_cost_trt_0 <- true_mon_cost_trt_0 <-
      est_mon_cost_ATE_trt_1 <- est_mon_cost_unweighted_trt_1 <- 	act_mon_cost_trt_1 <- true_mon_cost_trt_1 <-
      est_mon_cost_ATE_trt_2 <- est_mon_cost_unweighted_trt_2 <- 	act_mon_cost_trt_2 <- true_mon_cost_trt_2 <-	
      est_mon_cost_ATE_diff_1_0 <- est_mon_cost_unweighted_diff_1_0 <- 	act_mon_cost_diff_1_0 <- true_mon_cost_diff_1_0 <-
      est_mon_cost_ATE_diff_2_1 <- est_mon_cost_unweighted_diff_2_1 <- 	act_mon_cost_diff_2_1 <- 	true_mon_cost_diff_2_1 <-
      est_mon_cost_ATE_diff_2_0 <- est_mon_cost_unweighted_diff_2_0 <- 	act_mon_cost_diff_2_0 <- true_mon_cost_diff_2_0 <-    
    matrix(NA,nrow=months+1,ncol=M)
  
  #The Simulation loop begins
  #*********************************
  for (m in 1:M){ 
    
#******************************************************************************
#   STEP 1: SIMULATING THE DATA ####
#******************************************************************************
    #COVARIATES
    #*******************
    sim_X=data.frame(cbind(
      PATIENT_ID=c(1:n),
      X1=rbinom(n,1,.5)-0.5,
      X2=rbinom(n,1,.5)-0.5,
      X3=rnorm(n,0,1),
      X4=rnorm(n,0,1)
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
    CUM_TOT_DATA=data.frame(matrix(rep(NA, 9*n*months), nrow=n*months, 
                       ncol=9, dimnames=list(c(NULL),
                              c("PATIENT_ID", "MONTH","MON_TOT_CHRG", 
                                 "TRUE_MON_TOT_CHRG_0", "TRUE_MON_TOT_CHRG_1",
                                 "TRUE_MON_TOT_CHRG_2", "MON_TOT_ALIVE", "D", "TREATMENT"
                           ))))
    k=0
    for (i in 1:n){
      for (t in 1:months){
        k=k+1
        
        CUM_TOT_DATA[k,"PATIENT_ID"]=i
        CUM_TOT_DATA[k,"MONTH"]=t
        CUM_TOT_DATA[k,"D"]=sim_X$D[i]
        CUM_TOT_DATA[k,"TREATMENT"]=sim_X$TREATMENT[i]
        CUM_TOT_DATA[k,"MON_TOT_ALIVE"]=sim_X$MONTH_TOT[i]
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
        CUM_TOT_DATA[k,"MON_TOT_CHRG"]=ifelse(t<=time,
                                              ifelse(sim_X$TREATMENT[i]==0, Ft(t,time) + abs(s_error), #trt=0
                                              ifelse(sim_X$TREATMENT[i]==1, Delta_t_1(t)*( Ft(t, time) + abs(s_error)), #trt=1 
                                                     Delta_t_2(t)*(Ft(t, time) + abs(s_error)))), #trt=2 
                                             0)
   #calculating true monthly treatment under treatment
        CUM_TOT_DATA[k,"TRUE_MON_TOT_CHRG_0"]=ifelse(t<=time_cont_trt2,Ft(t,time_cont_trt2) ,0) #trt=0
        CUM_TOT_DATA[k,"TRUE_MON_TOT_CHRG_1"]=ifelse(t<=time_trt1,Delta_t_1(t)*Ft(t,time_trt1),0) #trt=1: +6 months of survival
        CUM_TOT_DATA[k,"TRUE_MON_TOT_CHRG_2"]=ifelse(t<=time_cont_trt2,Delta_t_2(t)*Ft(t,time_cont_trt2),0) #trt=2
      }
    }
    
   #Creating all months data
   #*********************************
    dummy_data=sim_X[rep(seq_len(nrow(sim_X)), each=months+1),]
    dummy_data$MONTH=rep(c(0:months),length(unique(dummy_data$PATIENT_ID)))
    
    sim_monthly_X=merge(dummy_data,subset(CUM_TOT_DATA, select = -c(TREATMENT)), by=c("PATIENT_ID","MONTH"), 
                        all.x = TRUE)
    sim_monthly_X[is.na(sim_monthly_X)] <- 0
    
   #Simulated Accumulative Cost
   #*********************************
    for (i in 1:n) {
      for (t in 1:(months + 1)) {
        l=(i-1)*(months + 1) + t #index for calculation
        
        sim_monthly_X$CUM_TOT_CHRG[l]=ifelse(sim_monthly_X$MONTH[l]==0,
                                        sim_monthly_X$MON_TOT_CHRG[l],
                                        sim_monthly_X$MON_TOT_CHRG[l] + sim_monthly_X$CUM_TOT_CHRG[l-1])
        
        sim_monthly_X$TRUE_CUM_TOT_CHRG_0[l]=ifelse(sim_monthly_X$MONTH[l]==0,
                                               sim_monthly_X$TRUE_MON_TOT_CHRG_0[l],
                                               sim_monthly_X$TRUE_MON_TOT_CHRG_0[l] + sim_monthly_X$TRUE_CUM_TOT_CHRG_0[l-1])
        
        sim_monthly_X$TRUE_CUM_TOT_CHRG_1[l]=ifelse(sim_monthly_X$MONTH[l]==0,
                                                sim_monthly_X$TRUE_MON_TOT_CHRG_1[l],
                                                sim_monthly_X$TRUE_MON_TOT_CHRG_1[l] + sim_monthly_X$TRUE_CUM_TOT_CHRG_1[l-1])
        
        sim_monthly_X$TRUE_CUM_TOT_CHRG_2[l]=ifelse(sim_monthly_X$MONTH[l]==0,
                                                sim_monthly_X$TRUE_MON_TOT_CHRG_2[l],
                                                sim_monthly_X$TRUE_MON_TOT_CHRG_2[l] + sim_monthly_X$TRUE_CUM_TOT_CHRG_2[l-1])
      }
    }

#******************************************************************************
#. STEP 2: DEFINING GAM MODEL####
#******************************************************************************
    gam_func= function(cost_var=CUM_TOT_CHRG, weigts=ps_weights) {
      GAM_MOD <- gam(cost_var ~ Beta1+	Beta2+	Beta3+	
                       Beta4+	Beta5+	Beta6+	Beta7+	Beta8+	
                       Beta9+	Beta10+	Beta11+	Beta12+	Beta13+	
                       Beta14+	Beta15+	Beta16+	Beta17+	Beta18+	
                       Beta19+	Beta20+	Beta21+	Beta22+	Beta23+	
                       Beta24+	Beta25+	Beta26+	Beta27+     
                              Gamma_1_1+	Gamma_1_2+	Gamma_1_3+	
                              Gamma_1_4+	Gamma_1_5+	Gamma_1_6+	Gamma_1_7+	Gamma_1_8+	
                              Gamma_1_9+	Gamma_1_10+	Gamma_1_11+	Gamma_1_12+	Gamma_1_13+	
                              Gamma_1_14+	Gamma_1_15+	Gamma_1_16+	Gamma_1_17+	Gamma_1_18+	
                              Gamma_1_19+	Gamma_1_20+	Gamma_1_21+	Gamma_1_22+	Gamma_1_23+	
                              Gamma_1_24+	Gamma_1_25+	Gamma_1_26+	Gamma_1_27+
                              Gamma_2_1+	Gamma_2_2+	
                              Gamma_2_3+	Gamma_2_4+	Gamma_2_5+	Gamma_2_6+	Gamma_2_7+	
                              Gamma_2_8+	Gamma_2_9+	Gamma_2_10+	Gamma_2_11+	Gamma_2_12+	
                              Gamma_2_13+	Gamma_2_14+	Gamma_2_15+	Gamma_2_16+	Gamma_2_17+	
                              Gamma_2_18+	Gamma_2_19+	Gamma_2_20+	Gamma_2_21+	Gamma_2_22+	
                              Gamma_2_23+	Gamma_2_24+	Gamma_2_25+	Gamma_2_26+Gamma_2_27 , 
                            method = "REML",  weights = weigts ,
                            correlation = corARMA(form = ~ MONTH|PATIENT_ID, 
                                                  p = 1), data=data_gamma)
    }

#******************************************************************************
#   STEP 3: I-SPLINE AND CUMULATIVE COST ESTIMATION ####
#******************************************************************************
    #Input Data
    knot_start=min(sim_monthly_X$MONTH)+1
    knot_end=max(sim_monthly_X$MONTH)-1
    
    #Beta I-pline
    I_spline_beta=iSpline(x=sim_monthly_X$MONTH, df = NULL, knots = c(knot_start:knot_end), 
                          degree = 3L,intercept = T)
    colnames(I_spline_beta)=paste0("Beta",colnames(I_spline_beta))
    
    #Adding Gamma_i I-splines

    #TREATMENT=1
    TRT_1=ifelse(sim_monthly_X$TREATMENT==1, 1, 0)
    I_spline_gamma_1=TRT_1*I_spline_beta
    colnames(I_spline_gamma_1)=paste0("Gamma_1_",1:ncol(I_spline_gamma_1))
    
    #TREATMENT=2
    TRT_2=ifelse(sim_monthly_X$TREATMENT==2, 1, 0)
    I_spline_gamma_2=TRT_2*I_spline_beta
    colnames(I_spline_gamma_2)=paste0("Gamma_2_",1:ncol(I_spline_gamma_2))
    
    #ALL Splines combined
    I_spline=as.matrix(cbind(I_spline_beta,  I_spline_gamma_1, I_spline_gamma_2))
    data_gamma=cbind(sim_monthly_X,I_spline_beta, I_spline_gamma_1, I_spline_gamma_2)

    #IPTW
    #******************************************************************************
    G2_casestudy_ATE_IS=gam_func(cost_var=data_gamma$CUM_TOT_CHRG, weigts=data_gamma$ps_weights)
   
    #Without IPTW
    #******************************************************************************
    G2_casestudy_unweighted_IS=gam_func(cost_var=data_gamma$CUM_TOT_CHRG, weigts=data_gamma$no_weight)
    
    
    #******************************************************************************
    #Estimations based on I-spline####
    #******************************************************************************
    #Estimating based on I-spline for t=1:24
    
    is=iSpline(x=c(0:months), df = NULL, knots = c(knot_start:knot_end), 
               degree = 3L,intercept = T)
    
    #TRT 0 
    is_trt_0=data.frame(cbind(is, 
                              matrix(data=0,nrow=nrow(is), ncol=ncol(is)), 
                              matrix(data=0,nrow=nrow(is), ncol=ncol(is))))
    colnames(is_trt_0)=colnames(I_spline)
    
    #TRT 1 
    is_trt_1=data.frame(cbind(is,
                              is, 
                              matrix(data=0,nrow=nrow(is), ncol=ncol(is))))
    colnames(is_trt_1)=colnames(I_spline)
    
    #TRT 2 
    is_trt_2=data.frame(cbind(is,
                              matrix(data=0,nrow=nrow(is), ncol=ncol(is)),
                              is))
    colnames(is_trt_2)=colnames(I_spline)
    
    #Diff TRT 0 vs TRT 1 : cost for chemo after surgery (TRT 1 -TRT 0)
    is_diff_1_0=data.frame(cbind(matrix(data=0,nrow=nrow(is), ncol=ncol(is)), 
                                is, 
                                matrix(data=0,nrow=nrow(is), ncol=ncol(is))))
    colnames(is_diff_1_0)=colnames(I_spline)
    
    
    #Diff TRT 2 vs TRT 1 : cost of chemo before surgery (TRT 2-TRT 1)
    is_diff_2_1=data.frame(cbind(matrix(data=0,nrow=nrow(is), ncol=ncol(is)), 
                                -is,
                                is))
    colnames(is_diff_2_1)=colnames(I_spline)
    
    #Diff TRT 2 vs TRT 0 : cost of chemo before & after surgery (TRT 2-TRT 0)
    is_diff_2_0=data.frame(cbind(matrix(data=0,nrow=nrow(is), ncol=ncol(is)), 
                                matrix(data=0,nrow=nrow(is), ncol=ncol(is)),
                                is))
    colnames(is_diff_2_0)=colnames(I_spline)
    
    ##TRT 0 
    #******************************************************************
    
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_trt_0[,m]= predict(G2_casestudy_ATE_IS,newdata=is_trt_0, se.fit = T)$fit

    #Unweighted Est
    #**************
    est_cum_cost_unweighted_trt_0[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_trt_0, se.fit = T)$fit
    
    #Actual Cost
    #**************
    act_cum_cost_trt_0[,m]=aggregate(subset(data_gamma,TREATMENT==0)$CUM_TOT_CHRG, 
                             by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    true_cum_cost_trt_0[,m]=aggregate(data_gamma$TRUE_CUM_TOT_CHRG_0, 
                              by=list(data_gamma$MONTH), FUN=mean)$x
    ##TRT 1 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_trt_1[,m]= predict(G2_casestudy_ATE_IS,newdata=is_trt_1, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_cum_cost_unweighted_trt_1[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_trt_1, se.fit = T)$fit
    #Actual Cost
    #**************
    act_cum_cost_trt_1[,m]=aggregate(subset(data_gamma,TREATMENT==1)$CUM_TOT_CHRG, 
                             by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_cum_cost_trt_1[,m]=aggregate(data_gamma$TRUE_CUM_TOT_CHRG_1, 
                              by=list(data_gamma$MONTH), FUN=mean)$x

    ##TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_trt_2[,m]= predict(G2_casestudy_ATE_IS,newdata=is_trt_2, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_cum_cost_unweighted_trt_2[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_trt_2, se.fit = T)$fit
    #Actual Cost
    #**************
    act_cum_cost_trt_2[,m]=aggregate(subset(data_gamma,TREATMENT==2)$CUM_TOT_CHRG, 
                             by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_cum_cost_trt_2[,m]=aggregate(data_gamma$TRUE_CUM_TOT_CHRG_2, 
                              by=list(data_gamma$MONTH), FUN=mean)$x
    
    ##Diff TRT 0 vs TRT 1 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_diff_1_0[,m]= predict(G2_casestudy_ATE_IS,newdata=is_diff_1_0, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_cum_cost_unweighted_diff_1_0[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_diff_1_0, se.fit = T)$fit
    #Actual Cost
    #**************
    act_cum_cost_diff_1_0[,m]=aggregate(subset(data_gamma,TREATMENT==1)$CUM_TOT_CHRG, 
                               by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==0)$CUM_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_cum_cost_diff_1_0[,m]= aggregate(data_gamma$TRUE_CUM_TOT_CHRG_1, 
                                     by=list(data_gamma$MONTH), FUN=mean)$x-
                            aggregate(data_gamma$TRUE_CUM_TOT_CHRG_0, 
                                      by=list(data_gamma$MONTH), FUN=mean)$x
    ##Diff TRT 1 vs TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_diff_2_1[,m]= predict(G2_casestudy_ATE_IS,newdata=is_diff_2_1, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_cum_cost_unweighted_diff_2_1[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_diff_2_1, se.fit = T)$fit
    #Actual Cost
    #**************
    act_cum_cost_diff_2_1[,m]=aggregate(subset(data_gamma,TREATMENT==2)$CUM_TOT_CHRG, 
                               by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==1)$CUM_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_cum_cost_diff_2_1[,m]= aggregate(data_gamma$TRUE_CUM_TOT_CHRG_2, 
                                     by=list(data_gamma$MONTH), FUN=mean)$x-
                          aggregate(data_gamma$TRUE_CUM_TOT_CHRG_1, 
                                    by=list(data_gamma$MONTH), FUN=mean)$x    

    ##Diff TRT 0 vs TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_cum_cost_ATE_diff_2_0[,m]= predict(G2_casestudy_ATE_IS,newdata=is_diff_2_0, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_cum_cost_unweighted_diff_2_0[,m]=predict(G2_casestudy_unweighted_IS,newdata=is_diff_2_0, se.fit = T)$fit
    #Actual Cost
    #**************
    act_cum_cost_diff_2_0[,m]=aggregate(subset(data_gamma,TREATMENT==2)$CUM_TOT_CHRG, 
                               by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==0)$CUM_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_cum_cost_diff_2_0[,m]= aggregate(data_gamma$TRUE_CUM_TOT_CHRG_2, 
                                     by=list(data_gamma$MONTH), FUN=mean)$x-
                          aggregate(data_gamma$TRUE_CUM_TOT_CHRG_0, 
                                    by=list(data_gamma$MONTH), FUN=mean)$x    
    

#******************************************************************************
#   STEP 4: M-SPLINE AND MONTHLY COST ESTIMATION ####
#******************************************************************************
    #Input Data
    knot_start=min(sim_monthly_X$MONTH)+1
    knot_end=max(sim_monthly_X$MONTH)-1
    
    #Beta I-pline
    M_spline_beta=mSpline(x=sim_monthly_X$MONTH, df = NULL, knots = c(knot_start:knot_end), 
                          degree = 3L,intercept = T)
    colnames(M_spline_beta)=paste0("Beta",colnames(M_spline_beta))
    
    #Adding Gamma_i M-splines

    #TREATMENT=1
    TRT_1=ifelse(sim_monthly_X$TREATMENT==1, 1, 0)
    M_spline_gamma_1=TRT_1*M_spline_beta
    colnames(M_spline_gamma_1)=paste0("Gamma_1_",1:ncol(M_spline_gamma_1))
    
    #TREATMENT=2
    TRT_2=ifelse(sim_monthly_X$TREATMENT==2, 1, 0)
    M_spline_gamma_2=TRT_2*M_spline_beta
    colnames(M_spline_gamma_2)=paste0("Gamma_2_",1:ncol(M_spline_gamma_2))
    
    #ALL Splines combined
    M_spline=as.matrix(cbind(M_spline_beta,  M_spline_gamma_1, M_spline_gamma_2))
    data_gamma=cbind(sim_monthly_X,M_spline_beta, M_spline_gamma_1, M_spline_gamma_2)
    
    #IPTW
    #******************************************************************************
    G2_casestudy_ATE_MS=gam_func(cost_var=data_gamma$MON_TOT_CHRG, weigts=data_gamma$ps_weights)
    
    #Without IPTW
    #******************************************************************************
    G2_casestudy_unweighted_MS=gam_func(cost_var=data_gamma$MON_TOT_CHRG, weigts=data_gamma$no_weight)
    
    
    #******************************************************************************
    #Estimations based on M-spline####
    #******************************************************************************
    #Estimating based on M-spline for t=1:24
    
    ms=mSpline(x=c(0:months), df = NULL, knots = c(knot_start:knot_end), 
               degree = 3L,intercept = T)
    
    #TRT 0 
    ms_trt_0=data.frame(cbind(ms, 
                              matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)), 
                              matrix(data=0,nrow=nrow(ms), ncol=ncol(ms))))
    colnames(ms_trt_0)=colnames(M_spline)
    
    #TRT 1 
    ms_trt_1=data.frame(cbind(ms,
                              ms, 
                              matrix(data=0,nrow=nrow(ms), ncol=ncol(ms))))
    colnames(ms_trt_1)=colnames(M_spline)
    
    #TRT 2 
    ms_trt_2=data.frame(cbind(ms,
                              matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)),
                              ms))
    colnames(ms_trt_2)=colnames(M_spline)
    
    #Diff TRT 0 vs TRT 1 : cost for chemo after surgery (TRT 1 -TRT 0)
    ms_diff_1_0=data.frame(cbind(matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)), 
                                 ms, 
                                 matrix(data=0,nrow=nrow(ms), ncol=ncol(ms))))
    colnames(ms_diff_1_0)=colnames(M_spline)
    
    
    #Diff TRT 2 vs TRT 1 : cost of chemo before surgery (TRT 2-TRT 1)
    ms_diff_2_1=data.frame(cbind(matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)), 
                                 -ms,
                                 ms))
    colnames(ms_diff_2_1)=colnames(M_spline)
    
    #Diff TRT 2 vs TRT 0 : cost of chemo before & after surgery (TRT 2-TRT 0)
    ms_diff_2_0=data.frame(cbind(matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)), 
                                 matrix(data=0,nrow=nrow(ms), ncol=ncol(ms)),
                                 ms))
    colnames(ms_diff_2_0)=colnames(M_spline)
    
    ##TRT 0 
    #******************************************************************
    
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_trt_0[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_trt_0, se.fit = T)$fit
    
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_trt_0[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_trt_0, se.fit = T)$fit
    
    #Actual Cost
    #**************
    act_mon_cost_trt_0[,m]=aggregate(subset(data_gamma,TREATMENT==0)$MON_TOT_CHRG, 
                                     by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    true_mon_cost_trt_0[,m]=aggregate(data_gamma$TRUE_MON_TOT_CHRG_0, 
                                  by=list(data_gamma$MONTH), FUN=mean)$x
    ##TRT 1 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_trt_1[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_trt_1, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_trt_1[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_trt_1, se.fit = T)$fit
    #Actual Cost
    #**************
    act_mon_cost_trt_1[,m]=aggregate(subset(data_gamma,TREATMENT==1)$MON_TOT_CHRG, 
                                     by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_mon_cost_trt_1[,m]=aggregate(data_gamma$TRUE_MON_TOT_CHRG_1, 
                                  by=list(data_gamma$MONTH), FUN=mean)$x
    
    ##TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_trt_2[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_trt_2, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_trt_2[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_trt_2, se.fit = T)$fit
    #Actual Cost
    #**************
    act_mon_cost_trt_2[,m]=aggregate(subset(data_gamma,TREATMENT==2)$MON_TOT_CHRG, 
                                     by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_mon_cost_trt_2[,m]=aggregate(data_gamma$TRUE_MON_TOT_CHRG_2, 
                                  by=list(data_gamma$MONTH), FUN=mean)$x
    
    ##Diff TRT 0 vs TRT 1 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_diff_1_0[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_diff_1_0, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_diff_1_0[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_diff_1_0, se.fit = T)$fit
    #Actual Cost
    #**************
    act_mon_cost_diff_1_0[,m]=aggregate(subset(data_gamma,TREATMENT==1)$MON_TOT_CHRG, 
                                        by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==0)$MON_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_mon_cost_diff_1_0[,m]= aggregate(data_gamma$TRUE_MON_TOT_CHRG_1, 
                                      by=list(data_gamma$MONTH), FUN=mean)$x-
      aggregate(data_gamma$TRUE_MON_TOT_CHRG_0, 
                by=list(data_gamma$MONTH), FUN=mean)$x
    
    ##Diff TRT 1 vs TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_diff_2_1[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_diff_2_1, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_diff_2_1[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_diff_2_1, se.fit = T)$fit
    #Actual Cost
    #**************
    act_mon_cost_diff_2_1[,m]=aggregate(subset(data_gamma,TREATMENT==2)$MON_TOT_CHRG, 
                                        by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==1)$MON_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==1)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_mon_cost_diff_2_1[,m]= aggregate(data_gamma$TRUE_MON_TOT_CHRG_2, 
                                      by=list(data_gamma$MONTH), FUN=mean)$x-
      aggregate(data_gamma$TRUE_MON_TOT_CHRG_1, 
                by=list(data_gamma$MONTH), FUN=mean)$x    
    
    ##Diff TRT 0 vs TRT 2 
    #******************************************************************
    #Weighted ATE 
    #**************
    est_mon_cost_ATE_diff_2_0[,m]= predict(G2_casestudy_ATE_MS,newdata=ms_diff_2_0, se.fit = T)$fit
    #Unweighted Est
    #**************
    est_mon_cost_unweighted_diff_2_0[,m]=predict(G2_casestudy_unweighted_MS,newdata=ms_diff_2_0, se.fit = T)$fit
    #Actual Cost
    #**************
    act_mon_cost_diff_2_0[,m]=aggregate(subset(data_gamma,TREATMENT==2)$MON_TOT_CHRG, 
                                        by=list(subset(data_gamma,TREATMENT==2)$MONTH), FUN=mean)$x-
      aggregate(subset(data_gamma,TREATMENT==0)$MON_TOT_CHRG, 
                by=list(subset(data_gamma,TREATMENT==0)$MONTH), FUN=mean)$x
    #True Cost
    #**************
    true_mon_cost_diff_2_0[,m]= aggregate(data_gamma$TRUE_MON_TOT_CHRG_2, 
                                      by=list(data_gamma$MONTH), FUN=mean)$x-
      aggregate(data_gamma$TRUE_MON_TOT_CHRG_0, 
                by=list(data_gamma$MONTH), FUN=mean)$x    
    
} #end of simulation loop (M loop)
  
  return(list(  #Accumulative
                est_cum_cost_ATE_trt_0= est_cum_cost_ATE_trt_0,
                est_cum_cost_unweighted_trt_0= est_cum_cost_unweighted_trt_0,
                act_cum_cost_trt_0= act_cum_cost_trt_0,
                true_cum_cost_trt_0=true_cum_cost_trt_0,
                
                est_cum_cost_ATE_trt_1= est_cum_cost_ATE_trt_1,
                est_cum_cost_unweighted_trt_1= est_cum_cost_unweighted_trt_1,
                act_cum_cost_trt_1= act_cum_cost_trt_1,
                true_cum_cost_trt_1=true_cum_cost_trt_1,
                
                est_cum_cost_ATE_trt_2= est_cum_cost_ATE_trt_2,
                est_cum_cost_unweighted_trt_2= est_cum_cost_unweighted_trt_2,
                act_cum_cost_trt_2= act_cum_cost_trt_2,
                true_cum_cost_trt_2=true_cum_cost_trt_2,
                
                est_cum_cost_ATE_diff_1_0= est_cum_cost_ATE_diff_1_0,
                est_cum_cost_unweighted_diff_1_0= est_cum_cost_unweighted_diff_1_0,
                act_cum_cost_diff_1_0= act_cum_cost_diff_1_0,
                true_cum_cost_diff_1_0= true_cum_cost_diff_1_0,
                
                est_cum_cost_ATE_diff_2_1= est_cum_cost_ATE_diff_2_1,
                est_cum_cost_unweighted_diff_2_1= est_cum_cost_unweighted_diff_2_1,
                act_cum_cost_diff_2_1= act_cum_cost_diff_2_1,
                true_cum_cost_diff_2_1= true_cum_cost_diff_2_1,
                
                est_cum_cost_ATE_diff_2_0= est_cum_cost_ATE_diff_2_0,
                est_cum_cost_unweighted_diff_2_0= est_cum_cost_unweighted_diff_2_0,
                act_cum_cost_diff_2_0= act_cum_cost_diff_2_0,
                true_cum_cost_diff_2_0= true_cum_cost_diff_2_0,
                
                #Monthly
                est_mon_cost_ATE_trt_0= est_mon_cost_ATE_trt_0,
                est_mon_cost_unweighted_trt_0= est_mon_cost_unweighted_trt_0,
                act_mon_cost_trt_0= act_mon_cost_trt_0,
                true_mon_cost_trt_0=true_mon_cost_trt_0,
                
                est_mon_cost_ATE_trt_1= est_mon_cost_ATE_trt_1,
                est_mon_cost_unweighted_trt_1= est_mon_cost_unweighted_trt_1,
                act_mon_cost_trt_1= act_mon_cost_trt_1,
                true_mon_cost_trt_1=true_mon_cost_trt_1,
                
                est_mon_cost_ATE_trt_2= est_mon_cost_ATE_trt_2,
                est_mon_cost_unweighted_trt_2= est_mon_cost_unweighted_trt_2,
                act_mon_cost_trt_2= act_mon_cost_trt_2,
                true_mon_cost_trt_2=true_mon_cost_trt_2,
                
                est_mon_cost_ATE_diff_1_0= est_mon_cost_ATE_diff_1_0,
                est_mon_cost_unweighted_diff_1_0= est_mon_cost_unweighted_diff_1_0,
                act_mon_cost_diff_1_0= act_mon_cost_diff_1_0,
                true_mon_cost_diff_1_0= true_mon_cost_diff_1_0,
                
                est_mon_cost_ATE_diff_2_1= est_mon_cost_ATE_diff_2_1,
                est_mon_cost_unweighted_diff_2_1= est_mon_cost_unweighted_diff_2_1,
                act_mon_cost_diff_2_1= act_mon_cost_diff_2_1,
                true_mon_cost_diff_2_1= true_mon_cost_diff_2_1,
                
                est_mon_cost_ATE_diff_2_0= est_mon_cost_ATE_diff_2_0,
                est_mon_cost_unweighted_diff_2_0= est_mon_cost_unweighted_diff_2_0,
                act_mon_cost_diff_2_0= act_mon_cost_diff_2_0,
                true_mon_cost_diff_2_0= true_mon_cost_diff_2_0,
                
                sim_mon_cost_data=data_gamma
                )
         )
  
} #end of simulation function (sim_fun)


#******************************************************************************
#                         SIMULATION RUN                                   ####
#******************************************************************************
set.seed(0707)
data_F1_S1_Sig10=sim_fun(n=2000, #sample size in a single simulation
                         months=24, #study duration in months
                         delta_e=c(0,1,-1,0), #for survival time
                         delta_1=c(2,1,2,1), #for Trt 1 selection
                         delta_2=c(1,2,1,2), #for Trt 2 selection
                         delta_b0= c(100,100,100,100), #for patient level cost variation
                         sigma=25, #for error
                         changpoint = c(3,8), #population level changepoints
                         Delta_t_1=function(t){1.5}, #difference between Trt 1 and Control
                         Delta_t_2=function(t){0.5}, #difference between Trt 1 and Control
                         M=1000 # No. of simulation run
                         )

#Saving on Simulation dataset for reference
mon_data=data_F1_S1_Sig10$sim_mon_cost_data