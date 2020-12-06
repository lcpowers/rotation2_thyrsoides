# Run through code through figure 6, changing growth rates

# Clear environment
rm(list=ls())

# Load packages
library(tidyverse)
library(nlme)
library(MASS)
library(cowplot)

setwd("~/Desktop/IQ/Rotation2/rotation2_thyrsoides")

## Part I. Fitting Models
#### Read first in dataset
# ===================================================
dataf = data.frame(read.table("ct.ipm.txt",header=T))

# Filter full dataframe for plants with 1 rosette
dataf_r1 = filter(dataf,ros==1)
attach(dataf_r1)

dataf_r1$size.t.all = log(nl.t*ll.t) ## plant sizes in year t
dataf_r1$size.t1.all = log(nl.t1*ll.t1) ## plant sizes in year t+1

## The next three are unnecessary since they are just duplicating data in another column, but keeping for consistency at this point
dataf_r1$flow.all = flow ## Logical vector if flowering: 1 = yes
dataf_r1$surv.all = surv ## Logical vector of survival: 1 = yes
dataf_r1$site.all = site ## Site vector

site.code.l=c("FU","SP") ## Used later for plots
pch.code=c(19,1)

# Detach old version, attach new version
detach(dataf_r1)
attach(dataf_r1)

all.sizes=c(size.t.all[flow.all==0],size.t1.all[year.t==2005]) ## all plant sizes
all.site=c(site.all[flow.all==0],site.all[year.t==2005]) ## All sites -- how is this one used?
detach(dataf_r1)

# Create growth_df
growth_df <- dataf_r1 %>% 
  filter(flow.all==0) %>% 
  dplyr::select(size.t = size.t.all,
                size.t1 = size.t1.all,
                site.s = site.all) %>% 
  filter(complete.cases(.))

sites = c("FU","SP")

growth_df$growth <- growth_df$size.t1 - growth_df$size.t

sets = 2
perc_change = seq(-0.75,0.75,0.01)

master_df = NULL

for(set in 1:sets){
  
  timestamp()
  for(ii in 1:length(perc_change)){
    
    # Print progress 
    print(noquote(paste0("rep ", set, ":", ii)))
    
    # If change in growth is negative, subtract to decrease growth, if zero or positive add to increase                         
     if(perc_change[ii] < 0) {
       growth_df$size.t1_new <- growth_df$size.t + growth_df$growth*(1 - abs(perc_change[ii]))
     } else{ growth_df$size.t1_new <- growth_df$size.t + growth_df$growth*(1 + abs(perc_change[ii]))}

    # start building tmp dataframe  
    tmp_df = growth_df %>% 
      group_by(site.s) %>% 
      summarise(mean.size.t = mean(size.t),
                mean.size.t1 = mean(size.t1_new),
                .groups = "keep")
    
    tmp_df$perc.change = perc_change[ii]

    ######### refit model with size and site main effects, and site specific decreasing variance #########
    fit.grow.gls<-gls(size.t1_new ~ site.s+size.t-1, ## Why -1?
                      na.action=na.omit,
                      weight=varExp(form=~fitted(.)|site.s),
                      method="ML",
                      data=growth_df)
    
    tmp_df$g.intercepts = fit.grow.gls$coef[1:2]
    tmp_df$g.slopes = rep(fit.grow.gls$coef[3],2)
    
    g.intercepts=fit.grow.gls$coef[1:2] ## Growth intercepts for each site
    g.slopes=rep(fit.grow.gls$coef[3],2) ## Growth slopes for each site
    var.exp.coef=fit.grow.gls$modelStruct$varStruct ## coef for variance structure
    sigma.g=fit.grow.gls$sigma ## Residual standard error
    
    #### Flowering
    flower_df <- dataf_r1 %>% 
      dplyr::select(flow.s = flow.all,
                    site.s = site.all,
                    size.t = size.t.all) %>% 
      filter(complete.cases(.))
    
    attach(flower_df)
    # table(site.s)
    
    store.size.flow=size.t[flow.s==1] # store size of plants that flowered
    store.site.flow=site.s[flow.s==1] # store sites of plants that flowered
    
    fit.flow=glm(flow.s~site.s/size.t-1,family=binomial) ## Why '/'  and why -1?
    
    f.intercepts=fit.flow$coef[1:2] ## Model intercepts
    f.slopes=c(fit.flow$coef[3:4]) ## Model slopes
    
    site.flow.SE=summary(fit.flow)$coef[5:6]
    
    
    #### plot flowering
    FUcol <- "black"
    SPcol <- "blue"
    
    flower_df$probability <- predict(fit.flow, site.s=flower_df$site.s, size.t = flower_df$size.t, type = "response")
    predict_df <- data.frame(site.s=c(rep("FU",1000),rep("SP",1000)),
                             size.t=rep(seq(2,10,length.out = 1000),2))
    predict_df <- cbind(predict_df,predict(fit.flow, newdata = predict_df, type = "response", se.fit = TRUE))
    
    ggplot(data = predict_df, aes(x = size.t, y = fit)) +
      geom_line(aes(color = site.s, linetype=site.s)) +
      scale_color_manual(name = "", values=c(SPcol, FUcol),labels = c("FU: n=760","SP: n=398"))+
      scale_linetype_discrete(name = "",labels = c("FU: n=760","SP: n=398"))+
      scale_x_continuous(limits = c(2, 10)) +
      # geom_point(data=flower_df,aes(x=size.t,y=probability))+
      theme_classic()+
      labs(x="Plant size year t", y = "Probability of flowering (year t+1)")+
      theme(legend.position = c(0.9, 0.2),
            legend.text = element_text(size=14))
    
    ### Don't understand where their points around the lines come from. Is that related to binning?
    detach(flower_df)
    
    #### Calculation: Survival
    
    survival_df <- dataf_r1 %>% 
      dplyr::select(surv.s = surv.all,
                    site.s = site.all,
                    size.t = size.t.all) %>% 
      filter(complete.cases(.))
    
    attach(survival_df)
    # table(site.s)
    
    fit.surv = glm(surv.s~site.s/size.t-1, family=binomial) ## Same as above -- why '/' and why -1?
    
    s.intercepts = fit.surv$coef[1:2]
    s.slopes = c(fit.surv$coef[3:4])
    
    #### Plot: Survival
    
    survival_df$probability <- predict(fit.surv, site.s=site.s, size.t = size.t, type = "response")
    
    predict_df <- data.frame(site.s=c(rep("FU",1000),rep("SP",1000)),
                             size.t=rep(seq(2,10,length.out = 1000),2))
    predict_df <- cbind(predict_df,predict(fit.surv, newdata = predict_df, type = "response", se.fit = TRUE))
    
    ggplot(data = predict_df, aes(x = size.t, y = fit)) +
      geom_line(aes(color = site.s, linetype=site.s)) +
      scale_color_manual(name = "", values=c(SPcol, FUcol),labels = c("FU: n=809","SP: n=417"))+
      scale_linetype_discrete(name = "",labels = c("FU: n=809","SP: n=417"))+
      scale_x_continuous(limits = c(2, 10)) +
      scale_y_continuous(limits = c(0,1))+
      # geom_point(data=flower_df,aes(x=size.t,y=probability))+
      theme_classic()+
      labs(x="Plant size year t", y = "Probability of survival to year t+1")+
      theme(legend.position = c(0.9, 0.2),
            legend.text = element_text(size=14))
    
    ### Don't understand where their points around the lines come from. Is that related to binning?
    detach(survival_df)
    
    #### Calculation: Fecundity
    
    IPM.fecundity=data.frame(read.table("CT.IPM.fecundity.txt", header=T))
    # names(IPM.fecundity)
    
    # Unbrowsed individuals
    fecundity_df <- IPM.fecundity %>% 
      filter(brows.t1==0) %>% 
      dplyr::select(size.t.f=size.t,
                    site.t.f=site,
                    si.t1.f=si.t1,
                    nl.t.f=nl.t)
    
    attach(fecundity_df)
    
    fit.fec = lm(log(si.t1.f)~ ## dep var = log of seeds per individual (unbrowsed) in year t+1
                   site.t.f + ## site of unbrowsed individuals
                   size.t.f-1) ## size of unbrowsed Why -1? #
    
    #### Plot: Fecundity
    
    ggplot() + 
      geom_abline(intercept = fit.fec$coefficients[1],slope = fit.fec$coefficients[3], color=FUcol)+ # FU line
      geom_abline(intercept = fit.fec$coefficients[2],slope = fit.fec$coefficients[3], linetype='longdash', color = SPcol)+ # SP line
      geom_point(data=fecundity_df,aes(x=size.t.f, y=log(si.t1.f), shape=site.t.f, color=site.t.f)) + # Size points for each site
      scale_shape_manual(name = "",values = c(16,1), labels = c("FU: n=735","SP: n=357")) +
      scale_color_manual(name = "",values = c(SPcol, FUcol), labels = c("FU: n=735","SP: n=357"))+
      theme_classic()+
      labs(x="Plant size year t", y="Viable seeds (log scale) year t+1")+
      xlim(2,10)+
      ylim(2,10) +
      theme(legend.position = c(0.8, 0.2),
            legend.text = element_text(size=14))
    
    detach(fecundity_df)
    
    #### Seedling sizes
    
    ####### Seedlings data #######
    IPM.seedlings=data.frame(read.table("CT.IPM.seedlings.txt",header=T)) %>% 
      filter(site != "JU") %>% 
      dplyr::select(seedlings.size.t=size.t,
                    seedlings.site=site)
    
    fit.seedlings=lm(seedlings.size.t~seedlings.site-1, data = IPM.seedlings) ## That -1! What is that?
    # summary(fit.seedlings)
    
    ## Add seedlings to all.sizes vector
    size.all.plus.seedlings=c(all.sizes, IPM.seedlings$size.t[IPM.seedlings$year.t!=2003])
    ## Add site codes for seedlings to all sites vector
    site.all.plus.seedlings=c(all.sizes, IPM.seedlings$site[IPM.seedlings$year.t!=2003])
    
    ####### Establishment data #######
    IPM.establishment=data.frame(read.table("IPM.establishment.data.txt",header=T)) %>% 
      dplyr::select(p.est.site=site,
                    p.est.seeds.t=seeds.t,
                    p.est.seedlings=sdl.t1)
    
    attach(IPM.establishment)
    
    ## Site specific seedling establishment rates
    est.p.est = IPM.establishment %>% 
      group_by(p.est.site) %>% 
      summarise(est_rate=sum(p.est.seedlings)/sum(p.est.seeds.t),
                .groups = "keep")
    est.p.est = as.numeric(est.p.est$est_rate)
    
    detach(IPM.establishment)
    
    #### Collect parameters
    
    ## Set upper and lower plant size bounds
    minsize<-1
    maxsize<-12
    
    # Global variables for midpoint rule approximation 
    # n.big.matrix is the number of mesh points for size, n.age is the number of age classes. 
    n.big.matrix = 250; n.age = 50; n=n.big.matrix ## 
    
    L= minsize; U= maxsize;
    
    # boundary points b and mesh points y
    b = L+c(0:n)*(U-L)/n
    y = 0.5*(b[1:n]+b[2:(n+1)])
    
    # step size for midpoint rule, see equations 4 and 5
    h = y[2]-y[1]
    
    ## initialize empty parameter matrix
    store.p.vec = matrix(NA,12,2)
    colnames(store.p.vec) <- c("FU","SP")
    p.vec.names <- c("1st survival param","2nd survival param",
                     "1st flow param    ","2nd flow param    ",
                     "ag                ","bg                ",
                     "sigma2 growth     ","intercept seeds   ",
                     "slope seeds       ","mean kids size   ",
                     "sigma2 kids size ","growth variance parameter ")
    
    ## initialize empty parameter matrix
    store.p.vec = matrix(NA,12,2)
    colnames(store.p.vec) <- c("FU","SP")
    
    ## Store parameters in store.p.vec 12x2 array
    for(i in 1:2){
      
      store.p.vec[1,i]<- as.numeric(s.intercepts[i])		# ; p.vec.names[1]<-"1st survival param";
      store.p.vec[2,i]<- s.slopes[i]				# ; p.vec.names[2]<-"2nd survival param";
      store.p.vec[3,i]<- f.intercepts[i]		# 	; p.vec.names[3]<-"1st flow param    ";
      store.p.vec[4,i]<- f.slopes[i]				# ; p.vec.names[4]<-"2nd flow param    ";
      store.p.vec[5,i]<- g.intercepts[i]		# 	; p.vec.names[5]<-"ag                ";
      store.p.vec[6,i]<- g.slopes[i]				# ; p.vec.names[6]<-"bg                ";
      store.p.vec[7,i]<- sigma.g^2				# ; p.vec.names[7]<-"sigma2 growth     ";
      store.p.vec[8,i]<- fit.fec$coef[i]   #			; p.vec.names[8]<-"intercept seeds   ";
      store.p.vec[9,i]<- fit.fec$coef[3]	 # 		; p.vec.names[9]<-"slope seeds       ";
      store.p.vec[10,i]<- fit.seedlings$coef[i]	# 	; p.vec.names[10]<-"mean kids size   ";
      store.p.vec[11,i]<- summary(fit.seedlings)$sigma^2	# ; p.vec.names[11]<-"sigma2 kids size ";
      store.p.vec[12,i]<- var.exp.coef[i]			# ; p.vec.names[12]<-"growth variance parameter ";
      
    }
    # summary(store.p.vec)
    ## Part II. Compute the kernel component functions from the fitted models
    
    #### Basic demographic functions
    ## Survival model
    
    sx <- function(x,params) {
      u <- exp(params[2]*x+params[1]); ## Logit model => odds
      return(u/(1+u)); ## odds to prob
    }
    
    ## Fecundity model
    fx<-function(x,params) {
      u<-exp(params[3]+params[4]*x) ## Logit model => odds
      return(u/(1+u))  ## odds to prob
    }
    
    ## Growth model
    gxy<-function(x,y,params) {
      mux<-params[5]+params[6]*x; ## growth slope and intercepts
      sigmax2<-(params[7])*exp(2*(params[12]*mux)) ## Variance around growth curve
      sigmax<-sqrt(sigmax2); ## Standard deviation around growth curve
      fac1<-sqrt(2*pi)*sigmax; ## ??
      fac2<-((y-mux)^2)/(2*sigmax2); ## ??
      return(exp(-fac2)/fac1); ## ??
    }
    
    ## Survival-growth function
    pxy<-function(x,y,params) {
      return(sx(x,params)*(1-fx(x,params))*gxy(x,y,params))}
    
    ## Fecundity function
    fxy<-function(x,y,params) {
      nkids<-p.est*exp(params[8]+params[9]*x) ## 8: seeds intercept; 9: seeds slope
      kidsize.mean<- params[10] ## Mean seedling size
      kidsize.var<- params[11] ## Sigma seedling size
      fac1<-sqrt(2*pi)*sqrt(kidsize.var) ## Same Q as 483 above
      fac2<-((y-kidsize.mean)^2)/(2*kidsize.var) ## 484 above
      f<-sx(x,params)*fx(x,params)*nkids*exp(-fac2)/fac1 ## End is 485 above
      return(f);
    }
    
    #### The 'big matrix' M of size n x n
    
    bigmatrix<-function(n,params) {
      
      # upper and lower integration limits
      L<-minsize; U<-maxsize;
      
      # boundary points b and mesh points y
      b<-L+c(0:n)*(U-L)/n; # Boundaries
      y<-0.5*(b[1:n]+b[2:(n+1)]); # Midpoints
      
      # construct the matrix
      I <- diag(n); ## Identity matrix
      
      ## Survival-growth matrix for all midpoint sizes
      P<-t(outer(y,y,pxy,params=params)) # P: survival-growth transpose outer product. Why transpose?
      
      ## fecundity matrix for all midpoint sizes
      B<-t(outer(y,y,fxy,params=params)) # B: fecundity
      
      ## Empty matrix of zeroes
      M=array(0,dim=c(2*n,2*n))
      
      ## Insert P matrix
      M[1:n,1:n]=P*(U-L)/n
      
      ## Insert B matrix to the right of the P matrix
      M[1:n,(n+1):(2*n)]=B*(U-L)/n
      
      ## Insert the identity matrix
      M[(n+1):(2*n),1:n]=diag(n)
      K<-M ## Call K to match paper?
      P<-(U-L)*P/n
      B<-(U-L)*B/n
      return(list(matrix=M,
                  kernel=K,
                  meshpts=y,
                  Pmatrix=P,
                  Bmatrix=B,
                  Imatrix=I)); 
    }
    
    
    R0.calc<-function(n,params){
      
      M<-bigmatrix(n,params) ## Construct matrix
      
      ## If any NAs in matrix, return NA for estimates
      if (any(is.na(M$matrix))){
        ave.R0=NA
        lam=NA
        T=NA } 
      
      else{
        N <- solve(M$Imatrix-M$Pmatrix) ## ?? Some linear algebra operation
        R <- M$Bmatrix %*% N ## Matrix multiplication
        ave.R0<-Re(eigen(R)$values[1]) ## Re: real; eigen:returns vector of eigenvalues. values[1] is the first eigenvalue in returned vector 
        lam<-Re(eigen(M$matrix)$values[1]); ## Re: real; eigen:returns vector of eigenvalues. values[1] is the first eigenvalue in returned vector 
        T=log(ave.R0)/log(lam) ## ?? understand why this is true
      }
      
      return(list(lam=lam,ave.R0=ave.R0,T=T))
    }
    
    R0.betas<-function(x){
      
      p.vec[3] <- x; ## p.vec[3] = 1st flow pram. Assign value of x to that
      nR0 <- R0.calc(n.big.matrix, p.vec) 
      return(nR0$ave.R0)
      
    }
    
    
    gen.time=rep(NA,2)
    
    tmp_df$lam = -9999
    tmp_df$R0 = -9999
    tmp_df$gentime = -9999
    tmp_df$ESS = -9999
    
    for(i in 1:2){
      #if(i==1) p.est= 8.604605e-05 else p.est=0.0001655622  # assuming dd-reg
      if(i==1) p.est= est.p.est[1] else p.est=est.p.est[2]   # actual ## site specific parameters
      p.vec=store.p.vec[,i] ## Site-specific parameters
      tmp=R0.calc(n.big.matrix,p.vec) 
      gen.time[i]=tmp$T ## get T from list of value returned from R0.calc
      tmp_df$lam[i] = tmp$lam
      tmp_df$R0[i] = tmp$ave.R0
      tmp_df$gentime[i] = tmp$T
      tmp_df$ESS[i] = optimize(R0.betas, c(-100,10), maximum=T, tol=0.01)$maximum
      # cat("Site ",i," lambda =",tmp$lam," R0 =",tmp$ave.R0," Generation time =",tmp$T,"\n") ## Print results
      # cat("ESS intercept ", optimize(R0.betas, c(-100,10), maximum=T, tol=0.01)$maximum,"\n") ## Print results
    }
    
    master_df = rbind(master_df,tmp_df)
    rm(tmp_df)
    
  }
}

master_df.store <- master_df

library(reshape2)
plot_df <- melt(data = master_df,id.vars = c(1,4),measure.vars = c(2:3,5:10)) %>% 
  filter(!(variable %in% c('g.slopes')))

ggplot(plot_df, aes(x = mean.size.t, y=value,color=site.s))+
  #geom_point()+
  geom_smooth()+
  facet_wrap(~variable,scales = 'free')+
  theme_classic()
