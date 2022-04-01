##### Script Version December 17 2021 - Bruno Levecke                                                                                                         ##########  
##### This script is part of the study entitled 'A general framework to make more evidence-based choices on fecal egg count methods and study designs to      ########## 
##### inform large-scale STH deworming programs – monitoring the therapeutic drug efficacy as a case study' by Johnny Vlaminck and colleagues.                ##########
##### It describes the methodology for both power analysis and the total operational costs to monitor drug efficacy.                                          ##########
########################################################################################################################################################################

print(Sys.time())

wd <- '/Users/brunolevecke/Dropbox/Grants/Bill Gates/2014/Ghent/Surveillance of AR/Publications/2020_WP1_timing&Cost/Analysis/'
wd <- getwd()

### 1. Libraries 
library(plyr)
library(MethComp)
library(prevalence)
library(data.table)
library(ggplot2)
library(readr)
library(plyr)
library(ggplot2)
library(ggpubr)
#library(xlsx)
library(writexl)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9999CC","#999999")

### 2. Data sets to be used in functions 
## CV between (cv_pop) and within subjects (day-to-day variation; cv_day), and between slides made from the same sample (slide-to-slide variation; cv_slide)
cv_pop <- c(3,1.1,3)                                                              # CV for mean FECs between individuals from the same population - based on Denwood/Coffeng
cv_day <- c(1.25,0.75,1.25)                                                       # CV for mean FECs within individuals across days - based on Denwood/Coffeng
cv_slide <- c(1,0.7,1)                                                            # CV for mean FECs within slides made from the same sample - based on Denwood/Coffeng
sth <- c('Ascaris','Trichuris','Hookworm')                                        # different STH species
data_cv <- data.frame(sth=sth,cv_pop=cv_pop,cv_day=cv_day,cv_slide=cv_slide)      # create summary table

## Output of analysis conducted by Luc E. Coffeng, Erasmus PC, Rotterdam (see Supplementary Info X)
sth <- c(rep('Ascaris',3),rep('Trichuris',3),rep('Hookworm',3))                                                           # different STH species
tech <- rep(c('KK','MF','FP'),3)                                                                                          # different methods: Kato-Katz (KK); mini-FLOTAC (MF), FECPAKG2 (FP)
n_slide <- rep(1,9)                                                                                                       # number of slides/devices prepared per sample; single (1); duplicate (2)
f <- c(1,0.645229010466808,0.248104909312327,1, 1.005206869784,0.151794514457339,1,0.801133520852417,0.56924106904857)  # egg counts (in EPG) relative to those of duplicate KK
k <- c(NA,0.578550858176212,0.519916243000555,NA,3.022465486,0.706335595,NA,1.465012279,0.574262714)                      # overdispersion parameter
m <- rep(c(1/24,1/10,1/34),3)                                                                                             # mass of stool examined (in gram) 
data_rel <- data.frame(sth=sth,n_slide=n_slide,m=m,tech=tech,f=f,k=k)                                                 # create summary table

## Time to process an individual stool sample for three different diagnostic methods (current manuscript)
tech <- c('KK','KK','MF','FP','MF','FP')                                                                        # different methods: Kato-Katz (KK); mini-FLOTAC (MF), FECPAKG2 (FP)
n_slide <- c(2,1,1,1,2,2)                                                                                       # number of slides/devices prepared per sample; single (1); duplicate (2)
prep <- c(133,67,128,590,192,1004)                                                                              # time (in sec) to prepare one sample 
demo <- c(12,12,12,31,12,31)                                                                                    # time (in sec) to enter demographic of one subject
count <- c(16,8,8,0,16,0)                                                                                       # time (in sec) to enter egg count data of one subject
int <- c(2.38961691,2.38961691,2.51542336,1.8348640,2.51542336,1.8348640)                                       # intercept of linear regression model log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
coef <- c(0.06612497,0.06612497,0.06609474,0.1730919,0.06609474,1.8348640)                                      # coefficient of linear regression model log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
m <- c(1/12,1/24,1/10,1/34,1/5,1/17)                                                                            # mass of stool (in gram) examined
n_lab <- rep(3,6)                                                                                               # number of laboratory technicians
h_lab <- rep(4,6)                                                                                               # number of hours the laboratory technicians can process samples
data_time <- data.frame(tech=tech,m=m,prep=prep,demo=demo,count=count,n_slide=n_slide,                          # create summary table
                        int=int,coef=coef, n_lab=n_lab,h_lab=h_lab)                                             

print(Sys.time())

## Cost for car rental, per diem of personnel (Let et al., 2018) and consumables for the three diagnostic methods (current manuscript)
driver <- rep(90,6)                                                                                             # cost for car rental (includes per diem driver)
nurse <- rep(22.5,6)                                                                                            # per diem nurse
lab <- rep(22.5,6)                                                                                              # per diem laboratory technician
tech <- c('KK','KK','MF','FP','MF','FP')                                                                        # different methods: Kato-Katz (KK); mini-FLOTAC (MF), FECPAKG2 (FP)
n_slide <- c(2,1,1,1,2,2)                                                                                       # number of slides/devices prepared per sample; single (1); duplicate (2)
sample <- rep(c(0.57),6)                                                                                        # cost for consumables for collecting one sample
test <- c(1.51,1.37,1.51,1.69,1.87,2.73)                                                                        # cost for consumables for testing one sample 
data_cost <- data.frame(tech=tech,test=test,n_slide=n_slide,sample=sample,driver=driver, nurse=nurse,lab=lab)   # create summary table

### 3. Functions
## Function to generate the cv when the STH specific true underlying FEC in that sample equals x. The values are based on the re-analysis of data kindly provided by Knopp et al., 2010 en Steinmann et al., 2011 
cv <- function(data, parasite) {                                                                                # data: data_cv; parasite: 'Ascaris', 'Trichuris' or 'Hookworm' 
  y <- subset(data, sth == as.character(parasite))                                                              # subset summary table so that CV values correspond with those of a given STH
  list(cv_pop = y$cv_pop,cv_day = y$cv_day, cv_slide=y$cv_slide)                                                # create summary list
  } 
## check a few things
cv(data = data_cv, parasite = "Ascaris")
   
## Function to generate true mean underlying FEC over different days given a true mean individual FEC (mu_ind)
tfec <- function(mu_ind, parasite) {                                                                            # mu_ind: a vector of true mean individual FECs; parasite: 'Ascaris', 'Trichuris' or 'Hookworm' 
  b <- data.frame(x=mu_ind, cv_day = cv(data_cv,parasite=as.character(parasite))$cv_day)                        # determine cv_day for corresponding STH species  
  k_day <- 1/b$cv_day^2                                                                                         # determine k_day for corresponding STH species based on cv_day
  b$day <- rgamma(length(b$x),shape = k_day, rate = k_day/b$x)                                                  # mean individual FEC for a day for all individuals in the vector mu_ind
  list(day = b$day)
}
## Check a few things
tfec(mu_ind=seq(1,1000,100), parasite = "Ascaris")

print(Sys.time())

## Function to generate STH specific egg counts accounting for the slide-to-slide variation and the differences in egg counts across methods (see paper by Cools et al., 2019)
epg2 <- function(x,data,method,parasite) {                                                                                    # x: vector of individual true underling FEC on a day (in EPG), data: data_rel; method: 'KK', 'MF' or 'FP'; parasite: 'Ascaris', 'Trichuris' or 'Hookworm' 
  v <- subset(data, tech==as.character(method) & sth==as.character(parasite))                                                 # subset summary table so that f and m correspond with those of a given method and STH
  b <- data.frame(x=x, cv_slide = cv(data_cv,parasite=as.character(parasite))$cv_slide, m=v$m,k=v$k,f=v$f)                    # determine cv_slide, amount of stool, dispersion parameter and factor expressing FEC relative to KK2 for corresponding STH species  
  k_slide <- 1/b$cv_slide^2                                                                                                   # determine k_slide based on cv_slide  
  #b$slide <- rgamma(length(b$x),shape = k_slide, rate = k_slide/b$x)                                                         # mean true individual FEC in a slide day for all individuals in the vector x
  b$slide <- ifelse(b$m==1/24,x*b$m,b$m*b$f*rgamma(length(b$x),shape = k_slide, rate = k_slide/b$x))                          # adapted based on Luc his comments                                
  #w <- ifelse(b$m==1/24,rpois(length(b$slide),b$slide*b$m),rnbinom(length(b$slide),mu = b$slide*b$m*b$f, size = b$k))        # raw egg counts in a slide made from a sample for all individuals in the vector x
  w <- rpois(length(b$slide),b$slide)                                                                                         # adapted based on Luc his comments 
   return(w)}
## Checking things
epg2(seq(0,1000,1),data_rel,method='MF',parasite='Ascaris')
mean(epg2(seq(0,1000,1),data_rel,method='KK',parasite='Ascaris'))*24
mean(epg2(seq(0,1000,1),data_rel,method='MF',parasite='Ascaris'))*10
mean(epg2(seq(0,1000,1),data_rel,method='FP',parasite='Ascaris'))*34

print(Sys.time())

## Alternative Function to generate ERR data when n subjects are recruited
data_err2 <- function(n, mu_pop, method, parasite, eff, width) {                                                      ## n: number of subjects recruited; mu_pop = mean true underlying FEC (in EPG) within study population;  method: 'KK', 'MF' or 'FP'; parasite: 'Ascaris', 'Trichuris' or 'Hookworm'; eff: mean individual drug response; width: 97.5 percentile of individual drug response - 2.5th percentile individual drug response  
  k_pop <- 1/cv(data_cv,parasite=as.character(parasite))$cv_pop^2                                                     ## determine k_pop
  mu_b <- rgamma(n,shape = k_pop, rate = k_pop/mu_pop)                                                                ## vector of mean true FEC (in EPG) at baseline - assuming a gamma distribution
  mu_b_d1 <- tfec(mu_ind=mu_b, parasite = as.character(parasite))$day                                                 ## vector of mean true FEC (in EPG) on day 1
  mu_b_d2 <- tfec(mu_ind=mu_b, parasite = as.character(parasite))$day                                                 ## vector of mean true FEC (in EPG) on day 2 
  
  egg_b_d1_r1 <- epg2(mu_b_d1,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for first reading/slide on the sample collected on day 1          
  egg_b_d1_r2 <- epg2(mu_b_d1,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for second reading/slide on the sample collected on day 1 
  egg_b_d1 <- egg_b_d1_r1 + egg_b_d1_r2                                                                            ## vector of observed FEC (not in EPG) based on 2 readings/slides on the sample collected on day 1
  
  egg_b_d2_r1 <- epg2(mu_b_d2,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for first reading/slide on the sample collected on day 2          
  egg_b_d2_r2 <- epg2(mu_b_d2,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for second reading/slide on the sample collected on day 2 
  egg_b_d2 <- egg_b_d2_r1 + egg_b_d2_r2                                                                            ## vector of observed FEC (not in EPG) based on 2 readings/slides on the sample collected on day 2
  
  ## B. Generate individual drug response (based on a beta distribution with mean eff and 95% of observations between eff-width/2 and eff+width/2)
  alpha <- betaExpert(best=eff,lower=max(0,eff-width/2),min(1,upper=eff+width/2,p=0.95), method='mean')$alpha         ## value for alpha 
  beta <- betaExpert(best=eff,lower=max(c(0,eff-width/2)),upper=min(c(0.999,eff+width/2)),p=0.95, method='mean')$beta ## value for beta       
  tde <- rbeta(n,shape1=alpha,shape2=beta)                                                                            ## vector of individual drug response
  
  ## C. Generate egg count data a follow-up (in EPG)
  mu_f <- (1-tde)*mu_b                                                                                                ## vector of true mean FEC at follow-up (in EPG)                       
  mu_f_d3 <- tfec(mu_ind=mu_f, parasite = as.character(parasite))$day                                                 ## vector of true mean FEC (in EPG) on day 3
  mu_f_d4 <- tfec(mu_ind=mu_f, parasite = as.character(parasite))$day                                                 ## vector of true mean FEC (in EPG) on day 4
  
  egg_f_d3_r1 <- epg2(mu_f_d3,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for first reading/slide on the sample collected on day 3          
  egg_f_d3_r2 <- epg2(mu_f_d3,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for second reading/slide on the sample collected on day 3 
  egg_f_d3 <- egg_f_d3_r1 + egg_f_d3_r2                                                                            ## vector of observed FEC (not in EPG) based on 2 readings/slides on the sample collected on day 3 
  
  egg_f_d4_r1 <- epg2(mu_f_d4,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for first reading/slide on sample collected on day 4          
  egg_f_d4_r2 <- epg2(mu_f_d4,data=data_rel,method=as.character(method),parasite=as.character(parasite))           ## vector of observed FEC (not in EPG) for second reading/slide on sample collected on day 4 
  egg_f_d4 <- egg_f_d4_r1 + egg_f_d4_r2                                                                            ## vector of observed FEC (not in EPG) based on 2 readings/slides on sample collected on day 4
  
  # D. ERR
  set <- data.frame(mu_b, mu_f, tde, mu_b_d1,egg_b_d1_r1, egg_b_d1_r2, egg_b_d1, mu_b_d2,egg_b_d2_r1, egg_b_d2_r2, egg_b_d2, mu_f_d3,egg_f_d3_r1, egg_f_d3_r2,egg_f_d3, mu_f_d4,egg_f_d4_r1, egg_f_d4_r2,egg_f_d4) ## make data set
  return(set)
}
### Checking a few things
a <- data_err2(n = 1000, mu_pop = 125,method = 'KK', parasite = 'Hookworm', eff = 0.75, width = 0.80)
1-sum(a$egg_f_d3)/sum(a$egg_b_d1)
epg2(a$mu_b_d2,data_rel,method="KK",parasite='Hookworm')

print(Sys.time())


## Function to calculate the ERR, number of days and cost for a certain study design.
calc3 <- function(n, mu_pop, method, parasite, eff, width,data1,data2,simul,thresh) {
  mat <- matrix(NA,nrow=simul,ncol= 7*6)
  for(q in 1:simul) { 
    set <- data_err2(n, mu_pop, method, parasite, eff, width)
    s <- subset(data1,tech==as.character(method))
    o <- subset(data2,tech==as.character(method)) 
    ### Design SS1x1/1x1      
    mat[q,1] <- ifelse(sum(set$egg_b_d1_r1)==0,99, 1 - mean(set$egg_f_d3_r1[set$egg_b_d1_r1>0])/mean(set$egg_b_d1_r1[set$egg_b_d1_r1>0])) ## ERR
    mat[q,2] <- length(set$egg_b_d1_r1[set$egg_b_d1_r1>0]) # Number of positive subjects that tested positive 
    mat[q,3] <- ceiling((n*(s$demo[s$n_slide==1] + s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                         + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)) ## Number of days of work at baseline                  
    mat[q,4] <- ceiling(ifelse(sum(set$egg_b_d1_r1)==0,0,(mat[q,2]*(s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                                                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1[set$egg_b_d1_r1>0])^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60))) ## Number of days of work at follow-up
    mat[q,5] <- (mat[q,3] + mat[q,4])*(o$driver[o$n_slide==1] + o$nurse[o$n_slide==1] + o$lab[o$n_slide==1]) ## cost for salary 
    mat[q,6] <- (n + mat[q,2])*(o$sample[o$n_slide==1] + o$test[o$n_slide==1]) ## cost for consumables                    
    mat[q,7] <- mat[q,5] + mat[q,6] ## total cost
    
    ### Design SS1x2/1x2
    mat[q,8] <- ifelse(sum(set$egg_b_d1)==0 ,99, 1 - mean(set$egg_f_d3[set$egg_b_d1>0])/mean(set$egg_b_d1[set$egg_b_d1>0]))
    mat[q,9] <- length(set$egg_b_d1[set$egg_b_d1>0])
    mat[q,10] <- ceiling((n*(s$demo[s$n_slide==2] + s$prep[s$n_slide==2] + s$count[s$n_slide==2]) 
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2))
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r2)^2)))/(s$n_lab[s$n_slide==2]*s$h_lab[s$n_slide==2]*60*60))                     
    mat[q,11] <- ceiling(ifelse(sum(set$egg_b_d1)==0,0,(mat[q,9]*(s$prep[s$n_slide==2] + s$count[s$n_slide==2])
                                                        + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1[set$egg_b_d1>0])^2)) 
                                                        + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r2[set$egg_b_d1>0])^2)))/(s$n_lab[s$n_slide==2]*s$h_lab[s$n_slide==2]*60*60)))
    mat[q,12] <- (mat[q,10] + mat[q,11])*(o$driver[o$n_slide==2] + o$nurse[o$n_slide==2] + o$lab[o$n_slide==2])
    mat[q,13] <- (n + mat[q,9])*(o$sample[o$n_slide==2] + o$test[o$n_slide==2])                   
    mat[q,14] <- mat[q,12] + mat[q,13]
    
    ### Design SSR1x1/1x1         
    mat[q,15] <- ifelse(sum(set$egg_b_d2_r1[set$egg_b_d1_r1>0])==0,99, 1 - mean(set$egg_f_d3_r1[set$egg_b_d1_r1>0])/mean(set$egg_b_d2_r1[set$egg_b_d1_r1>0]))
    mat[q,16] <- length(set$egg_b_d1_r1[set$egg_b_d1_r1>0])
    mat[q,17] <- ceiling((n*(s$demo[s$n_slide==1] + s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)) 
    + ceiling((mat[q,16]*(s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
               + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d2_r1[set$egg_b_d1_r1>0])^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60))                     
    mat[q,18] <- ceiling(ifelse(sum(set$egg_b_d2_r1[set$egg_b_d1_r1>0])==0,0,(mat[q,16]*(s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1[set$egg_b_d1_r1>0])^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)))
    mat[q,19] <- (mat[q,17] + mat[q,18])*(o$driver[o$n_slide==1] + o$nurse[o$n_slide==1] + o$lab[o$n_slide==1])  
    mat[q,20] <- (n + 2*mat[q,16])*(o$sample[o$n_slide==1] + o$test[o$n_slide==1])                    
    mat[q,21] <- mat[q,19] + mat[q,20]
    
    ### Design SSR1x1/1x2
    mat[q,22] <- ifelse(sum(set$egg_b_d2_r1[set$egg_b_d1_r1>0])==0, 99, 1 - mean(0.5*set$egg_f_d3[set$egg_b_d1_r1>0])/mean(set$egg_b_d2_r1[set$egg_b_d1_r1>0]))
    mat[q,23] <- length(set$egg_b_d1_r1[set$egg_b_d1_r1>0])
    mat[q,24] <- ceiling((n*(s$demo[s$n_slide==1] + s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)) 
    + ceiling((mat[q,23]*(s$prep[s$n_slide==1] + s$count[s$n_slide==1]) + 
                 sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d2_r1[set$egg_b_d1_r1>0])^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60))                     
    mat[q,25] <- ceiling(ifelse(sum(set$egg_b_d2_r1[set$egg_b_d1_r1>0])==0,0,(mat[q,23]*(s$prep[s$n_slide==2] + s$count[s$n_slide==2]) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1[set$egg_b_d1_r1>0])^2)) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r2[set$egg_b_d1_r1>0])^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)))
    mat[q,26] <- (mat[q,24] + mat[q,25])*(o$driver[o$n_slide==1] + o$nurse[o$n_slide==1] + o$lab[o$n_slide==1])  
    mat[q,27] <- (n + mat[q,23])*(o$sample[o$n_slide==1] + o$test[o$n_slide==1]) + mat[q,23]*(o$sample[o$n_slide==2] + o$test[o$n_slide==2])                     
    mat[q,28] <- mat[q,26] + mat[q,27]
    
    ### Design NS1x1/1x1
    mat[q,29] <- ifelse(sum(set$egg_b_d1_r1)==0,99, 1 - mean(set$egg_f_d3_r1)/mean(set$egg_b_d1_r1))
    mat[q,30] <- length(set$egg_b_d1_r1[set$egg_b_d1_r1>0])
    mat[q,31] <- ceiling((n*(s$demo[s$n_slide==1] + s$prep[s$n_slide==1] + s$count[s$n_slide==1])
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)) 
    mat[q,32] <- ceiling(ifelse(sum(set$egg_b_d1_r1)==0,0,(n*(s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)))
    mat[q,33] <- (mat[q,31] + mat[q,32])*(o$driver[o$n_slide==1] + o$nurse[o$n_slide==1] + o$lab[o$n_slide==1])  
    mat[q,34] <- 2*n*(o$sample[o$n_slide==1] + o$test[o$n_slide==1])                  
    mat[q,35] <- mat[q,33] + mat[q,34]
    
    ### Design NS1x1/1x2
    mat[q,36] <- ifelse(sum(set$egg_b_d1_r1)==0,99, 1 - mean(0.5*set$egg_f_d3)/mean(set$egg_b_d1_r1))
    mat[q,37] <- length(set$egg_b_d1_r1[set$egg_b_d1_r1>0])
    mat[q,38] <- ceiling((n*(s$demo[s$n_slide==1] + s$prep[s$n_slide==1] + s$count[s$n_slide==1]) 
                          + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_b_d1_r1)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)) 
    mat[q,39] <- ceiling(ifelse(sum(set$egg_b_d1_r1)==0,0,(n*(s$prep[s$n_slide==2] + s$count[s$n_slide==2]) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r1)^2)) 
                                                           + sum(10^(s$int[1] + s$coef[1]*log10(1+set$egg_f_d3_r2)^2)))/(s$n_lab[s$n_slide==1]*s$h_lab[s$n_slide==1]*60*60)))
    mat[q,40] <- (mat[q,38] + mat[q,39])*(o$driver[o$n_slide==1] + o$nurse[o$n_slide==1] + o$lab[o$n_slide==1])  
    mat[q,41] <- n*(o$sample[o$n_slide==1] + o$test[o$n_slide==1]) + n*(o$sample[o$n_slide==2] + o$test[o$n_slide==2])                     
    mat[q,42] <- mat[q,40] + mat[q,41]
  }
  
  powerSS11 <- sum(ifelse(mat[,1]<thresh,1,0))/simul
  powerSS12 <- sum(ifelse(mat[,8]<thresh,1,0))/simul
  powerSSR11 <- sum(ifelse(mat[,15]<thresh,1,0))/simul
  powerSSR12 <- sum(ifelse(mat[,22]<thresh,1,0))/simul
  powerNS11 <- sum(ifelse(mat[,29]<thresh,1,0))/simul
  powerNS12 <- sum(ifelse(mat[,36]<thresh,1,0))/simul
  
  missSS11 <- sum(ifelse(mat[,1]==99,1,0))/simul
  missSS12 <- sum(ifelse(mat[,8]==99,1,0))/simul
  missSSR11 <- sum(ifelse(mat[,15]==99,1,0))/simul
  missSSR12 <- sum(ifelse(mat[,22]==99,1,0))/simul
  missNS11 <- sum(ifelse(mat[,29]==99,1,0))/simul
  missNS12 <- sum(ifelse(mat[,36]==99,1,0))/simul
  
  costSS11 <- median(mat[,7])
  costSS12 <- median(mat[,14])
  costSSR11 <- median(mat[,21])
  costSSR12 <- median(mat[,28])
  costNS11 <- median(mat[,35])
  costNS12 <- median(mat[,42])
  
  costSS11_Q025 <- quantile(mat[,7], probs=c(0.025))
  costSS12_Q025 <- quantile(mat[,14], probs=c(0.025))
  costSSR11_Q025 <- quantile(mat[,21], probs=c(0.025))
  costSSR12_Q025 <- quantile(mat[,28], probs=c(0.025))
  costNS11_Q025 <- quantile(mat[,35], probs=c(0.025))
  costNS12_Q025 <- quantile(mat[,42], probs=c(0.025))
  
  costSS11_Q975 <- quantile(mat[,7], probs=c(0.975))
  costSS12_Q975 <- quantile(mat[,14], probs=c(0.975))
  costSSR11_Q975 <- quantile(mat[,21], probs=c(0.975))
  costSSR12_Q975 <- quantile(mat[,28], probs=c(0.975))
  costNS11_Q975 <- quantile(mat[,35], probs=c(0.975))
  costNS12_Q975 <- quantile(mat[,42], probs=c(0.975))
  #return(mat)
  list(powerSS11=powerSS11,powerSS12=powerSS12,powerSSR11=powerSSR11,powerSSR12=powerSSR12,powerNS11=powerNS11,powerNS12=powerNS12,
  costSS11=costSS11,costSS12=costSS12,costSSR11=costSSR11,costSSR12=costSSR12,costNS11=costNS11,costNS12=costNS12,
  costSS11_Q025=costSS11_Q025,costSS12_Q025=costSS12_Q025,costSSR11_Q025=costSSR11_Q025,costSSR12_Q025=costSSR12_Q025,costNS11_Q025=costNS11_Q025,costNS12_Q025=costNS12_Q025,
  costSS11_Q975=costSS11_Q975,costSS12_Q975=costSS12_Q975,costSSR11_Q975=costSSR11_Q975,costSSR12_Q975=costSSR12_Q975,costNS11_Q975=costNS11_Q975,costNS12_Q975=costNS12_Q975,
  missSS11=missSS11,missSS12=missSS12,missSSR11=missSSR11,missSSR12=missSSR12,missNS11=missNS11,missNS12=missNS12)
}
## Check a few things
calc3(n = 1000, mu_pop = 125,method = 'MF', parasite = 'Trichuris', eff = 0.35, width = 0.2,data_time,data_cost,simul=1000, thresh=0.40)
calc3(n = 1000, mu_pop = 125,method = 'FP', parasite = 'Hookworm', eff = 0.75, width = 0.2,data_time,data_cost,simul=1000, thresh=0.80)

print(Sys.time())

#### 4. Generation of data
## Trichuris
mu_pop_grid_tt <- c(2.8000, 12.9499, 49.7000, 124.7000)
n_grid_tt <- seq(100,1000,10)
method_grid_tt <-c('KK','MF','FP')

par_grid_tt <- as.data.table(expand.grid(n = n_grid_tt, mu_pop = mu_pop_grid_tt,method = method_grid_tt))
dim(par_grid_tt)
setkey(par_grid_tt, n, mu_pop, method)
sim_data <- par_grid_tt[, calc3(n = n, mu_pop = mu_pop,method = method, parasite = 'Trichuris', eff = 0.35, width = 0.2,data_time,data_cost,simul=1000,thresh=0.40)
,by = .(n, mu_pop, method)]
write.csv(sim_data, file.path(wd, 'data_tri_17DEC2021_POIS.csv'))

## Ascaris
mu_pop_grid_al <- c(9.6, 85.2, 360.0, 2195.5)
n_grid_al <- seq(100,1000,10)
method_grid_al <-c('KK','MF','FP')
eff <- 0.75
width = 0.2
par_grid_al <- as.data.table(expand.grid(n = n_grid_al, mu_pop = mu_pop_grid_al,method = method_grid_al))
dim(par_grid_al)
setkey(par_grid_al, n, mu_pop, method)
sim_data <- par_grid_al[, calc3(n = n, mu_pop = mu_pop,method = method, parasite = 'Ascaris', eff = 0.80, width = 0.2,data_time,data_cost,simul=1000,thresh=0.85)
                        ,by = .(n, mu_pop, method)]
sim_data
write.csv(sim_data, file.path(wd, 'data_al_17DEC2021_POIS.csv'))

## Hookworm
mu_pop_grid_hw <- c(3.7, 23.7, 61.7, 210.3)
n_grid_hw <- seq(100,1000,10)
method_grid_hw <-c('KK','MF','FP')

par_grid_hw <- as.data.table(expand.grid(n = n_grid_hw, mu_pop = mu_pop_grid_hw,method = method_grid_hw))
dim(par_grid_hw)
setkey(par_grid_hw, n, mu_pop, method)
sim_data <- par_grid_hw[, calc3(n = n, mu_pop = mu_pop,method = method, parasite = 'Hookworm', eff = 0.75, width = 0.2,data_time,data_cost,simul=1000,thresh=0.80)
                        ,by = .(n, mu_pop, method)]
sim_data
write.csv(sim_data, file.path(wd, 'data_hw_17DEC2021_POIS.csv'))

save(list=ls(), file="results.rda")

print(Sys.time())


### IGNORE FROM HERE ONWARDS

### 5. Analyzing data
data_al <- read.csv(file.path(wd, 'data_al_17DEC2021_POIS.csv'))
data_tri <- read.csv(file.path(wd, 'data_tri_17DEC2021_POIS.csv'))
data_hw <- read.csv(file.path(wd, 'data_hw_17DEC2021_POIS.csv'))


### Figure 4
power <-c(data_al$powerSS11,data_al$powerSS12,
          data_al$powerSSR11,data_al$powerSSR12,
          data_al$powerNS11,data_al$powerNS12)
method <-rep(data_al$method,6) 
n <-rep(data_al$n,6) 
mu_pop <-rep(data_al$mu_pop,6) 
cost <- c(data_al$costSS11,data_al$costSS12,
          data_al$costSSR11,data_al$costSSR12,
          data_al$costNS11,data_al$costNS12)
design <- c(rep("SS1x1/1x1",length(data_al$powerSS11)),
            rep("SS1x2/1x2",length(data_al$powerSS12)),
            rep("SSR1x1/1x1",length(data_al$powerSS11)),
            rep("SSR1x1/1x2",length(data_al$powerSS11)),
            rep("NS1x1/1x1",length(data_al$powerSS11)),
            rep("NS1x1/1x2",length(data_al$powerSS11)))
data_al2 <- data.frame(power=power,method=method,cost=cost,design=design,n=n,mu_pop=mu_pop)
data_al2$cost_power <- data_al2$cost/(100*data_al2$power)
data_al2$tel <- ifelse(100*data_al2$power>=80,1,0)
data_al22 <- data_al2[order(mu_pop,method,design,n),]
data_al22$tel1 <- shift(data_al22$tel,-1)
data_al22$tel2 <- shift(data_al22$tel1,-1)
data_al22$tel3 <- shift(data_al22$tel2,-1)
data_al22$tel4 <- shift(data_al22$tel3,-1)
data_al22$tel5 <- shift(data_al22$tel4,-1)
data_al22$tel6 <- shift(data_al22$tel5,-1)
data_al22$tel7 <- shift(data_al22$tel6,-1)
data_al22$tel8 <- shift(data_al22$tel7,-1)
data_al22$tel9 <- shift(data_al22$tel8,-1)
data_al22$tel10 <- shift(data_al22$tel9,-1)
length(data_al22$tel5)
data_al22$tel5[6542:6552]
data_al22[c(6542:6552),]

for (i in 1:length(data_al22$tel5)) {
  data_al22$sum3[i] <- sum(data_al22$tel[i],data_al22$tel1[i],data_al22$tel2[i],data_al22$tel3[i])
  data_al22$sum5[i] <- sum(data_al22$tel[i],data_al22$tel1[i],data_al22$tel2[i],data_al22$tel3[i],data_al22$tel4[i],data_al22$tel5[i])
  data_al22$sum10[i] <- sum(data_al22$tel[i],data_al22$tel1[i],data_al22$tel2[i],data_al22$tel3[i],data_al22$tel4[i],data_al22$tel5[i]
                            ,data_al22$tel6[i],data_al22$tel7[i],data_al22$tel8[i],data_al22$tel9[i],data_al22$tel10[i])
  
  } 
hist(data_al22$sum5)
hist(data_al22$sum10)
fig4 <- ggplot(data=data_al2, aes(x=n, y=100*power,group=design)) +
  geom_line(aes(color=as.factor(design)))+
  scale_color_manual(values=cbbPalette,name = "design")+
  theme(legend.position="bottom")+
  xlab('Number of subjects')+
  ylab('Power (%)')+
  ylim(c(0, 100)) +
  geom_hline(yintercept=80,linetype='dashed')+
  facet_grid(mu_pop ~ method)
fig4

### Table XX
data_al3 <- subset(data_al22,data_al22$sum3==4) 
#data_al3 <- subset(data_al2,data_al2$power >= 0.8) 
set <- ddply(data_al3, .(mu_pop,design,method), summarize, n = min(n))
data_al4 <- subset(data_al3,data_al3$method=='KK')
fig5 <- ggplot(data=data_al4, aes(x=100*power, y=cost,color=design)) +
  geom_point()  +
  geom_smooth(method="lm", formula = y ~ x,aes(fill=design))+
 # scale_color_manual(values=cbbPalette,name = "Study design")+
  theme(legend.position="bottom")+
  xlab('Power (%)')+
  ylab('Cost (US$)')+
  xlim(c(80, 90)) +
  ylim(c(3000, 10000)) +
  facet_grid(mu_pop ~ method)
fig5
data_al4$power2 <- 100*data_al4$power
data_SSR_1 <- subset(data_al4, data_al4$design=="SSR1x1/1x2" & data_al4$mu_pop==2195.5)
sum <- lm(cost ~ power2,data=data_SSR_1)
summary(sum)
80*276.32-18595.18
87.5*276.32-18595.18

data_SSR_2 <- subset(data_al4, data_al4$design=="SSR1x1/1x2" & data_al4$mu_pop==360.0)
sum2 <- lm(cost ~ power2,data=data_SSR_2)
summary(sum2)
81.5*304.37-21264.75

data_SSR_3 <- subset(data_al4, data_al4$design=="SSR1x1/1x2" & data_al4$mu_pop==85.2)
sum3 <- lm(cost ~ power2,data=data_SSR_3)
summary(sum3)
80*59.83-972.04

data_NS_1 <- subset(data_al4, data_al4$design=="NS1x1/1x2" & data_al4$mu_pop==2195.5)
sum4 <- lm(cost ~ power2,data=data_NS_1)
summary(sum4)
80*283.10-17672.29
87.5*283.10-17672.29
80*375.25 -26047.63 

### Trichuris
power <-c(data_tri$powerSS11,data_tri$powerSS12,
          data_tri$powerSSR11,data_tri$powerSSR12,
          data_tri$powerNS11,data_tri$powerNS12)
method <-rep(data_tri$method,6) 
n <-rep(data_tri$n,6) 
mu_pop <-rep(data_tri$mu_pop,6) 
cost <- c(data_tri$costSS11,data_tri$costSS12,
          data_tri$costSSR11,data_tri$costSSR12,
          data_tri$costNS11,data_tri$costNS12)
design <- c(rep("SS1x1/1x1",length(data_tri$powerSS11)),
            rep("SS1x2/1x2",length(data_tri$powerSS12)),
            rep("SSR1x1/1x1",length(data_tri$powerSS11)),
            rep("SSR1x1/1x2",length(data_tri$powerSS11)),
            rep("NS1x1/1x1",length(data_tri$powerSS11)),
            rep("NS1x1/1x2",length(data_tri$powerSS11)))
data_tri2 <- data.frame(power=power,method=method,cost=cost,design=design,n=n,mu_pop=mu_pop)
data_tri2$cost_power <- data_tri2$cost/(100*data_tri2$power)
data_tri2 <- data.frame(power=power,method=method,cost=cost,design=design,n=n,mu_pop=mu_pop)
data_tri2$cost_power <- data_tri2$cost/(100*data_tri2$power)
data_tri2$tel <- ifelse(100*data_tri2$power>=80,1,0)
data_tri22 <- data_tri2[order(mu_pop,method,design,n),]
data_tri22$tel1 <- shift(data_tri22$tel,-1)
data_tri22$tel2 <- shift(data_tri22$tel1,-1)
data_tri22$tel3 <- shift(data_tri22$tel2,-1)
data_tri22$tel4 <- shift(data_tri22$tel3,-1)
data_tri22$tel5 <- shift(data_tri22$tel4,-1)
data_tri22$tel6 <- shift(data_tri22$tel5,-1)
data_tri22$tel7 <- shift(data_tri22$tel6,-1)
data_tri22$tel8 <- shift(data_tri22$tel7,-1)
data_tri22$tel9 <- shift(data_tri22$tel8,-1)
data_tri22$tel10 <- shift(data_tri22$tel9,-1)
length(data_tri22$tel5)
data_tri22$tel5[6542:6552]
data_tri22[c(6542:6552),]

for (i in 1:length(data_tri22$tel5)) {
  data_tri22$sum3[i] <- sum(data_tri22$tel[i],data_tri22$tel1[i],data_tri22$tel2[i],data_tri22$tel3[i])
  data_tri22$sum5[i] <- sum(data_tri22$tel[i],data_tri22$tel1[i],data_tri22$tel2[i],data_tri22$tel3[i],data_tri22$tel4[i],data_tri22$tel5[i])
  data_tri22$sum10[i] <- sum(data_tri22$tel[i],data_tri22$tel1[i],data_tri22$tel2[i],data_tri22$tel3[i],data_tri22$tel4[i],data_tri22$tel5[i]
                            ,data_tri22$tel6[i],data_tri22$tel7[i],data_tri22$tel8[i],data_tri22$tel9[i],data_tri22$tel10[i])
  
} 
hist(data_tri22$sum3)
hist(data_tri22$sum5)
hist(data_tri22$sum10)
fig6 <- ggplot(data=data_tri2, aes(x=n, y=100*power,group=design)) +
  geom_line(aes(color=as.factor(design)))+
  scale_color_manual(values=cbbPalette,name = "Study design")+
  theme(legend.position="bottom")+
  xlab('Number of subjects')+
  ylab('Power (%)')+
  ylim(c(0, 100)) +
  geom_hline(yintercept=80,linetype='dashed')+
  facet_grid(mu_pop ~ method)
fig6


data_tri3 <- subset(data_tri22,data_tri22$sum3==4)
set <- ddply(data_tri3, .(mu_pop,design,method), summarize, n = min(n))
data_tri30 <- subset(data_tri3,!data_tri3$mu_pop==12.9499)
data_tri30$check <- ifelse(data_tri30$mu_pop==49.7000 & data_tri30$method=="MF",1,0)
data_tri31 <- subset(data_tri30,data_tri30$check==0)
set <- ddply(data_tri31, .(mu_pop,design,method), summarize, min = min(n),max=(max(n)))
data_tri33 <- merge(data_tri22,set,by=c('mu_pop','method','design'))
data_4 <- subset(data_tri33, data_tri33$power>0.8)
data_4$cost
data_4$power
fig7 <- ggplot(data=data_4, aes(x=100*power, y=cost,color=design)) +
  geom_point()  +
  geom_smooth(method="lm",formula = y ~ x, aes(fill=design))+
  # scale_color_manual(values=cbbPalette,name = "Study design")+
  theme(legend.position="bottom")+
  xlab('Power (%)')+
  ylab('Cost (US$)')+
  xlim(c(80, 90)) +
  ylim(c(3000, 10000)) +
  facet_grid(mu_pop ~ method)
fig7

data_4$power2 <- 100*data_4$power
data_NS_12 <- subset(data_4, data_4$design=="NS1x1/1x2" & data_4$mu_pop==49.7 & data_4$method=='KK')
sum <- lm(cost ~ power2,data=data_NS_12)
summary(sum)
80*299.36-18910.00

data_4$power2 <- 100*data_4$power
data_NS_12 <- subset(data_4, data_4$design=="NS1x1/1x2" & data_4$mu_pop==124.7 & data_4$method=='MF')
sum <- lm(cost ~ power2,data=data_NS_12)
summary(sum)
80*454.34-28801.10

data_NS_11 <- subset(data_4, data_4$design=="NS1x1/1x1" & data_4$mu_pop==49.7 & data_4$method=='KK')
sum <- lm(cost ~ power2,data=data_NS_11)
summary(sum)
80*221.09 -12737.93



data_SSR_2 <- subset(data_tri4, data_tri4$design=="SSR1x1/1x2" & data_tri4$mu_pop==360.0)
sum2 <- lm(cost ~ power2,data=data_SSR_2)
summary(sum2)
81.5*304.37-21264.75

data_SSR_3 <- subset(data_tri4, data_tri4$design=="SSR1x1/1x2" & data_tri4$mu_pop==85.2)
sum3 <- lm(cost ~ power2,data=data_SSR_3)
summary(sum3)
80*59.83-972.04

data_NS_1 <- subset(data_tri4, data_tri4$design=="NS1x1/1x2" & data_tri4$mu_pop==9.6)
sum4 <- lm(cost ~ power2,data=data_NS_1)
summary(sum4)
80*283.10-17672.29
87.5*283.10-17672.29

### Hookworm
power <-c(data_hw$powerSS11,data_hw$powerSS12,
          data_hw$powerSSR11,data_hw$powerSSR12,
          data_hw$powerNS11,data_hw$powerNS12)
method <-rep(data_hw$method,6) 
n <-rep(data_hw$n,6) 
mu_pop <-rep(data_hw$mu_pop,6) 
cost <- c(data_hw$costSS11,data_hw$costSS12,
          data_hw$costSSR11,data_hw$costSSR12,
          data_hw$costNS11,data_hw$costNS12)
design <- c(rep("SS1x1/1x1",length(data_hw$powerSS11)),
            rep("SS1x2/1x2",length(data_hw$powerSS12)),
            rep("SSR1x1/1x1",length(data_hw$powerSS11)),
            rep("SSR1x1/1x2",length(data_hw$powerSS11)),
            rep("NS1x1/1x1",length(data_hw$powerSS11)),
            rep("NS1x1/1x2",length(data_hw$powerSS11)))
data_hw2 <- data.frame(power=power,method=method,cost=cost,design=design,n=n,mu_pop=mu_pop)
data_hw2$cost_power <- data_hw2$cost/(100*data_hw2$power)
data_hw2$tel <- ifelse(100*data_hw2$power>=80,1,0)
data_hw22 <- data_hw2[order(mu_pop,method,design,n),]
data_hw22$tel1 <- shift(data_hw22$tel,-1)
data_hw22$tel2 <- shift(data_hw22$tel1,-1)
data_hw22$tel3 <- shift(data_hw22$tel2,-1)
data_hw22$tel4 <- shift(data_hw22$tel3,-1)
data_hw22$tel5 <- shift(data_hw22$tel4,-1)
data_hw22$tel6 <- shift(data_hw22$tel5,-1)
data_hw22$tel7 <- shift(data_hw22$tel6,-1)
data_hw22$tel8 <- shift(data_hw22$tel7,-1)
data_hw22$tel9 <- shift(data_hw22$tel8,-1)
data_hw22$tel10 <- shift(data_hw22$tel9,-1)
length(data_hw22$tel5)
data_hw22$tel5[6542:6552]
data_hw22[c(6542:6552),]

for (i in 1:length(data_hw22$tel5)) {
  data_hw22$sum3[i] <- sum(data_hw22$tel[i],data_hw22$tel1[i],data_hw22$tel2[i],data_hw22$tel3[i])
  data_hw22$sum5[i] <- sum(data_hw22$tel[i],data_hw22$tel1[i],data_hw22$tel2[i],data_hw22$tel3[i],data_hw22$tel4[i],data_hw22$tel5[i])
  data_hw22$sum10[i] <- sum(data_hw22$tel[i],data_hw22$tel1[i],data_hw22$tel2[i],data_hw22$tel3[i],data_hw22$tel4[i],data_hw22$tel5[i]
                             ,data_hw22$tel6[i],data_hw22$tel7[i],data_hw22$tel8[i],data_hw22$tel9[i],data_hw22$tel10[i])
  
} 
hist(data_hw22$sum3)
hist(data_hw22$sum5)
hist(data_hw22$sum10)

data_hw3 <- subset(data_hw22,data_hw22$sum3==4)
set <- ddply(data_hw3, .(mu_pop,design,method), summarize, n = min(n))
data_hw31 <- subset(data_hw3,!data_hw3$mu_pop==23.7 & data_hw3$method=="KK")
set <- ddply(data_hw31, .(mu_pop,design,method), summarize, min = min(n),max=(max(n)))
data_hw33 <- merge(data_hw22,set,by=c('mu_pop','method','design'))
data_4 <- subset(data_hw33, data_hw33$power>0.8)

fig6 <- ggplot(data=data_hw2, aes(x=n, y=100*power,group=design)) +
  geom_line(aes(color=as.factor(design)))+
  scale_color_manual(values=cbbPalette,name = "Study design")+
  theme(legend.position="bottom")+
  xlab('Number of subjects')+
  ylab('Power (%)')+
  ylim(c(0, 100)) +
  geom_hline(yintercept=80,linetype='dashed')+
  facet_grid(mu_pop ~ method)
fig6


data_4$power2 <- 100*data_4$power
data_NS_12 <- subset(data_4, data_4$design=="NS1x1/1x2" & data_4$mu_pop==61.7 & data_4$method=='KK')
sum <- lm(cost ~ power2,data=data_NS_12)
summary(sum)
80*553.3 -38400.4

data_4$power2 <- 100*data_4$power
data_NS_12 <- subset(data_4, data_4$design=="NS1x1/1x2" & data_4$mu_pop==210.3 & data_4$method=='KK')
sum <- lm(cost ~ power2,data=data_NS_12)
summary(sum)
80*195.7-9429.7 

fig8 <- ggplot(data=data_4, aes(x=100*power, y=cost,color=design)) +
  geom_point()  +
  geom_smooth(method="lm",formula = y ~ x, aes(fill=design))+
  # scale_color_manual(values=cbbPalette,name = "Study design")+
  theme(legend.position="bottom")+
  xlab('Power (%)')+
  ylab('Cost (US$)')+
  xlim(c(80, 90)) +
  ylim(c(3000, 10000)) +
  facet_grid(mu_pop ~ method)
fig8

figure <- ggarrange(fig7, fig8,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure

save(list=ls(), file="results_and_plots.rda")

print(Sys.time())
