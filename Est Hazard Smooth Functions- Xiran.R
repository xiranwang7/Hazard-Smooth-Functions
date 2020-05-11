##############################################################################################################
####            
####                        Rongjie Huang, Xiran Wang, Anja Zgodic
####                              University of South Carolina
####   Smoothing Hazard Function by using Uniform kernel Epachnikov kernel and  Biweight kernel
####                                       3/2020
##############################################################################################################



#Reference
#  Klein, J.P. and Moeschberger, M.L., Survival Analysis: Techniques for Censored 
#  and Truncated Data, Pg165-198, Springer, Berlin, 2nd edition, 2003.



#############################   Read Example Data   #################################
rm(list=ls())
library(survival)
library(KMsurv)
#############################################################################

#setwd( "/Users/effie/Documents/BIOS812")
########### Use 2 Dataset as Examples to Test Out The Functions ############

############################# Read in Bone Data #############################
bone <- read.csv("bone.csv") 
bone$time <- bone$t2
bone$status <- bone$dfree
############################# Read in Kidney Data ###########################
kid <- read.csv("kidney.csv")

kid$time   <- kid$time/365
kid$status <- kid$death

#############################################################################



#############################   Kernel functions   #################################

#### Uniform kernel
k.unif <- function(t,ti,td,b){
  x<-(t-ti)/b
  #Symmetric uniform kernel
  if (t>=b & t<=td-b & x>=-1 & x<=1){
    k<-1/2
    return(k) 
  }
  #Asymmetric uniform kernel
  if (t<b & x>=-1 & x<=t/b){
    q<-t/b
    k<-4*(1+q^3)/(1+q)^4 + 6*(1-q)*x/(1+q)^3
    return(k) 
  }
  #Asymmetric uniform kernel(negative x)
  if (t>td-b & t<td & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    k<-4*(1+q^3)/(1+q)^4 + 6*(1-q)*x/(1+q)^3 
    return(k)
  }  
  #Asymmetric uniform kernel(negative x) x>=td
  if (t>=td & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    k<-4*(1+q^3)/(1+q)^4 + 6*(1-q)*x/(1+q)^3 
    return(k)  
  }
  #x is out of bound  
  else
  {
    k<-0
    return(k)
  }
}


#### Epachnikov kernel
k.ep <- function(t,ti,td,b){
  x<-(t-ti)/b
  #Symmetric Epanechnikov kernel
  if (t>=b & t<=td-b & x>=-1 & x<=1){
    k<-0.75*(1-x^2)
    return(k)
  }
  #Asymmetric Epanechnikov kernel  
  if (t<b & x>=-1 & x<=t/b){ 
    q<-t/b
    alpha<-64*(2-4*q+6*q^2-3*q^3)/((1+q)^4*(19-18*q+3*q^2))
    beta<-240*(1-q)^2/((1+q)^4*(19-18*q+3*q^2))
    k<-0.75*(1-x^2)*(alpha+beta*x)
    return(k)
  }
  #Asymmetric Epanechnikov kernel (negative x)
  if (t>td-b & t<td & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    alpha<-64*(2-4*q+6*q^2-3*q^3)/((1+q)^4*(19-18*q+3*q^2))
    beta<-240*(1-q)^2/((1+q)^4*(19-18*q+3*q^2))
    k<-0.75*(1-x^2)*(alpha+beta*x)
    return(k)
  }
  #Asymmetric Epanechnikov kernel (negative x) when t >=td
  if (t>=td  & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    alpha<-64*(2-4*q+6*q^2-3*q^3)/((1+q)^4*(19-18*q+3*q^2))
    beta<-240*(1-q)^2/((1+q)^4*(19-18*q+3*q^2))
    k<-0.75*(1-x^2)*(alpha+beta*x)
    return(k)
  }
  
  
  else
  {
    k<-0
    return(k)
  }
}

#### Biweight kernel
k.biw <- function(t,ti,td,b){
  x<-(t-ti)/b
  #Symmetric biweight kernel
  if(t>=b & t<=td-b & x>=-1 & x<=1){
    k<-15/16*(1-x^2)^2
    return(round(k,4))
  }
  #Asymmetric biweight kernel
  if (t<b & x>=-1 & x<=1){
    q<-t/b
    alpha<-64*(8-24*q+48*q^2-45*q^3+15*q^4)*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    beta<-1120*(1-q)^3*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    k<-15/16*(1-x^2)^2*(alpha+beta*x)
    return(round(k,4))
  }
  #Asymmetric biweight kernel (negative x)
  if (t>td-b & t<td & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    alpha<-64*(8-24*q+48*q^2-45*q^3+15*q^4)*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    beta<-1120*(1-q)^3*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    k<-15/16*(1-x^2)^2*(alpha+beta*x)
    return(round(k,4))
  }
  
  #Asymmetric biweight kernel (negative x) t>=td
  if ( t>=td & -x>=-1 & -x<=(td-t)/b){
    q<-(td-t)/b
    x<--x  
    alpha<-64*(8-24*q+48*q^2-45*q^3+15*q^4)*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    beta<-1120*(1-q)^3*(1+q)^(-5)*(81-168*q+126*q^2-40*q^3+5*q^4)^(-1)
    k<-15/16*(1-x^2)^2*(alpha+beta*x)
    return(round(k,4))
  }
  #x is out of bound
  else
  {
    k<-0
    return(round(k,4))
  }
}

#############################################################################








##############################################################

smoothing_calculation_ep <- function(data,b,user_t,group=''){ 
  
  if(group==''){
    fit <- survfit(Surv(time, status) ~ 1,data=data, ctype=1)
    t_i <- summary(fit)$time
    cumhaz <- summary(fit)$cumhaz
    delta_cumhaz <- c(cumhaz[1], diff(cumhaz))
    varcumhar <- (summary(fit)$std.chaz)^2
    delta_varcumhar <- c(varcumhar[1], diff(varcumhar))
    
    x <- matrix(NA,nrow=length(t_i),ncol = length(user_t))
    kernels<-matrix(NA,nrow=length(t_i),ncol = length(user_t))
    td <- max(t_i)
    for (t in 1:length (user_t)){
      for (i in 1:length(t_i)){
        x[i,t] <- (user_t[t]-t_i[i])/b
        kernels[i,t]<- k.ep(t=user_t[t],ti=t_i[i],td,b)
      }
    }
    
    out <- cbind(t_i, delta_cumhaz, delta_varcumhar, x, kernels)
    out <- data.frame(out)
    colnames(out) <- c("t_i", "delta_cumhaz", "delta_varcumhar", paste0(user_t, "-t_i/", b),  paste0("K(",user_t, "-t_i/", b,")"))
    return(out)
  }
  
  
  if(group!=''){
    fit <- survfit(Surv(time, status) ~ data[,group], data=data, ctype=1)
    g <- vector("list", length = length(summary(fit)$n))
    t_i <- vector("list", length = length(summary(fit)$n))
    cumhaz <- vector("list", length = length(summary(fit)$n))
    delta_cumhaz <- vector("list", length = length(summary(fit)$n))
    varcumhar <- vector("list", length = length(summary(fit)$n))
    delta_varcumhar <- vector("list", length = length(summary(fit)$n))
    out <- vector("list", length = length(summary(fit)$n))
    
    breaks <- c(which(diff(summary(fit)$time)<0), length(summary(fit)$time))
    t_i[[1]] <- summary(fit)$time[1:breaks[1]]
    g[[1]] <- rep(unique(data[, group])[1], length(t_i[[1]]))
    cumhaz[[1]] <- summary(fit)$cumhaz[1:breaks[1]]
    varcumhar[[1]] <- (summary(fit)$std.chaz[1:breaks[1]])^2
    delta_cumhaz[[1]] <- c(cumhaz[[1]][1], diff(cumhaz[[1]]))
    delta_varcumhar[[1]] <- c(varcumhar[[1]][1], diff(varcumhar[[1]]))
    
    for(l in 2:length(breaks)){
      t_i[[l]] <- summary(fit)$time[(breaks[l-1]+1):breaks[l]]
      g[[l]] <- rep(unique(data[, group])[l], length(t_i[[l]]))
      cumhaz[[l]] <- summary(fit)$cumhaz[(breaks[l-1]+1):breaks[l]]
      varcumhar[[l]] <- (summary(fit)$std.chaz[(breaks[l-1]+1):breaks[l]])^2
      delta_cumhaz[[l]] <- c(cumhaz[[l]][1], diff(cumhaz[[l]])) 
      delta_varcumhar[[l]] <- c(varcumhar[[l]][1], diff(varcumhar[[l]])) 
    }
    
    x = kernels = td <- c()
    for(l in 1:length(summary(fit)$n)){
      x[[l]] <- matrix(NA,nrow=length(t_i[[l]]),ncol = length(user_t))
      kernels[[l]] <-matrix(NA,nrow=length(t_i[[l]]),ncol = length(user_t))
      td[[l]] <- max(t_i[[l]])
    }
    
    for(l in 1:length(summary(fit)$n)){
      for (t in 1:length(user_t)) {
        for (i in 1:length(t_i[[l]])) {
          x[[l]][i,t] <- (user_t[t]-t_i[[l]][i])/b
          kernels[[l]][i,t]<- k.ep(t=user_t[t],ti=t_i[[l]][i],td[[l]],b) #if statement for different kernels, I can add later
        }
      }
    }
    
    for(l in 1:length(summary(fit)$n)){
      out[[l]]  <- cbind(g[[l]], t_i[[l]], delta_cumhaz[[l]], delta_varcumhar[[l]], x[[l]], kernels[[l]])
    }  
    
    out <- do.call("rbind", out)
    out <- data.frame(out)
    colnames(out) <- c("group", "t_i", "delta_cumhaz", "delta_varcumhar", paste0(user_t, "-t_i/", b),  paste0("K(",user_t, "-t_i/", b,")"))
    
    return(out)
  }
}

user_t <- c(150, 50, 600)
timegrid <- seq(1,662,1)
b <- 100

table1 <- smoothing_calculation_ep(data=bone, b=b, user_t=user_t, group='g') # : output group='g'by group 
table1 <- smoothing_calculation_ep(data=bone, b=b, user_t=user_t, group='')  # : null
table1 <- smoothing_calculation_ep(data=bone[bone$g==1,], b=b, user_t=user_t, group='')

#############################################################################











###########################     est_haz function  ####################################

est_haz <-function(data,b,timegrid,group='',kernel='ep'){
  
  if(group==''){
    fit <- survfit(Surv(time, status) ~ 1,data=data, ctype=1)
    t_i <- summary(fit)$time
    cumhaz <- summary(fit)$cumhaz
    delta_cumhaz <- c(cumhaz[1], diff(cumhaz))
    td <- max(t_i)
    
    kernels<-matrix(NA,nrow=length(timegrid),ncol=length(t_i))
    
    if(kernel=="ep"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- k.ep(t=timegrid[i],ti=t_i[j],td,b)*delta_cumhaz[j]
        }
      }
    }
    
    if(kernel=="unif"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- k.unif(t=timegrid[i],ti=t_i[j],td,b)*delta_cumhaz[j]
        }
      }
    }
    
    if(kernel=="biw"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- k.biw(t=timegrid[i],ti=t_i[j],td,b)*delta_cumhaz[j]
        }
      }
    }
    
    haz <- rowSums(kernels, na.rm = T)/b
    haz <- data.frame(timegrid, haz)
    return(haz)
  }
  
  
  if(group!=''){
    fit <- survfit(Surv(time, status) ~ data[,group], data=data, ctype=1)
    g <- vector("list", length = length(summary(fit)$n))
    t_i <- vector("list", length = length(summary(fit)$n))
    cumhaz <- vector("list", length = length(summary(fit)$n))
    delta_cumhaz <- vector("list", length = length(summary(fit)$n))
    
    breaks <- c(which(diff(summary(fit)$time)<0), length(summary(fit)$time))
    t_i[[1]] <- summary(fit)$time[1:breaks[1]]
    g[[1]] <- rep(unique(data[, group])[1], length(t_i[[1]]))
    cumhaz[[1]] <- summary(fit)$cumhaz[1:breaks[1]]
    delta_cumhaz[[1]] <- c(cumhaz[[1]][1], diff(cumhaz[[1]]))
    
    for(l in 2:length(breaks)){
      t_i[[l]] <- summary(fit)$time[(breaks[l-1]+1):breaks[l]]
      g[[l]] <- rep(unique(data[, group])[l], length(t_i[[l]]))
      cumhaz[[l]] <- summary(fit)$cumhaz[(breaks[l-1]+1):breaks[l]]
      delta_cumhaz[[l]] <- c(cumhaz[[l]][1], diff(cumhaz[[l]])) 
    }
    
    x = kernels = td <- c()
    for(l in 1:length(summary(fit)$n)){
      kernels[[l]] <-matrix(NA,nrow=length(timegrid), ncol = length(t_i[[l]]))
      td[[l]] <- max(t_i[[l]])
    }
    
    
    if(kernel=="ep"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            kernels[[l]][i,j] <- k.ep(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b)*delta_cumhaz[[l]][j]
          }
        }
      }
    }
    
    if(kernel=="unif"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            kernels[[l]][i,j] <- k.unif(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b)*delta_cumhaz[[l]][j]
          }
        }
      }
    }
    
    if(kernel=="biw"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            kernels[[l]][i,j] <- k.biw(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b)*delta_cumhaz[[l]][j]
          }
        }
      }
    }
    
    haz_tmp <- c()
    haz_tmp <- lapply(kernels, function(x) (rowSums(x, na.rm = T)/b))
    gr <- rep(unique(data[, group]), each=length((haz_tmp[[1]])))
    haz <- c()
    haz <- data.frame(rep(timegrid,length(unique(data[, group]))), unlist(haz_tmp), gr)
    colnames(haz) <- c("timegrid", "est_haz", "group")
    return(haz)
  }
}



timegrid <- seq(1,730,1)
b <- 100
kernel <- "ep"

haz <- est_haz(data=bone, b=b, timegrid=timegrid, group='', kernel=kernel)
haz.group <- est_haz(data=bone, b=b, timegrid=timegrid, group='g', kernel=kernel)
est_haz(data=bone, b=b, timegrid=50, group='g', kernel=kernel)

#############################################################################







###########################    Est. hazard Plot  Function     ####################################

est_haz_plot <- function(data=data, b=b, timegrid=timegrid, group='', max_haz = max_haz, kernel = kernel){
  
  if(group==''){
    out <- est_haz(data=data,b=b,timegrid=timegrid,group=group,kernel=kernel)
    plot(out, type="l", ylim=c(0,max_haz), lty=1,xlab="Timegrid", ylab="Estimated Hazard Rate", xaxt='n')
  }
  
  if(group!=''){
    out <- est_haz(data=data,b=b,timegrid=timegrid,group=group,kernel=kernel)
    gr <- unique(data[, group])
    plot(out[out$group==gr[1],"est_haz"], type="l",lty=1, ylim=c(0,max_haz), xlab="Timegrid", ylab="Estimated Hazard Rate",  xaxt='n')
    for(i in 2:length(gr)){
      lines(out[out$group==gr[i],"est_haz"], type="l",lty=(5-i))
    }
  }
}

max_haz <- 0.02
b=100
timegrid <- seq(1,690,1)
kernel='ep'

est_haz_plot(data=bone,b=b,timegrid=timegrid,group='g',max_haz=max_haz, kernel=kernel)

#############################################################################







###########################     Cum. Hazard Plot      ####################################

cum_haz_plot <- function(data){ 
  fit <- survfit(Surv(time, status) ~ 1, data=data, ctype=1)
  t_i <- summary(fit)$time
  cumhaz <- summary(fit)$cumhaz
  plot(t_i,cumhaz, ylab="Cumulative Hazard", xlab = "Time", type="S", col="blue")
  return(data.frame(t_i, cumhaz))
}

test <- cum_haz_plot(data=kid)
#########################################################################################






###########################     Est. hazard plot by kernels    ####################################

est_haz_plot_by_kernel <- function(data=data, b=b, timegrid=timegrid, max_haz = max_haz){
  
  ep <- est_haz(data=data, b=b, timegrid=timegrid, group='', kernel='ep')
  unif <- est_haz(data=data, b=b, timegrid=timegrid, group='', kernel='unif')
  biw <- est_haz(data=data, b=b, timegrid=timegrid, group='', kernel='biw')
  
  plot(x=ep$timegrid, y=ep$haz, type="l", ylim=c(0, max_haz), lty=3, ylab="Estimated Hazard Rate", xlab="Timegrid",  xaxt='n')
  lines(x=ep$timegrid, y=unif$haz, lty=1)
  lines(x=ep$timegrid, y=biw$haz, lty=2)
  legend( x="topright", 
          legend=c("Uniform kernel","Epachnikov kernel","Biweight kernel"), lwd=1, lty=c(3,1,2),cex=0.7)
  
  out <- data.frame(timegrid, ep$haz, unif$haz, biw$haz)
  return(out)
  
}

max_haz <- 0.2
b=1
timegrid <- seq(0,8,0.01)
test <- est_haz_plot_by_kernel(data=kid, b=b, timegrid=timegrid, max_haz=max_haz)


##########################################################################






###########################     Est. Hazard Plot Using Bandwidth Method    ####################################

est_haz_plot_by_bandwidth <- function(data=data, b=b, timegrid=timegrid, max_haz = max_haz, kernel='ep'){
  
  ep_b_0.5 <- est_haz(data=data, b=0.5, timegrid=timegrid, group='', kernel=kernel)
  ep_b_1 <- est_haz(data=data, b=1, timegrid=timegrid, group='', kernel=kernel)
  ep_b_1.5 <- est_haz(data=data, b=1.5, timegrid=timegrid, group='', kernel=kernel)
  ep_b_2 <- est_haz(data=data, b=2, timegrid=timegrid, group='', kernel=kernel)
  
  plot(x=timegrid, y=ep_b_0.5$haz, type="l", ylim=c(0, max_haz), lty=1, ylab="Estimated Hazard Rate", xlab="Timegrid", xaxt='n')
  lines(x=timegrid, y=ep_b_1$haz, lty=2)
  lines(x=timegrid, y=ep_b_1.5$haz, lty=5)
  lines(x=timegrid, y=ep_b_2$haz, lty=4)
  
  out <- data.frame(timegrid, ep_b_0.5$haz, ep_b_1$haz, ep_b_1.5$haz, ep_b_2$haz)
  return(out)
} 

max_haz <- 0.17
timegrid <- seq(0,8,0.01)
kid <- read.csv("kidney.csv")
kid$time <- kid$time/365
kid$status <- kid$death
data=kid

test <- est_haz_plot_by_bandwidth(data=kid, b=b, timegrid=timegrid, max_haz=max_haz, kernel='ep')


##########################################################################










###########################     Variance of Est. Hazard function  ####################################

var_haz <- function(data, b, timegrid, group='', kernel='biw'){
  if(group==''){
    fit <- survfit(Surv(time, status) ~ 1,data=data, ctype=1)
    t_i <- summary(fit)$time
    cumhaz <- summary(fit)$cumhaz
    delta_cumhaz <- c(cumhaz[1], diff(cumhaz))
    varcumhar <- (summary(fit)$std.chaz)^2
    delta_varcumhar <- c(varcumhar[1], diff(varcumhar))
    td <- max(t_i)
    
    kernels<-matrix(NA,nrow=length(timegrid),ncol=length(t_i))
    
    if(kernel=="ep"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- ((k.ep(t=timegrid[i],ti=t_i[j],td,b))^2)*(delta_varcumhar[j])
          #print((k.ep(t=timegrid[i],ti=t_i[j],td,b)))
        }
      }
    }
    
    if(kernel=="unif"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- ((k.unif(t=timegrid[i],ti=t_i[j],td,b))^2)*(delta_varcumhar[j])
        }
      }
    }
    
    if(kernel=="biw"){
      for (i in 1:length(timegrid)){
        for (j in 1:length(t_i)){
          kernels[i,j] <- ((k.biw(t=timegrid[i],ti=t_i[j],td,b))^2)*(delta_varcumhar[j])
        }
      }
    }
    
    varhaz <- (rowSums(kernels, na.rm = T))/(b^2)
    varhaz <- data.frame(timegrid, varhaz)
    return(varhaz)
  }
  
  
  if(group!=''){
    
    fit <- survfit(Surv(time, status) ~ data[,group], data=data, ctype=1)
    g <- vector("list", length = length(summary(fit)$n))
    t_i <- vector("list", length = length(summary(fit)$n))
    cumhaz <- vector("list", length = length(summary(fit)$n))
    delta_cumhaz <- vector("list", length = length(summary(fit)$n))
    varcumhar <- vector("list", length = length(summary(fit)$n))
    delta_varcumhar <- vector("list", length = length(summary(fit)$n))
    
    breaks <- c(which(diff(summary(fit)$time)<0), length(summary(fit)$time))
    t_i[[1]] <- summary(fit)$time[1:breaks[1]]
    g[[1]] <- rep(unique(data[, group])[1], length(t_i[[1]]))
    cumhaz[[1]] <- summary(fit)$cumhaz[1:breaks[1]]
    delta_cumhaz[[1]] <- c(cumhaz[[1]][1], diff(cumhaz[[1]]))
    varcumhar[[1]] <- (summary(fit)$std.chaz[1:breaks[1]])^2
    delta_varcumhar[[1]] <- c(varcumhar[[1]][1], diff(varcumhar[[1]]))
    
    for(l in 2:length(breaks)){
      t_i[[l]] <- summary(fit)$time[(breaks[l-1]+1):breaks[l]]
      g[[l]] <- rep(unique(data[, group])[l], length(t_i[[l]]))
      cumhaz[[l]] <- summary(fit)$cumhaz[(breaks[l-1]+1):breaks[l]]
      delta_cumhaz[[l]] <- c(cumhaz[[l]][1], diff(cumhaz[[l]]))
      varcumhar[[l]] <- (summary(fit)$std.chaz[(breaks[l-1]+1):breaks[l]])^2
      delta_varcumhar[[l]] <- c(varcumhar[[l]][1], diff(varcumhar[[l]])) 
    }
    
    x = kernels = td <- c()
    for(l in 1:length(summary(fit)$n)){
      kernels[[l]] <-matrix(NA,nrow=length(timegrid), ncol = length(t_i[[l]]))
      td[[l]] <- max(t_i[[l]])
    }
    
    
    if(kernel=="ep"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            #if(l==1 & i %in% c(50, 150, 600)) {print(c("l=", l, "i=", i, "j=", j));print(((k.ep(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b))^2)); print((delta_varcumhar[[l]][j]))}
            kernels[[l]][i,j] <- ((k.ep(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b))^2)*(delta_varcumhar[[l]][j])
            #if(l==1 & i %in% c(50, 150, 600)) { print(kernels[[l]][i,j])}
          }
        }
      }
    }
    
    if(kernel=="unif"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            kernels[[l]][i,j] <- ((k.unif(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b))^2)*(delta_varcumhar[[l]][j])
          }
        }
      }
    }
    
    if(kernel=="biw"){
      for(l in 1:length(summary(fit)$n)){
        for (i in 1:length(timegrid)){
          for (j in 1:length(t_i[[l]])){
            kernels[[l]][i,j] <- ((k.biw(t=timegrid[i],ti=t_i[[l]][j],td[[l]],b))^2)*(delta_varcumhar[[l]][j])
          }
        }
      }
    }
    
    varhaz_tmp <- c()
    varhaz_tmp <- lapply(kernels, function(x) (rowSums(x, na.rm = T))/(b^2))
    gr <- rep(unique(data[, group]), each=length((varhaz_tmp[[1]])))
    varhaz <- c()
    varhaz <- data.frame(rep(timegrid,length(unique(data[, group]))), unlist(varhaz_tmp), gr)
    colnames(varhaz) <- c("timegrid", "var_haz", "group")
    return(varhaz)
  }
}

##########################################################################










################################# g(b)  Optimal Function ##########################################
####g(b) function
opti_gb <- function(b,u,data,kernel){
  
  fit_optim <- survfit(Surv(time, status) ~ 1, data=data, ctype=1)
  t_i_optim <- summary(fit_optim)$time
  cumhaz_optim <- summary(fit_optim)$cumhaz
  diffcumhz_optim <-c(cumhaz_optim[1], diff(cumhaz_optim) )
  
  ### gb1
  gb1<- vector()
  for (i in 1: ( length(u)-1) ) {
    
    est_hi <- est_haz(data=data, b=b, timegrid=u[i],group='',kernel= kernel)
    est_hi1 <- est_haz(data=data, b=b, timegrid=u[i+1],group='',kernel=kernel)
    gb1[i] <- ( u[i+1]-u[i] ) /2 *  (  ( est_hi$haz )^2 +  ( est_hi1$haz )^2   ) 
  }
  
  ### gb2
  gb2 <- matrix(NA,nrow= length(t_i_optim),ncol = length(t_i_optim))
  
  if(kernel=="unif"){ 
    
    for (i in 1:length(t_i_optim) ){
      for (j in 1:length(t_i_optim) ){
        gb2[i,j]<- - (2/b) * k.unif(t=t_i_optim[i],ti=t_i_optim[j],td=max(t_i_optim),b=b) * diffcumhz_optim[i]* diffcumhz_optim[j]
      }
    }
  }      
  
  if(kernel=="ep") { 
    
    for (i in 1:length(t_i_optim) ){
      for (j in 1:length(t_i_optim) ){
        gb2[i,j]<- - (2/b) * k.ep(t=t_i_optim[i],ti=t_i_optim[j],td=max(t_i_optim),b=b) * diffcumhz_optim[i]* diffcumhz_optim[j]
      } 
    }
  }
  
  
  if(kernel=="biw") { 
    
    for (i in 1:length(t_i_optim) ){
      for (j in 1:length(t_i_optim) ){
        gb2[i,j]<- - (2/b) * k.biw(t=t_i_optim[i],ti=t_i_optim[j],td=max(t_i_optim),b=b) * diffcumhz_optim[i]* diffcumhz_optim[j]
      }
    }
  }
  
  
  diag(gb2) <-0 # dont add i= j case so make diag(gb2 matrix )
  gb <- sum(gb1) + sum(gb2)
  
  return(gb)
  
}

################################# Optimal b plot ##########################################



#####Fig 6.5  function
optimal_b_plot <- function (u, bgrid,data ){
  
  pts_ep<-vector()
  pts_unif<-vector()
  pts_biw<-vector()
  
  for (p in 1:length(bgrid)) {
    pts_ep[p]<-opti_gb(b=bgrid[p],u=u,data=data, kernel='ep')
    pts_unif[p]<-opti_gb(b=bgrid[p],u=u,data=data, kernel='unif')
    pts_biw[p]<-opti_gb(b=bgrid[p],u=u,data=data, kernel='biw') 
  }
  
  ### fig6.5
  plot(bgrid,pts_ep,type = 'l',lty=2,xlab = "b", ylab = "g(b)", ylim=c(-0.02, 0.02))
  lines(bgrid,pts_unif,lty=1 )
  lines(bgrid, pts_biw,lty=5)
  legend("topleft", legend=c("Epanechnikov kernel", "Uniform kernel","Biweight kernel"),lty=c(2,1,5), cex=0.8)
  ### optim bgrid
  ep_optim <- bgrid[(pts_ep==min(pts_ep))]
  unif_optim <-bgrid[(pts_unif==min(pts_unif))]
  biw_optim <- bgrid[(pts_biw==min(pts_biw))]
  
  out <- data.frame(ep_optim,unif_optim,biw_optim)
  
  return(out)
  
}


kid <- read.csv("kidney.csv")
kid$status <- kid$death
kid <- kid[order(kid$time), ]
kid$time <- kid$time/365
u <-  seq(0,6,1)
bgrid <- seq(0.01,1,0.01) 
test <- optimal_b_plot(u=u, bgrid=bgrid, data=kid)
test
############################################








###########################     Est. Hazard with optimak b plot   ################################


est_haz_optim_plot <- function(data=data, timegrid=timegrid, max_haz = max_haz, kernel='biw', optimal_b){
  opti_haz <- est_haz(data = data, timegrid = timegrid, b = optimal_b, group = '',kernel = kernel)
  var.haz <- var_haz(data = data, timegrid = timegrid, b = optimal_b, group = '', kernel = kernel)
  
  haz_lo <- (opti_haz$haz)*exp((-1.96*sqrt(var.haz$varhaz))/(opti_haz$haz))
  haz_hi <- (opti_haz$haz)*exp((1.96*sqrt(var.haz$varhaz))/(opti_haz$haz))
  
  plot(timegrid,opti_haz$haz,type = 'l' ,ylim = c(0,max_haz),ylab=c('Estimated Hazard Rate'),
       xlab=c('Years Post Transplant'), xaxt='n')
  points(timegrid,haz_lo,type = 'l', lty=2)
  points(timegrid,haz_hi,type = 'l', lty=2)
  
  out <- data.frame(timegrid, opti_haz$haz, var.haz$varhaz, haz_lo, haz_hi)
  names(out) <- c("timegrid", "est.haz", "var.haz", "haz.low", "haz.high")
  return(out)
} 


max_haz = 0.35
timegrid=seq(0,8,0.001)
kernel='biw'
#optimal_b <- test$ep_optim
optimal_b=0.23 
kid <- read.csv("/cloud/project/kidney.csv")
kid$time <- kid$time/365
kid$status <- kid$death

test <- est_haz_optim_plot(data=data, timegrid=timegrid, max_haz = max_haz, kernel=kernel, optimal_b)

##########################################################################






#############################################
####g(b) function by group 
opti_gb_group <- function(b,u,data,kernel,group='g'){
  
  fit_optim <- survfit(Surv(time, status) ~ data[,group], data=data, ctype=1)
  g_optim <- vector("list", length = length(summary(fit_optim)$n))
  t_i_optim <- vector("list", length = length(summary(fit_optim)$n))
  cumhaz_optim <- vector("list", length = length(summary(fit_optim)$n))
  diffcumhz_optim <- vector("list", length = length(summary(fit_optim)$n))
  
  breaks_optim <- c(which(diff(summary(fit_optim)$time)<0), length(summary(fit_optim)$time))
  t_i_optim[[1]] <- summary(fit_optim)$time[1:breaks_optim[1]]
  g_optim[[1]] <- rep(unique(data[, group])[1], length(t_i_optim[[1]]))
  cumhaz_optim[[1]] <- summary(fit_optim)$cumhaz[1:breaks_optim[1]]
  diffcumhz_optim[[1]] <- c(cumhaz_optim[[1]][1], diff(cumhaz_optim[[1]]))
  
  for(l in 2:length(breaks_optim)){
    t_i_optim[[l]] <- summary(fit_optim)$time[(breaks_optim[l-1]+1):breaks_optim[l]]
    g_optim[[l]] <- rep(unique(data[, group])[l], length(t_i_optim[[l]]))
    cumhaz_optim[[l]] <- summary(fit_optim)$cumhaz[(breaks_optim[l-1]+1):breaks_optim[l]]
    diffcumhz_optim[[l]] <- c(cumhaz_optim[[l]][1], diff(cumhaz_optim[[l]])) 
  }
  
  
  ### gb1
  gb11=gb12=gb13<- vector()
  
  for (i in 1: ( length(u)-1) ) {
    
    est_hi <- est_haz(data=data, b=b, timegrid=u[i],group='g',kernel= kernel)
    est_hi1 <- est_haz(data=data, b=b, timegrid=u[i+1], group='g',kernel=kernel)
    
    group1<- est_hi[est_hi$group ==1,][2]
    group11<-est_hi1[est_hi1$group ==1,][2]
    
    group2<- est_hi[est_hi$group ==2,][2]
    group21<-est_hi1[est_hi1$group ==2,][2]
    
    group3 <- est_hi[est_hi$group ==3,][2]
    group31<-est_hi1[est_hi1$group ==3,][2]
    
    
    gb11[i] <-  ( u[i+1]-u[i] ) /2 *  (  ( group1 )^2 +  ( group11 )^2   )  
    gb12[i] <- ( u[i+1]-u[i] ) /2 *  (  ( group2 )^2 +  ( group21 )^2   ) 
    gb13[i] <- ( u[i+1]-u[i] ) /2 *  (  ( group3 )^2 +  ( group31 )^2   ) 
    
  }
  
  ### gb2
  gb21 <- matrix(NA,nrow= length(t_i_optim[[1]]),ncol = length(t_i_optim[[1]]))
  gb22 <- matrix(NA,nrow= length(t_i_optim[[2]]),ncol = length(t_i_optim[[2]]))
  gb23 <- matrix(NA,nrow= length(t_i_optim[[3]]),ncol = length(t_i_optim[[3]]))
  
  if(!(dim(gb21)[1] %in% 0 & dim(gb21)[1] %in% 0)){ 
    for (i in 1:length(t_i_optim[[1]]) ){
      for (j in 1:length(t_i_optim[[1]]) ){
        gb21[i,j]<- - (2/b) * k.ep(t=t_i_optim[[1]][i],ti=t_i_optim[[1]][j],td=max(t_i_optim[[1]]),b=b) * diffcumhz_optim[[1]][i]* diffcumhz_optim[[1]][j]
      } 
    }
  }
  
  if(!(dim(gb22)[1] %in% 0 & dim(gb22)[1] %in% 0)){ 
    for (i in 1:length(t_i_optim[[2]]) ){
      for (j in 1:length(t_i_optim[[2]]) ){
        gb22[i,j]<- - (2/b) * k.ep(t=t_i_optim[[2]][i],ti=t_i_optim[[2]][j],td=max(t_i_optim[[2]]),b=b) * diffcumhz_optim[[2]][i]* diffcumhz_optim[[2]][j]
      } 
    }
  }
  
  if(!(dim(gb23)[1] %in% 0 & dim(gb23)[1] %in% 0)){ 
    for (i in 1:length(t_i_optim[[3]]) ){
      for (j in 1:length(t_i_optim[[3]]) ){
        gb23[i,j]<- - (2/b) * k.ep(t=t_i_optim[[3]][i],ti=t_i_optim[[3]][j],td=max(t_i_optim[[3]]),b=b) * diffcumhz_optim[[3]][i]* diffcumhz_optim[[3]][j]
      } 
    }
  }
  
  if(dim(gb21)[1]==0 | dim(gb21)[2]==0){ gb21 <- matrix(0, ncol=1, nrow=1) }
  if(dim(gb22)[1]==0 | dim(gb22)[2]==0){ gb22 <- matrix(0, ncol=1, nrow=1) }
  if(dim(gb23)[1]==0 | dim(gb23)[2]==0){ gb23 <- matrix(0, ncol=1, nrow=1) }
  
  
  diag(gb21) <-0 # dont add i= j case so make diag(gb2 matrix )
  diag(gb22) <-0
  diag(gb23) <-0
  
  
  gb1 <- sum(unlist( gb11) ) + sum(gb21)
  gb2 <- sum(unlist( gb12) ) + sum(gb22)
  gb3 <- sum(unlist( gb13) ) + sum(gb23)
  
  out <- data.frame(gb1, gb2,gb3)
  names(out) <- c("gb_grp1", "gb_grp2", "gb_grp3")
  return(out)
  
}
##############################################





##############################################
####fig  
'optimal_b_plot_group <- function( u,bgrid,data,kernel){
pts_ep1<-pts_ep2<-pts_ep3<-vector()
for (p in 1:length(bgrid)) {
optbs<-opti_gb_group(b=bgrid[p],u=u,data=data, kernel=kernel)
pts_ep1[p]<-optbs[1,1]
pts_ep2[p]<-optbs[1,2]
pts_ep3[p]<-optbs[1,3]
}

plot(bgrid,pts_ep1,type = "l",lty=2,ylab="g(b)", xlab="b", xaxt="n" )
lines(bgrid,pts_ep2 )
lines(bgrid, pts_ep3,lty=5)

### optim bs
optim1 <- bgrid[(pts_ep1==min(pts_ep1))]
optim2 <-bgrid[(pts_ep2==min(pts_ep2))]
optim3 <- bgrid[(pts_ep3==min(pts_ep3))]

out <- data.frame(optim1, optim2, optim3)
names(out) <- c("b_grp1", "b_grp2", "b_grp3")
return( out )

}


u <-  seq(1, 662, 10)
bgrid <- seq(1,300,10) 

test <- optimal_b_plot_group(u=u, bgrid=bgrid, data=bone, kernel="ep")
'
##############################################


################################## est_haz_optim_group_plot  ####################################


est_haz_optim_group_plot <- function(data, grp='g', timegrid, max_haz = max_haz, kernel=kernel, optimal_b){
  
  haz_1 <- est_haz(data=data[data[,grp]==1,], b=optimal_b[1], timegrid=timegrid, group='', kernel=kernel)$haz
  haz_2 <- est_haz(data=data[data[,grp]==2,], b=optimal_b[2], timegrid=timegrid, group='', kernel=kernel)$haz
  haz_3 <- est_haz(data=data[data[,grp]==3,], b=optimal_b[3], timegrid=timegrid, group='', kernel=kernel)$haz
  
  plot(timegrid, haz_1,type = 'l', ylim=c(0,max_haz ), xlim=c(0,max(timegrid)), ylab = 'Estimate Hazard Rate', xlab = "Timegrid", xaxt='n' )
  lines(timegrid,  haz_2,type = 'l',lty=2)
  lines(timegrid,  haz_3,type = 'l',lty=3)
  
  out <- data.frame(timegrid, haz_1, haz_2, haz_3)
  names(out) <- c("timegrid", "est_haz_grp1", "est_haz_grp2", "est_haz_grp3")
  return(out)
}

max_haz=0.01
optimal_b <- c(161, 50, 112)
#optimal_b <- optimal_b_plot_group(u=u, bgrid=bgrid, data=bone, kernel='ep')
timegrid<-seq(1,730,1)
grp='g'
data=bone
kernel='ep'

test <- est_haz_optim_group_plot(data=data, grp=grp, timegrid=timegrid, max_haz = max_haz, kernel=kernel, optimal_b)

###############################################










