library(MASS)

anglealpha <- function(data,alpha,start){
  data <- data[((data[,1]>start[1])&(data[,2]>start[2])),]
  
  outnum <- floor(dim(data)[1]*(1-alpha))
  
  ab <- function(point,start){
    a <- (point[2]-start[2])/(point[1]-start[1])
    b <- (start[2]*point[1]-point[2]*start[1])/(point[1]-start[1])
    return(c(a,b))
  }
  
  abangle <- function(a,b,point,start){
    pointtrans <- point-start
    
    if(a>=0){
      angle <- atan(a)/pi*180
    }else if((a<0)&(pointtrans[1]<0)){
      angle <- atan(a)/pi*180+180
    }else{
      angle <- atan(a)/pi*180
    }
    return(angle)
  }
  
  Tan <- c();Angle <- c();Coef <- c()
  for(i in 1:dim(data)[1]){
    point <- data[i,]
    coef <- ab(point,start)
    Tan <- c(Tan,coef[1])
    Coef <- rbind(Coef,coef)
    
    angle <- abangle(coef[1],coef[2],point,start)
    Angle <- c(Angle,angle)
  }
  
  comple <- function(outnum,Angle,orderup,orderdown){
    updown <- c();anglevs <- c()
    for(i in 1:(outnum+1)){
      up <- orderup[i]
      down <- orderdown[outnum+2-i]
      updown <- rbind(updown,c(up,down))
      anglevs <- c(anglevs,abs(Angle[up]-Angle[down]))
    }
    
    return(list(updown,anglevs))
  }
  
  orderup <- order(-Angle)
  orderdown <- order(Angle)
  if(outnum==0){
    whereup <- orderup[1]
    wheredown <- orderdown[1]
    anglemin <- abs(Angle[whereup]-Angle[wheredown])
    tanup <- Tan[whereup];tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    pointup <- data[whereup,];coefup <- Coef[whereup,]
    pointdown <- data[wheredown,];coefdown <- Coef[wheredown,]
    
    angleminlist <- rbind(c(anglemin,tan2),pointup,pointdown,coefup,coefdown)
    angleequallist <- angleminlist;anglemeanlist <- angleminlist;angleneighlist <- angleminlist;anglemedianlist <- angleminlist
  }else{
    updownanglevs <- comple(outnum,Angle,orderup,orderdown)
    updown <- updownanglevs[[1]]
    anglevs <- updownanglevs[[2]]
    
    where <- order(anglevs)[1]
    whereup <- updown[where,1]
    wheredown <- updown[where,2]
    anglemin <- anglevs[where]
    tanup <- Tan[whereup];tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    pointup <- data[whereup,];coefup <- Coef[whereup,]
    pointdown <- data[wheredown,];coefdown <- Coef[wheredown,]
    angleminlist <- rbind(c(anglemin,tan2),pointup,pointdown,coefup,coefdown)
    
    if(outnum==1){
      angleequallist <- angleminlist
    }else if(outnum%%2==0){
      where <- ceiling(outnum/2)
      whereup <- updown[where,1]
      wheredown <- updown[where,2]
      angleequal <- anglevs[where]
      tanup <- Tan[whereup];tandown <- Tan[wheredown]
      tan2 <- (tanup-tandown)/(1+tanup*tandown)
      pointup <- data[whereup,];coefup <- Coef[whereup,]
      pointdown <- data[wheredown,];coefdown <- Coef[wheredown,]
      angleequallist <- rbind(c(angleequal,tan2),pointup,pointdown,coefup,coefdown)
    }else{
      where <- floor(outnum/2)
      whereup <- updown[where,1];whereupup <- updown[where+1,1]
      wheredown <- updown[where,2];wheredowndown <- updown[where+1,2]
      pointup <- apply(data[c(whereup,whereupup),],2,mean);coefup <- ab(pointup,start)
      pointdown <- apply(data[c(wheredown,wheredowndown),],2,mean);coefdown <- ab(pointdown,start)
      angleequal <- abs(mean(Angle[c(whereup,whereupup)])-mean(Angle[c(wheredown,wheredowndown)]));tan2 <- tan(angleequal)
      angleequallist <- rbind(c(angleequal,tan2),pointup,pointdown,coefup,coefdown)
    }
    
    anglemean <- mean(anglevs);tan2 <- tan(anglemean)
    pointup <- apply(data[updown[,1],],2,mean);coefup <- ab(pointup,start)
    pointdown <- apply(data[updown[,2],],2,mean);coefdown <- ab(pointdown,start)
    anglemeanlist <- rbind(c(anglemean,tan2),pointup,pointdown,coefup,coefdown)
    
    anglemedian <- median(anglevs);tan2 <- tan(anglemedian)
    pointup <- apply(data[updown[,1],],2,median);coefup <- ab(pointup,start)
    pointdown <- apply(data[updown[,2],],2,median);coefdown <- ab(pointdown,start)
    anglemedianlist <- rbind(c(anglemedian,tan2),pointup,pointdown,coefup,coefdown)
    
    where <- order(anglevs)[1]
    whereup <- updown[where,1];upnum <- which(orderup==whereup);whereupup <- orderup[ifelse(upnum==1,upnum,upnum-1)]
    wheredown <- updown[where,2];downnum <- which(orderdown==wheredown);wheredowndown <- orderdown[ifelse(downnum==1,downnum,downnum-1)]
    pointup <- apply(data[c(whereup,whereupup),],2,mean);coefup <- ab(pointup,start)
    pointdown <- apply(data[c(wheredown,wheredowndown),],2,mean);coefdown <- ab(pointdown,start)
    angleneigh <- abs(mean(Angle[c(whereup,whereupup)])-mean(Angle[c(wheredown,wheredowndown)]));tan2 <- tan(angleneigh)
    angleneighlist <- rbind(c(angleneigh,tan2),pointup,pointdown,coefup,coefdown)
  }
  
  out <- list(anglemin=angleminlist,anglemedian=anglemedianlist,angleequal=angleequallist,anglemean=anglemeanlist,angleneigh=angleneighlist)
  
  return(out)
}

sigmamatrix <- function(sigma,rho,i){
  Sigma <- matrix(c(sigma[1],rho*(1-(i-1)/5)*sqrt(sigma[1]*sigma[2]),rho*(1-(i-1)/5)*sqrt(sigma[1]*sigma[2]),sigma[2]),2,2)
  return(Sigma)
}

split <- function(mu,sigma,rho,i,N,Tij){
  data <- mvrnorm(Tij*N^i,mu,sigmamatrix(sigma,rho,i))
  return(data)
}

treedataresidual <- function(mu,sigma,rho,T,N,Tij){
  treeresidual1 <- vector('list',T)
  treeresidual2 <- vector('list',T)
  for(i in 1:T){
    residual1 <- split(mu,sigma,rho,i,N,Tij)
    treeresidual1[[i]] <- residual1[,1]
    treeresidual2[[i]] <- residual1[,2]
  }
  return(treeresidualdata=list(treeresidual1=treeresidual1,treeresidual2=treeresidual2))
}

musigma <- function(data){
  mu <- apply(data,2,mean)
  sigma <- matrix(rep(0,(length(mu))^2),length(mu),length(mu))
  for(i in 1:(dim(data)[1])){
    x0 <- matrix(data[i,]-mu,length(mu),1)
    sigma <- sigma+(x0)%*%t(x0)
  }
  sigma <- sigma/(dim(data)[1])
  
  return(c(mu,diag(sigma),sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])))
}

datauni <- function(data){
  T<- length(data[[1]])
  Sigma <- sqrt(musigma(cbind(data[[1]][[T]],data[[2]][[T]]))[3:4])
  for(i in 1:T){
    data[[1]][[i]] <- data[[1]][[i]]/Sigma[1]*sqrt(0.0001)
    data[[2]][[i]] <- data[[2]][[i]]/Sigma[2]*sqrt(0.0001)
    data[[1]][[i]] <- data[[1]][[i]]-mean(data[[1]][[i]])+sqrt(-2*log(0.05)*0.0001 + sum(1/(1:i))*0.1)
    data[[2]][[i]] <- data[[2]][[i]]-mean(data[[2]][[i]])+sqrt(-2*log(0.05)*0.0001 + sum(1/(1:i))*0.1)
  }
  data <- cbind(unlist(data[[1]]),unlist(data[[2]]))
  return(data)
}

deltaanglefunc <- function(rho,mu,sigma,T,N,Tij,alpha,start){
  data <- treedataresidual(mu,sigma,rho,T,N,Tij)
  data <- list(without=cbind(unlist(data[[1]]),unlist(data[[2]])),with=datauni(data))
  deltaangle <- c()
  for(i in 1:2){
    out <- anglealpha(data[[i]],alpha,start)
    deltaangle <- rbind(deltaangle,unlist(lapply(out,function(out0){return(out0[1,1])})))
  }
  row.names(deltaangle) <- c('without_normalization','with_normalization')
  colnames(deltaangle) <- c('min','median','equal','mean','neigh')
  return(deltaangle)
}