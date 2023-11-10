#' Identifier of the first short-term outcome
#' @param x the first short-term outcome
#' @returns returns the input if this function is used directly. See `details` for more explanation
#' @export
#' @details
#' This is a helperfunction used by \code{multipredict()} to identify the first
#' short-term outcome. If used directly it will only return the input. Therefore it should
#' not be called directly and only be used in the `formula` argument in \code{multipredict()}
#' to identify the variable that is the first short-term outcome. See `Examples` section of
#' \code{multipredict()} for more details.
s1 <- function(x) x


#' Identifier of the second short-term outcome
#' @param x the second short-term outcome
#' @returns returns the input if this function is used directly. See `details` for more explanation
#' @export
#' @details
#' This is a helperfunction used by \code{multipredict()} to identify the second
#' short-term outcome. It should not be called directly. If used directly it will only return the input. Therefore it should
#' not be called directly and only be used in the `formula` argument in \code{multipredict()}
#' to identify the variable that is the second short-term outcome. See `example` section of
#' \code{multipredict()} for more details.
s2 <- function(x) x



#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
g.logit <- function(xx){exp(xx)/(1+exp(xx))}


#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
Kern.FUN.log  <- function(S1i,ti,hi) {
  out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
  norm.k = dnorm(out)/hi
  norm.k
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
# calculate kernel matrix using bivariate normal
Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
{ out = (log(rbind(S1i,S2i))- log(ti))
norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
norm.k = det(Hi)^(-1/2) * norm.k
norm.k
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr2 <- function(da2,t0,L,h,folds,reps,s.seq)
{
Y     = da2$Y
XS1   = da2$XS1
XS2   = da2$XS2
delta = da2$delta
X     = data.matrix(subset(da2, select = -c(Y, XS1, XS2, delta, what, group)))
Wi    = da2$what

# logit link function
g.logit <- function(xx){exp(xx)/(1+exp(xx))}

# calculate kernel matrix (a log is added. Need to discuss)
Kern.FUN.log  <- function(S1i,ti,hi) {
  out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
  norm.k = dnorm(out)/hi
  norm.k
}
loc.fun.ex <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
n = dim(da2)[1]; nv = floor(n/folds)
BW.vec = vector(length=folds)
replic = matrix(nrow = reps, ncol =folds)
for(lay in 1:reps) {
  for(k in 1:folds) {
    ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
    # subset.val = ind.val[dS1[ind.val]==1 & XS1[ind.val]<= t0 & Y[ind.val] > t0];
    # subset.tra = ind.tra[dS1[ind.tra]==1 & XS1[ind.tra]<= t0 & Y[ind.tra] > t0];
    # calculate the est.coef for P in the MSE
    P.md= loc.fun.ex(s.seq = s.seq[1], data=da2[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
    p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
    BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(Wi[ind.val]))
    BW.vec[k] = BW
  }
  replic[lay,] = BW.vec
}
mean(replic)
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}


#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr3 <- function(da3,t0,L,h,folds,reps, s.seq)
{
  Y     = da3$Y
  XS1   = da3$XS1
  XS2   = da3$XS2
  delta = da3$delta
  X     = data.matrix(subset(da3, select = -c(Y, XS1, XS2, delta, what, group)))
  Wi    = da3$what


  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)
  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.gr3 <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
    Wi    = data$what

    fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS2 <= t0)
    K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da3)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS2[ind.val]==1 & XS2[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS2[ind.tra]==1 & XS2[ind.tra]<= t0 & Y[ind.tra] > t0];
      # s.seq = s.seq
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr3(s.seq = s.seq[1], data=da3[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}


#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex.gr3 <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS2 <= t0)
  K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}


#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr4 <- function(da4,t0,L,h,folds,reps, s.seq)
{
  H=rbind(c(h[1],0), c(0,h[2]))
  Y     = da4$Y
  XS1   = da4$XS1
  XS2   = da4$XS2
  delta = da4$delta
  X     = data.matrix(subset(da4, select = -c(Y, XS1, XS2, delta, what, group)))
  Wi    = da4$what

  g.logit <- function(xx){exp(xx)/(1+exp(xx))}
  Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
  { out = (log(rbind(S1i,S2i))- log(ti))
  norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
  norm.k = det(Hi)^(-1/2) * norm.k
  norm.k
  }

  loc.fun.ex.gr4 <- function(s.seq, data, t0,  L, H, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    dS1   = data$deltaS1
    dS2   = data$deltaS2
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
    Wi    = data$what

    fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
    #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
    K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
    K=t(as.matrix(K))
    est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(t(as.matrix(s.seq))[,1])) {
      m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }


  n = dim(da4)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)

      #calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr4(s.seq = s.seq, data=da4[ind.tra,], t0=t0, L=L, H=H, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}


#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex.gr4 <- function(s.seq, data, t0,  L, H, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  dS1   = data$deltaS1
  dS2   = data$deltaS2
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
  #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
  K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
  K=t(as.matrix(K))
  est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(t(as.matrix(s.seq))[,1])) {
    m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}

###############
#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr2.SE <- function(da2,t0,L,h,folds,reps,s.seq)
{
  Y     = da2$Y
  XS1   = da2$XS1
  XS2   = da2$XS2
  delta = da2$delta
  X     = data.matrix(subset(da2, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da2$what

  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)
  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS1 <=t0)
    K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da2)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS1[ind.val]==1 & XS1[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS1[ind.tra]==1 & XS1[ind.tra]<= t0 & Y[ind.tra] > t0];
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.SE(s.seq = s.seq[1], data=da2[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr3.SE <- function(da3,t0,L,h,folds,reps, s.seq)
{
  Y     = da3$Y
  XS1   = da3$XS1
  XS2   = da3$XS2
  delta = da3$delta
  X     = data.matrix(subset(da3, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da3$what

  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)

  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.gr3.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS2 <= t0)
    K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da3)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS2[ind.val]==1 & XS2[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS2[ind.tra]==1 & XS2[ind.tra]<= t0 & Y[ind.tra] > t0];
      # s.seq = s.seq
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr3.SE(s.seq = s.seq[1], data=da3[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
min.BW.cens.ex.gr4.SE <- function(da4,t0,L,h,folds,reps, s.seq)
{
  H=rbind(c(h[1],0), c(0,h[2]))
  Y     = da4$Y
  XS1   = da4$XS1
  XS2   = da4$XS2
  delta = da4$delta
  X     = data.matrix(subset(da4, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da4$what

  g.logit <- function(xx){exp(xx)/(1+exp(xx))}
  Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
  { out = (log(rbind(S1i,S2i))- log(ti))
  norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
  norm.k = det(Hi)^(-1/2) * norm.k
  norm.k
  }
  loc.fun.ex.gr4.SE <- function(s.seq, data, t0,  L, H, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    dS1   = data$deltaS1
    dS2   = data$deltaS2
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
    #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
    K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
    K=t(as.matrix(K))
    est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(t(as.matrix(s.seq))[,1])) {
      m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da4)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)

      #calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr4.SE(s.seq = s.seq, data=da4[ind.tra,], t0=t0, L=L, H=H, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex.gr3.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS2 <= t0)
  K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}

#' helperfunction that should not be run by user directly
#' @noRd
#' @keywords internal
loc.fun.ex.gr4.SE <- function(s.seq, data, t0,  L, H, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  dS1   = data$deltaS1
  dS2   = data$deltaS2
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
  #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
  K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
  K=t(as.matrix(K))
  est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(t(as.matrix(s.seq))[,1])) {
    m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}


