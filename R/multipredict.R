#' Landmark prediction with multiple short-term events
#'
#' @include HelperFunctions.R
#' @param data Input dataset
#' @param formula a \code{formula} object, with a \code{Surv()} object, such as \code{Surv(time, event)},
#' on the left of a \code{~} operator, and the terms on the right. On the right-hand-side, the
#' time to the occurrence of short-term event 1 and short-term event 2 should be called by
#' statement \code{s1()} and \code{s2()}, respectively. The details of model specification
#' are given under ‘Details’
#' @param t0 Landmark time
#' @param L Length of time into the future (starting from the landmark time) for
#' which we want to make a risk prediction. This is called the `prediction horizon`
#' in the dynamic prediction literature
#' @param s1_beta1 A scalar or a vector. Time to the occurrence of short-term event 1 for the estimation
#' of the regression coefficient beta1 in group 2. If a \code{Null}
#' is given, then the coefficients for group 2 will NOT be estimated
#' @param s2_beta2 A scalar or a vector. Time to the occurrence of short-term event 2 for the estimation
#' of the regression coefficient beta2 in group 3.  If a \code{Null} is given, then the
#' coefficients for group 3 will NOT be estimated
#' @param s1s2_beta3 A matrix or a dataframe with two columns. The first column should be s1
#' and the second should be s2. Time to the occurrence of short-term event 1 & 2 for the estimation
#' of the regression coefficient beta3 in group 4. If a \code{Null}
#' is given, then the coefficients for group 4 will NOT be estimated.
#' @param SE Logical. `True` if user wants to estimate SE for the coefficient using
#' the perturbation-resampling method
#' @param SE.gs Logical. `True` if user wants to conduct grid search for the bandwidth in each
#' perturbation. It is expected to give more accurate results but will consume longer time.
#' `False` if user wants to use the same bandwidth found in the point estimation for
#' all perturbations
#' @param grid1 A prespecified grid for bandwidth search for group2
#' @param grid2 A prespecified grid for bandwidth search for group3
#' @param grid3 A list with prespecified grids for bandwidth search for group4
#' @param folds.grid The number of folds in cross-validation
#' @param reps.grid The number of repetitions of cross-validation
#' @param c01 A constant to shrink the bandwidth for group2
#' @param c02 A constant to shrink the bandwidth for group3
#' @param c03 A constant to shrink the bandwidth for group4
#' @param B Number of perturbations for estimating SE
#' @param gs.method Method used by gridsearch. Default is `loop`. Use `snow` will implement
#' parallel computing and will speed up the calculation
#' @param gs.cl Default is \code{Null}. Number of clusters used in parallel computing
#' in gridsearch. Specify when gs.method = `snow`
#' @param gs.seed An integer to set the seed for parallel computing to ensure reproducible
#' outcome, or `NULL` if not to set reproducible outcome
#'
#' @export
#'
#' @import survival
#' @import emdbook
#' @import NMOF
#' @import landpred
#' @import snow
#
#'
#' @return returns estimated coefficients for each short-term outcome and the long-term outcome:
#' \item{coefficients}{A named vector of the estimated regression coefficients}
#' \item{SE}{The standard error of coefficients estimated by perturbation resampling}
#'
#' @details The \code{multipredict} function fits time-fixed model and univariate/bivariate
#' varying-coefficient models using the data from subgroups formed based on the
#' information on the short-term outcomes (such as HF hospitalization and CHD hospitalization)
#' before landmark time t0, among those who haven't experienced the long-term outcome (such as death) at t0.
#' In this way the short-term outcome information are incorporated into the prediction
#' of long-term survival outcomes, and the risk prediction can vary based on the
#' event times of the short-term outcomes.
#'
#' The \code{+s1()} statement specified the column that determines the occurrence time of the first short-term outcome.
#' The \code{+s2()} statement specified the column that determines the occurrence time of the second short-term outcome.
#'
#' User may set the statement \code{gs.method} = `True`.
#'
#' By default the regression coefficients for group 1 is calculated in each run of this function.
#'
#' Currently, parameter estimates from parallel computing are slightly different in each run because of the
#' different (uncontrolled) random numbers used in the estimation. This will be solved in the near future.
#'
#'
#' @examples
#' library(survival)
#' library(emdbook)
#' library(NMOF)
#' library(landpred)
#' library(snow)
#' set.seed(1234)
#' res <- multipredict(data = simulation, formula = Surv(time, outcome) ~ age + s1(st1) + s2(st2),
#'                 t0 = 5, L = 20, SE = FALSE,
#'                 gs.method = "loop", gs.cl = 2, SE.gs = FALSE, B = 200, gs.seed = 100,
#'                 s1_beta1 = 1.5, grid1 = seq(0.01, 5, length.out=20),
#'                 s2_beta2 = 1.5, grid2 = seq(0.01, 5, length.out=20),
#'                 s1s2_beta3 = NULL, grid3=list(seq(0.01, 5, length.out=20),
#'                                                 seq(0.01, 5, length.out=20)))
#' print(res)
#'
#'
#' @author Wen Li, Qian Wang
#'
#' @references Li, Wen. (2023), "Landmarking Using A Flexible Varying Coefficient Model to Improve Prediction Accuracy of Long-term Survival Following Multiple Short-term Events An Application to the Atherosclerosis Risk in Communities (ARIC) Study",
#' \emph{Statistics in Medicine} \strong{90}(7) 1-29. doi:10.18637/jss.v090.i07
#' @references Parast, Layla, Su-Chun Cheng, and Tianxi Cai. (2012),
#' "Landmark Prediction of Long Term Survival Incorporating Short Term Event Time Information",
#' \emph{J Am Stat Assoc} \strong{107}(500) 1492-1501. doi: 10.1080/01621459.2012.721281
#' @references "Incorporating short-term outcome information to predict long-term survival with discrete markers".
#' \emph{Biometrical Journal} \strong{53.2} (2011): 294-307. doi: 10.1080/01621459.2012.721281






multipredict <- function(data, formula, t0, L,
                         SE=FALSE, SE.gs = FALSE,
                         s1_beta1=NULL, s2_beta2=NULL, s1s2_beta3=NULL,
                         grid1 = seq(0.01, 5, length.out=20),
                         grid2 = seq(0.01, 5, length.out=20),
                         grid3 = list(seq(0.01, 5, length.out=20),
                                      seq(0.01, 5, length.out=20)),
                         folds.grid = 8, reps.grid = 3,
                         c01 = 0.1,
                         c02 = 0.1,
                         c03 = 0.05,
                         B = 500,
                         gs.method = "loop",
                         gs.cl = NULL,
                         gs.seed = NULL
) {

  # -------------------------------------------------------#
  # Qian's new code for data manipulation
  mf <- model.frame(formula, data)

  pos_s1 <- grep("s1", names(mf))
  pos_s2 <- grep("s2", names(mf))
  # pos_delta <- grep("delta", names(mf))

  XS1 <- mf[[pos_s1]]
  XS2 <- mf[[pos_s2]]

  delta <- mf[[1]][,2]

  Y <- mf[[1]][,1]
  X <- mf[-c(1, pos_s1, pos_s2)]

  Xnames <- names(mf)[-c(1, pos_s1, pos_s2)]

  mydata <- as.data.frame(cbind(Y, delta, XS1, XS2, X))

  mydata$what=Wi.FUN(data=mydata, t0=t0, tau=L, weight.given=NULL)
  mydata$group=0
  mydata$group[(mydata$Y>t0)&(mydata$XS1>t0)&(mydata$XS2>t0)]=1
  mydata$group[(mydata$Y>t0)&(mydata$XS1<=t0)&(mydata$XS2>t0)]=2
  mydata$group[(mydata$Y>t0)&(mydata$XS1>t0)&(mydata$XS2<=t0)]=3
  mydata$group[(mydata$Y>t0)&(mydata$XS1<=t0)&(mydata$XS2<=t0)]=4

  da1 <- subset(mydata, group==1)
  da2 <- subset(mydata, group==2)
  da3 <- subset(mydata, group==3)
  da4 <- subset(mydata, group==4)
  n.f.b=nrow(mydata)

  #-------------------- calculate beta_t0 for group 1 -------------------------#

  fml1 <- as.formula(paste("1*(Y <= t0+L) ~", paste(Xnames, collapse = "+")))
  G1_model = glm(fml1,  data=da1, family = "binomial", weights = what)
  G1_coef = G1_model$coeff
  G1_con= G1_model$conv
  res1 <- list("group1_coef" = G1_coef,
               "group1_convergence" = G1_con)

  res <- list("group1" = res1)

  if (!is.null(s1_beta1)) {

    #-------------------- calculate beta_s1 for group2 ---------------------------#
    # find the bandwidth that optimize the MSE
    fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(Xnames, collapse = "+")))

    h2.min1 = gridSearch(fun = min.BW.cens.ex.gr2, da2=da2,t0=t0,L=L, s.seq = s1_beta1,
                         folds= folds.grid, reps= reps.grid, npar = 1, levels = list(a=grid1),
                         method = gs.method, cl=gs.cl)

    # shrink bandwidth to get the final bandwidth
    h2.final1= h2.min1$minlevels/(n.f.b^c01)

    # apply final bandwidth to the estimating equation
    G2_model1 = loc.fun.ex(s.seq=s1_beta1, data=da2, t0=t0,  L=L, h=h2.final1, weight = da2$what)
    G2_coef1 = G2_model1$est.mat
    G2_con1= G2_model1$conv


    res2 <- list("group2_coef" = G2_coef1,
                 "group2_convergence" = G2_con1)
    colnames(res2$group2_coef) <- c("Intercept", Xnames)
    res <- append(res, res2)
  }

  #-------------------- calculate beta_s2 for group3 ---------------------------#
  # find the bandwidth that optimize the MSE

  if (!is.null(s2_beta2)) {

    h3.min1 = gridSearch(fun = min.BW.cens.ex.gr3, da3=da3,t0=t0, L=L,
                         folds= folds.grid, reps= reps.grid, npar = 1, s.seq=s2_beta2, levels = list(a=grid2),
                         method = gs.method, cl=gs.cl)

    # shrink bandwidth to get the final bandwidth
    h3.final1= h3.min1$minlevels/(n.f.b^c02)

    # apply final bandwidth to the estimating equation
    G3_model1 = loc.fun.ex.gr3(s.seq=s2_beta2, data=da3, t0=t0,  L=L, h=h3.final1, weight = da3$what)
    G3_coef1 = G3_model1$est.mat
    G3_con1= G3_model1$conv

    res3 <- list("group3_coef" = G3_coef1,
                 "group3_convergence" = G3_con1)
    colnames(res3$group3_coef) <- c("Intercept", Xnames)
    res <- append(res, res3)
  }

  #-------------------- calculate beta_s3 for group4 ---------------------------#
  # find the bandwidth that optimize the MSE

  if (!is.null(s1s2_beta3)) {
    res4 <- c()
    res4c <- c()

    # convert s1s2_beta3 into required format
    if (class(s1s2_beta3)[1] == "numeric") {
      s1s2_beta3 <- matrix(s1s2_beta3, nrow = 1)} else {
        s1s2_beta3 <- as.matrix(s1s2_beta3)
      }

    if (gs.method == "snow" & !is.null(gs.cl)) {
      cl <- snow::makeCluster(c(rep("localhost", gs.cl)))
      clusterSetupRNGstream(cl, gs.seed)

      for (i in 1:nrow(s1s2_beta3)) {
        status4_11=gridSearch(fun=min.BW.cens.ex.gr4,da4=da4,t0=t0,L=L,
                              s.seq = s1s2_beta3[i,],
                              npar= 2, folds=folds.grid, reps=reps.grid,
                              levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                              method = gs.method, cl=gs.cl)
        h4.band11=status4_11$minlevels
        h4.final11=(sqrt(h4.band11)/(n.f.b^c03))^2
        G4_model11 = loc.fun.ex.gr4(s.seq=s1s2_beta3[i,], data=da4, t0=t0,
                                    L=L, H=rbind(c(h4.final11[1],0), c(0,h4.final11[2])), weight = da4$what)
        G4_coef11 = G4_model11$est.mat
        G4_con11 = G4_model11$conv[1,1]
        res4i <- G4_coef11
        res4ic <- G4_con11
        res4 <- rbind(res4, res4i)
        res4c <-rbind(res4c, res4ic)
        # colnames(res4$group4_coef) <- c("Intercept", Xnames)
      }

      stopCluster(cl)
    } else {
      for (i in 1:nrow(s1s2_beta3)) {
        status4_11=gridSearch(fun=min.BW.cens.ex.gr4,da4=da4,t0=t0,L=L,
                            s.seq = s1s2_beta3[i,],
                            npar= 2, folds=folds.grid, reps=reps.grid,
                            levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                            method = gs.method, cl=gs.cl)
        h4.band11=status4_11$minlevels
        h4.final11=(sqrt(h4.band11)/(n.f.b^c03))^2
        G4_model11 = loc.fun.ex.gr4(s.seq=s1s2_beta3[i,], data=da4, t0=t0,
                                  L=L, H=rbind(c(h4.final11[1],0), c(0,h4.final11[2])), weight = da4$what)
        G4_coef11 = G4_model11$est.mat
        G4_con11 = G4_model11$conv[1,1]
        res4i <- G4_coef11
        res4ic <- G4_con11
        res4 <- rbind(res4, res4i)
        res4c <-rbind(res4c, res4ic)
        # colnames(res4$group4_coef) <- c("Intercept", Xnames)
      }
    }

    res4 <- list("group4_coef" = res4, "group4_convergence" = res4c)
    colnames(res4$group4_coef) <- c("Intercept", Xnames)
    res <- append(res, res4)

  }



  #------------ Calculate empirical variance if `SE == T` ---------------------#
  if(SE) {

    if(SE.gs == T){
      resap.SE.vec <- c()
      for (i in 1:B) {
        mydata.o=mydata
        mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
        mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
        mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
        mydata.o$group=0
        mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2>t0)]=1
        da1.sap=mydata.o[mydata.o$group==1,]

        G1_model.sap = glm(fml1,  data=da1.sap, family = "binomial", weights = da1.sap$W.resap)
        G1_coef.sap = G1_model.sap$coeff
        resap.SE.vec=rbind(resap.SE.vec, G1_coef.sap)
      }
      resap.SE.vec=as.data.frame(resap.SE.vec)
      group1_SE=apply(resap.SE.vec, 2, sd)
      group1_SE <- as.data.frame(group1_SE)
      res <- append(res, group1_SE)
      names(res$group1_SE) = c("Intercept", Xnames)

      #------------------ SE for group2 ------------------#
      if (!is.null(s1_beta1)) {
        se.list2 <- list()
        for (i in 1:length(s1_beta1)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2>t0)]=2
            da2.sap=mydata.o[mydata.o$group==2,]


            if (gs.method == "snow" & !is.null(gs.cl)) {
              cl <- snow::makeCluster(c(rep("localhost", gs.cl)))
              clusterSetupRNGstream(cl, gs.seed)

              h2.min1 = gridSearch(fun = min.BW.cens.ex.gr2.SE, da2=da2.sap,t0=t0, L=L,
                                   folds= folds.grid, reps= reps.grid, npar = 1,
                                   s.seq=s1_beta1, levels = list(a=grid1),
                                   method = gs.method, cl= cl)
              stopCluster(cl)
            } else {
              h2.min1 = gridSearch(fun = min.BW.cens.ex.gr2.SE, da2=da2.sap,t0=t0, L=L,
                                 folds= folds.grid, reps= reps.grid, npar = 1,
                                 s.seq=s1_beta1, levels = list(a=grid1),
                                 method = gs.method, cl= cl)
              }

            # shrink bandwidth to get the final bandwidth
            h2.final1.sap= h2.min1$minlevels/(nrow(da2.sap)^c02)
            G2_model1.sap = loc.fun.ex.SE(s.seq=s1_beta1[i], data=da2.sap, t0=t0,  L=L, h=h2.final1.sap,
                                          weight = da2.sap$W.resap)
            G2_coef1.sap = G2_model1.sap$est.mat
            resap.SE.vec=rbind(resap.SE.vec, G2_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group2.row",i)
          se.list2[[df_name]] <- SE

          #resap.SE.seq = append(resap.SE.seq, G2_coef1.sap)
        }
        res <- append(res, se.list2)
      }
      #------------------ SE for group3 ------------------#
      if (!is.null(s2_beta2)) {

        se.list3 <- list()

        for (i in 1:length(s2_beta2)){
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2<=t0)]=3
            da3.sap=mydata.o[mydata.o$group==3,]


            if (gs.method == "snow" & !is.null(gs.cl)) {
              cl <- snow::makeCluster(c(rep("localhost", gs.cl)))
              clusterSetupRNGstream(cl, gs.seed)
              h3.min1 = gridSearch(fun = min.BW.cens.ex.gr3.SE, da3=da3.sap,t0=t0, L=L,
                                   folds= folds.grid, reps= reps.grid, npar = 1,
                                   s.seq=s2_beta2, levels = list(a=grid2),
                                   method = gs.method, cl=gs.cl)
              stopCluster(cl)
            } else {
              h3.min1 = gridSearch(fun = min.BW.cens.ex.gr3.SE, da3=da3.sap,t0=t0, L=L,
                                   folds= folds.grid, reps= reps.grid, npar = 1,
                                   s.seq=s2_beta2, levels = list(a=grid2),
                                   method = gs.method, cl=gs.cl)
            }

            # shrink bandwidth to get the final bandwidth
            h3.final1.sap= h3.min1$minlevels/(nrow(da3.sap)^c02)
            G3_model1.sap=loc.fun.ex.gr3.SE(s.seq=s2_beta2[i], data=da3.sap, t0=t0,  L=L, h=h3.final1.sap,
                                            weight = da3.sap$W.resap)
            G3_coef1.sap = G3_model1.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G3_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group3.row",i)
          se.list3[[df_name]] <- SE
        }
        res <- append(res, se.list3)
      }
      #------------------ SE for group4 ------------------#
      if (!is.null(s1s2_beta3)) {

        se.list4 <- list()

        for (i in 1:nrow(s1s2_beta3)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2<=t0)]=4
            da4.sap=mydata.o[mydata.o$group==4,]


            if (gs.method == "snow" & !is.null(gs.cl)) {
              cl <- snow::makeCluster(c(rep("localhost", gs.cl)))
              clusterSetupRNGstream(cl, gs.seed)
              status4_11=gridSearch(fun=min.BW.cens.ex.gr4.SE,da4=da4.sap,t0=t0,L=L,
                                    s.seq = s1s2_beta3[i,],
                                    npar= 2, folds=folds.grid, reps=reps.grid,
                                    levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                                    method = gs.method, cl=gs.cl)
              stopCluster(cl)
            } else {
              status4_11=gridSearch(fun=min.BW.cens.ex.gr4.SE,da4=da4.sap,t0=t0,L=L,
                                    s.seq = s1s2_beta3[i,],
                                    npar= 2, folds=folds.grid, reps=reps.grid,
                                    levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                                    method = gs.method, cl=gs.cl)
            }
            h4.band11_sap=status4_11$minlevels
            h4.final11_sap=(sqrt(h4.band11_sap)/(nrow(da4.sap)^.05))^2
            G4_model11.sap = loc.fun.ex.gr4.SE(s.seq=as.matrix(s1s2_beta3)[1,], data=da4.sap, t0=t0,
                                               L=L, H=rbind(c(h4.final11_sap[1],0), c(0,h4.final11_sap[2])), weight = da4.sap$W.resap)
            G4_coef11.sap = G4_model11.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G4_coef11.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group4.row",i)
          se.list4[[df_name]] <- SE
        }
        res <- append(res, se.list4)
      }

    } else if(SE.gs == F)
    {
      se.list1 <- list()
      resap.SE.vec <- c()
      for (i in 1:B) {
        # print(i)
        mydata.o=mydata
        mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
        mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
        mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
        mydata.o$group=0
        mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2>t0)]=1
        da1.sap=mydata.o[mydata.o$group==1,]

        G1_model.sap = glm(fml1,  data=da1.sap, family = "binomial", weights = da1.sap$W.resap)
        G1_coef.sap = G1_model.sap$coeff
        resap.SE.vec=rbind(resap.SE.vec, G1_coef.sap)
      }
      resap.SE.vec=as.data.frame(resap.SE.vec)
      group1_SE=apply(resap.SE.vec, 2, sd)
      group1_SE <- as.data.frame(group1_SE)
      res <- append(res, group1_SE)
      names(res$group1_SE) = c("Intercept", Xnames)

      #------------------ SE for group2 ------------------#
      if (!is.null(s1_beta1)) {

        se.list2 <- list()

        for (i in 1:length(s1_beta1)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2>t0)]=2
            da2.sap=mydata.o[mydata.o$group==2,]

            G2_model1.sap = loc.fun.ex.SE(s.seq=s1_beta1[i], data=da2.sap, t0=t0,  L=L, h=h2.final1,
                                          weight = da2.sap$W.resap)
            G2_coef1.sap = G2_model1.sap$est.mat
            resap.SE.vec=rbind(resap.SE.vec, G2_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group2.row",i)
          se.list2[[df_name]] <- SE

        }
        res <- append(res, se.list2)
      }
      #------------------ SE for group3 ------------------#
      if (!is.null(s2_beta2)) {

        se.list3 <- list()

        for (i in 1:length(s2_beta2)){
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2<=t0)]=3
            da3.sap=mydata.o[mydata.o$group==3,]
            G3_model1.sap = loc.fun.ex.gr3.SE(s.seq=s2_beta2, data=da3.sap, t0=t0,  L=L, h=h3.final1,
                                              weight = da3.sap$W.resap)
            G3_coef1.sap = G3_model1.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G3_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group3.row",i)
          se.list3[[df_name]] <- SE
        }
        res <- append(res, se.list3)
      }
      #------------------ SE for group4 ------------------#
      if (!is.null(s1s2_beta3)) {

        se.list4 <- list()

        for (i in 1:nrow(s1s2_beta3)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2<=t0)]=4
            da4.sap=mydata.o[mydata.o$group==4,]

            # status4_11=gridSearch(fun=min.BW.cens.ex.gr4.SE,da4=da4.sap,t0=t0,L=L,
            #                       s.seq = s1s2_beta3[i,],
            #                       npar= 2, folds=folds.grid, reps=reps.grid,
            #                       levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
            #                       method = gs.method, cl=gs.cl)
            h4.band11_sap=status4_11$minlevels
            h4.final11_sap=(sqrt(h4.band11_sap)/(nrow(da4.sap)^.05))^2
            G4_model11.sap = loc.fun.ex.gr4.SE(s.seq=as.matrix(s1s2_beta3)[1,], data=da4.sap, t0=t0,
                                               L=L, H=rbind(c(h4.final11_sap[1],0), c(0,h4.final11_sap[2])), weight = da4.sap$W.resap)
            G4_coef11.sap = G4_model11.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G4_coef11.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group4.row",i)
          se.list4[[df_name]] <- SE
        }
        res <- append(res, se.list4)
      }
    } else {stop("`SE.gs` must be logical")}

  }
  return(res)
}


