
###################################################################
# rdwinselect: window selection for randomization inference in RD
# !version 0.3 13-Mar-2018
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Window selection for RD designs under local randomization
#'
#' \code{rdwinselect} implements the window-selection procedure
#'  based on balance tests for RD designs under local randomization.
#'  Specifically, it constructs a sequence of nested windows around the RD cutoff
#'  and reports binomial tests for the running variable runvar and covariate balance
#'  tests for covariates covariates (if specified). The recommended window is the
#'  largest window around the cutoff such that the minimum p-value of the balance test
#'  is larger than a prespecified level for all nested (smaller) windows. By default,
#'  the p-values are calculated using randomization inference methods.
#'
#' @author
#' Matias Cattaneo, University of Michigan. \email{cattaneo@umich.edu}
#'
#' Rocio Titiunik, University of Michigan. \email{titiunik@umich.edu}
#'
#' Gonzalo Vazquez-Bare, University of Michigan. \email{gvazquez@umich.edu}
#'
#' @references
#' M.D. Cattaneo, B. Frandsen and R. Titiunik. (2015).  \href{https://sites.google.com/site/rdpackages/rdlocrand/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate}. \emph{Journal of Causal Inference} 3(1): 1-24.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2016). \href{https://sites.google.com/site/rdpackages/rdlocrand/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf}{Inference in Regression Discontinuity Designs under Local Randomization}. \emph{Stata Journal} 16(2): 331-367.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2017). \href{https://sites.google.com/site/rdpackages/rdlocrand/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf}{Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality}. \emph{Journal of Policy Analysis and Management} 36(3): 643-681.
#'
#'
#' @param R a vector containing the values of the running variable.
#' @param X the matrix of covariates to be used in the balancing tests. The matrix is optional but the recommended window is only provided when at least one covariate is specified. This should be a matrix of size n x k where n is the total sample size and $k$ is the number of covariates.
#' @param cutoff the RD cutoff (default is 0).
#' @param obsmin the minimum number of observations above and below the cutoff in the smallest window. Default is 10.
#' @param wmin the smallest window to be used.
#' @param wobs the number of observations to be added at each side of the cutoff at each step.
#' @param wstep the increment in window length.
#' @param nwindows the number of windows to be used. Default is 10.
#' @param statistic the statistic to be used in the balance tests. Allowed options are \code{diffmeans} (difference in means statistic), \code{ksmirnov} (Kolmogorov-Smirnov statistic), \code{ranksum} (Wilcoxon-Mann-Whitney standardized statistic) and \code{hotelling} (Hotelling's T-squared statistic). Default option is \code{diffmeans}. The statistic \code{ttest} is equivalent to \code{diffmeans} and included for backward compatibility.
#' @param approx forces the command to conduct the covariate balance tests using a large-sample approximation instead of finite-sample exact randomization inference methods.
#' @param p the order of the polynomial for outcome adjustment model (for covariates). Default is 0.
#' @param evalat specifies the point at which the adjusted variable is evaluated. Allowed options are \code{cutoff} and \code{means}. Default is \code{cutoff}.
#' @param kernel specifies the type of kernel to use as weighting scheme. Allowed kernel types are \code{uniform} (uniform kernel), \code{triangular} (triangular kernel) and \code{epan} (Epanechnikov kernel). Default is \code{uniform}.
#' @param reps number of replications. Default is 1000.
#' @param seed the seed to be used for the randomization tests.
#' @param level the minimum accepted value of the p-value from the covariate balance tests. Default is .15.
#' @param plot draws a scatter plot of the minimum p-value from the covariate balance test against window length.
#' @param quietly suppress output
#' @param obsstep the minimum number of observations to be added on each side of the cutoff for the sequence of fixed-increment nested windows. Default is 2. This option is deprecated and only included for backward compatibility.
#'
#' @return
#' \item{window}{recommended window (NA is covariates are not specified)}
#' \item{results}{table including window lengths, minimum p-value in each window, corresponding number of the variable with minimum p-value (i.e. column of covariate matrix), Binomial test p-value and sample sizes to the left and right of the cutoff in each window.}
#' \item{summary}{summary statistics.}
#'
#' @examples
#' # Toy dataset
#' X <- array(rnorm(200),dim=c(100,2))
#' R <- X[1,] + X[2,] + rnorm(100)
#' # Window selection adding 5 observations at each step
#' # Note: low number of replications to speed up process.
#' tmp <- rdwinselect(R,X,obsmin=10,wobs=5,reps=500)
#' # Window selection setting initial window and step
#' # The user should increase the number of replications.
#' tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,reps=500)
#' # Window selection with approximate (large sample) inference and p-value plot
#' tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,approx=TRUE,nwin=80,quietly=TRUE,plot=TRUE)
#'
#'
#' @export


rdwinselect = function(R, X,
                       cutoff=0,
                       obsmin=NULL,
                       wmin = NULL,
                       wobs = NULL,
                       wstep = NULL,
                       nwindows = 10,
                       statistic = 'diffmeans',
                       approx = FALSE,
                       p = 0,
                       evalat = 'cutoff',
                       kernel = 'uniform',
                       reps = 1000,
                       seed = NULL,
                       level = .15,
                       plot = FALSE,
                       quietly = FALSE,
                       obsstep = NULL) {


  #################################################################
  # Parameters and error checking
  #################################################################

  if (cutoff<=min(R,na.rm=TRUE) | cutoff>=max(R,na.rm=TRUE)){stop('Cutoff must be within the range of the running variable')}
  if (p<0){stop('p must be a positive integer')}
  if (p>0 & approx==TRUE & statistic!='ttest' & statistic!='diffmeans'){stop('approximate and p>1 can only be combined with ttest')}
  if (statistic!='diffmeans' & statistic!='ttest' & statistic!='ksmirnov' & statistic!='ranksum' & statistic!='hotelling'){stop(paste(statistic,'not a valid statistic'))}
  if (evalat!='cutoff' & evalat!='means'){stop('evalat only admits means or cutoff')}
  if (kernel!='uniform' & kernel!='triangular' & kernel!='epan'){stop(paste(kernel,'not a valid kernel'))}
  if (kernel!='uniform' & evalat!='cutoff'){stop('kernel can only be combined with evalat(cutoff)')}
  if (kernel!='uniform' & statistic!='ttest' & statistic!='diffmeans'){stop('kernel only allowed for diffmeans')}

  Rc = R - cutoff
  D = Rc >= 0
  
  if (!missing(X)){
    data = data.frame(Rc,D,X)
    data = data[order(data$Rc),]
    Rc = data[,1]
    D = data[,2]
    X = data[,c(-1,-2)]
  } else {
    data = data.frame(Rc,D)
    data = data[order(data$Rc),]
    Rc = data[,1]
    D = data[,2]
  }

  n = length(Rc[!is.na(Rc)])
  n1 = sum(D,na.rm=TRUE)
  n0 = n - n1
  dups = merge(Rc,table(Rc),by=1)
  dups = dups[,2]

  if (!missing(X)){X = as.matrix(X)}
  if (!is.null(seed)){set.seed(seed)}

  if (approx==FALSE){testing.method='rdrandinf'}else{testing.method='approximate'}

  ## Display upper-right panel

  if (quietly==FALSE){
    cat('\n')
    cat('\nWindow selection for RD under local randomization \n')
    cat('\n')
    cat(paste0(format('Number of obs  =', width=25), toString(n))); cat("\n")
    cat(paste0(format('Order of poly  =', width=25), p)); cat("\n")
    cat(paste0(format('Kernel type    =', width=25), kernel)); cat("\n")
    cat(paste0(format('Reps           =', width=25), reps)); cat("\n")
    cat(paste0(format('Testing method =', width=25), testing.method)); cat("\n")
    cat(paste0(format('Balance test   =', width=25), statistic))
    cat('\n\n')
  }


  #################################################################
  # Initial window, steps and window list
  #################################################################

  ## Define smallest window

  if (!is.null(wmin)){
    if (!is.null(obsmin)){
      stop('Cannot specify both obsmin and wmin')
    } else{
      obsmin = 10
    }
  } else{
    if (is.null(obsmin)){
      obsmin = 10
    }
    posl = sum(Rc<0,na.rm=TRUE)
    posr = posl + 1
    wmin = findwobs(obsmin,1,posl,posr,Rc,dups)
  }

  ## Define step
  
  if (is.null(wstep) & is.null(wobs)){
    if (is.null(obsstep)){
      obsstep = 2
    }
    wstep = findstep(Rc,D,obsmin,obsstep,10)
    window.list = seq(from=wmin,by=wstep,length.out=nwindows)
  } else if (!is.null(wstep)){
    if (!is.null(wobs)){
      stop('Cannot specify both wobs and wstep')
    }
    if (!is.null(obsstep)){
      stop('Cannot specify both obsstep and wstep')
    }
    window.list = seq(from=wmin,by=wstep,length.out=nwindows)
  } else if (!is.null(wobs)){
    pos0 = sum(Rc<0,na.rm=TRUE)
    posl = pos0 - sum(Rc<0 & Rc>=-wmin,na.rm=TRUE)
    posr = pos0 + 1 + sum(Rc>=0 & Rc<wmin,na.rm=TRUE)
    window.list = findwobs(wobs,nwindows-1,posl,posr,Rc[!is.na(Rc)],dups[!is.na(Rc)])
    window.list = c(wmin,window.list)
  }


  #################################################################
  # Summary statistics
  #################################################################

  table.sumstats = array(NA,dim=c(5,2))
  table.sumstats[1,] = c(n0,n1)

  qq0 = round(quantile(abs(Rc[D==0]),probs = c(.01,.05,.1,.2),type=1,na.rm=TRUE),5)
  qq1 = round(quantile(Rc[D==1],probs = c(.01,.05,.1,.2),type=1,na.rm=TRUE),5)

  n0.q1 = sum(Rc>=-qq0[1]& Rc<0,na.rm=TRUE)
  n0.q2 = sum(Rc>=-qq0[2]& Rc<0,na.rm=TRUE)
  n0.q3 = sum(Rc>=-qq0[3]& Rc<0,na.rm=TRUE)
  n0.q4 = sum(Rc>=-qq0[4]& Rc<0,na.rm=TRUE)
  n1.q1 = sum(Rc<=qq1[1]& Rc>=0,na.rm=TRUE)
  n1.q2 = sum(Rc<=qq1[2]& Rc>=0,na.rm=TRUE)
  n1.q3 = sum(Rc<=qq1[3]& Rc>=0,na.rm=TRUE)
  n1.q4 = sum(Rc<=qq1[4]& Rc>=0,na.rm=TRUE)

  table.sumstats[2,] = c(n0.q1,n1.q1)
  table.sumstats[3,] = c(n0.q2,n1.q2)
  table.sumstats[4,] = c(n0.q3,n1.q3)
  table.sumstats[5,] = c(n0.q4,n1.q4)

  ## Display upper left panel

  if (quietly==FALSE){
    cat(paste0(format(paste0("Cutoff c = ", toString(round(cutoff, 3))), width=25), format("Left of c", width=20), format("Right of c", width=20))); cat("\n")
    cat(paste0(format("Number of obs",   width=25), format(toString(n0), width=20), format(toString(n1), width=20))); cat("\n")
    cat(paste0(format("1st percentile",  width=25), format(toString(n0.q1), width=20), format(toString(n1.q1), width=20))); cat("\n")
    cat(paste0(format("5th percentile",  width=25), format(toString(n0.q2), width=20), format(toString(n1.q2), width=20))); cat("\n")
    cat(paste0(format("10th percentile", width=25), format(toString(n0.q3), width=20), format(toString(n1.q3), width=20))); cat("\n")
    cat(paste0(format("20th percentile", width=25), format(toString(n0.q4), width=20), format(toString(n1.q4), width=20)))
    cat("\n\n")
  }


  #################################################################
  # Balance tests
  #################################################################

  table.rdw = array(NA,dim=c(length(window.list),5))
  col = 1

  ## Being main panel display

  if (quietly==FALSE){
    cat(paste(format('Window length / 2', width=23),
              format('p-value ',          width=12),
              format('Var. name',         width=15),
              format('Bin.test',          width=12),
              format('Obs<c',             width=6),
              format('Obs>=c',            width=5)))
    cat('\n\n')
  }

  for (w in window.list){

    ww = (Rc >= -w) & (Rc <= w)

    Dw = D[ww]
    Rw = Rc[ww]

    ## Drop NA values

    if (!missing(X)){
      Xw = X[ww,]
      data = cbind(Rw,Dw,Xw)
      data = data[complete.cases(data),]
      Rw = data[,1]
      Dw = data[,2]
      Xw = data[,c(-1,-2)]
    } else {
      data = cbind(Rw,Dw)
      data = data[complete.cases(data),]
      Rw = data[,1]
      Dw = data[,2]
    }

    ## Sample sizes

    n0.w = sum(Rw<0)
    n1.w = sum(Rw>=0)
    n.w = n0.w+n1.w
    table.rdw[col,4] = n0.w
    table.rdw[col,5] = n1.w

    ## Binomial test

    bitest = binom.test(sum(Dw),length(Dw),p=0.5)
    table.rdw[col,3] = bitest$p.value

    if (!missing(X)){

      ## Weights

      kweights = rep(1,n.w)

      if (kernel=='triangular'){
        kweights = (1-abs(Rw/w))*(abs(Rw/w)<1)
      }
      if (kernel=='epan'){
        kweights = .75*(1-(Rw/w)^2)*(abs(Rw/w)<1)
      }

      ## Model adjustment

      if (p>0){
        X.adj = NULL
        Xk.adj = NULL
        if (evalat==''|evalat=='cutoff'){
          evall = cutoff
          evalr = cutoff
        } else {
          evall = mean(Rw[D==0])
          evalr = mean(Rw[D==1])
        }
        R.adj = Rw - Dw*evalr - (1-Dw)*evall

        for (k in 1:ncol(X)){
          Rpoly = poly(R.adj,order=p,raw=TRUE)
          lfit.t = lm(Xw[Dw==1,k] ~ Rpoly[Dw==1,],weights=kweights[Dw==1])
          Xk.adj[Dw==1] = lfit.t$residuals + lfit.t$coefficients[1]
          lfit.c = lm(Xw[Dw==0,k] ~ Rpoly[Dw==0,],weights=kweights[Dw==0])
          Xk.adj[Dw==0] = lfit.c$residuals + lfit.c$coefficients[1]
          X.adj = cbind(X.adj,Xk.adj)
        }
        Xw = X.adj
      }


      ## Statistics and p-values

      if (statistic=='hotelling'){
        aux = hotelT2(Xw,Dw)
        obs.stat = as.numeric(aux$statistic)
        if (approx==FALSE){
          stat.distr = array(NA,dim=c(reps,1))
          for (i in 1:reps){
            D.sample = sample(Dw,replace=FALSE)
            aux.sample = hotelT2(Xw,D.sample)
            obs.stat.sample = as.numeric(aux.sample$statistic)
            stat.distr[i] = obs.stat.sample
          }
          p.value = mean(abs(stat.distr) >= abs(obs.stat))
        } else {
          p.value = as.numeric(aux$p.value)
        }
        table.rdw[col,1] = p.value
        varname = NA

      } else {
        aux = rdrandinf.model(Xw,Dw,statistic=statistic,kweights=kweights,pvalue=TRUE)
        obs.stat = as.numeric(aux$statistic)
        if (approx==FALSE){
          stat.distr = array(NA,dim=c(reps,ncol(X)))
          for (i in 1:reps){
            D.sample = sample(Dw,replace=FALSE)
            aux.sample = rdrandinf.model(Xw,D.sample,statistic=statistic,kweights=kweights)
            obs.stat.sample = as.numeric(aux.sample$statistic)
            stat.distr[i,] = obs.stat.sample
          }
          p.value = rowMeans(t(abs(stat.distr)) >= abs(obs.stat))
        } else {
          if (p==0){
            p.value = as.numeric(aux$p.value)
          } else {
            p.value = NULL
            for (k in 1:ncol(X)){
              lfit = lm(Xw[,k] ~ Dw + Rpoly + Dw*Rpoly,kweights=kweights)
              tstat = lfit$coefficients['Dw']/sqrt(vcov(lfit)['Dw','Dw'])
              p.value = c(p.value,2*pnorm(-abs(tstat)))
            }
          }
        }

        table.rdw[col,1] = min(p.value)
        tmp = which.min(p.value)
        table.rdw[col,2] = tmp

        if (!is.null(colnames(X)[tmp])){
          if (colnames(X)[tmp]!='') {
            varname = substring(colnames(X)[tmp],1,15)
          }
          else{
            varname = tmp
          }
        } else {
          varname = tmp
        }
      }

    } else {
      table.rdw[col,1] = NA
      table.rdw[col,2] = NA
      varname = NA
    }

    if (quietly==FALSE){
      cat(paste0(format(toString(round(w,3)),width=25),
                 format(toString(round(table.rdw[col,1],3)),width=10),
                 format(toString(varname),width=20),
                 format(toString(round(table.rdw[col,3],3)),width=13),
                 format(toString(table.rdw[col,4]),width=7),
                 format(toString(table.rdw[col,5]),width=5)))
      cat('\n')
    }
    col = col + 1

  }


  #################################################################
  # Find recommended window
  #################################################################

  if (!missing(X)){
    Pvals = table.rdw[,1]
    if (all(Pvals>level)){
      tmp = length(Pvals)
    } else {
      tmp = min(which(Pvals<level))
      if (tmp==1){
        warning('Smallest window does not pass covariate test. Decrease smallest window or reduce level')
      } else {
        tmp = tmp - 1
      }
    }

    rec.length = window.list[tmp]
    rec.window = c(cutoff-rec.length,cutoff+rec.length)

    if (quietly==FALSE & tmp!=1){
      cat('\n\n')
      cat(paste0('Recommended window is [',round(rec.window[1],3),';',round(rec.window[2],3),'] with ',table.rdw[tmp,4]+table.rdw[tmp,5],' observations (',
                 table.rdw[tmp,4],' below, ',table.rdw[tmp,5],' above).'))
      cat('\n\n')
    }

  } else {
    if(quietly==FALSE){warning('Need to specify covariates to find recommended length')}
    rec.window = NA
  }



  #################################################################
  # Plot p-values
  #################################################################

  if (plot==TRUE){
    if (!missing(X)){
      plot(window.list,Pvals)
    } else {
      stop('Cannot draw plot without covariates')
    }
  }


  #################################################################
  # Output
  #################################################################

  colnames(table.sumstats) = c('Left of c','Right of c')
  rownames(table.sumstats) = c('Number of obs','1th percentile','5th percentile','10th percentile','20th percentile')

  table.rdw = cbind(window.list,table.rdw)
  colnames(table.rdw)= c('W.length','p-value','Variable','Bi.test','Obs<c','Obs>=c')

  output = list(window = rec.window,
                results = table.rdw,
                summary = table.sumstats)

  return(output)

}
