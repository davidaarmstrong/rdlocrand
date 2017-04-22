
###################################################################
# rdsensitivity: sensitivity analysis for randomization inference in RD
# !version 0.2 22-Apr-2017
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Sensitivity analysis for RD designs under local randomization
#'
#' \code{rdsensitivity} analyze the sensitivity of randomization p-values
#' and confidence intervals to different window lengths.
#'
#'
#' @author
#' Matias Cattaneo, University of Michigan. \email{cattaneo@umich.edu}
#'
#' Rocio Titiunik, University of Michigan. \email{titiunik@umich.edu}
#'
#' Gonzalo Vazquez-Bare, University of Michigan. \email{gvazquez@umich.edu}
#'
#' @references
#' M.D. Cattaneo, B. Frandsen and R. Titiunik. (2015).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate}. \emph{Journal of Causal Inference} 3(1): 1-24.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2016). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf}{Inference in Regression Discontinuity Designs under Local Randomization}. \emph{Stata Journal} 16(2): 331-367.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2017). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf}{Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality}. \emph{Journal of Policy Analysis and Management}, forthcoming.
#'
#'
#' @param Y a vector containing the values of the outcome variable.
#' @param R a vector containing the values of the running variable.
#' @param cutoff the RD cutoff (default is 0).
#' @param wlist the list of window lengths to be evaluated. By default the program constructs 10 windows around the cutoff, the first one including 10 treated and control observations and adding 5 observations to each group in subsequent windows.
#' @param tlist the list of values of the treatment effect under the null to be evaluated. By default the program employs ten evenly spaced points within the asymptotic confidence interval for a constant treatment effect in the smallest window to be used.
#' @param ci returns the confidence interval corresponding to the indicated window length. \code{ci} has to be a scalar or a two-dimensional vector, where the first value needs to be one of the values in \code{wlist}. The second value, if specified, indicates the level of the confidence interval. Default level is .05 (95\% level CI).
#' @param statistic the statistic to be used in the balance tests. Allowed options are \code{diffmeans} (difference in means statistic), \code{ksmirnov} (Kolmogorov-Smirnov statistic) and \code{ranksum} (Wilcoxon-Mann-Whitney standardized statistic). Default option is \code{diffmeans}. The statistic \code{ttest} is equivalent to \code{diffmeans} and included for backward compatibility.
#' @param nodraw suppresses contour plot.
#' @param p the order of the polynomial for outcome adjustment model. Default is 0.
#' @param evalat specifies the point at which the adjusted variable is evaluated. Allowed options are \code{cutoff} and \code{means}. Default is \code{cutoff}.
#' @param kernel specifies the type of kernel to use as weighting scheme. Allowed kernel types are \code{uniform} (uniform kernel), \code{triangular} (triangular kernel) and \code{epan} (Epanechnikov kernel). Default is \code{uniform}.
#' @param reps number of replications. Default is 1000.
#' @param seed the seed to be used for the randomization tests.
#' @param fuzzy indicates that the RD design is fuzzy. \code{fuzzy} can be specified as a vector containing the values of the endogenous treatment variable, or as a list where the first element is the vector of endogenous treatment values and the second element is a string containing the name of the statistic to be used. Allowed statistics are \code{ar} (Anderson-Rubin statistic) and \code{tsls} (2SLS statistic). Default statistic is \code{ar}. The \code{tsls} statistic relies on large-sample approximation.
#' @param quietly suppresses the output table.
#'
#' @return
#' \item{tlist}{treatment effects grid}
#' \item{wlist}{window grid}
#' \item{results}{table with corresponding p-values for each window and treatment effect pair.}
#' \item{ci}{confidence interval (if \code{ci} is specified).}
#'
#' @examples
#' # Toy dataset
#' R <- runif(100,-1,1)
#' Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(100)
#' # Sensitivity analysis
#' # Note: low number of replications to speed up process.
#' # The user should increase the number of replications.
#' tmp <- rdsensitivity(Y,R,wlist=seq(.75,2,by=.25),tlist=seq(0,5,by=1),reps=500)
#'
#'
#' @export


rdsensitivity = function(Y,R,
                         cutoff = 0,
                         wlist,
                         tlist,
                         ci,
                         statistic = 'diffmeans',
                         nodraw = FALSE,
                         p = 0,
                         evalat = 'cutoff',
                         kernel = 'uniform',
                         reps = 1000,
                         seed = '',
                         fuzzy = '',
                         quietly = FALSE){



  #################################################################
  # Parameters and error checking
  #################################################################

  if (cutoff<=min(R,na.rm=TRUE) | cutoff>=max(R,na.rm=TRUE)){stop('Cutoff must be within the range of the running variable')}
  if (statistic!='diffmeans' & statistic!='ttest' & statistic!='ksmirnov' & statistic!='ranksum'){stop(paste(statistic,'not a valid statistic'))}
  if (evalat!='cutoff' & evalat!='means'){stop('evalat only admits means or cutoff')}


  data = cbind(Y,R)
  data = data[complete.cases(data),]
  Y = data[,1]
  R = data[,2]

  Rc = R - cutoff


  #################################################################
  # Default window list
  #################################################################

  if (missing(wlist)){
    aux = rdwinselect(Rc,obsstep=5,quietly=TRUE)
    wlist = round(aux$results[,1],2)
  }


  #################################################################
  # Default tau list
  #################################################################

  if (missing(tlist)){
    D = as.numeric(Rc >= 0)
    aux = lm(Y ~ D)
    ci.ub = round(aux$coefficients['D']+1.96*sqrt(vcov(aux)['D','D']),2)
    ci.lb = round(aux$coefficients['D']-1.96*sqrt(vcov(aux)['D','D']),2)
    wstep = round((ci.ub-ci.lb)/10,2)
    tlist = seq(ci.lb,ci.ub,by=wstep)
  }


  #################################################################
  # Sensitivity analysis
  #################################################################

  results = array(NA,dim=c(length(tlist),length(wlist)))
  if (quietly==FALSE) {cat('\nRunning sensitivity analysis...\n')}
  row = 1
  for (t in tlist){
    col = 1
    for (w in wlist){
      if (evalat=='means'){
        ww = (Rc >= -w) & (Rc <= w)
        Rw = R[ww]
        Dw = D[ww]
        evall = mean(Rw[Dw==0])
        evalr = mean(Rw[Dw==1])
      } else{
        evall=''
        evalr=''
      }

      aux = rdrandinf(Y,Rc,wl=-w,wr=w,p=p,reps=reps,nulltau=t,
                      statistic=statistic,kernel=kernel,evall=evall,evalr=evalr,
                      fuzzy=fuzzy,seed=seed,quietly=TRUE)
      results[row,col] = aux$p.value
      col = col + 1
    }
    row = row + 1
  }
  if (quietly==FALSE) {cat('\nSensitivity analysis complete.\n')}


  #################################################################
  # Confidence interval
  #################################################################

  if (!missing(ci)){
    ci.window = ci[1]
    if (length(ci)>1){
      ci.level = ci[2]
      if (ci.level<=0 | ci.level>=1){stop('ci level has to be between 0 and 1')}
    }else {
      ci.level = .05
    }
    if (is.element(ci.window,wlist)==TRUE){
      col = which(wlist==ci.window)
      aux = results[,col]
      if (all(aux>ci.level)){
        index.lb = 1
        index.ub = length(aux)
        ci.lb = tlist[index.lb]
        ci.ub = tlist[index.ub]
      } else if (all(aux<ci.level)){
        warning('no valid confidence interval in specified window grid')
        ci.lb = NA
        ci.ub = NA
      } else {
        index.lb = min(which(aux>=.05))
        index.ub = max(which(aux>=.05))
        ci.lb = tlist[index.lb]
        ci.ub = tlist[index.ub]
      }


      conf.int = c(ci.lb,ci.ub)

    } else{
      stop('window specified in ci not in wlist')
    }
  }


  #################################################################
  # Output
  #################################################################

  if (missing(ci)){
    output = list(tlist = tlist, wlist = wlist, results = results)
  } else {
    output = list(tlist = tlist, wlist = wlist, results = results, ci = conf.int)
  }


  #################################################################
  # Plot
  #################################################################

  if (nodraw==FALSE){
    if (dim(results)[2]==1){
      warning('need a window grid to draw plot')
    } else if (dim(results)[1]==1){
      warning('need a tau grid to draw plot')
    } else {
      filled.contour(wlist,tlist,t(results),
                     xlab='window',ylab='treatment effect',
                     key.title=title(main = 'p-value',cex.main=.8),
                     levels=seq(0,1,by=.01),col=gray.colors(100,1,0))

    }
  }

  return(output)

}


