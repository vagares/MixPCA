
################################################################################

doOne <- function(K
{
  
  data = data(n=c(100,100,100),K ,q = 4,p = 10,nx=4,
              pi = c(0.2,0.35,0.45),
              mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3),
              beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4),
              SNR1 = 100,
              SNR2 = 3)
  
  est = estimates(data,K,maxits)
  result <-c(pi = est$pi,mu=est$mu,beta=est$beta,Q=est$Q,theta2=est$theta2,sigma2=est$sigma2)
  return(result)
}


vlis1 <- function(nsim=1000) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    k = list(value = c(2,3))
  )
  return(vList)
}


runSims1 <- function(vList=vlis1(), doOne ,parallel=TRUE,
                     seedList=NULL, # set_seed(vList=vList)
                     sfile = NULL # .rds save simulation
) {

  if(parallel){
    res <- simsalapar::doMclapply(vList,  doOne = doOne, seed=NULL,
                                  sfile = sfile, #"res42_n1000_lapply_LEc.rds"
                                  monitor = interactive())
  } else{
    res <- doLapply(vList, doOne = doOne,
                    sfile = sfile, 
                    monitor = interactive())
  }
  return(res)
}
#########################################################################
