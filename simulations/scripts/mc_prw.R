#	Loading required libraries. The dependencies should be installed as well. 
require('parsimonious')
require('reshape2')
require('ggplot2')
require('xtable')
require('parallel')
require('knitr')

#Sourcing the subs
source('subs/mc_parsimonious.R')


ncores 	<- 16
alasso 	<- TRUE
mc 	<- 10000
npath	<- 5

set.seed(123)

r <- 1
veps = 0.1
# sample length
nT   <- 100
varpar <- c(0.1,1,10)	
alphaT <- c(0.01,0.031,0.1)

mcxp   <- list()
#mcpath <- list()

for(vpar in varpar){
	for(at in alphaT){
		
		stime <- proc.time()
			
		# Estimation
		#mcpar <- lapply(1:mc,mcptv,X,ppath,s2,alasso=alasso)
		mcpar <- mclapply(1:mc,mcptv,r=r,nT=nT,alpha=at,veps=veps,vpar=vpar,alasso=alasso,peninit=FALSE,mc.cores=ncores)
		mcxp[[paste0('Var(eta)==',vpar,' alpha[T]==',at)]] <- mcpar
	#	mcpath[[paste0('sigma_nu =',vpar,' alpha=',at)]] 	<- ppath
		
		
		tm <- (proc.time()-stime)[3]
		cat('\nDone with: ',paste0('sigma_nu=',vpar,', alpha=',at,' in ',round(tm,1),'sec. \n'))
		cat('Average Time per iteration: ',tm/(mc*ncores),'\n')
	}
}

xpn <- 'prw'

#save(file=paste0('mcsaves/mc_path_',xpn),mcpath)
save(file=paste0('mcsaves/mc_',xpn),mcxp)
