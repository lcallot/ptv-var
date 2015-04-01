#	Loading required libraries. The dependencies should be installed as well. 
require('parsimonious')
require('reshape2')
require('ggplot2')
require('xtable')
require('parallel')

#Sourcing the subs
source('subs/mc_parsimonious.R')


set.seed(123)
ncores 	<- 16
alasso 	<- TRUE
mc 	<- 10000
npath	<- 5
r <- 1


#storage
mcxp   <- list()
mcX	   <- list()
mcpath <- list()

# sample length 
vT <- c(100,1000)
# Varying residual variance
n2s <- c(0.1,1,10)

for(nT in vT){
	for(ns in n2s){
		
		stime <- proc.time()

		# Parameter generation:
		path <- matrix(1,ncol=r,nrow=nT)
				
		# Estimation
		#mcpar <- lapply(1:mc,mcptv,x=X,path=path,alasso=alasso,peninit=FALSE)
		mcpar <- mclapply(1:mc,mcptv,r=r,nT=nT,path=path,veps=ns,alasso=alasso,peninit=FALSE,mc.cores=ncores)
		
		#saving
		mcxp[[paste0('T=',nT,' n2s=',ns)]] <- mcpar
		#mcX[[paste0('T=',nT,' n2s=',ns)]]		<- mcpar$x
		mcpath[[paste0('T=',nT,' n2s=',ns)]] 	<- path
		
		tm <- (proc.time()-stime)[3]
		cat('\nDone with: ',paste0('T=',nT,' n2s=',ns,' in ',round(tm,1),'sec. \n'))
		cat('Average Time per iteration: ',tm/(mc*ncores),'\n')
	}
}

xpn <- 'n2s'

save(file=paste0('mcsaves/mc_X',xpn),mcX)
save(file=paste0('mcsaves/mc_path_',xpn),mcpath)
save(file=paste0('mcsaves/mc_',xpn),mcxp)
