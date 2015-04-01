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
# residual standard deviation
veps <- 0.1

# sample length
nT <- 100

#storage
mcxp   <- list()
mcpath <- list()

count  <- c(1,2,3)
for(a in count){
	for(b in count){
		
		stime <- proc.time()

		# Parameter generation:
		# Break location
		if(a==1){	loc <- c(0.1,0.5,0.9)
					path <- matrix(1,ncol=r,nrow=nT)
				  	path[1:(loc[b]*nT),] <- 0}
		# break size
		if(a==2) {	sz <- c(0.1,1,10) 
					path <- sz[b]*matrix(rep(c(0,1),each=nT/2),ncol=r,nrow=nT)}
		# nbr breaks
		if(a==3){	if(b==1)path <- matrix(rep(c(0,1,0,1),each=nT/4),ncol=r,nrow=nT)
					if(b==2)path <- matrix(rep(c(0,1,0,1,0,1,0,1,0,1),each=nT/10),ncol=r,nrow=nT)
					if(b==3)path <- matrix(rep(c(0,1,0,0,0,0,0,0,1,0),each=nT/10),ncol=r,nrow=nT)}
				
		# Estimation
		#mcpar <- lapply(1:mc,mcptv,r=r,nT=nT,path,alasso=alasso,peninit=FALSE,conslas=TRUE)
		mcpar <- mclapply(1:mc,mcptv,r=r,nT=nT,path=path,veps=veps,alasso=alasso,peninit=FALSE,conslas=TRUE,mc.cores=ncores)
		
		#saving
		mcxp[[paste0('a=',a,' b=',b)]] <- mcpar
		mcpath[[paste0('a=',a,' b=',b)]] 	<- path
		
		tm <- (proc.time()-stime)[3]
		cat('\nDone with: ',paste0('a=',a,' b=',b,' in ',round(tm,1),'sec. \n'))
		cat('Average Time per iteration: ',tm/(mc*ncores),'\n')
	}
}

xpn <- 'brk'

save(file=paste0('mcsaves/mc_path_',xpn),mcpath)
save(file=paste0('mcsaves/mc_',xpn),mcxp)
