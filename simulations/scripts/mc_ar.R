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

# sample length
nT <- 100
pmax <- nT/2#sqrt(nT)
mcxp   <- list()

count  <- c(1,2,3)
for(a in count){
	for(b in count[-3]){
	
		aar <- c(0,0.5,0.9)
		bar <- c(0.5,0.9)
		
		path <- matrix(bar[b],ncol=r,nrow=nT)
		path[1:(0.5*nT),] <- aar[a]
	
		stime <- proc.time()
			
		# Estimation
		#mcpar <- lapply(1:mc,mcptv,r=r,nT=nT,path=path,veps=0.1,alasso=alasso,peninit=FALSE,ar=TRUE)
		mcpar <- mclapply(1:mc,mcptv,r=r,nT=nT,path=path,veps=0.1,alasso=alasso,peninit=FALSE,ar=TRUE,pmax=pmax,mc.cores=ncores)
		mcxp[[paste0('gamma^0==',aar[a],' gamma^1==',bar[b])]] <- mcpar
		
		
		tm <- (proc.time()-stime)[3]
		cat('\nDone with: ',paste0('aar=',aar[a],', bar=',bar[b],' in ',round(tm,1),'sec. \n'))
		cat('Average Time per iteration: ',tm/(mc*ncores),'\n')
	}
}

xpn <- 'ar'

#save(file=paste0('mcsaves/mc_path_',xpn),mcpath)
save(file=paste0('mcsaves/mc_',xpn),mcxp)
