

# Constructs that stat table after a ptv estimation
ptv.xptab <- function(mcxp,statn,en,mc)
{
	xptbl <- array(0,dim=c(length(mcxp),length(statn),1+length(en))
				   ,dimnames=list('Experiment'=names(mcxp),'Stat'=statn,'Estimator'=c(en,'DGP')))
	
	# Looping over the experiments
	for(xp in names(mcxp)){
		xpl <- mcxp[[xp]]
		
		# storage init
		x <- rep(0,length(en))
		names(x) <- en
		nbrk <- eerr <- perr <- lmbd <- fneg <- fpos <- tpos <- rmse <- x
		
		nzada <- 0
		dgpbrk <- 0
		
		for(m in xpl){# Looping over the iterations
			xpX <- m[['x']]	
			dgpbrk <- dgpbrk + sum((m[['DGP']][-1] - m[['DGP']][-length(m[['DGP']])])!=0)
			for(n in en){ # Looping over the estimators
				if(is.null(m[[n]])) m[[n]] <- matrix(0*m[['DGP']])
				if(n=='aLasso') nzada <- nzada + 1*(m[['lambda']][n]!=0)
				go <- TRUE
				if(n=='aLasso') go <- m[['lambda']][n]!=0
				nbrk[n] <- nbrk[n] + sum((m[[n]][-1]-m[[n]][-length(m[[n]])])!=0) # nbr breaks
				fpos[n] <- fpos[n] + sum(((m[[n]][-1]-m[[n]][-length(m[[n]])])!=0) * ((m[['DGP']][-1] - m[['DGP']][-length(m[['DGP']])])==0)) # nbr of false positive
				tpos[n] <- tpos[n] + sum(((m[[n]][-1]-m[[n]][-length(m[[n]])])!=0) * ((m[['DGP']][-1] - m[['DGP']][-length(m[['DGP']])])!=0)) # nbr of false positive
					fneg[n] <- fneg[n] + sum(((m[[n]][-1]-m[[n]][-length(m[[n]])])==0) * ((m[['DGP']][-1] - m[['DGP']][-length(m[['DGP']])])!=0)) # nbr of false negative
				if(go){
					eerr[n] <- eerr[n] + sqrt(mean(abs(m[[n]]-m[['DGP']])))# l1 Estimation error
					perr[n] <- perr[n] + sqrt(mean( (m[[n]]*xpX - m[['DGP']]*xpX)^2) )# l2 Prediction error
					rmse[n] <- rmse[n] + m[['rmse']][n] # rmse! 
					lmbd[n] <- lmbd[n] + m[['lambda']][n] # Penalty parameter
				}
			}		
		}
		xptbl[xp,1,'DGP'] <- dgpbrk/mc
		
		sc <- c(mc,nzada,ifelse(length(en)==4,c(mc,mc),mc))
		
		xptbl[xp,1,en] <- nbrk/mc
		xptbl[xp,2,en] <- fpos/mc
		xptbl[xp,3,en] <- tpos/mc
		xptbl[xp,4,en] <- fneg/mc
		xptbl[xp,5,en] <- eerr/sc
		xptbl[xp,6,en] <- perr/sc
		xptbl[xp,7,en] <- rmse/sc
		xptbl[xp,8,en] <- lmbd/sc
	}
	
	return(xptbl)
}



# Estimator wrapper function to be called by mclapply / other parallel routines
mcptv <- function(m,r,nT,path=NULL,alphaT=NULL,veps=1,vpar=1,alasso=TRUE,peninit=FALSE,conslasso=FALSE,ar=FALSE,pmax=NULL)
{	
	# Some checks	
	if(is.null(alphaT)&&is.null(path))stop('Either a path XOR alphaT in (0,1) must be given, not nothing')
	if((!is.null(alphaT))&&(!is.null(path)))stop('Either a path XOR alphaT in (0,1) must be given, not both')
	
	if(is.null(path)){
		if((!is.numeric(alphaT))||(alphaT<0)||(alphaT>1)) stop('alpha must be a numeric between 0 and 1.')
		if((!is.numeric(vpar))||(vpar<0)) stop('the variance of the parameters vpar must be numeric and be >0. ')
		
		path <- matrix(sqrt(vpar)*rnorm(r*nT),ncol=r)
		censr <- matrix(runif(r*nT),ncol=r)
		path <- (censr<alphaT)*path
		path <- apply(path,2,cumsum)
	}
	if(!is.null(path)) { if(sum(c(nT,r)-dim(path))!=0) stop('The dimension of x must be the same as that of path.')}
	
	
	# DGP setting, variance fixed to 1 since there are enough other things to vary.
	if(!ar){
		# Generate the covariates
		x <- matrix(rnorm(r*nT),ncol=r,nrow=nT)
		colnames(x) <- paste('x',1:r,sep='')
		# Generate the innovations
		eps <- sqrt(veps) * rnorm(nT)
		# Dep variable
		Y <- (x*path)%*%matrix(1,nrow=r,ncol=1) + eps
	}
	if(ar){
		#init 
		x <- Y <- matrix(0,nrow=nT,ncol=1)	
		#making some noise and initial value		
		y0 <- x[1] <- rnorm(1)
		eps  <- sqrt(veps) * rnorm(nT)
		# generating the ar process
		for(i in 1:nT){
			Y[i] <- path[i] * x[i] + eps[i]
			if(i<nT)x[i+1] <- Y[i]
		}
	}
	
	# Estimation
	ptv <- list('Lasso'=NULL,'aLasso'=NULL,'DGP'=NULL,'Post'=NULL,'x'=x)
	# Computing the Lasso
	las		<- ptvfit(Y~x,peninit=peninit,pmax=pmax)
	
	# Adaptive Lasso
	if((sum(las$rawcoef[-c(1:r)]!=0)>0) && alasso){
		alas 	<- ptvfit(Y~x,wada=1/abs(las$rawcoef),peninit=peninit)
		ptv[['aLasso']] <- alas$coefficients
	}
	else ptv[['aLasso']] <- las$coefficients
	
	# Conservative Lasso
	if(conslasso){
		clas 	<- ptvfit(Y~x,wada=1/abs(las$rawcoef),peninit=peninit,conslas=TRUE)
		ptv[['cLasso']] <- clas$coefficients
	}
	else ptv[['cLasso']] <- las$coefficients
		
	# Storing the parameters
	ptv[['Lasso']]	<- las$coefficients
	ptv[['Post']]   <- las$postcoef
	ptv[['DGP']]	<- path
	
	# Storing the rmse
	rmse <- c(sqrt(mean(residuals(las)^2)),
			  ifelse((sum(las$rawcoef[-c(1:r)]!=0)>0) && alasso,sqrt(mean(residuals(alas)^2)),0),
			  ifelse(conslasso,sqrt(mean(residuals(clas)^2)),0),
			  ifelse(!is.null(las$postres),sqrt(mean((las$postres)^2)),0))
	names(rmse) <- c('Lasso','aLasso','cLasso','Post')
	ptv[['rmse']] <- rmse
	
	
	#Storing the penalty
	lambda <- c(las$lambda,
				ifelse((sum(las$rawcoef[-c(1:r)]!=0)>0)&&alasso,alas$lambda,0),
				ifelse(conslasso,clas$lambda,0))
	names(lambda) <- c('Lasso','aLasso','cLasso')
	ptv[['lambda']] <- lambda
	
	return(ptv)
}

mcptv.plot <- function(mcpar,npath=10,alasso){
	
	mcpar <- mcpar
	
	# Initialization	
	r<-ncol(mcpar[[1]]$DGP)	
	vn <- paste0('Variable ',1:r)

	# final prep of the data
	mparam$iter <- factor(mparam$iter)
	colnames(mparam)[1] <- 'Time'
	
	# plotting
	gpar <- ggplot(mparam,aes(x=Time,y=value,colour=Estimator,shape=Estimator))   + scale_colour_manual(values=c('DGP'='grey','adaptive Lasso'='steelblue','Lasso'='firebrick'))
	
	gpar <- gpar + geom_line(data=subset(mparam,Estimator!='DGP'),alpha=0.7,size=0.3,aes(group=interaction(iter,Estimator))) 
	gpar <- gpar + geom_point(data=subset(mparam,Estimator!='DGP'),size=0.7,aes(group=interaction(iter,Estimator)))
	
	gpar <- gpar + geom_line(data=subset(mparam,Estimator=='DGP'),size=1,aes(group=interaction(iter,Estimator)))
	gpar <- gpar + geom_point(data=subset(mparam,Estimator=='DGP'),size=1,aes(group=interaction(iter,Estimator)))
	
	gpar <- gpar + facet_wrap(~Var2,ncol=1) + theme_minimal()
	return(gpar)
}



mcplt <- function(mxp,ncol=1,post=FALSE){
	
	# index for the dots:
	thin <- seq(1,100,5)
	thin2 <- seq(3,100,5)
	
	gpar <- ggplot(subset(mxp,Estimator%in%c('Lasso','aLasso','DGP')),aes(x=Time,y=value,shape=Estimator,colour=Estimator))  
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='DGP'),size=0.5,aes(group=interaction(iter,Estimator)))
	gpar <- gpar + geom_point(data=subset(mxp,(Estimator=='DGP')&(Time%in%thin)),aes(group=interaction(iter,Estimator),shape='DGP'))
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='aLasso'),alpha=0.7,size=0.4,aes(group=interaction(iter,Estimator))) 
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='Lasso'),alpha=0.7,size=0.4,aes(group=interaction(iter,Estimator))) 
	gpar <- gpar + geom_point(data=subset(mxp,(Estimator=='Lasso')&(Time%in%thin2)),aes(group=interaction(iter,Estimator),shape='Lasso'))
	if(post)gpar <- gpar + geom_line(data=subset(mxp,Estimator=='Post'),alpha=0.7,size=0.4,aes(group=interaction(iter,Estimator))) 

	
	gpar <- gpar + facet_wrap(~Experiment,scale='free',ncol=ncol) +	scale_y_continuous(name=expression(beta))
	gpar <- gpar + theme_minimal() +
		theme(legend.text=element_text(size=8),legend.box = "horizontal",legend.position="bottom") +
		scale_colour_manual("Estimator", labels=c( "Lasso", "Adaptive Lasso",'DGP'),
							breaks=c("Lasso", "aLasso",'DGP'),
							values=c( "firebrick",'black','navyblue' )) +
		scale_shape_manual("Estimator", labels=c("Lasso", "Adaptive Lasso",'DGP'),
						   breaks=c( "Lasso", "aLasso",'DGP'),
						   values=c( 46,3,16 )) +
		ylab("") + xlab("") 
	
	return(gpar)
}

# trying to get fancy facet labels with facet_wrap
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
	#works with R 3.0.1 and ggplot2 0.9.3.1
	require(gridExtra)
	
	g <- ggplotGrob(gg.plot)
	gg <- g$grobs      
	strips <- grep("strip_t", names(gg))
	
	for(ii in seq_along(labels))  {
		modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
						   grep=TRUE, global=TRUE)
		gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[[ii]])
	}
	
	g$grobs <- gg
	class(g) = c("arrange", "ggplot",class(g)) 
	g
}


mcpltgrid <- function(mxp,ar=FALSE){

	thin <- seq(1,100,5)
	thin2 <- seq(3,100,5)	
	
	gpar <- ggplot(subset(mxp,Estimator%in%c('Lasso','aLasso','cLasso','DGP')),aes(x=Time,y=value,shape=Estimator,colour=Estimator))
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='DGP'),size=0.5,aes(group=interaction(iter,Estimator)))
	gpar <- gpar + geom_point(data=subset(mxp,(Estimator=='DGP')&(Time%in%thin)),aes(group=interaction(iter,Estimator),shape='DGP'))
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='aLasso'),alpha=0.7,size=0.4,aes(group=interaction(iter,Estimator))) 
	gpar <- gpar + geom_line(data=subset(mxp,Estimator=='Lasso'),alpha=0.7,size=0.4,aes(group=interaction(iter,Estimator))) 
	gpar <- gpar + geom_point(data=subset(mxp,(Estimator=='Lasso')&(Time%in%thin2)),aes(group=interaction(iter,Estimator),shape='Lasso'))

	if(!ar)gpar <- gpar + facet_wrap(alpha~sigma,scales='free') +	scale_y_continuous(name=expression(beta)) + geom_blank() 
	if(ar) gpar <- gpar + facet_wrap(theta1~theta2,scales='free') +	scale_y_continuous(name=expression(beta)) + geom_blank() 
	gpar <- gpar + theme_minimal() +
		theme(legend.text=element_text(size=8),legend.box = "horizontal",legend.position="bottom") +
		scale_colour_manual("Estimator", labels=c( "Lasso", "Adaptive Lasso",'DGP'),
							breaks=c("Lasso", "aLasso",'DGP'),
							values=c( "firebrick",'black','navyblue' )) +
		scale_shape_manual("Estimator", labels=c("Lasso", "Adaptive Lasso",'DGP'),
						   breaks=c( "Lasso", "aLasso",'DGP'),
						   values=c( 46,3,16 )) +
		ylab("") + xlab("") 	
	labexp <- NULL
	
	# Hard coded labels (boooo) for the facets of the prw xp
	if(!ar){
		av <- c(0.01,0.031,0.1)
		nv <- c(0.1,1,10)
		for(i in 1:length(unique(mxp$alpha))){ 
			for(j in 1:length(unique(mxp$sigma))){ 
				labexp <- c(labexp,as.expression(bquote(Var(eta)==.(nv[j])~alpha[T]==.(av[i])))[[1]])
			}
		}	
	}
	if(ar){
		# Hard coded labels (boooo) for the facets of the prw AR
		aar <- c(0,0.5,0.9)
		for(i in 1:length(unique(mxp$gamma1))){ 
			for(j in 1:length(unique(mxp$gamma2))){ 
				labexp <- c(labexp,as.expression(bquote(gamma^0==.(aar[j+1])~gamma^1==.(aar[i])))[[1]])
			}
		}
	}
	gpar <- facet_wrap_labeller(gpar,labels = labexp)
	return(gpar)
}

catmctab <- function(tbl,statn,clabel,tlabel,nbrk,clabel2=NULL){
	
	cnt <- NULL 
	cat('\\begin{table}\n')
	cat('\\tiny\n')
	cat('\\begin{tabular}{','p{2.5cm}',rep('l',1),rep('r',ncol(tbl)), '}\n',sep=' ')
	cat('\\toprule\n')
	
	#column label 
	cat(' & ')
	#for(n in alphaT){cat('& & \\multicolumn{1}{c}{$\\alpha_T=',n,'$}')}
	for(o in clabel){cat('&',o)}
	cat('\\\\\n')
	if(!is.null(clabel2)){
		cat(' & ')
		for(o in clabel2){cat('&',o)}
	    cat('\\\\\n')
	}
	
	for(n in statn){
		cat('\\cmidrule{3-',ncol(tbl)+2,'}\n',sep='')
		rind <- which(n == statn)
		
		if(rind==1)
		{
			cat('\\multirow{3}{*}{',n,'} & ')
			cat('DGP & ')
			cat(round(nbrk,3),sep=' & ')
			cat('\\\\\n')
			#Lasso		
			cat(' & Lasso & ')
			cat(round(tbl[3*(rind-1)+1,],3),sep=' & ')
			cat('\\\\\n')
		}
		
		if(rind>1)
		{
			cat('\\multirow{2}{*}{',n,'} & ')
			#Lasso
			cat('  Lasso & ')
			cat(round(tbl[3*(rind-1)+1,],3),sep=' & ')
			cat('\\\\\n')
		}
		
		#adaptive Lasso
		cat('& aLasso & ')
		cat(round(tbl[3*(rind-1)+2,],3),sep=' & ')
		cat('\\\\\n')
	
		if(!(rind%in%c(1,2,3,4,8))){
			#post Lasso
			cat('& Post & ')
			cat(round(tbl[3*rind,],3),sep=' & ')
			cat('\\\\\n')
		}
	#if(i%in%mdru.ind)cat('\\midrule\n')
	
	}
	cat('\n')
	cat('\\bottomrule\n')			
	cat('\\end{tabular}')
	cat('\n')
	cat('\\label{tab:',tlabel,'}')
	cat('\n')
	cat(paste('\\caption{Constant parameter, varying sample size: ',mc,' iterations.}'))
	cat('\n')
	cat('\\end{table}')
}

