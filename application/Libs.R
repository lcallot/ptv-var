GetData <- function()
{
	if (!exists("FREDData.RData"))
	{
		require(quantmod)
		require(xts)
		effr <- apply.quarterly(getSymbols("FEDFUNDS", src="FRED", auto.assign=FALSE), FUN="mean")
		R <- ts(100*log(1+effr/100), freq=4, start=c(1954,3)) # OK!
		gdpdef <- ts(getSymbols("GDPDEF", src="FRED", auto.assign=FALSE), freq=4, start=c(1947,1))
		inf <- ts(400*diff(log(gdpdef)), freq=4, start=c(1947,1)) # OK!
		gdp <- ts(getSymbols("GDPC1", src="FRED", auto.assign=FALSE), freq=4, start=c(1947,1))
		y <- ts(DeTrend1(log(gdp), 0.000625)[,1]*100, freq=4, start=c(1947,1))
		save(R, inf, y, file="FREDData.RData")
	}
	else load("FREDData.RData")
	return(list(R=R, inf=inf, y=y)) 
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

mf_labeller <- function(var, value){
	if (all(value == c(1,2,3))) return(expression(alpha[1], alpha[2], alpha[1]+alpha[2]))
	if (all(value == c(2,3,4))) return(expression(alpha[1], alpha[2], alpha[1]+alpha[2]))
	if (all(value == c(3,4,5))) return(expression(alpha[1], alpha[2], alpha[1]+alpha[2]))
	if (all(value == c(4,5,6))) return(expression(beta[1], beta[2], beta[1]+beta[2]))
	if (all(value == c(5,6,7))) return(expression(beta[1], beta[2], beta[1]+beta[2]))
	if (all(value == c(6,7,8))) return(expression(beta[1], beta[2], beta[1]+beta[2]))
	if (all(value==c("R", "inf", "y"))) return(expression(R[t], pi[t], y[t]))
}

DeTrend1 <- function(y, q)
{
#/* -- Procedure for producing one-sided detrended versions of a 
#      series, using a white-noise + I(2) model.  This produces a 
#      very smooth version of the trend component.  The two-sided 
#      version of the model produces the HP filter.
#
#      A nice description of what this does is given in 
#      Harvey and Jaeger, JAE, July-Sep. 1993, pp. 231-248
#
#      Inputs:
#      y = series to be detrended
#      q = relative variance of I(2) component
#          Note: HP quarterly uses q=.000625 (Kydland Prescott)
#                for monthly data a value of q=.00000075
#                matches quarterly gain at 50%, 80% and 90% periods
#*/
#	decl i, vague, f, x, p, xf, h, e, k, xc;
	
	vague <- 1e+4; # // Vague Prior on I(2) component //
	
#	// -- Initialize System Matrices -- //
	f     <- matrix(0, 2, 2)#  = zeros(2,2);
	f[1,1] <- 2;#f[0][0] = 2;
	f[1,2] <- -1;# f[0][1] = -1;
	f[2,1] <- 1;#f[1][0] = 1;
	x   <- matrix(0, 2, 1)    #= zeros(2,1);
	p   <- matrix(vague, 2, 2);#    = vague*ones(2,2);
	p[1,1] <- vague+q#[0][0] = vague+q;
	xf      <- matrix(NA, length(y), 1)#zeros(rows(y),1)*.NaN;
	
	for (i in 1:length(y))
	{
		x <- f%*%x
		p <- f%*%p%*%t(f)
		p[1,1] <- p[1,1] + q
		if (!is.na(y[i]))
		{
			h <- p[1,1] + 1
			e <- y[i]-x[1]
			k <- p[,1]/h
			x <- x+k*e
			p <- p - (k%*%p[1,,drop=FALSE])
		}
		xf[i] <- x[1]
	}
	xc <- y-xf
	return(cbind(xc,xf));
}


plot_par <- function(plotdata,plot_title){
	
	thin <- plotdata$date[seq(1,length(plotdata$date),5)]

	mplt <- melt(plotdata,id.vars = c(1,8),measure.vars = c(2,4,5))
	
	g <- ggplot(mplt, aes(x=date)) + 
		facet_grid(var~., scales="free", labeller = mf_labeller) +
		geom_rect(aes(xmin=1957.5, xmax=1958.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2007.75, xmax=2009.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_ribbon(data=plotdata,aes(ymin=olscilower, ymax=olsciupper,x=date), fill="gray90", alpha=0.5)+
		geom_line(data=subset(mplt,variable%in%c('ols','ptv','hptv')),aes(y=value,colour=variable,shape=variable)) + 
		geom_point(data=subset(mplt,(date%in%thin)&(variable%in%c('ptv'))),aes(y=value,colour=variable,shape=variable)) + 
		theme_minimal() +
		theme(panel.grid.major = element_line(colour = "grey50"), legend.box = "horizontal",
			  legend.position="bottom", plot.title = element_text(vjust=2, size=12)) +
		scale_x_continuous(breaks = seq(round(min(plotdata$date)), round(max(plotdata$date)), by = 5)) +
		scale_colour_manual("Estimator", breaks=c("ols", "ptv",'hptv'),
							labels=c("OLS",'Lasso (max 16 breaks)', "Lasso"),
							values=c( "firebrick",'navyblue',"grey25" )) +
		scale_shape_manual("Estimator", breaks=c("ols", "ptv",'hptv'),
							labels=c("OLS",'Lasso (max 16 breaks)', "Lasso"),
							values=c(16,46,46)) +	
		ylab("") + xlab("") + ggtitle(plot_title)
	return(g)
}

plot_data <- function(plotdata){
	g <- ggplot(plotdata, aes(date)) + 
		facet_grid(variable~., scales="free", labeller = mf_labeller) +
		geom_rect(aes(xmin=1957.5, xmax=1958.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2007.75, xmax=2009.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_line(aes(y=value)) + 
		theme_minimal() +
		theme(panel.grid.major = element_line(colour = "grey50"), legend.box = "horizontal",
			  legend.position="bottom", plot.title = element_text(vjust=2, size=12)) +
		scale_x_continuous(breaks = seq(1955, 2010, by = 5))+
		ylab("") + xlab("") + ggtitle("Interest Rate, Inflation, and Output Gap.")
	
	return(g)
}
