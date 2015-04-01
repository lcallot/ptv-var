library(parsimonious)
library(sandwich)
library(xtable)


source("Libs.R")

# Load data
data <- GetData()
Y <- window(data$R, start=c(1954,4), end=c(2014,3))
X <- window(cbind(lag(data$R, -1), lag(data$inf,-1), lag(data$y,-1)), start=c(1954,4), end=c(2014,3))
colnames(X) <- c("FFR_1", "INF_1", "OUT_1")

# Standardization of the penalized covariates:
Xstand <- FALSE
if(Xstand){
	Xvars <-matrix( c(1,apply(X[,-1],2,var)),nrow=nrow(X),ncol=ncol(X),byrow = TRUE)
	X <- X/sqrt(Xvars)
}

# Lasso
w_rho <- rbind(rep(0, ncol(X)), cbind(matrix(1, nrow(X)-1,1), matrix(1, nrow(X)-1, ncol(X)-1)))
w <- rbind(rep(0, ncol(X)), cbind(matrix(Inf, nrow(X)-1,1), matrix(1, nrow(X)-1, ncol(X)-1)))
w_cste_rho <- rbind(c(0,rep(0, ncol(X))), cbind(matrix(1, nrow(X)-1,1), matrix(1, nrow(X)-1,1), matrix(1, nrow(X)-1, ncol(X)-1)))
w_cste <- rbind(c(0,rep(0, ncol(X))), cbind(matrix(1, nrow(X)-1,1), matrix(Inf, nrow(X)-1,1), matrix(1, nrow(X)-1, ncol(X)-1)))
res <- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)))

# Lasso with ptv rho
ptv_rho	<- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w_rho)),ptvcste=FALSE)
if(Xstand) hres$coef <- hres$coef/sqrt(Xvars)
ptv_rho_pmax <- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w_rho)),ptvcste=FALSE,pmax=20)
if(Xstand) hres$coef <- hres$coef/sqrt(Xvars)

# Lasso with fixed rho
ptv	<- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)),ptvcste=FALSE)
if(Xstand) lres$coef <- lres$coef/sqrt(Xvars)
ptv_pmax <- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)),ptvcste=FALSE,pmax=20)
if(Xstand) ares$coef <- ares$coef/sqrt(Xvars)


mk_coefs <- function(co,a0){cbind(a0,co[,1], co[,2]/(1-co[,1]), co[,3]/(1-co[,1]))} 
mk_coefs_ptv_cste <- function(co){cbind(co[,1],co[,2], co[,3]/(1-co[,2]), co[,4]/(1-co[,2]))}

# Storag	
co <- ts(ptv$coef, freq=4, start=c(1954,4))
coefs1 <- mk_coefs(co,ptv$intercept)
aco <- ts(ptv_pmax$coef, freq=4, start=c(1954,4))
coefs2 <- mk_coefs(aco,ptv_pmax$intercept)
hco <- ts(ptv_rho$coef, freq=4, start=c(1954,4))
coefs3 <- cbind(ptv_rho$intercept,hco)# mk_coefs(hco,hres$intercept)
hco <- ts(ptv_rho_pmax$coef, freq=4, start=c(1954,4))
coefs4 <- cbind(ptv_rho_pmax$intercept,hco)# mk_coefs(hco,hres$intercept)


# OLS
olsres <- lm(Y~X)
if(Xstand)olsres$coef <- olsres$coef/c(1,sqrt(Xvars[1,]))
olscoefs <- c(olsres$coef[2], 
			  olsres$coef[3]/(1-olsres$coef[2]),
			  olsres$coef[4]/(1-olsres$coef[2]),
			  olsres$coef[1]/(1-olsres$coef[2]))
olsse <- deltamethod(list(~x2,
						  ~x3/(1-x2),
						  ~x4/(1-x2),
						  ~x1/(1-x2)), olsres$coef, vcovHAC(olsres))


plotdata <- c()
n <- nrow(ptvcoefs)
for (i in c(2,3,4))
{
	plotdata <- rbind(plotdata, 
					  cbind(time(co), coefs1[,i], coefs2[,i], coefs3[,i], coefs4[,i],
					  	  ifelse(i==1,0,rep(olscoefs[i-1], n)),
					  	  ifelse(i==1,0,rep(olscoefs[i-1]-qnorm(0.95)*olsse[i-1], n)), 
					  	  ifelse(i==1,0,rep(olscoefs[i-1]+qnorm(0.95)*olsse[i-1], n)),
						  rep(ifelse(i==1,1,i-1), n)))	
}
colnames(plotdata) <- c("date", "ptv1", "ptv2",'ptv3','ptv4', "ols", "olscilower", "olsciupper", "var")
plotdata <- data.frame(plotdata)
#plotdata$date <- as.double(plotdata$date)


lab_func <- function(var, value)c(expression(rho, alpha, beta))

#mpltd <- melt(plotdata[,-c(4],id.vars = c(1,))
thin <- plotdata$date[seq(1,length(plotdata$date),5)]
mplt <- melt(plotdata,id.vars=c(1,9),measure.vars = c(2,3,6)) 

g <- ggplot(mplt, aes(date)) + 
		facet_grid(var~., scales="free",labeller = lab_func) +
		theme_minimal() +
		theme(panel.grid.major = element_line(colour = "grey50"), legend.box = "horizontal",
			  legend.position="bottom", plot.title = element_text(vjust=2, size=12)) +
		scale_x_continuous(breaks = seq(round(min(plotdata$date)), round(max(plotdata$date)), by = 5)) +
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
		geom_line(data=subset(mplt,variable%in%c('ols','ptv1','ptv2')),aes(y=value,colour=variable,shape=variable)) + 
		geom_point(data=subset(mplt,(date%in%thin)&(variable%in%c('ptv2'))),aes(y=value,colour=variable,shape=variable)) + 
		scale_colour_manual("Estimator", breaks=c("ols", "ptv1", "ptv2"),
							labels=c("OLS", "Lasso", "Lasso (max 16 breaks)"),
							values=c( 'navyblue',"firebrick","grey25" )) +
		scale_shape_manual("Estimator", breaks=c("ols", "ptv1", "ptv2"),
						   labels=c("OLS", "Lasso", "Lasso (max 16 breaks)"),
						   values=c( 46,16,46 )) +
		ylab("") + xlab("") 

print(g)
ggsave(filename = '1lag.pdf',plot=g,height=6,width=10)


lab_func_rho <- function(var, value)c(expression(rho[t], (1-rho[t])*alpha[t], (1-rho[t])*beta[t]))

mplt2 <- melt(plotdata,id.vars=c(1,9),measure.vars = c(4,5)) 
g_rho <- ggplot(mplt2, aes(date)) + 
		facet_grid(var~., scales="free",labeller = lab_func_rho) +
		geom_rect(aes(xmin=1957.5, xmax=1958.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_rect(aes(xmin=2007.75, xmax=2009.25, ymin=-Inf, ymax=+Inf), fill='gray90', alpha=0.5) +
		geom_line(data=subset(mplt2,variable%in%c('ptv3','ptv4')),aes(y=value,colour=variable,shape=variable)) + 
		geom_point(data=subset(mplt2,(date%in%thin)&(variable%in%c('ptv4'))),aes(y=value,colour=variable,shape=variable)) + 
		theme_minimal() +
		theme(panel.grid.major = element_line(colour = "grey50"), legend.box = "horizontal",
			  legend.position="bottom", plot.title = element_text(vjust=2, size=12)) +
		scale_x_continuous(breaks = seq(round(min(plotdata$date)), round(max(plotdata$date)), by = 5)) +
			scale_colour_manual("Estimator", breaks=c( "ptv3", "ptv4"),
							labels=c( "Lasso", "Lasso (max 16 breaks)"),
							values=c( 'navyblue',"firebrick" )) +
		scale_shape_manual("Estimator", breaks=c( "ptv3", "ptv4"),
						   labels=c( "Lasso", "Lasso (max 16 breaks)"),
						   values=c( 46,16 )) +
		ylab("") + xlab("") 
print(g_rho)
ggsave(filename = 'ptv_rho.pdf',plot = g_rho,height=6,width=10)
