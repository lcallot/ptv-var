library(parsimonious)
library(sandwich)
library(xtable)
source("Libs.R")

# Load data
data <- GetData()
Y <- window(data$R, start=c(1954,4), end=c(2014,3))
X <- window(cbind(lag(data$R, -1), lag(data$inf, -1), lag(data$inf, -2), 
				  lag(data$y, -1), lag(data$y, -2)), start=c(1954,4), end=c(2014,3))
colnames(X) <- c("FFR_1", "INF_1", "INF_2", "OUT_1", "OUT_2")

# Standardization of the penalized covariates:
Xstand <- FALSE
if(Xstand){
	Xvars <-matrix( c(1,apply(X[,-1],2,var)),nrow=nrow(X),ncol=ncol(X),byrow = TRUE)
	X <- X/sqrt(Xvars)
}

# Lasso
w <- rbind(rep(0, 5), cbind(matrix(Inf, nrow(X)-1,1), matrix(1, nrow(X)-1, 1), matrix(1, nrow(X)-1,1), matrix(1, nrow(X)-1,2)))
w_cste <- rbind(rep(0, ncol(X)+1), cbind(matrix(1, nrow(X)-1,1), matrix(Inf, nrow(X)-1,1), matrix(1, nrow(X)-1, ncol(X)-1)))
res <- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)))


# Lasso 
hres	<- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)))
if(Xstand) hres$coef <- hres$coef/sqrt(Xvars)

#Lasso with constraint
lres	<- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=c(t(w)),ptvcste = FALSE,pmax=22)
if(Xstand) lres$coef <- lres$coef/sqrt(Xvars)


# Adaptive Lasso
aw <- 1/(abs((lres$rawcoef)))
ares <- ptvfit(Y~X, peninit=FALSE, postest=FALSE, wada=aw,ptvcste = FALSE)
if(Xstand) ares$coef <- ares$coef/sqrt(Xvars)


# Storage	

mk_co <- function(co,ptvcste,a0=0){
	if(!ptvcste) ptvco <- cbind(a0,co[,1],co[,2]/(1-co[,1]), co[,3]/(1-co[,1]),(co[,2]+co[,3])/(1-co[,1]),
								co[,4]/(1-co[,1]), co[,5]/(1-co[,1]), (co[,4]+co[,5])/(1-co[,1]))
	if(ptvcste) ptvco <- cbind(co[,1],co[,2], co[,3]/(1-co[,2]), co[,4]/(1-co[,2]),(co[,3]+co[,4])/(1-co[,2]),
							   co[,5]/(1-co[,2]), co[,6]/(1-co[,2]), (co[,5]+co[,6])/(1-co[,2]))
	
	return(ptvco)	
}

co <- ts(lres$coef, freq=4, start=c(1954,4))
ptvcoefs <- mk_co(co,FALSE,lres$intercept)
aco <- ts(ares$coef, freq=4, start=c(1954,4))
aptvcoefs <- mk_co(aco,FALSE,ares$intercept)
hco <- ts(hres$coef, freq=4, start=c(1954,4))
hptvcoefs <- mk_co(hco,FALSE,hres$intercept)

# (a)Lasso table
rc <- matrix(hres$rawcoef, ncol=5, byrow=TRUE)[-1,-1]
arc <- matrix(ares$rawcoef, ncol=5, byrow=TRUE)[-1,-1]
ae <- ptvcoefs[,c(2,3,5,6)]
aae <- aptvcoefs[,c(2,3,5,6)]

meannz <- function(x) { mean(x[x!=0]) }


# OLS
olsres <- lm(Y~X)
if(Xstand)olsres$coef <- olsres$coef/c(1,sqrt(Xvars[1,]))
olscoefs <- c(olsres$coef[2], 
			  olsres$coef[3]/(1-olsres$coef[2]),
			  olsres$coef[4]/(1-olsres$coef[2]),
			  (olsres$coef[3]+olsres$coef[4])/(1-olsres$coef[2]),
			  olsres$coef[5]/(1-olsres$coef[2]),
			  olsres$coef[6]/(1-olsres$coef[2]),
			  (olsres$coef[5]+olsres$coef[6])/(1-olsres$coef[2]),
			  olsres$coef[1]/(1-olsres$coef[2]))
olsse <- deltamethod(list(~x2,
						  ~x3/(1-x2),
						  ~x4/(1-x2),
						  ~(x3+x4)/(1-x2),
						  ~x5/(1-x2),
						  ~x6/(1-x2),
						  ~(x5+x6)/(1-x2),
						  ~x1/(1-x2)), olsres$coef, vcovHAC(olsres))


# Figure 1: Inflation
plotdata <- c()
n <- nrow(ptvcoefs)
for (i in 3:5)
{
	plotdata <- rbind(plotdata, 
					  cbind(time(co), ptvcoefs[,i], aptvcoefs[,i], hptvcoefs[,i], rep(olscoefs[i-1], n),
							rep(olscoefs[i-1]-qnorm(0.95)*olsse[i-1], n), 
							rep(olscoefs[i-1]+qnorm(0.95)*olsse[i-1], n),
							rep(i, n)))	
}
colnames(plotdata) <- c("date", "ptv", "aptv",'hptv', "ols", "olscilower", "olsciupper", "var")
plotdata <- data.frame(plotdata)

plot_title <- "Inflation Response"
ginf <- plot_par(plotdata,plot_title)
print(ginf)
ggsave('inf.pdf',ginf,height=6,width=9)

# Figure 2: Output
plotdata <- c()
n <- nrow(ptvcoefs)
for (i in 6:8)
{
	plotdata <- rbind(plotdata, 
					  cbind(time(co), ptvcoefs[,i], aptvcoefs[,i],hptvcoefs[,i], rep(olscoefs[i-1], n),
							rep(olscoefs[i-1]-qnorm(0.95)*olsse[i-1], n), 
							rep(olscoefs[i-1]+qnorm(0.95)*olsse[i-1], n),
							rep(i, n)))
}
colnames(plotdata) <- c("date", "ptv", "aptv", 'hptv' , "ols", "olscilower", "olsciupper", "var")
plotdata <- data.frame(plotdata)

plot_title <- "Output Gap Response"
gout <- plot_par(plotdata,plot_title)
print(gout)
ggsave('out.pdf',gout,height=6,width=9)


# Data figure
plotdata <- cbind(time(data$R), data$R, data$inf, data$y)
colnames(plotdata) <- c("date", "R", "inf", "y")
plotdata <- data.frame(plotdata)
plotdata <- melt(plotdata, id.vars="date")
gdata <- plot_data(plotdata)




