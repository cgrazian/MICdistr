##############################################################################
##############################################################################
### Censored Gaussian Mixture Models
### by Clara Grazian
### 
### Code to reproduce the analysis of the paper
### "Estimating MIC distributions and cutoffs through mixture models: 
### an application to establish \textit{M. Tuberculosis} resistance"
##############################################################################
##############################################################################

# 1) Data loading

dst_data <- read.table("~/Dropbox/2018_2019_Oxford/CRyPTIC/PNAS/mixtures/CRyPTIC_MIC.txt",
					sep="\t",header=T,stringsAsFactor=F)
str(dst_data)

table(dst_data$key_id)
table(dst_data$drug)

# 2) Estimating a censored GMM with a reversible jump algorithm
#	 (Richardson & Green, JRSSB, 59, 731-792 (1997)).

library(coda); library(MCMCpack); library(lattice)

# Hyperparameter

alpha=1    #cluster weights Dirichlet(alpha,...,alpha)

# Likelihood function

log.lkd.da <- function(y,w,mu,sigma,S)
{
	n <- length(y)
	K <- length(w)
	tot <- 0
	for(i in 1:n){
		for(k in 1:K)
		{
			tot <- tot + (S[i]==k) * ( log(w[k]) + log(dnorm(y[i],mu[k],sigma[k]) ) )
		}		
	}
	return(tot)
}

# Prior distribution

log.prior<-function(w,mu,sigma) {                    
  mdp=dpois(nc<-length(w),lambda,log=TRUE)           
  sp=sum(dgamma(sigma,shape=1.5,rate=0.5,log=TRUE))  
  wp=log(ddirichlet(w,rep(alpha,length(w))))         
  mup=sum(dnorm(mu,mean=20,sd=10,log=TRUE))         
  return(sum(mup+sp+wp+mdp))
}

# Relabeling to deal with the label switching

change.label <- function(x,old.lab,new.lab)
{
	G <- length(old.lab)
	x1 <- letters[x]
	for(g in 1:G)
	{
		x1[x1==letters[g]] <- new.lab[g]
	}	
	x1 <- as.numeric(x1)
	return(x1)
}

censored_GMM <- function(y,lambda,K=10^6,SS=100,alpha=1)
{
	#	y		:	MIC data in log2(MIC)
	#	K		:	number of simulations
	#	SS		:	thin
	
	n=length(y)
		
	#Initialise MCMC with one component
	mu=mean(y)
	sigma=sd(y)
	w=1
	Sdraw <- rep(1,n)
	G <- 1
	nh <- n

	u 		<- c()
	ystar 	<- rep(Inf,n)
	while(sum(abs(ystar)==Inf)!=0)
	{
		for(l in 1:n)
		{
			treshb   =  pnorm(y[l]+1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
			tresha   =  pnorm(y[l]-1, mu[Sdraw[l]],sigma[Sdraw[l]]) 

			if(tresha != treshb)
			{
				u[l] = runif(1,tresha , treshb)	
				ystar[l] <- qnorm(u[l],mu[Sdraw[l]], sigma[Sdraw[l]])
			} else {
				ystar[l] <- y[l]
			}

		}	# end of the observations loop
	}
	
	oll=log.lkd.da(ystar,w,mu,sigma)
	olp=log.prior(w,mu,sigma)

	#MCMC
	Nsamp=K/SS;  
	PL=matrix(NA,Nsamp,3); colnames(PL)<-c("d","olp","oll"); TH=list();

	for (k in 1:K) {

		### Step 1: Update w, mu, sigma
			
		OK=TRUE
		d=length(w)
		move=sample(1:5,1) #choose a MCMC move (2 birth death, 3 fixed dimension)		
		
		if (move==1) {
			i=sample(1:d,1)		#add a cluster
			wn=runif(1,min=0,max=w[i]); wp=c(w,wn); wp[i]=w[i]-wn
			mun=runif(1,min(y)-2,max(y)+2); mup=c(mu,mun)
			sn=rgamma(1,shape=1.5,rate=0.5); sigmap=c(sigma,sn)
			qd=dunif(mun,min(y)-2,max(y)+2,log=TRUE)+dgamma(sn,shape=1.5,rate=0.5,log=TRUE)-
						log(w[i])
			qn=0
			rhod=-log(d)
			rhon=-log((d+1)*d)
		}

		if (move==2) {		# kill a cluster
			if (d>1) {
				ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
				wp=w; wp[j]=w[j]+w[i]; wp=wp[-i]
				mup=mu[-i]; sigmap=sigma[-i]
				qd=0
				qn=dunif(mun,min(y)-2,max(y)+2,log=TRUE)+
							dgamma(sigma[i],shape=1.5,rate=0.5,log=TRUE)-log(w[i]+w[j])
				rhod=-log(d*(d-1))
				rhon=-log(d-1)
			} else {
				OK=FALSE
			}
		}

		if (move==3) {	# fixed dimension: move mu
			i=sample(1:d,1); 
			mup=mu; mup[i]=rnorm(1,mean=mu[i],sd=1);
			sigmap=sigma; wp=w;
			qd=qn=rhod=rhon=0      
		}

		if (move==4) {		# fixed dimension: move sigma
			i=sample(1:d,1); 
			sigmap=sigma; u=runif(1,0.5,2); sigmap[i]=sigma[i]*u
			mup=mu; wp=w;
			qn=-log(u); qd=rhod=rhon=0      
		}

		if (move==5) {		# fixed dimension: redistribute a weight
			if (d>1) {
				ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
				wp=w; wp[i]=runif(1,min=0,max=w[i]+w[j]); wp[j]=w[i]+w[j]-wp[i]
				mup=mu; sigmap=sigma;
				qn=qd=rhod=rhon=0
			} else {
				OK=FALSE
			}
		}

		if (OK) {
			nll=log.lkd.da(ystar,wp,mup,sigmap)
			nlp=log.prior(wp,mup,sigmap)  
			MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
			if (log(runif(1))<MHR) {
				w=wp; mu=mup; sigma=sigmap
				oll=nll; olp=nlp
			}
		}
 
		### Step 2: Multinomial sampler

		G <- length(w)
		Sdraw <- c()
		for(l in 1:n)
		{
			probs <- c()
			for(g in 1:G)
			{
				probs[g] <-  dnorm(y[l], mu[g], sigma[g])
			}
			probs[probs<0] <- 0
			Sdraw[l] <- sample(1:G,size=1,prob=probs)
		}

 		### Step 3: Data augmentation

		u 		<- c()
		ystar 	<- rep(Inf,n)
		while(sum(abs(ystar)==Inf)!=0)
		{
			for(l in 1:n)
			{
				treshb   =  pnorm(y[l]+1, mu[Sdraw[l]],sigma[Sdraw[l]]) 
				tresha   =  pnorm(y[l], mu[Sdraw[l]],sigma[Sdraw[l]]) 
				if(tresha != treshb)
				{
					u[l] = runif(1,tresha , treshb)
					ystar[l] <- qnorm(u[l],mu[Sdraw[l]], sigma[Sdraw[l]])
				} else {
					ystar[l] <- y[l]
				}
			}	# end of the observations loop
		}
			
		# Relabeling
	
		ordering <- order(mu)
		w <- w[ordering]
		mu <- mu[ordering]
		sigma <- sigma[ordering]
		old.lab <- sort(unique(Sdraw))
		Sdraw <- change.label(Sdraw,old.lab,ordering)

		G <- length(w)
		  
		# Saving
		if (k%%SS==0) {
			TH[[k/SS]]=list(w=w,mu=mu,sigma=sigma,Sall=Sdraw)
			PL[k/SS,]=c(length(w),olp,oll)
		}
		
	}
	
	return(list(TH=TH,PL=PL))
	
}	# end of function

drugs <- unique(dst_data$drug)

# Select the drug: an example with Amikacin
S <- 1
drug <- drugs[S]	

#Load the data
y=dst_data$miclog2[dst_data$drug==drug]
y=y[is.na(y)==F]

#Hyperparameters
lambda=length(unique(dst_data$miclog2[dst_data$drug==drug]))  
	
GS_simulation <- censored_GMM(y=y,lambda=lambda)


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

mix.dens <- function(x,theta)
{
	k <- length(theta)/3
	ww <- theta[1:k]
	mu <- theta[(k+1):(2*k)]
	sigma <- theta[(2*k+1):(3*k)]
	dens <- c()
	for(i in 1:length(x)){
		dens[i] <- sum( ww * dnorm(x[i], mu , sigma) ) 		
	}
	return(dens)
}

dmax <- getmode(GS_simulation$PL[,1])

wmat <- matrix(NA,ncol=dmax,nrow=Nsamp)
mumat <- matrix(NA,ncol=dmax,nrow=Nsamp)
sigmamat <- matrix(NA,ncol=dmax,nrow=Nsamp)

Nsamp <- 10^6/100  

for(k in 1:Nsamp)
{
	if(length(GS_simulation$TH[k][[1]]$w)==dmax)
	{
		wmat[k,] <- GS_simulation$TH[k][[1]]$w	
		mumat[k,] <- GS_simulation$TH[k][[1]]$mu
		sigmamat[k,] <- GS_simulation$TH[k][[1]]$sigma	
	}
}

wmean <- apply(wmat,2,mean,na.rm=T)
mumean <- apply(mumat,2,mean,na.rm=T)
sigmamean <- apply(sigmamat,2,mean,na.rm=T)

xseq <- seq(min(y),max(y),0.01)
breakpoint_dens <- c()
for(i in 1:length(xseq))
{
	breakpoint_dens[i] <- ( wmean[1]*dnorm(xseq[i],mumean[1],sigma[1]) +
							wmean[2]*dnorm(xseq[i],mumean[2],sigma[2]) )/ 	 	
							mix.dens(xseq[i],c(wmean,mumean,sigmamean))	
	breakpoint_dens[i] <- wmean[1]*dnorm(xseq[i],mumean[1],sigma[1]) / 	
							mix.dens(xseq[i],c(wmean,mumean,sigmamean))	
}
	
max(xseq[breakpoint_dens>0.5])
