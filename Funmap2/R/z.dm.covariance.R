## This common library for Funmap and fGWAS
## The newest source is in fGWAS library
## !!!! ALL modifies in fGWAS !!!!

##-----------------------------------------------------------
##
##
##-----------------------------------------------------------

setClass("fg.covariance.base",
	representation(
		type           = "character",    #covariance_type
		description    = "character"
  )
)

setMethod("show", signature(object="fg.covariance.base"), function(object){
   cat("     Class :", class(object), "\n");
   cat("Covar.Type :", object@type, "\n");
   
   info <- get_param_info(object, NULL );
   cat("Parameters :", info$names, "\n");
});


##-----------------------------------------------------------
## First-order Autoregressive Moving [type="ARMA(1,1)"]
##
##-----------------------------------------------------------

ARMA1_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

ARMA1_get_matrix <- function(object, par, times, options=list())
{
	rho<- par[1];
	phi<- abs(par[2]);
	s2 <- abs(par[3]);
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	sigma <-  abs(s2) *
	          ( (matrix(1, nrow=n, ncol=n) - diag(n)) * phi + diag(n)) *
	          rho^abs( matrix( 1:n, nrow=n, ncol=n, byrow=T) -  matrix(1:n, nrow=n, ncol=n, byrow = F) ) 
	return(sigma);

}

ARMA1_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("rho", "phi", "sigma2")));
}

ARMA1_check_param<-function(object, par, times, options=list())
{
	if ( par[1]>1 || par[1]<0 || par[2]>1 || par[2]<0)
		return(FALSE)
	else
		return(TRUE);

}

ARMA1_get_simu_param<-function(object, times, options=list())
{
	return( c(0.75, 0.9, 1.2) );
}

ARMA1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	s2  <- var( pheY[,1], na.rm=T );
	rho <- cor( pheY[!is.na(pheY[,3]) & !is.na(pheY[,2]),3], pheY[!is.na(pheY[,3]) & !is.na(pheY[,2]),2])
	phi <- cor( pheY[!is.na(pheY[,2]) & !is.na(pheY[,1]),2], pheY[!is.na(pheY[,2]) & !is.na(pheY[,1]),1])/rho

	return( c(rho, phi, s2*runif(1, 0.9, 1.1) ) );
}

##-----------------------------------------------------------
## S4 Class "fg.covariance.ARMA1"
##
##-----------------------------------------------------------

setClass("fg.covariance.ARMA1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature(object="fg.covariance.ARMA1"), ARMA1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.ARMA1" ), ARMA1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.ARMA1" ), ARMA1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.ARMA1" ), ARMA1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.ARMA1"), ARMA1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.ARMA1"), ARMA1_est_init_param)


##-----------------------------------------------------------
## Heterogeneous Autoregressive [type="ARH(1)"]
##
##-----------------------------------------------------------

ARH1_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

ARH1_get_matrix <- function(object, par, times, options=list())
{
	rho<- par[1];
	s2 <- abs(par[-1]);
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	sigma <-   rho^abs( matrix(1:n, nrow=n, ncol=n, byrow=T) -
		                matrix(1:n, nrow=n, ncol=n, byrow = F) )*
	           sqrt( matrix( s2, nrow=n, ncol=n, byrow = T) *
		             matrix( s2, nrow=n, ncol=n, byrow = F) )

	return(sigma);

}

ARH1_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	return(list(count=1+n, names=c("rho", paste("sigma2", 1:n, sep="_") ) ) );
}

ARH1_check_param<-function(object, par, times, options=list())
{
	if (par[1]>1 || par[1]<0)
		return(FALSE)
	else
		return(TRUE);

}

ARH1_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( c(0.75, seq(1.2, 1+n*0.2, 0.2) ) );
}

ARH1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	n <- ifelse ( is.vector(pheT), length(pheT), NCOL(pheT) );
	rho = mean(unlist(lapply( 1:(NCOL(pheY)-1), function(i){
			sel = is.na(pheY[,i]) | is.na(pheY[,i+1]);
			return( cor(pheY[!sel,i], pheY[!sel, i+1]) );
			})));
			
	s2  = colSds(pheY)^2;
	return( c(rho, s2*runif(n, 0.9, 1.1) ) );
}

##-----------------------------------------------------------
## S4 Class "fg.covariance.ARH1"
##
##-----------------------------------------------------------

setClass("fg.covariance.ARH1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature(object="fg.covariance.ARH1"), ARH1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.ARH1" ), ARH1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.ARH1" ), ARH1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.ARH1" ), ARH1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.ARH1"), ARH1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.ARH1"), ARH1_est_init_param)



##-----------------------------------------------------------
## First-order Autoregressive [type="AR(1)"]
##
##-----------------------------------------------------------

AR1_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

AR1_get_matrix <- function(object, par, times, options=list())
{
	rho<- par[1];
	s2 <- abs(par[2]);
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	sigma <-  abs(s2) * rho^abs( matrix(rep(1:n, n), nrow=n, byrow=T) -
		                         matrix(rep(1:n, n), nrow=n, byrow = F) )
	return(sigma);

}

AR1_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("rho", "sigma2")));
}

AR1_check_param<-function(object, par, times, options=list())
{
	if (par[1]>1 || par[1]<0)
		return(FALSE)
	else
		return(TRUE);

}

AR1_get_simu_param<-function(object, times, options=list())
{
	return( c(0.75, 1.2) );
}

AR1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	rho = mean(unlist(lapply( 1:(NCOL(pheY)-1), function(i){
			sel = is.na(pheY[,i]) | is.na(pheY[,i+1]);
			return( cor(pheY[!sel,i], pheY[!sel, i+1]) );
			})));

	s2  = sd(as.double(unlist(pheY)), na.rm=T)^2;
	return( c(rho, s2*runif(1, 0.9, 1.1) ) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.AR1"
##
##-----------------------------------------------------------

setClass("fg.covariance.AR1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.AR1"), AR1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.AR1" ), AR1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.AR1" ), AR1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.AR1" ), AR1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.AR1"), AR1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.AR1"), AR1_est_init_param)


##-----------------------------------------------------------
## Compound Symmetry [type="CS"]
##
##-----------------------------------------------------------

CS_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

CS_get_matrix <- function(object, par, times, options=list())
{
	rho<- par[1];
	s2 <- abs(par[2]);
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	sigma <-  abs(s2) * rho^ abs (matrix(1, nrow=n , ncol=n) -diag(n));
	return(sigma);

}

CS_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("rho","sigma2")));
}

CS_check_param<-function(object, par, times, options=list())
{
	if (par[1]>1 || par[1]<0)
		return(FALSE)
	else
		return(TRUE);

}

CS_get_simu_param<-function(object, times, options=list())
{
	return( c(0.75, 1.2) );
}

CS_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	sel <- is.na(pheY[,1]) | is.na(pheY[,2]);
	rho <- cor(pheY[!sel,1], pheY[!sel, 2]);
	s2  <- sd(as.double(unlist(pheY)), na.rm=T)^2;

	return( c(rho, s2*runif(1, 0.9, 1.1) ) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.CS"
##
##-----------------------------------------------------------

setClass("fg.covariance.CS",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.CS"), CS_show)

setMethod("get_matrix",  signature( object = "fg.covariance.CS" ), CS_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.CS" ), CS_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.CS" ), CS_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.CS"), CS_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.CS"), CS_est_init_param)


##-----------------------------------------------------------
## Heterogeneous Compound Symmetry [type="CSH"]
##
##-----------------------------------------------------------

CSH_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

CSH_get_matrix <- function(object, par, times, options=list())
{
	rho<- par[1];
	s2 <- abs(par[-1]);
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	sigma <-  rho^ abs (matrix(1, nrow=n , ncol=n) -diag(n)) *
				     sqrt( matrix(s2, nrow=n, ncol=n, byrow = T) *
		                   matrix(s2, nrow=n, ncol=n, byrow = F) )
	return(sigma);

}

CSH_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(list(count=1+n, names=c("rho", paste("sigma2", 1:n, sep="_"))));
}

CSH_check_param<-function(object, par, times, options=list())
{
	if (par[1]>1 || par[1]<0)
		return(FALSE)
	else
		return(TRUE);

}

CSH_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( c(0.75, seq(1.2, 1+n*0.2, 0.2) ) );
}

CSH_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	n <- ifelse ( is.vector(pheT), length(pheT), NCOL(pheT) );

	sel <- is.na(pheY[,1]) | is.na(pheY[,2]);
	rho <- cor(pheY[!sel,1], pheY[!sel, 2]);
	s2  <- colSds( pheY, na.rm=T )^2;
	
	return( c(rho, s2*runif(n, 0.9, 1.1) ) );
}

##-----------------------------------------------------------
## S4 Class "fg.covariance.CSH"
##
##-----------------------------------------------------------

setClass("fg.covariance.CSH",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.CSH"), CSH_show)

setMethod("get_matrix",  signature( object = "fg.covariance.CSH" ), CSH_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.CSH" ), CSH_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.CSH" ), CSH_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.CSH"), CSH_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.CSH"), CSH_est_init_param)



##-----------------------------------------------------------
## Variance Components [type="VS"]
##
##-----------------------------------------------------------

VS_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

VS_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	s2 <- abs(par);

	sigma <-  diag(s2);
	return(sigma);

}

VS_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(list(count=n, names=paste("sigma2", 1:n, sep="_")));
}

VS_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

VS_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( seq(1.2, 1.0+ n*0.2, 0.2)  );
}

VS_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return( colSds(pheY)^2 );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.VS"
##
##-----------------------------------------------------------

setClass("fg.covariance.VS",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.VS"), VS_show)

setMethod("get_matrix",  signature( object = "fg.covariance.VS" ), VS_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.VS" ), VS_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.VS" ), VS_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.VS"), VS_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.VS"), VS_est_init_param)


##-----------------------------------------------------------
## Factor Analytic _ First order  [type="FA1"]
##
##-----------------------------------------------------------

FA1_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

FA1_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	d <- abs(par[1]);
	s2 <- abs(par[-1]);

	sigma <-  d*diag(n) + sqrt( matrix(s2, nrow=n, ncol=n, byrow = T) *
		                   matrix(s2, nrow=n, ncol=n, byrow = F) )
	return(sigma);

}

FA1_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(list(count=1+n, names=c("d",paste("sigma2", 1:n, sep="_"))));
}

FA1_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

FA1_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( c(0.5, seq(1.2, 1.0+ n*0.2, 0.2) ) );
}

FA1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return(c( sd(pheY, na.rm=T)^2/100, colSds(pheY)^2 * runif(1, 0.8, 1.2)) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.FA1"
##
##-----------------------------------------------------------

setClass("fg.covariance.FA1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.FA1"), FA1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.FA1" ), FA1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.FA1" ), FA1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.FA1" ), FA1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.FA1"), FA1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.FA1"), FA1_est_init_param)



##-----------------------------------------------------------
## Hetergenous Factor Analytic - First order  [type="FAH1"]
##
##-----------------------------------------------------------

FAH1_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

FAH1_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	d <- abs(par[1:n]);
	s2 <- abs(par[-(1:n)]);

	sigma <-  diag(d) + sqrt( matrix(s2, nrow=n, ncol=n, byrow = T) *
		                   matrix(s2, nrow=n, ncol=n, byrow = F) )
	return(sigma);

}

FAH1_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(list(count=2*n, names=c(paste("d", 1:n, sep="_"), paste("sigma2", 1:n, sep="_"))));
}

FAH1_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

FAH1_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( c( runif(n, 0.5, 1.0), seq(1.2, 1.0+ n*0.2, 0.2) ) );
}

FAH1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return( c( colSds(pheY)^2*runif(1, 0.05, 0.1), colSds(pheY)^2 ) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.FAH1"
##
##-----------------------------------------------------------

setClass("fg.covariance.FAH1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.FAH1"), FAH1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.FAH1" ), FAH1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.FAH1" ), FAH1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.FAH1" ), FAH1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.FAH1"), FAH1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.FAH1"), FAH1_est_init_param)


##-----------------------------------------------------------
## Huynh-Feldt  [type="HF"]
##
##-----------------------------------------------------------

HF_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

HF_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	lambda <- abs(par[1]);
	s2 <- abs(par[-1]);

	sigma <-  ( matrix(s2, nrow=n, ncol=n, byrow = T)  + matrix(s2, nrow=n, ncol=n, byrow = F) )/2 - 
	          lambda*( matrix(1, ncol=n, nrow=n) - diag(n)  );
if(sum(sigma<0)>0) show(sigma);
	return(sigma);

}

HF_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(list(count=1+n, names=c("lambda", paste("sigma2", 1:n, sep="_"))));
}

HF_check_param<-function(object, par, times, options=list())
{
	mat <- HF_get_matrix(object, par, times, options);
	return(ifelse(any(which(mat<=0)), FALSE, TRUE));

}

HF_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( c( runif(n, 0.5, 0.8), seq(1.2, 1.0+ n*0.2, 0.2) ) );
}

HF_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return( c( min(colSds(pheY)^2)*runif(1, 1/32, 1/16), (colSds(pheY)^2/2) * runif(NCOL(pheY), 0.8, 1.2 ) ) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.HF"
##
##-----------------------------------------------------------

setClass("fg.covariance.HF",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.HF"), HF_show)

setMethod("get_matrix",  signature( object = "fg.covariance.HF" ), HF_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.HF" ), HF_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.HF" ), HF_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.HF"), HF_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.HF"), HF_est_init_param)


##-----------------------------------------------------------
## Toeplitz [ type="TOEP" ]
##
##-----------------------------------------------------------

TOEP_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

TOEP_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	s2 <- abs( par );
	mat <- matrix(1:n, nrow=n, ncol = n);
	sigma <-  matrix( s2[ abs( mat - t(mat) ) + 1 ], nrow=n, ncol=n);

	return(sigma);
}

TOEP_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( list(count=n, names=paste("sigma", 1:n, sep="_") ) );
}

TOEP_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

TOEP_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( seq(1.2, 1+n*0.2, 0.2) );
}

TOEP_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	s2 <- colSds(pheY)^2;
	return( sort(s2*runif(length(s2), 0.8, 1.2), decreasing=T) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.TOEP"
##
##-----------------------------------------------------------

setClass("fg.covariance.TOEP",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.TOEP"), TOEP_show)

setMethod("get_matrix",  signature( object = "fg.covariance.TOEP" ), TOEP_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.TOEP" ), TOEP_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.TOEP" ), TOEP_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.TOEP"), TOEP_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.TOEP"), TOEP_est_init_param)



##-----------------------------------------------------------
## Heterogeneous Toplitz [type="TOEPH" ]
##
##-----------------------------------------------------------

TOEPH_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

TOEPH_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	rho <- c(1, abs(par[1:(n-1)]) )
	s2 <- abs( par[-c(1:(n-1))] );
	mat <- matrix(1:n, nrow=n, ncol = n);
	s2.mat <- matrix(s2, nrow=n, ncol = n);
	sigma <- sqrt(s2.mat*t(s2.mat)) * matrix( rho[ abs( mat - t(mat) ) + 1 ], nrow=n, ncol=n);

	return(sigma);
}

TOEPH_get_param_info<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return( list(count=n + n-1, names=c(paste("rho", 1:(n-1), sep="_"), paste("sigma2", 1:n, sep="_")) ) );
}

TOEPH_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

TOEPH_get_simu_param<-function(object, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	c( seq(1, 1-n*0.1, -0.1)[-1],  seq(1.2, 1+n*0.2, 0.2) )
	return( c( seq(1, 1-n*0.1, -0.1)[-c(1:2)],  seq(1.2, 1+n*0.2, 0.2) ) );
}

TOEPH_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	n <- ifelse ( is.vector(pheT), length(pheT), NCOL(pheT) );
	rho = unlist( lapply( 2:NCOL(pheY), function(i){
			sel = is.na(pheY[,1]) | is.na(pheY[,i]);
			return( cor(pheY[!sel,1], pheY[!sel, i]) );
			}));

	return( c( rho, colSds(pheY)^2) );
}

##-----------------------------------------------------------
## S4 Class "fg.covariance.TOEPH"
##
##-----------------------------------------------------------

setClass("fg.covariance.TOEPH",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.TOEPH"), TOEPH_show)

setMethod("get_matrix",  signature( object = "fg.covariance.TOEPH" ), TOEPH_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.TOEPH" ), TOEPH_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.TOEPH" ), TOEPH_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.TOEPH"), TOEPH_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.TOEPH"), TOEPH_est_init_param)


##-----------------------------------------------------------
## Scaled Identity [type="SI" ]
##
##-----------------------------------------------------------

SI_show <- function(object)
{
	cat("Covariance type=", object@type, "\n");
}

SI_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );
	return(abs(par[1])*diag(n));
}

SI_get_param_info<-function(object, times, options=list())
{
	return( list(count=1, names=c("sigma2") ) );
}

SI_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

SI_get_simu_param<-function(object, times, options=list())
{
	return( runif(1, 1, 2 ) );
}

SI_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return( mean(colSds(pheY)^2) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.SI"
##
##-----------------------------------------------------------

setClass("fg.covariance.SI",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature( object = "fg.covariance.SI"), SI_show)

setMethod("get_matrix",  signature( object = "fg.covariance.SI" ), SI_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.SI" ), SI_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.SI" ), SI_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.SI"), SI_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.SI"), SI_est_init_param)



##-----------------------------------------------------------
## first-order structured antedependence  model[SAD(1) or SAD1 ]
## specified as in Jaffr´ezic et al. [2003]
##
##-----------------------------------------------------------

SAD1_show<-function(object)
{
	cat("Covariance type=", object@type, "\n");
}

SAD1_get_matrix <- function(object, par, times, options=list())
{
	n <- ifelse ( is.vector(times), length(times), NCOL(times) );

	phi<- par[1];
	v2 <- par[2];
	tmp <- (1-phi^2);

	sigma <- array(1, dim=c(n,n));
	for(i in 1:n)
	{
		sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp;
		sigma[i:n,i] <- sigma[i,i:n];
	}

	sigma <- sigma*abs(v2);
	return(sigma);
}

SAD1_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("phe", "v2")));
}

SAD1_check_param<-function(object, par, times, options=list())
{
	return(ifelse(par[1]==1, FALSE, TRUE) );

}

SAD1_get_simu_param<-function(object, times, options=list())
{
	return(c(0.75, 1.35));
}

SAD1_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	rho = mean(unlist(lapply( 1:(NCOL(pheY)-1), function(i){
			sel = is.na(pheY[,i]) | is.na(pheY[,i+1]);
			return( cor(pheY[!sel,i], pheY[!sel, i+1]) );
			})));

	s2  = sd(as.double(unlist(pheY)), na.rm=T)^2;
	return( c(rho, s2*runif(1, 0.9, 1.1) ) );
}


##-----------------------------------------------------------
## S4 Class "fg.covariance.SAD1"
##
##-----------------------------------------------------------

setClass("fg.covariance.SAD1",
	representation(
		par_num= "integer"
	), contains = "fg.covariance.base"
)

#setMethod("show", signature(object="fg.covariance.SAD1"), SAD1_show)

setMethod("get_matrix",  signature( object = "fg.covariance.SAD1" ), SAD1_get_matrix )

setMethod("get_param_info",  signature( object = "fg.covariance.SAD1" ), SAD1_get_param_info )

setMethod("check_param",  signature( object = "fg.covariance.SAD1" ), SAD1_check_param )

setMethod("get_simu_param", signature( object = "fg.covariance.SAD1"), SAD1_get_simu_param)

setMethod("est_init_param", signature( object = "fg.covariance.SAD1"), SAD1_est_init_param)



#------------------------------------------------------------------------------------------
# getCovariance
#------------------------------------------------------------------------------------------

fg.getCovariance<-function(type)
{
	obj.covar <- NULL;
	if(is.character(type))
	{
		for(obj in fg.allCovariances())
			if(toupper(obj@type)==toupper(type))
				obj.covar <- obj;
	}

	if(is.numeric(type))
		obj.covar <- fg.allCovariances()[[type]];

	return(obj.covar);
}

fg.allCovariances<-function()
{
	list.covariance <- list()
	cn<-0
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.AR1",     type = "AR1",        description = "First-order Autoregressive");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.SAD1",    type = "SAD1",       description = "First-order Structured Antedependence");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.ARMA1",   type = "ARMA(1,1)",  description = "First-order Autoregressive Moving Average");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.ARH1",    type = "ARH1",       description = "Heterogeneous Autoregressive");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.CS",      type = "CS",         description = "Compound Symmetry");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.CSH",     type = "CSH",        description = "Heterogeneous Compound Symmetry");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.VS",      type = "VS",         description = "Variance Components");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.SI",      type = "SI",         description = "Scaled Identity");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.FA1",     type = "FA1",        description = "Factor Analytic - First-order");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.FAH1",    type = "FAH1",       description = "Heterogeneous Factor Analytic - First-order");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.TOEP",    type = "TOEP",       description = "Toeplitz");
	list.covariance[[cn<-cn+1]] <- new("fg.covariance.TOEPH",   type = "TOEPH",      description = "Heterogeneous Toplitz");

	list.covariance[[cn<-cn+1]] <- new("fg.covariance.HF",      type = "HF",         description = "Huynh-Feldt");
	return(list.covariance);
}

fg_get_covar_count<-function()
{
	return(length(fg.allCovariances()));
}


fg.addCovariance<-function()
{

}