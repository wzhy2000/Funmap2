## This common library for Funmap and fGWAS
## The newest source is in fGWAS library
## !!!! ALL modifies in fGWAS !!!!

##-----------------------------------------------------------
##
##
##-----------------------------------------------------------

setClass("fg.curve.base",
	representation(
		type           = "character",    #curve_type
		description    = "character"
  )
)

setMethod("show", signature(object="fg.curve.base"), function(object){
   cat("     Class :", class(object), "\n");
   cat("Curve Type :", object@type, "\n");
   
   info <- get_param_info(object, NULL );
   cat("Parameters :", info$names, "\n");
   cat("   Formula :", info$formula, "\n");
});

#------------------------------------------------------------------------------------------
# Curve function
#------------------------------------------------------------------------------------------
fg.registerCurve<-function(fun, curve.type, param.name, formula.string=NULL, simu.param=NULL, est.fun=NULL)
{
	idx.userdef <- -1;
	if( is.null(.RR("userdef.curves", NULL)) )
	{
		for(obj in fg.allBuiltinCurves())
			if(toupper(obj@type)==toupper(curve.type))
				stop("! The curve type exists in the built-in curves, please rename the curve type name.");
	}
	else
	{
		userdef <- .RR("userdef.curves", NULL)
		for(i in  1:length(userdef) )
		{
			obj <- userdef[[i]]
			if(toupper(obj@type)==toupper(curve.type))
			{
				warning("! The curve type exists as a user defined curves. The content will be updated.");
				idx.userdef <- i;
				break;
			}	
		}		
	}

	if(idx.userdef==-1)
	{
		.RW("userdef.curves", list());
		idx.userdef <- 1;
	}	
		
	userdef <- .RR("userdef.curves", NULL);
	
	default.est.fun<-function(object, pheY, pheX, pheT, options=list())
	{
		r.max <- max( pheY, na.rm=T );
		r.min <- min( pheY, na.rm=T );
	
		return( runif( length(object@param.name), r.min, r.max) );
	}

	userdef[[idx.userdef]] <- new("fg.curve.userdef", description=paste(curve.type, "curve", sep=" "), 
						type = as.character(curve.type), 
						param.name = as.character(param.name), 
						fun = fun,
						formula = as.character(formula.string), 
						simu.param = as.list(simu.param), 
						est.fun = if(!is.null(est.fun)) est.fun else default.est.fun );
	
	.RW("userdef.curves", userdef );
	x <- userdef[[idx.userdef]];
	
	invisible(x);
}

fg.getCurve<-function(type)
{
	userdef <- .RR("userdef.curves", NULL);
	if(!is.null(userdef))
	{
		for(obj in userdef)
			if(toupper(obj@type)==toupper(type))
				return(obj);
	}

	obj.curve <- NULL;
	if(is.character(type))
	{
		for(obj in fg.allBuiltinCurves())
			if(toupper(obj@type)==toupper(type))
				obj.curve <- obj;
	}

	if(is.numeric(type))
		obj.curve <- fg.allBuiltinCurves()[[type]];

	return(obj.curve);
}

fg_get_curve_count<-function()
{
	return(length(fg.allBuiltinCurves()));
}

fg.allBuiltinCurves<-function()
{
	list.curve <- list()
	cn<-0
	list.curve[[cn<-cn+1]] <- new("fg.curve.log",   type = "Logistic",   description = "logistic curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.log2",  type = "Bi-Logistic",  description = "Double logistic curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.abrk",  type = "ABRK",  description = "ABRK model");

	list.curve[[cn<-cn+1]] <- new("fg.curve.pharmacology",  type = "Pharmacology",  description = "Pharmacology curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.exponential",  type = "Exponential",  description = "Exponential curve");
	list.curve[[cn<-cn+1]] <- new("fg.curve.bi.exponential",  type = "Bi-Exponential",  description = "Bi-exponential curve");

	list.curve[[cn<-cn+1]] <- new("fg.curve.power",  type = "Power",  description = "power curve");

	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre2",  type = "Legendre2",  description = "Legendre Polynomial(2nd-order)");
	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre3",  type = "Legendre3",  description = "Legendre Polynomial(3rd-order)");
	list.curve[[cn<-cn+1]] <- new("fg.curve.legendre4",  type = "Legendre4",  description = "Legendre Polynomial(4th-order)");
	
	#list.curve[[cn<-cn+1]] <- new("fg.curve.ChapmanRichard",  type = "ChapmanRichard",  description = "Chapman-Richard");

	return(list.curve);
}

##-----------------------------------------------------------
## user-define curve
##
##    unknown
##
##-----------------------------------------------------------
userdef_get_cueve <- function(object, par, times, options=list())
{
	y <- object@fun( object, par, times, options);
	return(y);
}

userdef_get_param_info<-function(object, times, options=list())
{
	return(list(count=length(object@param.name), names=object@param.name, formula=object@formula ));
}

userdef_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

userdef_get_simu_param<-function(object, times, options=list())
{
	if(length(object@simu.param)>0)
		return( rbind( object@simu.param[[1]], object@simu.param[[2]], object@simu.param[[3]]) )
	else
		return(NULL);
}

userdef_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return( object@est.fun( object, pheY, pheX, pheT, options) )
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.userdef
##-----------------------------------------------------------

setClass("fg.curve.userdef",
	representation(
	param.name = "character", 
	fun = "function",
	formula = "character",
	simu.param = "list", 
	est.fun = "function" ), 
	contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.userdef" ), userdef_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.userdef" ), userdef_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.userdef" ), userdef_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.userdef"), userdef_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.userdef"), userdef_est_init_param)


##-----------------------------------------------------------
## Logistic curve
##
##    y = a/(1+b*exp(-r*t))
##
##-----------------------------------------------------------
log_get_cueve <- function(object, par, times, options=list())
{
	y <- par[1]/(1+par[2]*exp(-1 * par[3]*times) )

	return(y);
}

log_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("a", "b", "r"), formula="y = a/(1+b*exp(-r*t))" ));
}

log_check_param<-function(object, par, times, options=list())
{
	return(TRUE);

}

log_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99), c(17.08, 9.78, 0.97), c(15.95, 9.88, 0.98)	) );
}

log_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	r.max <- max( pheY, na.rm=T );
	r.min <- min( pheY, na.rm=T );

	return( runif( 3, r.min, r.max) );
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.log
##-----------------------------------------------------------

setClass("fg.curve.log",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.log" ), log_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.log" ), log_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.log" ), log_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.log"), log_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.log"), log_est_init_param)


##-----------------------------------------------------------
## Double Logistic
## 
##    y = a1/(1+b1*exp(-r1*t)) +  a2/(1+b2*exp(-r2*t))
##
##-----------------------------------------------------------

log2_show<-function(object)
{

}

log2_get_cueve <- function(object, par, times, options=list())
{
	y <- par[1]/(1+par[2]*exp(-par[3]*times) ) + par[4]/(1+par[5]*exp(-par[6]*times) )

	return(y);
}

log2_get_param_info<-function(object, times, options=list())
{
	return(list(count=6, 
	        names=c("a1", "b1", "r1","a2", "b2", "r2"),
	        formula="y = a1/(1+b1*exp(-r1*t)) +  a2/(1+b2*exp(-r2*t))" ));
}

log2_check_param<-function(object, times, options=list())
{
	return(TRUE);

}

log2_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99), c(17.08, 9.78, 0.97), c(15.95, 9.88, 0.98)	) );
}

log2_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	r.max <- max( pheY, na.rm=T );
	r.min <- min( pheY, na.rm=T );

	return( runif( 6, r.min, r.max) );
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.log2
##-----------------------------------------------------------

setClass("fg.curve.log2",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.log2" ), log2_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.log2" ), log2_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.log2" ), log2_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.log2"), log2_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.log2"), log2_est_init_param)


##-----------------------------------------------------------
## ABRK:
##
##    y = a*(1+b*exp(-r*t))^(1/(1-k))
##
## Reference:<no>
##-----------------------------------------------------------

abrk_show<-function(object)
{

}

abrk_get_cueve <- function(object, par, times, options=list())
{
	y<- par[1]*( 1 + par[2]*exp(-par[3]*times) )^(1/(1-par[4]) );
	return(y);
}

abrk_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("a","b","r","k"),
	       formula="y = a*(1+b*exp(-r*t))^(1/(1-k))" ) );
}

abrk_check_param<-function(object, times, options=list())
{
	return(TRUE);
}

abrk_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(18.18, 9.98, 0.99), c(17.08, 9.78, 0.97), c(15.95, 9.88, 0.98)	) );
}

abrk_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	r.max <- max( pheY, na.rm=T );
	r.min <- min( pheY, na.rm=T );

	return( runif( 4, r.min, r.max) );
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.abrk
##-----------------------------------------------------------

setClass("fg.curve.abrk",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.abrk" ), abrk_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.abrk" ), abrk_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.abrk" ), abrk_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.abrk"), abrk_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.abrk"), abrk_est_init_param)


#-----------------------------------------------------------------
# Pharmacology Curve
#
#    y = E0 + Emax*t/(E50+t)
#
#-----------------------------------------------------------------

pc_show <- function(object)
{

}

pc_get_cueve <- function(object, par, times, options=list())
{
	return(  par[1] + par[3]*times/(par[2] + times) );
}

pc_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("E0", "E50", "Emax"),
	       formula="y = E0 + Emax*t/(E50+t)" ));
}

pc_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

pc_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(10.9824, 15.909, 20.7768), c( 8.9824, 16.098, 20.7768), c(6.9507, 12.090,18.5737) ) );
}

pc_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	return(c(min(pheY,na.rm=T), mean(as.matrix(pheY), na.rm=T), max(pheY,na.rm=T)));
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.pharmacology
##-----------------------------------------------------------

setClass("fg.curve.pharmacology",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.pharmacology" ), pc_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.pharmacology" ), pc_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.pharmacology" ), pc_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.pharmacology"), pc_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.pharmacology"), pc_est_init_param)


#-----------------------------------------------------------------
# Exponential  Curve
#
#    y = a*exp(r*t)
#
#-----------------------------------------------------------------

exp_show <- function(object)
{
}

exp_get_cueve <- function(object, par, times, options=list())
{
	return(  par[[1]]*exp( par[[2]]*times ) );
}

exp_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("a", "r"), formula="y = a*exp(r*t)"));
}

exp_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

exp_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c( 2, 0.0128), c( 1.8, 0.02), c(1.6, 0.024) ) );
}

exp_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	r <-  mean( (log( pheY[, 2]) - log(pheY[, 1]) )/(pheT[,2]-pheT[,1]), na.rm=T)
	w0 <- mean( pheY[, 1]/exp(r*pheT[,1]), na.rm=T);
	
	return(c(w0,r ));
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.exponential
##-----------------------------------------------------------

setClass("fg.curve.exponential",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.exponential" ), exp_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.exponential" ), exp_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.exponential" ), exp_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.exponential"), exp_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.exponential"), exp_est_init_param)



#-----------------------------------------------------------------
# Bi-exponential Curve
#
#    y = a1*exp(-r1*t) + a2*exp(-r2*t)
#
#-----------------------------------------------------------------

biexp_show <- function(object)
{

}

biexp_get_cueve <- function(object, par, times, options=list())
{
	return(  par[1]*exp(-par[2]*times) + par[3]*exp(-par[4]*times) );
}

biexp_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("a1", "r1", "a2", "r2"), 
	          formula="y = a1*exp(-r1*t) + a2*exp(-r2*t)"));
}

biexp_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

biexp_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(19.9824, 0.4699, 8.7768,  1.4699), c( 17.9824, 0.0699, 9.7768, 1.0699), c(15.9507, 0.1836, 10.5737, 1.8836) ) );
}

biexp_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	r <-  -1* mean((log( pheY[, 2]) - log(pheY[, 1]))/(pheT[,2]-pheT[,1]), na.rm=T)
	a_double <- mean(pheY[, 1]/exp(-1*r*pheT[,1]), na.rm=T);

	return(c(a_double/2, r, a_double/2, r))
}

##-----------------------------------------------------------
## S4 class: 
##           fg.curve.bi.exponential
##-----------------------------------------------------------

setClass("fg.curve.bi.exponential",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.bi.exponential" ), biexp_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.bi.exponential" ), biexp_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.bi.exponential" ), biexp_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.bi.exponential"), biexp_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.bi.exponential"), biexp_est_init_param)



##-----------------------------------------------------------------
## Power Curve
##
##         y = a*t^b
##
##-----------------------------------------------------------------
power_show <- function(object)
{

}

power_get_cueve <- function(object, par, times, options=list())
{
	return(  par[1]*( times^par[2]) );
}

power_get_param_info<-function(object, times, options=list())
{
	return(list(count=2, names=c("a","b"), formula="y = a*t^b"));
}

power_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

power_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(simu_a = 11.049, simu_b = 1.151), c(simu_a = 9.049, simu_b = 1.251), c(simu_a = 7.148,  simu_b = 1.359) ) );
}

power_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	b <- mean( log(pheY[, NCOL(pheY)])-log(pheY[, 1]) / log(pheT[, NCOL(pheT)])-log(pheT[, 1]), na.rm=T);
	a <- exp( mean( log(pheY[, NCOL(pheY)]) - b*log(pheT[, NCOL(pheY)]), na.rm=T ))
	return(c(a, b));
}


##-----------------------------------------------------------
## S4 class: 
##           fg.curve.power
##-----------------------------------------------------------

setClass("fg.curve.power",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.power" ), power_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.power" ), power_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.power" ), power_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.power"), power_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.power"), power_est_init_param)


##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(2nd-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1)
##
##-----------------------------------------------------------------

Legendre2_show <- function(object)
{
}

Legendre2_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times - options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3]  );
}

Legendre2_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("u0","u1","u2"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1)"));
}

Legendre2_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre2_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:3),drop=F ] );
}

Legendre2_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	ti <- -1 + 2*(pheT-options$min.time)/( options$max.time - options$min.time );
    y1 <- pheY[,1];
    y2 <- pheY[,2];
    y3 <- pheY[,3];
    
	u2 <- mean( ( (y1-y2)/(ti[,1]-ti[,2]) - (y1-y3)/(ti[,1]-ti[,3]) ) / ( (ti[,1]+ti[,2])-(ti[,1]+ti[,3]) ) / 1.5, na.rm = T);
	u1 <- mean( (y1-y2) /(ti[,1]-ti[,2]) - u2 * 1.5 * (ti[,1]+ti[,2]), na.rm=T);
	u0 <- mean( y3 + u2*0.5 - u1*ti[,3] - u2 * 1.5 * ti[,3]^2, na.rm=T);
	
	return(c(u0, u1, u2));
}


##-----------------------------------------------------------
## S4 class: 
##           fg.curve.legendre2
##-----------------------------------------------------------

setClass("fg.curve.legendre2",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre2" ), Legendre2_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.legendre2" ), Legendre2_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre2" ), Legendre2_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre2"), Legendre2_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre2"), Legendre2_est_init_param)



##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(3rd-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t)
##
##-----------------------------------------------------------------

Legendre3_show <- function(object)
{
}

Legendre3_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3] +  0.5*(5*ti^3-3*ti)*par[4] );
}

Legendre3_get_param_info<-function(object, times, options=list())
{
	return(list(count=4, names=c("u0","u1","u2","u3"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) "));
}

Legendre3_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre3_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:4),drop=F] );
}

Legendre3_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	u3 <- Legendre2_est_init_param(object, pheY, pheX, pheT, options);
	return(c(u3, 0.0001));
}


##-----------------------------------------------------------
## S4 class: 
##           fg.curve.legendre3
##-----------------------------------------------------------

setClass("fg.curve.legendre3",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre3" ), Legendre3_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.legendre3" ), Legendre3_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre3" ), Legendre3_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre3"), Legendre3_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre3"), Legendre3_est_init_param)



##-----------------------------------------------------------------
## Legendre:  Legendre Polynomial(4th-order)
##
##        y = u0 + u1 *t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3)
##
##-----------------------------------------------------------------

Legendre4_show <- function(object)
{
}

Legendre4_get_cueve <- function(object, par, times, options=list())
{
	ti <- -1 + 2*(times-options$min.time)/( options$max.time - options$min.time );
	return( par[1] + ti*par[2] + 0.5*(3*ti*ti-1)* par[3] +  0.5*(5*ti^3-3*ti)*par[4] +  0.125*(35*ti^4-30*ti^2+3)* par[5] );
}

Legendre4_get_param_info<-function(object, times, options=list())
{
	return(list(count=5, names=c("u0","u1","u2","u3", "u4"), formula="y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3)"));
}

Legendre4_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

Legendre4_get_simu_param<-function(object, times, options=list())
{
	QQ2 = c(simu_u0 = 11.049, simu_u1 = 1.551, simu_u2 = -8.019, simu_u3 = 3.151, simu_u4 = 0.652, simu_u5 = -0.597, simu_u6 = 0.821);
	Qq1 = c( simu_u0 = 9.049, simu_u1 = 1.151, simu_u2 = -6.019, simu_u3 = 2.651, simu_u4 = 0.652, simu_u5 = -0.797, simu_u6 = 0.621);
	qq0 = c( simu_u0 = 7.148, simu_u1 = 1.379, simu_u2 = -4.489, simu_u3 = 2.004, simu_u4 = 0.662, simu_u5 = -0.836, simu_u6 = 0.432)

	return( rbind( QQ2, Qq1, qq0)[, c(1:5),drop=F] );
}

Legendre4_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	u3 <- Legendre2_est_init_param(object, pheY, pheX, pheT, options);
	return(c(u3, 0.0001, 0.0001));
}


##-----------------------------------------------------------
## S4 class: 
##           fg.curve.legendre4
##-----------------------------------------------------------

setClass("fg.curve.legendre4",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.legendre4" ), Legendre4_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.legendre4" ), Legendre4_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.legendre4" ), Legendre4_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.legendre4"), Legendre4_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.legendre4"), Legendre4_est_init_param)


##-----------------------------------------------------------
## Chapman-Richard
## 
##    y = a*(1-exp(-rt))^b
##
##-----------------------------------------------------------

cr_show <- function(object)
{
}

cr_get_cueve <- function(object, par, times, options=list())
{
	return ( par[1]*(1 - exp(-par[3]*times))^par[2] );
}

cr_get_param_info<-function(object, times, options=list())
{
	return(list(count=3, names=c("a","b","r"), formula="y = a*(1-exp(-rt))^b"));
}

cr_check_param<-function(object, par, times, options=list())
{
	return(TRUE);
}

cr_get_simu_param<-function(object, times, options=list())
{
	return( rbind( c(21.98, 0.47, 9.78), c(19.98, 0.47, 8.77), c(15.95, 0.48, 7.58)	) );
}

cr_est_init_param<-function(object, pheY, pheX, pheT, options=list())
{
	a <- max(pheY, na.rm=T);
	b <- 1
	r <- mean( log(1-pheY/a)/pheT/(-1), trim=0.2, na.rm=T)
	if( is.infinite(r)) r<-0;
	
	return(c(a, b, r));
}



##-----------------------------------------------------------
## S4 class: 
##           fg.curve.ChapmanRichard
##-----------------------------------------------------------

setClass("fg.curve.ChapmanRichard",
	representation(
	), contains = "fg.curve.base"
)

setMethod("get_curve",  signature( object = "fg.curve.ChapmanRichard" ), cr_get_cueve )

setMethod("get_param_info",  signature( object = "fg.curve.ChapmanRichard" ), cr_get_param_info )

setMethod("check_param",  signature( object = "fg.curve.ChapmanRichard" ), cr_check_param )

setMethod("get_simu_param", signature( object = "fg.curve.ChapmanRichard"), cr_get_simu_param)

setMethod("est_init_param", signature( object = "fg.curve.ChapmanRichard"), cr_est_init_param)


