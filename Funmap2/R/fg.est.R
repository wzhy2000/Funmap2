fg_dat_est<-function( obj.phe, curve.type="auto", covariance.type="auto", file.plot.pdf=NULL, options=list() )
{
	if( is.null( obj.phe$obj.curve ) )
	{
		cat("  No curve is specified, curve fitting is being performed.\n");
		r.test <- fg_fit_curve( obj.phe$pheY, obj.phe$pheX, obj.phe$pheT, curve.type=curve.type, file.plot.pdf=file.plot.pdf );

		if(r.test$error)
			stop("? No curve is fitted to this data.\n");

		obj.phe$obj.curve <- fg.getCurve( r.test$type );
		obj.phe$summary.curve <- r.test;
	}

	cat("  Curve Type==> ", obj.phe$obj.curve@type, " <==\n");

	r <- proc_est_curve( obj.phe$pheY, obj.phe$pheX, obj.phe$pheT, obj.phe$obj.curve )
	if( r$error )
		return(list(error=T, err.info="Can not estimate the parameter of mean vector according to the curve function" ) )

	parX.len <- NCOL(obj.phe$pheX)
	if(is.null(obj.phe$pheX)) parX.len <- 0;

	cat(" Parameter range estimation ...... \n")
	range <- proc_est_curve_range(obj.phe$pheY, obj.phe$pheX, obj.phe$pheT, obj.phe$obj.curve, par.init = r$par);

	obj.phe$est.curve <- list( type = obj.phe$obj.curve@type,
							   param = r$par[-(1:(1+parX.len))],
							   param.lower = range$lower[-(1:(1+parX.len))],
							   param.upper = range$upper[-(1:(1+parX.len))],
							   parx  = r$par[1:(1+parX.len)],
							   parX.lower = range$lower[(1:(1+parX.len))],
							   parX.upper = range$upper[(1:(1+parX.len))] );

	if( is.null(obj.phe$obj.covar) || ( toupper(obj.phe$obj.covar@type) != toupper(covariance.type) ) )
	{
		r.test <- fg_fit_covar( obj.phe$pheY, obj.phe$pheX, obj.phe$pheT, r$y.resd, obj.phe$obj.curve, covariance.type );

		if(r.test$error)
			stop("? No covariance is fitted to this data.\n")
		else
		{
			obj.phe$summary.covar <- r.test;
			obj.phe$obj.covar <- fg.getCovariance( r.test$type );
		}
	}

	cat("  Covariance Type==> ", obj.phe$obj.covar@type, " <==\n");

	r.est <- proc_est_covar( r$y.resd, NULL, obj.phe$pheT, obj.phe$obj.curve, obj.phe$obj.covar );
	if ( r.est$error )
		return(list(error=T, err.info="Can not estimate the parameter of mean vector according to the curve function" ) );

	obj.phe$est.covar<- list( type = obj.phe$obj.covar@type, param = r.est$par);
	obj.phe$error <- F;

	return(obj.phe);
}

fn_get_resd<-function(pheY, pheX, pheT, obj.curve, parin  )
{
	par_X <- c( parin[1] );
	if ( !is.null(pheX) )
		par_X <- c( par_X, parin[ 2:(NCOL(pheX)+1)]);

	par_c <- parin[ -1 ];
	if ( !is.null(pheX) )
		par_c <- parin[ -(1:(NCOL(pheX)+1)) ]

	mu_gen <- get_curve( obj.curve, par_c, pheT, options=list(max.time=max(pheT, na.rm=T), min.time=min(pheT, na.rm=T)) )
	if(all(is.na( mu_gen )))
		return(NULL);

	if( is.vector( mu_gen ) )
		y_resd <- t(t(pheY) - mu_gen) -  matrix( rep( as.matrix( cbind(rep(1, NROW(pheY)), pheX)) %*% par_X, NCOL(pheY)), ncol=NCOL(pheY), byrow=F)
	else
		y_resd <- pheY - mu_gen  - matrix( rep( as.matrix( cbind(rep(1, NROW(pheY)), pheX)) %*% par_X, NCOL(pheY)), ncol=NCOL(pheY), byrow=F)

	return( y_resd );
}

proc_est_curve<-function(  pheY, pheX, pheT, obj.curve, par.init=NULL, options=list(n.loop=10)  )
{
	options$max.time <- max(pheT, na.rm=T)
	options$min.time <- min(pheT, na.rm=T)

	get_init_curve_par<-function( pheY, pheX, pheT, f.obj )
	{
		par.curve <- est_init_param( f.obj, pheY, pheX, pheT, options=options );
		pheX.len  <- 1 + NCOL(pheX)
		if( is.null(pheX) )  pheX.len <- 1;;
		par.X <- rep( mean(pheY, na.rm=T), pheX.len );

		return( c(par.X, par.curve)  );
	}

	get_rand_curve_par<-function( pheY, pheX, pheT, f.obj, parin )
	{
		uc <- check_fittness( pheY, pheX, pheT, f.obj, parin )

		if(uc$fit)
		{
			return( parin*runif(length(parin), 0.9, 1.1));
		}
		else
		{
			if ( uc$over0.05 < NCOL(pheY)/4 )
				return( parin * runif(length(parin), 0.5, 1.5))
			else
			{
				par.curve <- get_init_curve_par( pheY, pheX, pheT, f.obj );
				#pheX.len  <- 1 + NCOL(pheX)
				#if( is.null(pheX) )  pheX.len <- 1;;
				#par <- c( rep( mean(pheY, na.rm=T), pheX.len ),  par.curve );
				return( par.curve * runif( length(par.curve) , 0.5, 1.5 ) );
			}
		}
	}

	# parin : par.X[1..n], par.curve.
	fn_mle_est<-function( parin, extra_par )
	{
		y_resd <- fn_get_resd( extra_par$pheY, extra_par$pheX, extra_par$pheT, extra_par$obj.curve, parin );
		if(is.null(y_resd))
			return(NaN);

		A <- sum( y_resd^2, na.rm=T );

		return(A);
	}

	if( is.null( par.init) )
	{
		par.init <-  c( 0 );
		if( !is.null(pheX)) par.init <- c( par.init, mean(colMeans(pheY, na.rm=T), na.rm=T)/colMeans(pheX, na.rm=T) );
		par.init <- c( par.init, est_init_param( obj.curve, pheY, pheX, pheT, options ) );
	}

	if(.RR("debug")) cat("par.init", par.init, "\n")

	h0 <- proc_mle_loop(  pheY, pheX, pheT,
			obj.curve,
			fn_mle_est,
			mle_extra_par = list( pheY=pheY, pheX=pheX, pheT=pheT, obj.curve=obj.curve ),
			parin = par.init,
			fn_init_par = get_init_curve_par,
			fn_rand_par = get_rand_curve_par,
			options=options )

	if(.RR("debug")) cat( "MU[F]", h0$value, h0$par, "\n");

	r.check <- check_fittness( pheY, pheX, pheT, obj.curve, h0$par );
	if( !r.check$fit || !is.finite( h0$value ) )
		return(list(error=T, par=NA, val=NA))
	else
	{
		y.resd <- fn_get_resd( pheY, pheX, pheT,  obj.curve, h0$par )
		## The sum of the squared differences between each observation and its predicted value.
		y.SSE <- sum(y.resd^2,na.rm=T) ;

		##Gives the average of the squares of the errors of each value.
		y.MSE <- y.SSE/length(which(!is.na(pheY)));

		##The square root of the MSE that estimates the standard deviation of the random error.
		y.RMSE <- sqrt(y.MSE)
		## pramater count
		K <- get_param_info(obj.curve, pheT)$count;
		n.sample <- NROW(pheY);

		AIC  <- 2*K +n.sample*log(y.SSE/n.sample)
		AICc <- log(y.SSE/n.sample) + (n.sample+K)/(n.sample-K-2)
		BIC  <- n.sample*log(y.SSE/n.sample) + K*log(n.sample)

		pheY_hat <- matrix(NA, nrow=NROW(pheY), ncol=NCOL(pheY) );
		pheY_hat_vec <- c(pheY_hat);

		pheT_vec <- c(pheT);
		pheY_vec <- c(pheY);
		for (ti in unique(pheT_vec) )
		{
			if (!is.na(ti))
				pheY_hat_vec[which(pheT_vec==ti)] <- mean(pheY_vec[which(pheT_vec==ti)], na.rm=T)
		}

		pheY_hat <- matrix(pheY_hat_vec, nrow=NROW(pheY), ncol=NCOL(pheY) );
		R2   <- sum(((pheY-y.resd)-pheY_hat)^2, na.rm=T)/sum((pheY-pheY_hat)^2, na.rm=T)

		return(list(error=F, par.count = K, AIC = AIC, AICc = AICc, BIC = BIC,
			              SSE = y.SSE, MSE = y.MSE, RMSE = y.RMSE, R2 = R2,
			              par=h0$par, y.resd=y.resd))
	}
}

proc_est_curve_range<-function( pheY, pheX, pheT, f.curve, par.init )
{
	n.obs <- NROW( pheY );
	mu.pars <- c();

	loop <- 0;
	while( loop < .RR("mu.range.loop")  )
	{
		y.sub <- sample(n.obs)[1:round( n.obs * runif(1,0.5,0.9) )]

		pheX0 <- NULL;
		pheY0 <- pheY[y.sub,,drop=F];
		if( !is.null(pheX) ) pheX0 <- pheX[y.sub,,drop=F];
		if( !is.null(dim(pheT)) )  pheT0 <- pheT[y.sub,,drop=F] else pheT0 <- pheT;

		r <- proc_est_curve( pheY0, pheX0, pheT0, f.curve, par.init, options=list(n.loop=2) );
		if (r$error) next;

		mu.pars <- rbind(mu.pars, r$par);
		loop <- loop + 1;
	}

	mu.lower <-c();
	mu.upper <-c();

	for(i in 1:NCOL(mu.pars) )
	{
		mu.lower <- c(mu.lower, min(mu.pars[,i]) )
		mu.upper <- c(mu.upper, max(mu.pars[,i]) )
	}

	return(list(lower=mu.lower, upper=mu.upper))
}


proc_est_covar<-function( Y.resd, pheX, pheT, obj.curve, obj.covar, par.init=NULL, options=list(n.loop=10)  )
{
	get_init_covar_par<-function( Y.resd, pheX, pheT, f.covar)
	{
		return( est_init_param( f.covar, Y.resd, pheX, pheT, options=options ) );
	}

	get_rand_covar_par<-function( Y.resd, pheX, pheT, f.covar, parin )
	{
		return( parin* runif(length(parin), 0.9, 1.1) );
	}

	#parin:
	# phi1, s1, phi2, s2
	fn_mle_est<-function( parin, extra_par)
	{
		y.resd  <- extra_par$Y.resd;
		pheX    <- extra_par$pheX;
		pheT    <- extra_par$pheT;
		f.covar <- extra_par$f.covar;

		cov.mat  <- get_matrix(f.covar, parin, pheT);
		if(any(is.na(cov.mat))) return(NaN);

		pv <- dmvnorm_fast( y.resd, rep(0, NCOL(y.resd)), cov.mat, log=T);
		if(any(is.infinite(pv)))
			return(NaN);

		A <- sum( pv );
		return( -A );
	}

	if( is.null( par.init) )
		par.init <- get_init_covar_par(  Y.resd, pheX, pheT, obj.covar );

	h0 <- proc_mle_loop( Y.resd, pheX, pheT,
			obj.covar,
			fn_mle_est,
			mle_extra_par=list(Y.resd=Y.resd, pheX=pheX, pheT=pheT, f.covar=obj.covar),
			parin = par.init,
			fn_init_par = get_init_covar_par,
			fn_rand_par = get_rand_covar_par,
			options=options )

	if(.RR("debug"))  cat( "COV(F)", h0$value, h0$par, "\n");

	if( is.finite( h0$value ) )
	{
		K <- get_param_info(obj.covar, pheT)$count;
		###Note: L=-log(L)
		AIC  <- 2*K + 2*h0$value;
		BIC  <- 2*h0$value + K*log(NROW(Y.resd))
		return(list(error=F, AIC=AIC, BIC=BIC, par=h0$par, logL = -h0$value ))
	}
	else
		return(list(error=T, par=NA, err.info="Failed to estimate covariance structure."));
}

proc_mle_loop<-function(  pheY, pheX, pheT, f.obj, fn_mle, mle_extra_par, parin, fn_init_par, fn_rand_par, options=list(n.loop=10)  )
{
	h0<-list( value=Inf, par=parin );
	control <- list(maxit = 50000, reltol=1e-8);

	init.loop <- 1
	while( is.infinite(h0$value) )
	{
		try( h0 <- optim( parin, fn_mle, extra_par=mle_extra_par,
				method = ifelse(runif(1)>0.75, "Nelder-Mead", "BFGS"),
				control= control ), .RR("try.silent", FALSE)  );

		if ( class(h0)=="try-error" || is.na(h0) || is.infinite(h0$value) || is.null(h0$convergence) || (h0$convergence !=0) )
		{
			if ( !is.na(h0) && is.list(h0) && !is.null(h0$convergence) )
			{
				if (h0$convergence == 1)
				{
					control$maxit <- control$maxit*2;
					if ( control$maxit > 500*4096 )
					{
						control$maxit <- 500*4096;
						control$reltol <- control$reltol*2;
					}
				}
				else
					if(.RR("debug")) cat("h0$convergence =", h0$convergence , "\n");
			}

			reset_seed();
			parin <- fn_init_par( pheY, pheX, pheT, f.obj );

			init.loop <- init.loop + 1;
			if(init.loop > 100) return( list(error=T, par=NA, value=NA) );

			next;
		}
	}

	if(.RR("debug")) cat( "X[0]", h0$value, h0$par, "\n");

	parin0 <- h0$par;
	loop <- 1;
	unfit <- 0;
	h2.better <-NULL;
	loop.optim <-1;

	min.fit<-Inf;

	while ( loop < options$n.loop )
	{
		h2 <- NA;
		parinx<- fn_rand_par( pheY, pheX, pheT, f.obj, parin0 );

		loop.optim <- loop.optim+1;
		try( h2 <- optim( parinx, fn_mle, extra_par=mle_extra_par,
				method = ifelse(loop.optim%%2==1, "Nelder-Mead", "BFGS"), control=control ), .RR("try.silent", FALSE)  );

		if (class(h2)=="try-error" || any(is.na(h2)) || (h2$convergence !=0)  )
		{
			if ( is.list(h2) && h2$convergence == 1)
			{
				control$maxit <- control$maxit*2;
				if ( control$maxit > 500*4096 )
				{
					control$maxit <- 500*4096;
					control$reltol <- control$reltol*2;
				}
			}

			if(loop.optim>100)
			{
				loop.optim <- 0;
				loop <- loop+1;

				if (is.infinite(h0$value) && n.loop - loop< 2 )
					n.loop <- n.loop + 1;
			}

			reset_seed();
			next;
		}
		else
		{
			uc <- check_fittness( pheY, pheX, pheT, f.obj, h2$par )
			if ( !uc$fit )
			{
				if (uc$over0.05<min.fit)
				{
					parin0    <- h2$par;
					min.fit <- uc$over0.05;
				}

				reset_seed();
				next;
			}

			if(.RR("debug")) cat( "MU(L", loop, ")", uc$over0.05, "/", h2$value, h2$par, "\n");

			if ( h2$value < h0$value && h2$value>0 ) h0 <- h2;
			if ( h2$value >0 && h0$value<0 ) h0 <- h2;

			parin0 <- h0$par;
			min.fit <- Inf;

			if (is.infinite(h0$value) && n.loop - loop< 2 )
				n.loop2 <- n.loop2 + 1;
		}

		loop.optim <- 0;
		loop <- loop+1;
	}

	if(.RR("debug")) cat( "X[F]", h0$value, h0$par, "\n");

	uc <- check_fittness( pheY, pheX, pheT, f.obj, h0$par );
	if( uc$fit && is.finite(h0$value) )
		return(list(error=F, par=h0$par, value=h0$value))
	else
		return(list(error=T, par=NA, value=NA));
}

check_fittness<-function( pheY, pheX, pheT, f.obj, parin )
{
	if( !inherits(f.obj, "fg.curve.base") )
		return(list(fit=TRUE));

	y_resd <- fn_get_resd( pheY, pheX, pheT, f.obj, parin );
    if( is.null(y_resd) )
   		return(list(fit=TRUE));

if(NCOL(y_resd)==1) browser();

	y_sd   <- apply( pheY, 2, function(x) sd(x) )

	py.prob  <- pnorm( colMeans(y_resd), mean=0, sd=y_sd, lower.tail=F );
	py.05 <- length( which(py.prob<0.05) );

	return(list(fit = ifelse( py.05 > 1, F, T), over0.05=py.05) );
}

fg_fit_curve<-function( pheY, pheX, pheT, curve.type="auto", file.plot.pdf=NULL )
{
	cat(" Searching curve type ......\n" );

	obj.curves <- list();
	if(curve.type=="auto" || curve.type=="" || is.null(curve.type) )
	{
		for ( i in 1:fg_get_curve_count()) obj.curves[[i]] <- fg.getCurve( i );
	}
	else
	{
		for(i in 1:length(curve.type))
			obj.curves[[i]] <- fg.getCurve(curve.type[i]);
	}

	if(!is.null(file.plot.pdf)) pdf(file.plot.pdf);

	est.curve <- list();
	index = 0
	for ( obj.curve in obj.curves)
	{
		cat( "* [", index<-index+1, "/", length(obj.curves), "] try curve: ", obj.curve@type, "\n" );
		r <- proc_est_curve( pheY, pheX, pheT, obj.curve )
		if( r$error )
			warning("Can not estimate the parameters of mean vector according to the curve function[ curve.type=", obj.curve$type, "]" )
		else
		{
			r$type <- obj.curve@type;
			est.curve[[index]] <- r;
		}

	    if(!is.null(file.plot.pdf))
	    {
            plot(1,1, type="n", xlim=c(min(pheT, na.rm=T), max(pheT, na.rm=T)), ylim=c(min(pheY, na.rm=T), max(pheY, na.rm=T)), main=r$type, xlab="Time", ylab="Phenotype");
            for(i in 1:NROW(pheY))
            	lines(pheT[i,], pheY[i,], col="gray", lwd=0.2);

            pheT_vec <- c(pheT)
            pheY_vec <- c(pheY)
	        ti <- unique(c(pheT));
	        ti <- sort( ti [ !is.na(ti) ] );
	        y.mu <- rep(NA, length(ti));
	        for(i in 1:length(ti) )
	           y.mu[i] <-  mean(pheY_vec[pheT_vec==ti[i]], na.rm=T);

	        lines(ti, y.mu, col="black", lty="22", lwd=1);

	        par_X <- c( r$par[1] );
	        if ( !is.null(pheX) )
	   	        par_X <- c( par_X, r$par[ 2:(NCOL(pheX)+1)]);
	        par_c <- r$par[ -1 ];
	        if ( !is.null(pheX) )
		       par_c <- r$par[ -(1:(NCOL(pheX)+1)) ]
		    mu_X <- mean( cbind(1,pheX) %*% par_X, na.rm=T );

	        ti <- seq( min(pheT, na.rm=T), max(pheT, na.rm=T), 1);
	        mu_gen <- get_curve( obj.curve, par_c, ti, options=list(max.time=max(pheT, na.rm=T), min.time=min(pheT, na.rm=T)))

	        lines(ti, mu_gen + mu_X, col="black", lwd=1.5);
	    }
	}

	if(!is.null(file.plot.pdf)) dev.off();


	est.summary = data.frame(type = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$type);     }) ),
							 parm = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$par.count);}) ),
							 AIC  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$AIC);      }) ),
							 AICc = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$AICc);     }) ),
							 BIC  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$BIC);      }) ),
							 SSE  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$SSE);      }) ),
							 MSE  = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$MSE);      }) ),
							 RMSE = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$RMSE);     }) ),
							 R2   = unlist( lapply(1:length(est.curve), function(i){return(est.curve[[i]]$R2);       }) ) );

	## use AIC to determine the best curve type.
	cat("  Curve Fitting Summary:\n");
	old.width <- options(width=132)
	show(est.summary);
	options(width=old.width$width)

	fit.idx <- which.min( est.summary[,3] );
	return(list(error=F, type=est.curve[[fit.idx]]$type, y.resd=est.curve[[fit.idx]]$y.resd, par = est.curve[[fit.idx]]$par, summary=est.summary));

}

fg_fit_covar<-function( pheY, pheX, pheT, y.resd, obj.curve, covariance.type="auto")
{
	cat(" Searching covariance matrix .......\n" );

	obj.covars <- list();
	if(covariance.type=="auto" || covariance.type=="" || is.null(covariance.type) )
	{
		for ( i in 1:fg_get_covar_count()) obj.covars[[i]] <- fg.getCovariance( i );
	}
	else
	{
		for(i in 1:length(covariance.type))
			obj.covars[[i]] <- fg.getCovariance(covariance.type[i]);
	}

	est.covar <- list();
	index = 1
	for ( obj.covar in obj.covars)
	{
		cat(" *[", index, "/", fg_get_covar_count(), "]: try covariance matrix: ", obj.covar@type, "\n" );

		r.est <- proc_est_covar( y.resd, NULL, pheT, obj.curve, obj.covar );
		if( r.est$error )
			warning("Can not estimate the parameters of covariance matrix according to the covariance function[ covariance.type=", obj.covar@type, "]" )
		else
		{
			r.est$type <- obj.covar@type;
			est.covar[[index]] <- r.est;
			index <- index+1;
		}
	}

	if(length(est.covar)>=1)
	{
		est.summary = data.frame(type = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$type);     }) ),
							 L    = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$logL[1]);      }) ),
							 AIC  = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$AIC);      }) ),
							 BIC  = unlist( lapply(1:length(est.covar), function(i){return(est.covar[[i]]$BIC);      }) ) );


		cat("  Covariance Estimation Summary:\n");
		old.width <- options(width=132)
		show(est.summary);
		options(width=old.width$width)

		fit.idx <- which.min(est.summary[,3]);
		return(list(error=F, type=as.character(est.covar[[fit.idx]]$type), par = est.covar[[fit.idx]]$par, summary=est.summary));
	}
	else
		return(list(error=T));
}
