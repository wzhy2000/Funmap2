LOG_ERR  <- 0;
LOG_WARN <- 1;
LOG_INFO <- 2;
LOG_DEBUG<- 3;

reset_seed<-function(options)
{
	.RW("n.seed", .RR("n.seed") + runif(1, 1, 5) );
	set.seed( .RR("n.seed" ) );
	if( .RR("debug")) cat(".");
}

.RR<-function(key, def=NULL)
{
	val <- try( get( key, FG_ENV$i_hash), silent = TRUE);
	if(class(val)=="try-error" || is.null(val) )
		return(def)
	else
		return(val);
}

.RW<-function(key, value)
{
	old_value <- NULL;
	old_value <-  try( get( key, FG_ENV$i_hash), silent = TRUE);
	if(class(old_value)=="try-error" || is.null(old_value) )
		old_value <- NA;

	assign( key, value, FG_ENV$i_hash );
	return (old_value);
}

get_con_param<-function(parm.id)
{
	for (e in commandArgs())
	{
		ta = strsplit(e,"=", fixed=TRUE);
		if(! is.na( ta[[1]][2]))
		{
			temp = ta[[1]][2];
			if( ta[[1]][1] == parm.id) {
				return (as.character(temp));
			}
		}
	}

	return(NA);
}

remove_extname<-function(filename)
{
	rs <- unlist(strsplit(filename, "\\."))
	if (length(rs)==1)
		return(filename);

	rs.last <- rs[length(rs)];

	if (grepl(pattern = "[\\/]", x = rs.last))
		return(filename);

	rs <- rs[-length(rs)];
	file.noext <- paste(rs, collapse=".");

	return(file.noext);
}


colSds<-function(tb, na.rm=T)
{
	unlist( lapply(1:NCOL(tb), function(i){sd(tb[,i], na.rm=na.rm)}) );
}

get_non_na_number<-function(y.resd)
{
	nna.mat <- array(1, dim=c(NROW(y.resd), NCOL(y.resd)))
	nna.mat[ which(is.na(y.resd)) ] <- 0;

	nna.vec <- rep(0, NROW(nna.mat));
	for(i in 1:NCOL(nna.mat))
		nna.vec <- nna.vec*2 + nna.mat[,i];

	return(nna.vec);
}

dmvnorm_fast<-function( y.resd, mu, cov.mat, nna.vec=NULL, log=T)
{
	if(is.null(nna.vec))
		nna.vec = get_non_na_number(y.resd);

	pv <- rep(NA, NROW(nna.vec));
	for( nna in unique(nna.vec) )
	{
		if(nna>0)
		{
			nna.idx <- which(nna.vec==nna);
			col.idx <- !is.na(y.resd[nna.idx[1],]);
			pv[nna.idx] <- dmvnorm( y.resd[nna.idx, col.idx, drop=F ],
						   mu[col.idx] ,
						   cov.mat[col.idx, col.idx, drop=F], log=log);
		}
	}

	return(pv);
}

optim_BFGS<-function ( par.x, par.curve, par.covar, proc_h0, ... )
{
    options <- list(...)$options;
	mle.control <- list(loop = 0, loop.optim = 0);
	parinx <- c(par.x, par.curve, par.covar);
	n.par <- length( parinx );
	control <- list(maxit = 500 + n.par * 200, reltol=1e-8);
	h0.best <- list(value = Inf, par = rep(NA, length(parinx)) );

	while ( mle.control$loop < options$max.loop )
	{
		h0 <- try( optim( parinx, proc_h0, ...,
					method  = ifelse(mle.control$loop.optim%%2==0, "BFGS", "Nelder-Mead" ),	control = control ),
					.RR("try.slient") );

		mle.control$loop.optim <- mle.control$loop.optim + 1;
		if (class(h0)=="try-error" || any(is.na(h0)) || h0$convergence!=0 )
		{
			if ( is.list(h0) )
			{
				if( h0$convergence == 1)
				{
					control$maxit <- control$maxit*2;

					if ( control$maxit > 500*4096 )
					{
						control$maxit <- 500*4096;
						control$reltol <- control$reltol*2;
					}
				}
				else
					cat("optim, convergence=", h0$convergence, "", h0$message, "\n")
			}

			reset_seed();
			parinx <- c(par.x*runif(length(par.x), 0.95, 1.05), par.curve*runif(length(par.curve), 0.95, 1.05), par.covar );
			next;
		}
		else
		{
			if ( h0$value < h0.best$value ) { h0.best <- h0; }
		}

		mle.control$loop <- mle.control$loop + 1;
	}

	return(h0.best);
}
