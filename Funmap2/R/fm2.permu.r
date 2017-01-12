#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    5) permutation()
#    6) summary_perm_ret()
#    7) plot_perm_ret()
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# permu.execute
#
# do a permutation to get a p-value vs cutoff table.
# Same as T10 test.
#
# Input
#    dat    : data object
# Output
#    Cut_off table
#--------------------------------------------------------------
permu.execute<-function( dat, grp.idx=NULL, permu.loop, filter.ratio=1, scan.step=1, n.cores=1 )
{
	if(filter.ratio<1)
		dat$clustering <- permu.cluster(dat);
	
	#p0 <- task_start("Execute the permutation, nCount=", permu.loop, "...\n");

	cpu.fun <- function(i)
	{
		dat0 <- dat;
		new_index <- sample( NROW( dat0$obj.phe$pheY ) );
		dat0$obj.phe$pheY <- dat0$obj.phe$pheY[new_index, ];
		dat0$obj.phe$pheX <- dat0$obj.phe$pheX[new_index, ];
		dat0$obj.phe$pheT <- dat0$obj.phe$pheT[new_index, ];
		dat0$clustering$Q <- dat0$clustering$Q[new_index, ];
		
		qtl_table <- NULL;
		if(filter.ratio<1)
			qtl_table <- permu.qtlfilter( dat0, grp.idx, filter.ratio, scan.step );

		r <- try( Qtlmle.qtlscan( dat0, qtl_table ), FALSE);
		if (class(r)=="try-error")
			return(c(i,NA, NA));

		lr2.range <- c( i, mean(r$full.res[,3], na.rm=TRUE ), max(r$full.res[,3], na.rm=TRUE) );
		cat("Permutation Test(", i,"): LR2=", lr2.range[3], "\n");
		return( lr2.range  )
	}

	r.cluster <- mclapply(1:permu.loop, cpu.fun, mc.cores=n.cores);

	x2 <- do.call("rbind", r.cluster);
	p_cut <- c( 0.9,  0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
			0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01);

	if (permu.loop>=1000)
		p_cut <- c(p_cut, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001 );

	if (permu.loop>=10000)
		p_cut <- c(p_cut, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001);

	if (permu.loop>=100000)
		p_cut <- c(p_cut, 0.00001);

	pv2<-array(Inf, dim=c(length(p_cut), 2) );
	pv2[,1] <- p_cut
	for(i in 1:length(p_cut))
	{
		pv2.order <- round (length(x2[,3]) * p_cut[i]);
		if ( pv2.order>0)
		{
			pv2[i,1] <- p_cut[i];
			pv2[i,2] <- sort(x2[,3], decreasing=TRUE )[ pv2.order ];
		}
	}

	ret<-list(
		name      = paste( dat$obj.phe$curve.type, ".permu",sep=""),
		pheno.csv = dat$obj.phe$pheno_csv,
		cross.type= dat$obj.phe$cross.type,
		curve.type= dat$obj.phe$curve.type,
		covar.type= dat$obj.phe$covar.type,
		permu.loop= permu.loop,
		pv.table  = pv2,
		full.res  = x2,
		param     = list(permu.loop=permu.loop, permu.filter.ratio=filter.ratio, scan.step=scan.step, n.cores=n.cores));

	return (ret);
}

#--------------------------------------------------------------
# permu.summary
#
# Summarize the result object of the permutation tests.
# Used by summary( XX.perm.ret object );
#--------------------------------------------------------------
permu.summary<-function( object )
{
	res <- object;

	str <- "";
	str <- paste( str, "Permutation result: \n", sep="" );
	str <- paste( str, "------------------------------------\n", sep="" );
	str0 <- sprintf("%15s: %s\n", 	 "Cross", 	  res$cross.type  );
	str1 <- sprintf("%15s: %s\n", 	 "Curve",	  res$curve.type );
	str2 <- sprintf("%15s: %s\n", 	 "Covariance",res$covar.type  );
	str3 <- sprintf("%15s: %d\n", 	 "Permutation", res$permu.loop);
	str <- paste( str, str0, str1, str2, str3, sep="" );
	str <- paste( str, "------------------------------------\n", sep="" );
	str <- paste( str, "\n", sep="" );

	str<- paste( str, "p-value\tCutoff\n", "\n", sep=" " );
	for (i in 1:NROW(res$pv.table) )
	{
		st0 <- sprintf("%7.5f\t%10.2f\n", res$pv.table[i,1], res$pv.table[i,2] )
		str <- paste( str, st0, sep="" );
	}

	str<- paste( str, "\n", sep=" " );

	return(str);
}

#--------------------------------------------------------------
# permu.plot
#
# Plot the result object of the permutation tests.
# Used by plot( XX.perm.ret object );
#--------------------------------------------------------------
permu.plot<-function( object, file.pdf =NULL )
{
	if(is.null(file.pdf)) X11() else pdf(file.pdf);

	fpt.plot_permutation ( object$pv.table );

	if(!is.null(file.pdf)) 
	{
		dev.off();
		cat( "* The permutation cutoff are shown to ", file.pdf, ".\n");
	}	

	invisible();
}

fin.find_cutoff <- function(perm, p=0.05)
{
	perm.cutoff <- NULL;
	if (!is.null(perm))
	{
		pv.table <- perm$full.res;
		idx <- round(NROW(pv.table)*p);
		if( idx > NROW(pv.table) )
		{
			warning("pvalue is too big and no exact cutoff in the permutation results.")
			idx <- NROW(pv.table);
		}
		if (idx < 1 )
		{
			warning("pvalue is too small and no exact cutoff in the permutation results.")
			idx <- 1;
		}
					
		perm.cutoff <- sort(pv.table[,3], decreasing=TRUE)[idx];
	}

	return( perm.cutoff );
}


#------------------------------------------------
#
#  filter
#
#
#------------------------------------------------
permu.qtlfilter <- function(dat, grp.idx=NULL, percent = 1, scan.step=2 )
{
	qtl.table <- fin.get_qtl( dat$obj.gen$marker.table, scan.step );
	if(!is.null(grp.idx))
		qtl.table <- qtl.table[which(qtl.table$grp %in% grp.idx), ]

	qtl.cor <- list();
	for ( nQtl in 1:NROW(qtl.table) )
	{
		par.cross <-list(lmarker = qtl.table[nQtl, 2],
					 rmarker = qtl.table[nQtl, 2]+1,
					 qtl.pos = qtl.table[nQtl, 4],
					 dist    = qtl.table[nQtl, 3])

		Q <- dat$clustering$Q;

		flank_snps <- dat$obj.gen$genos.matrix[, c(par.cross$lmarker, par.cross$rmarker)];
		missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );
		if ( length(missing) > 0)
		{
			flank_snps <- flank_snps[ -(missing), ];
			Q <- Q[-(missing), ];
		}

		P <- dat$obj.cross$get_qtl_prob( flank_snps[,1], flank_snps[,2], par.cross$dist , par.cross$qtl.pos );
		P[ which(P==0) ]<-0.005;
		P[ which(P==1) ]<-0.995;

		qtl.cor[[nQtl]] <- c( grp=qtl.table[nQtl, 1], idx=qtl.table[nQtl, 2], pos=qtl.table[nQtl, 5],
				L1=getpossibleLogL1( Q, P ),
				L2=getpossibleLogL2( Q, P ),
				L3=getpossibleLogL3( Q, P ),
				L4=getpossibleLogL4( Q, P ));
	}

	L <- do.call("rbind", qtl.cor);

	L1.order <- order(L[,"L1"], decreasing=TRUE)[c(1:round(NROW(L)*percent))];
	L2.order <- order(L[,"L2"], decreasing=TRUE)[c(1:round(NROW(L)*percent))];
	L3.order <- order(L[,"L3"], decreasing=TRUE)[c(1:round(NROW(L)*percent))];
	L4.order <- order(L[,"L4"], decreasing=TRUE)[c(1:round(NROW(L)*percent))];

	L.filted <- sort( unique(c(L1.order, L2.order, L3.order, L4.order)) );
	qtls.filter <- qtl.table[L.filted, ,drop=F];

	return(qtls.filter)
}

getpossibleLogL1<-function( Q, P)
{
	get_x2_num <- function(q,p)
	{
		vec.x <- rowSums(q*p);
		vec.x <- vec.x[which(vec.x!=0)]
		return(-sum(log(vec.x)))
	}

	if( NCOL(Q) ==3)
	{
		sig.x <- max( get_x2_num(Q[,c(1,2,3)],P),
		   	 get_x2_num(Q[,c(1,3,2)],P),
		   	 get_x2_num(Q[,c(2,1,3)],P),
		   	 get_x2_num(Q[,c(2,3,1)],P),
		   	 get_x2_num(Q[,c(3,1,2)],P),
	   		 get_x2_num(Q[,c(3,2,1)],P) )
	}

	if( NCOL(Q) ==2)
	{
		sig.x <- max( get_x2_num(Q[,c(1,2)],P), get_x2_num(Q[,c(2, 1)],P))
	}

	return(sig.x)
}

getpossibleLogL2<-function( Q, P)
{
	get_x2_num <- function(q,p)
	{
   		x<- lapply(1:NROW(q), function(i)
   		{
     		M <- round(cbind(q[i,], p[i,])*NROW(q));

			if (min(M)<=5)
				f <- fisher.test(M)
			else
				f <- chisq.test(M)

			return(f$p.value >0.01);
		});

		sum(unlist(x));
	}

	Q <- Q/rowSums(Q);
	if( NCOL(Q) ==3)
	{
		sig.x <- max( get_x2_num( Q[,c(1,2,3)],P),
		   	 get_x2_num( Q[,c(1,3,2)],P),
		   	 get_x2_num( Q[,c(2,1,3)],P),
		   	 get_x2_num( Q[,c(2,3,1)],P),
		   	 get_x2_num( Q[,c(3,1,2)],P),
	   		 get_x2_num( Q[,c(3,2,1)],P) )
	}

	if( NCOL(Q) ==2)
   	{
   		sig.x <- max( get_x2_num( Q[,c(1,2)],P ), get_x2_num( Q[,c(2, 1)],P) )
   	}

	return(sig.x)
}

getpossibleLogL3<-function( Q, P)
{
	get_x2_num <- function(q,p)
	{
    	x<- lapply(1:NROW(q), function(i){ cor(q[i,], p[i,]) } );
    	return(mean(unlist(x)));
	}


	if( NCOL(Q) ==3)
	{
		sig.x <- max( get_x2_num(Q[,c(1,2,3)],P),
		   	 get_x2_num(Q[,c(1,3,2)],P),
		   	 get_x2_num(Q[,c(2,1,3)],P),
		   	 get_x2_num(Q[,c(2,3,1)],P),
		   	 get_x2_num(Q[,c(3,1,2)],P),
   			 get_x2_num(Q[,c(3,2,1)],P) )
   	}

   	if( NCOL(Q) ==2)
   	{
   		sig.x <- max( get_x2_num(Q[,c(1,2)],P), get_x2_num(Q[,c(2, 1)],P))
   	}

	return(sig.x)
}

getpossibleLogL4<-function( Q, P)
{
	if( NCOL(Q) ==3)
	{
		sig.x <- max(cor(c(Q[,c(1,2,3)]), c(P)),
			cor(c(Q[,c(1,3,2)]), c(P)),
			cor(c(Q[,c(2,1,3)]), c(P)),
			cor(c(Q[,c(2,3,1)]), c(P)),
			cor(c(Q[,c(3,1,2)]), c(P)),
			cor(c(Q[,c(3,2,1)]), c(P)))
	}

	if(NCOL(Q)==2)
	{
		sig.x <- max(cor(c(Q[,c(1,2)]), c(P)), cor(c(Q[,c(2, 1)]), c(P)))
	}

	return(sig.x)
}


#--------------------------------------------------------------
# permu.cluster
#
# Curve Clustering method used for filtering QTL in the permutation
#
# 1) estimating the paramerters for the mixture model
# 2) calculating the probability matrix for the curve clusterring
#
#--------------------------------------------------------------
permu.cluster <- function( dat )
{
	obj.curve <- dat$obj.curve;
	obj.covar <- dat$obj.covar;
	obj.cross <- dat$obj.cross;

	n.grp     <- obj.cross$gen_num;
	n.par.cov <- get_param_info(obj.covar, dat$obj.phe$pheT)$count;
	n.par.curve <- get_param_info(obj.curve, dat$obj.phe$pheT)$count;

	mle <- function(par, W, pheY, pheT, pheX, obj.curve, obj.covar, obj.cross, optim=T)
	{
		par.covar <-  par[1:n.par.cov];
		if (!check_param( obj.covar, par.covar ))
			return(NaN);

		cov <- get_matrix(obj.covar, par.covar, pheT  );
		if (is.null(cov) ) return(NaN);

		par.m <- matrix(par[-c(1:n.par.cov)], nrow=n.grp, byrow=TRUE)

		i <- 1;
		Q <- matrix(NA, ncol=n.grp, nrow=NROW(pheY));
		L <- matrix(NA, ncol=n.grp, nrow=NROW(pheY));
		if (obj.cross$gen_QQ)
		{
			y.delt <- pheY- get_curve( obj.curve, par.m[i,], pheT  );
			Q[,i] <- dmvnorm_fast( y.delt, rep(0, n.par.cov), cov, log=F);
			L[,i] <- Q[,i]*W[i]
			i <- i+1;
		}

		if (obj.cross$gen_Qq)
		{
			y.delt <- pheY- get_curve( obj.curve, par.m[i,], pheT  );
			Q[,i] <- dmvnorm_fast( y.delt, rep(0, NCOL(y.delt)), cov, log=F);
			L[,i] <- Q[,i]*W[i]
			i <- i+1;
		}

		if (obj.cross$gen_qq)
		{
			y.delt <- pheY- get_curve( obj.curve, par.m[i,], pheT  );
			Q[,i] <- dmvnorm_fast( y.delt, rep(0, NCOL(y.delt)), cov, log=F);
			L[,i] <- Q[,i]*W[i]
			i <- i+1;
		}

		if(optim)
			return( -sum(log(rowSums(L))) )
		else
			return( list(L=L, Q=Q) );
	}

	W <- rep(1, n.grp)
	W.new <- rep(1/n.grp, n.grp);
	par.new <- c(dat$obj.phe$est.covar$param, rep( dat$obj.phe$est.curve$param, n.grp) )
	while( max(abs( c(W.new) - c(W) ) ) > 1e-5 )
	{
		W <- W.new;
		par <- par.new;
		r<- optim( par, mle, W=W, pheY=dat$obj.phe$pheY, pheT=dat$obj.phe$pheT, pheX=dat$obj.phe$pheX, obj.curve=dat$obj.curve, obj.covar= dat$obj.covar, obj.cross=dat$obj.cross, optim=T, control = list(maxit=50000));

		if(r$convergence==0)
		{
			par.new <- r$par
			r.est <- mle(par.new, W, dat$obj.phe$pheY, dat$obj.phe$pheT, dat$obj.phe$pheX, dat$obj.curve, dat$obj.covar, dat$obj.cross, optim=F);
			Q <- r.est$Q
			L <- r.est$L / rowSums(r.est$L);
			W.new <- colSums(L)/NROW(dat$obj.phe$pheY);
cat("*******w", W.new, "\n");

			if( any(which(W.new <0.05)))
			{
				W.idx <- which(W < 0.05);
				W.new[ W.idx ] <- runif(1,1,n.grp)/n.grp;
				W.new <- W.new/sum(W.new);

				par.m <- matrix(par.new[-c(1:n.par.cov)], nrow=n.grp, byrow=TRUE)
				par.m[ W.idx, ] <- colMeans( par.m[ -W.idx,,drop=F ] );
				par.new <- c( dat$obj.phe$est.covar$param, c(t(par.m)));
cat("!!!!!!!!", par.new, "\n");
			}
		}
		else
		{
			W.new <- (W.new+0.1)/sum( W.new+0.1 )
cat("++++++++++w.new", W.new, "\n")
		}
	}

	cat("full.par=", par.new, "\n");
	return(list(par=par.new, W=W, Q=Q))
}
