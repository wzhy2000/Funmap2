Qtlmle.qtlscan<-function( dat, qtl.table=NULL, grp_idx=NULL, options=list(scan.step=1))
{
	n.cores <- options$n.cores
	if(is.null(n.cores)	) n.cores <- 1;

	# search for all marker for current group(chromosome)
	if(is.null(qtl.table))
		qtl.table <- fin.get_qtl( dat$obj.gen$marker.table, scan.step=options$scan.step );

	if(!is.null(grp_idx))
		qtl.table <- qtl.table[which(qtl.table$grp %in% grp_idx), ]

	task0 <- task_start( "Execute the hypothesis test ", dat$obj.cross$hp_desc, "...\n" );

	grps <- unique(qtl.table[,1]);
	res.list <- mclapply( grps, function(i)
	{
		qtl.table.i <- 	qtl.table[which(qtl.table$grp==i),,drop=FALSE];

		grp.res <- list();
		h0 <- NULL;
		old_geno_pos <- -1;
		for(nQtl in 1:NROW( qtl.table.i))
		{
			par.cross <-list(lmarker = qtl.table.i[nQtl, 2],
					 rmarker = qtl.table.i[nQtl, 2]+1,
					 qtl.pos = qtl.table.i[nQtl, 4],
					 dist    = qtl.table.i[nQtl, 3])

			if ( par.cross$lmarker != old_geno_pos )
			{
				h0 = NULL;
				old_geno_pos <- par.cross$lmarker;
			}

			ltest <- try( Qtlmle.get_est_LR2(dat, par.cross, h0 ) );
			if (ltest$error)
				grp.res[[nQtl]] <- c( qtl.table.i[nQtl,1], qtl.table.i[nQtl,5], rep(NA, 5) )
			else
			{
				lr2 <- 2*( ltest$h0$value - ltest$h1$value );
				grp.res[[nQtl]] <- c( grp = qtl.table.i[nQtl,1],
							pos = qtl.table.i[nQtl,5],
							LR  = lr2,
							H1  = ltest$h1$value,
							H0  = ltest$h0$value,
							ltest$h1$par,
							ltest$h0$par );
				h0 <- ltest$h0;
			}

			if(.RR("debug", TRUE))
				cat(nQtl, "\t", grp.res[[nQtl]]$grp, "\t", grp.res[[nQtl]]$pos, "\t", grp.res[[nQtl]]$LR, "\n");

		}

		if(n.cores==1)
			task_elapsed(task0, finished = which(grps==i)[1]/length(grps), "Group:", which(grps==i)[1], "/", length(grps), ",", "$SYS_PROMPT$","\n" );

		return ( do.call("rbind", grp.res) );
	}, mc.cores=n.cores );

	task_stop("The hypothesis test is done\n");

	res <- do.call("rbind", res.list);

	ret<-list(
			name      = "Qtlmle.MLE",
			error     = FALSE,
			param     = list(grp_idx=grp_idx, scan.step=options$scan.step, n.cores=options$n.cores),
			cross.type= dat$obj.cross$type,
			covar.type= dat$obj.covar@type,
			curve.type= dat$obj.curve@type,
			full.res  = res);

	ret = fin.qtl_locate( ret, peak.count = .RR("peak.count", 5 ) );
	return( ret );
}

Qtlmle.get_est_LR2<-function( dat, par.cross, h0=NULL )
{
	pheY <- dat$obj.phe$pheY;
	pheT <- dat$obj.phe$pheT;
	pheX <- dat$obj.phe$pheX;

	# current two markers in current inteval
	flank_snps <- dat$obj.gen$genos.matrix[, c(par.cross$lmarker, par.cross$rmarker)];
	missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );

	if ( length(missing) > 0)
	{
		flank_snps <- flank_snps[ -(missing),,drop=F ];
		pheY <- pheY[ -(missing),,drop=F ];
		pheT <- pheT[ -(missing),,drop=F ];
		if(!is.null(pheX)) pheX <- pheX[ -(missing),,drop=F];
	}

	# probability table for QQ,Qq,qq at markers MiMjNiNj;
	# 0: homozygote(qq)
	# 1: heterozygote(Qq)
	# 2: homozygote+additive(QQ)

	allprob <- dat$obj.cross$get_qtl_prob( flank_snps[,1], flank_snps[,2], par.cross$dist , par.cross$qtl.pos );
	allprob[which(allprob==0)]<-0.005;
	allprob[which(allprob==1)]<-0.995;

	n.par.cov <- get_param_info(dat$obj.covar, pheT )$count;
	n.par.curve <- get_param_info(dat$obj.curve, pheT )$count;
	nna.vec <- get_non_na_number(pheY);
	options <- list(max.loop=4, min.time=min(dat$obj.phe$pheT, na.rm=T), max.time=max(dat$obj.phe$pheT, na.rm=T) );

	if ( is.null( h0) )
	{
		parin.x <- dat$obj.phe$est.curve$parX;
		parin.curve <- dat$obj.phe$est.curve$param;
		parin.covar <- dat$obj.phe$est.covar$param;
		names(parin.curve) <- paste("H0_", get_param_info(dat$obj.curve, pheT )$names, sep="");
		names(parin.covar) <- paste("H0_", get_param_info(dat$obj.covar, pheT )$names, sep="");
		if(!is.null(parin.x))
			names(parin.x) <- paste("H0_", colnames(pheX), sep="");

		h0 <-  optim_BFGS ( parin.x, parin.curve, parin.covar,
					  Qtlmle.H0, pheY = pheY, pheT = pheT, pheX = pheX,
					  obj.curve = dat$obj.curve, obj.covar = dat$obj.covar,
					  nna.vec=nna.vec, n.par.curve=n.par.curve, options=options );

		if ( is.null(h0) )
			return(list(error=TRUE));

		if(!is.null(h0) && !is.null(h0$par) )
			names(h0$par) <- c( names(parin.x), names(parin.curve), names(parin.covar) )
	}

	parin.x <- c();
	if(!is.null(pheX))
		parin.x <- h0$par[ 1:NCOL(pheX) ];

	## WARNING: length(parin.x) != NCOL(pheX);
	n.gentype <- NCOL(allprob);
	parin.curve <- rep(h0$par[ (length(parin.x)+1):(length(parin.x)+n.par.curve) ], n.gentype );
	parin.covar <- h0$par[ -(1:(length(parin.x)+n.par.curve)) ] ;

	names.curve <- c();
	for(i in 1:n.gentype)
		names.curve <- c(names.curve, paste( "H1", get_param_info(dat$obj.curve, pheT )$names, i-1, sep="_"));
	names(parin.covar) <- paste("H1", get_param_info(dat$obj.covar, pheT )$names, sep="_");
	if(!is.null(parin.x))
		names(parin.x) <- paste("H1_", colnames(pheX), sep="");


	h1 <-  optim_BFGS ( parin.x, parin.curve, parin.covar,
					  Qtlmle.H1,  qtl.prob = allprob, pheY = pheY, pheT = pheT, pheX = pheX,
					  obj.curve = dat$obj.curve, obj.covar = dat$obj.covar,
					  nna.vec=nna.vec, n.par.curve=n.par.curve, options=options );

	if(!is.null(h1) && !is.null(h1$par) )
		names(h1$par) <- c( names(parin.x), names.curve, names(parin.covar) )

	return(list(h0=h0, h1=h1, error=ifelse ( !is.null(h1), FALSE, TRUE )));
}

Qtlmle.H0<-function( parin, pheY, pheT, pheX, obj.curve, obj.covar , nna.vec, n.par.curve, options )
{
	if (!is.null(pheX))
	{
		parin.X <- parin[1:NCOL(pheX)];
		parin <- parin[-c(1:NCOL(pheX))];
	}

	parin.curve <- parin[1:n.par.curve];
	parin.covar <- parin[-c(1:n.par.curve)];

	mat.cov <- get_matrix( obj.covar, parin.covar, pheT );

	if(!is.null(pheX))
		X  <- ( pheX %*% parin.X ) %*% t(rep(1,NCOL(pheY)))
	else
		X  <- 0;

	Y.delt <- pheY - get_curve(obj.curve, parin.curve, pheT, options=options ) - X;
	pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec=nna.vec, log=T), .RR("try.silent", FALSE) );
	if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

	A <- -sum( pv );
	return (A);
}

Qtlmle.H1<-function( parin, qtl.prob, pheY, pheT, pheX, obj.curve, obj.covar, nna.vec, n.par.curve, options)
{
	if (!is.null(pheX))
	{
		parin.X <- parin[1:NCOL(pheX)];
		parin <- parin[-c(1:NCOL(pheX))];
	}

	n.gentype <- NCOL(qtl.prob);
	parin.curve <- parin[1:(n.par.curve * n.gentype)];
	parin.covar <- parin[-c(1:(n.par.curve*n.gentype))];

	mat.cov <- get_matrix( obj.covar, parin.covar, pheT );

	fy <- matrix( NA, nrow=NROW(pheY), ncol=n.gentype );
	for(k in 1:n.gentype)
	{
		parin.curve0 <- parin.curve[1:n.par.curve];
		parin.curve  <- parin.curve[-c(1:n.par.curve)];

		if(!is.null(pheX))
			X  <-  ( pheX %*% parin.X ) %*% t(rep(1,NCOL(pheY)))
		else
			X  <- 0;

		Y.delt <- pheY - get_curve(obj.curve, parin.curve0, pheT, options=options ) - X;
		pv <- try(dmvnorm_fast( Y.delt, rep(0, NCOL(mat.cov)), mat.cov, nna.vec=nna.vec, log=F), .RR("try.silent", FALSE) );
		if ( class(pv)=="try-error" || is.na(pv) )
			return (NaN);

		fy[,k] <- pv;
	}

	pf <- fy * qtl.prob;
	A <- - sum( log( rowSums( pf ) ) );

	.RW( "mle.pdf", length(which( pf>1 | pf<0 ) ) );
	return (A);

}



#--------------------------------------------------------------
# private: fin.get_qtl_scanpoints_grp
#
# Get QTL intervals between two markers.
#
# input
#       mrk_dist: distance between two markers.
#           step: step
# output: the vector of the intervals.
#--------------------------------------------------------------
fin.get_qtl<-function( marker.table, scan.step=2)
{
	qtl_list <- lapply (1:(NROW(marker.table)-1), function(i){
		if( marker.table[i+1,3] != marker.table[i,3] )
			return(c());

		mrk_dist  <- marker.table[i+1,2] - marker.table[i,2];
		if ( mrk_dist > scan.step )
			qtls <- seq(0, mrk_dist, scan.step)
		else
			qtls <- c(0,mrk_dist);

		if(! mrk_dist %in% qtls) qtls <- c( qtls, mrk_dist);

		# remove last one!
		if( (i != NROW(marker.table)-1) && ( i+2<=NROW(marker.table) && marker.table[i+1,3] == marker.table[i+2,3] ) )
			qtls <- qtls[-length(qtls)]

		return( data.frame( grp=marker.table[i,3], idx=i, dist=mrk_dist, qtl=qtls, pos=qtls + marker.table[i,2] ) );
	});

	return ( do.call("rbind", qtl_list));
}

fin.qtl_locate <- function( res, pvalue=NULL, cutoff=NULL, peak.count=NULL)
{
	full.res <- res$full.res [ !is.na(rowSums(res$full.res)),,drop=F];
	if(NROW(full.res)<=0)
	{
		warning("No QTL result.");
		return( res );
	}

	res$full.res <- full.res;
	peaks.candidate <- fin.select_peaks_by_simple_way( full.res[,3] );
	if ( length(peaks.candidate) <= 0 )
	{
		warning("No peaks in QTL result.");
		return( res );
	}

	peaks.idx <- c();
	if ( !is.null(peak.count) )
	{
		chromo5 <-c();
		for (i in 1:length(peaks.candidate))
		{
			if (! (full.res[peaks.candidate[i],1] %in% chromo5) )
			{
				chromo5 <- c( chromo5, full.res[peaks.candidate[i],1] );
				peaks.idx <- c( peaks.idx, peaks.candidate[i] );
				if (length(chromo5) >= peak.count )
					break;
			}
		}

		res$threshold.type <- "count";
		res$threshold <- peak.count;
	}

	if(!is.null(cutoff))
	{
		res$threshold.type <- "LR";
		res$threshold <- cutoff;
	}

	if (!is.null(pvalue) && !is.null( res$obj.permu) )
	{
		pv.table <- res$obj.permu$full.res;
		idx <- round(NROW(pv.table)*pvalue);
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

		cutoff <- sort(pv.table[,3], decreasing=TRUE)[idx];
		res$threshold.type <- "pvalue";
		res$threshold <- pvalue;
	}

	if(!is.null(cutoff))
	{
		peaks.idx <- peaks.candidate [ full.res[peaks.candidate,3] > cutoff ]
	}

	res$qtl.peaks <- peaks.idx;
	return (res );
}


#--------------------------------------------------------------
# Qtlmle.summary
#
# Summarize the result object of the hypothesis test.
# Used by summary( Qtlmle.ret object );
#--------------------------------------------------------------
Qtlmle.summary <- function( res )
{
	if (res$error )
	{
		str<- paste( "Failed to do hypothesis test(10).\n\n", sep="" )
		return(str);
	}

	str <- "";
	str<- paste( str, "Hypothesis test: \n", sep="" );
	str<- paste( str, res$obj.cross$hp_desc, "\n", sep="" )

	str<- paste( str, "------------------------------------\n", sep="" );

	str1 <- sprintf("%15s: %s\n", 	 "Cross", 	   res$cross.type );
	str0 <- sprintf("%15s: %s\n", 	 "Curve",	   res$curve.type );
	str2 <- sprintf("%15s: %s\n", 	 "Covariance", res$covar.type );
	str3 <- sprintf("%15s: %s=%f\n", 	 "QTL criterion", res$threshold.type, res$threshold);
	str  <- paste(str, str0, str1, str2, "\n", str3, sep="");

	st0 <- sprintf("\nIndex \t Grp\t Pos.\t LR" )
	str <- paste(str, st0, sep="");
	if(length(res$qtl.peaks)>0)
	{
		for(i in res$qtl.peaks)
		{
			st0 <- sprintf("%d\t%d\t%4.1f\t%4.3f", i, res$full.res[i,1], res$full.res[i,2], res$full.res[i,3] ) ;
			str <- paste(str, st0, sep="\n");
		}
	}
	else
		str <- paste(str, "No QTL peaks are selected.", sep="\n");

	str<- paste( str, "\n------------------------------------\n\n", sep="" );

	return(str);
}

#--------------------------------------------------------------
# Qtlmle.plot
#
# Plot the figure for the hypothesis test.
# Used by plot( NP.ret object );
#--------------------------------------------------------------
Qtlmle.plot<-function( res, plot_type=NULL, pdf_file=NULL  )
{
	if (res$error )
		stop( "Failed to do hypothesis test(10).");

	if (is.null(pdf_file))  X11() else pdf(pdf_file);

	if( is.null(plot_type) || plot_type==1)
	{
	  	cutoff.05 <- NULL;
	  	cutoff.01 <- NULL;
	  	if(!is.null(res$obj.permu))
	  	{
	  		cutoff.05 <- fin.find_cutoff(res$obj.permu, 0.05);
	  		cutoff.01 <- fin.find_cutoff(res$obj.permu, 0.01);
	  	}

		grp.idx <- unique(res$full.res[,1]);
		marker.table <- res$obj.gen$marker.table[ res$obj.gen$marker.table$grp_idx %in% grp.idx,,drop=F]

	  	err.fig<-try ( fpt.plot_qtl_map( res$full.res, marker.table, cutoff.05=cutoff.05, cutoff.01=cutoff.01  ) );
		if (class(err.fig)!="try-error")
			title("The LR profile");
	}

	if( is.null(plot_type) || plot_type==2 )
	{
		if( length(res$qtl.peaks) <= 0)
			cat(" ! No QTL peaks in the result.\n")
		else
		for (i in 1:length(res$qtl.peaks))
		{
			##QTL profile of single chromosome
			qtl.info <- res$full.res[ res$qtl.peaks[i], ]
			err.fig <- try ( fpt.plot_qtl_pos( res$full.res, res$obj.gen$marker.table,  qtl.info[1], qtl.info[2], cutoff.05=cutoff.05, cutoff.01=cutoff.01 ) );
			if (class(err.fig)!="try-error")
				title(paste("The LR profile for QTL position(Group:", qtl.info[1], ")", sep="") );

			QQ_par <- qtl.info[paste( "H1", get_param_info(res$obj.curve, res$obj.phe$pheT )$names, 0, sep="_")];
			Qq_par <- qtl.info[paste( "H1", get_param_info(res$obj.curve, res$obj.phe$pheT )$names, 1, sep="_")];
			qq_par <- qtl.info[paste( "H1", get_param_info(res$obj.curve, res$obj.phe$pheT )$names, 2, sep="_")];
			if (any(is.na(QQ_par))) QQ_par <- NULL;
			if (any(is.na(Qq_par))) Qq_par <- NULL;
			if (any(is.na(qq_par))) qq_par <- NULL;

			##QTL curve
			err.fig <- try ( fpt.plot_com_curve ( max(res$obj.phe$pheT, na.rm=T),
									round(max(res$obj.phe$pheT, na.rm=T) *1.1),
									res$obj.curve,
									pheY   = res$obj.phe$pheY,
									pheT   = res$obj.phe$pheT,
									QQ_par = QQ_par,
									Qq_par = Qq_par,
									qq_par = qq_par ) );

			if (class(err.fig)!="try-error")
				title(paste("QTL Curves (Group:", qtl.info[1], ", Pos.:",qtl.info[2], ")", sep="") );
		}
	}

	if( is.null(plot_type) || plot_type==3)
	{
		if(!is.null(res$obj.permu))
			fpt.plot_permutation (res$obj.permu$pv.table )
		else
			cat("No permutation data for the plot.\n");
	}

	if (!is.null(pdf_file))
	{
		cat( "*The results of QTL mapping are saved to ", pdf_file, ".\n" );
		dev.off()
	}
}

#--------------------------------------------------------------
# public: fin.select_peaks_by_simple_way
#
#--------------------------------------------------------------
fin.select_peaks_by_simple_way<-function( curve )
{
	curv2 <- curve[2:length(curve)];
	diff  <- curv2 - curve[1:(length(curve)-1)];

	peaks1<-c();
	for (i in 2:length(diff) )
	{
		if ( diff[i-1]>0 && diff[i]<0 )
			peaks1 <- c(peaks1, i);
	}

	index <- sort(curve[peaks1], index.return=TRUE )$ix;
	peak_index <- c(index[1]);

	x <- curve[peaks1];
	o <- order(x, decreasing = TRUE);
	return ( peaks1[o] );
}

#--------------------------------------------------------------
# Qtlmle.set_lr2_cutoff
#
#
#--------------------------------------------------------------
Qtlmle.set_lr2_cutoff<-function( res, p05, p01)
{
	newres <- res;
	newres$qtl_p05<- p05;
	newres$qtl_p01<- p01;
	if (newres$qtl_pvalue<p05)
		newres$qtl_pvalue<- ">0.05"
	else
	if (newres$qtl_pvalue<p01)
		newres$qtl_pvalue<- "<0.05"
	else
		newres$qtl_pvalue<- "<0.01"

	return (newres);
}


ZZZ.regmodel<-function()
{
	FM2_ENV$models <-list(
		name 	       = "QTL Model",
		get_est_LR2    = Qtlmle.get_est_LR2,
		qtlscan        = Qtlmle.qtlscan,
		set_lr2_cutoff = Qtlmle.set_lr2_cutoff);
}
