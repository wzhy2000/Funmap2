#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) dat.get_simuparam()
#    2) dat.summary_par()
#    3) dat.simulate()
#    4) dat.load()
#    5) dat.summary()
#    6) dat.report()
#    7) dat.plot()
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# dat.summary_par
#
# Summarize the parameter object for the F2 simulation
# Used by summary( XX.F2.par object)
#--------------------------------------------------------------
dat.summary_par<-function( obj.par )
{
	strt <- sprintf("The parameters for simulation:\n");
	stru <- sprintf("------------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	   "Date", 		    Sys.time() );
	str1 <- sprintf("%15s: %s\n",	   "Cross",			obj.par$cross.type );
	str2 <- sprintf("%15s: %-10.0f\n", "Sample size", 	obj.par$simu.obs );
	str3 <- sprintf("%15s: %-10.0f\n", "Sample times",	length(obj.par$simu.times) );

	marker_s <- paste( cumsum( obj.par$simu.mrkdist ), collapse=",", sep="");
	str4 <- sprintf("%15s: %s\n", 	   "Marker pos", 	marker_s );
	str5 <- sprintf("%15s: %-10.0f\n", "QTL pos", 	  	obj.par$simu.qtlpos );

	str6 <- sprintf("%15s: %s\n",	   "Curve",			obj.par$curvr.type );

	str7.1 <- sprintf("%15s: %s\n", 	   "QQ", 		"" );
	strcomma <- paste( obj.par$par2, collapse=",", sep="");
	str7.2 <- sprintf("%15s: %s\n", 	   strcomma, 		"" );

	str7.3 <- sprintf("%15s: %s\n", 	   "Qq", 		"" );
	strcomma <- paste( obj.par$par1, collapse=",", sep="");
	str7.4 <- sprintf("%15s: %s\n", 	   strcomma, 		"" );

	str7.5 <- sprintf("%15s: %s\n", 	   "qq", 		"" );
	strcomma <- paste( obj.par$par0, collapse=",", sep="");
	str7.6 <- sprintf("%15s: %s\n", 	   strcomma, 		"" );

	str8 <- sprintf("%15s: %s\n",	   "Cobvariance",	obj.par$covar.type );
	strcomma <- paste( obj.par$par.covar, collapse=",", sep="");
	str8.1 <- sprintf("%15s: %s\n", 	   strcomma, 		"" );

	strd <- sprintf("------------------------------------\n\n");
	str <- paste(strt, stru, str0, str1, str2, str3, str4, str5, str6,
					str7.1, str7.2, str7.3, str7.4, str7.5, str7.6, str8, str8.1, strd, sep="" );

	return ( str );
}


#--------------------------------------------------------------
# dat.simulate
#
# Simulate a data set for backcros, including phenotype, genotype,
# marker table.
#
# Input
#       parameter object for Backcross
# output:
#       data object
#--------------------------------------------------------------

dat.simulate <- function( obj.cross, obj.curve, obj.covar, simu.mrkdist, simu.qtlpos, simu.obs, simu.times,
			par.X=NULL, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, phe.missing=0.01, marker.missing=0.01 )
{
	if( missing(par0) || is.null(par0) || missing(par1) || is.null(par1) || missing(par2) || is.null(par2) )
	{
		simu_parm <- get_simu_param( obj.curve, simu.times );
		par0 <- simu_parm[1,];
		par1 <- simu_parm[2,];
		par2 <- simu_parm[3,];
	}
	if(missing(par.covar) || is.null(par.covar) )
		par.covar <- get_simu_param( obj.covar, simu.times );

	if(length(simu.times)==1)
		simu.times <- 1:simu.times;

	obj.par <- list (
		covar.type     = obj.covar@type,
		curve.type     = obj.curve@type,
		cross.type     = obj.cross$type,
		simu.times     = simu.times,
		simu.obs       = simu.obs,
		simu.mrkdist   = simu.mrkdist,
		simu.qtlpos    = simu.qtlpos,
		cross.options  = NULL,
		par.X          = par.X,
		parQ           = par0,
		par1           = par1,
		par2           = par2,
		par.covar      = par.covar,
		phe.missing    = phe.missing,
		marker.missing = marker.missing );
	class(obj.par) <-  "FM2.par";

	geno_simu <- function()
	{
		genos.matrix <- obj.cross$get_simu_marker( obj.par$simu.obs, obj.par$simu.mrkdist, obj.cross );
		mk_name  <- c();
		mk_dist  <- c()
		mk_index <- c();
		mk_group <- c();
		mrkplace = cumsum(obj.par$simu.mrkdist);
		for (i in 1:length(obj.par$simu.mrkdist) )
		{
			mk_name  <- c( mk_name, paste('marker',i) );
			mk_dist  <- c( mk_dist, mrkplace[i])
			mk_index <- c( mk_index, 1);
			mk_group <- c( mk_group, 'G1' );
		}

		marker.table <- data.frame(Marker=mk_name, Dist=mk_dist,grp_idx=mk_index, Group=mk_group);
		#marker.obj   <- fin.get_markerobj( marker.table );

		return( list(   obj.cross    = obj.cross,
		 				marker.table = marker.table,
						#marker.obj   = marker.obj,
						genos.matrix = genos.matrix) );
	}

	obj.gen <- geno_simu();
	pheX <- NULL;
	if(length(par.X)>0)
	{
		pheX <- matrix(0, nrow=simu.obs, ncol=0 );
		rownames(pheX) <- paste("N", 1:simu.obs, sep="_");

		for(i in 1:length(par.X) )
		{
			if(i==1)
				pheX <- round( runif(simu.obs, 1, 2) )
			else
				pheX <- cbind( pheX, runif(simu.obs, -1, 1) );
		}
		colnames(pheX) <- paste("X", 1:length(par.X), sep="_");
	}


	pheY <- array( 0, dim=c( obj.par$simu.obs, length(obj.par$simu.times)));

	colnames( pheY ) <- paste("Y", obj.par$simu.times, sep="_");
	rownames( pheY ) <- paste("N", 1:obj.par$simu.obs, sep="_");

	#generate traits
	options <- list( max.time=max( obj.par$simu.times, na.rm=T), min.time=min( obj.par$simu.times, na.rm=T) );
	sim.mu   <-  get_curve( obj.curve, par0, obj.par$simu.times, options=options );
	sim.mu   <-  rbind(sim.mu, get_curve( obj.curve, par1, obj.par$simu.times, options=options ) );
	sim.mu   <-  rbind(sim.mu, get_curve( obj.curve, par2, obj.par$simu.times, options=options ) );

	get_gencode <- function(d.gen, d.snpinfo)
	{
		d.g <- array( 9, simu.obs );
		d.gen2 <- as.character(unlist(d.gen));
		snpB <- as.character(unlist(d.snpinfo[5]));
		snpA <- as.character(unlist(d.snpinfo[4]));

		QQ2<- paste(snpB, snpB, sep="");
		qq0<- paste(snpA, snpA, sep="");
		Qq1<- c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;

		d.g[which(d.gen2==QQ2)]<-2;
		d.g[which(d.gen2==qq0)]<-0;
		d.g[which(d.gen2==Qq1[1])]<-1;
		d.g[which(d.gen2==Qq1[2])]<-1;

		return(d.g);
	}

	mrkplace = cumsum(obj.par$simu.mrkdist);
	idx     <- max(which( mrkplace < obj.par$simu.qtlpos ));
	qtlmrk1 <- idx[1];
	qtlmrk2 <- idx[1]+1;
	qtl.gen <- obj.cross$get_simu_qtl( obj.par$simu.obs, obj.gen$genos.matrix[,qtlmrk1], obj.gen$genos.matrix[,qtlmrk2],
									  obj.par$simu.qtlpos, mrkplace[qtlmrk1], mrkplace[qtlmrk2], obj.par$cross.options );

	sim.covar <- get_matrix( obj.covar, par.covar, 1:length(obj.par$simu.times) );
	for (i in 1:obj.par$simu.obs)
	{
		 if (qtl.gen[i]==9) qtl.gen[i] <- round(runif(1, 0, 2));

		 pheY[i, ] <- rmvnorm(1, sim.mu[ qtl.gen[i] + 1, ], sim.covar );
		 if( !is.null(pheX) )
		 	pheY[i, ]  <- pheY[i, ] + sum( pheX[i,] * par.X );
	}
	
	obj.phe <- list();
	obj.phe$pheno_csv <- paste("Simu-", obj.curve@type, "-", obj.covar@type, ".csv", sep="");
	obj.phe$obj.curve <- obj.curve;
	obj.phe$obj.covar <- obj.covar;
	obj.phe$ids       <- rownames( pheY );
	obj.phe$pheY      <- pheY;
	obj.phe$pheX      <- pheX;
	obj.phe$pheT      <- matrix( rep( obj.par$simu.times, obj.par$simu.obs), nrow=obj.par$simu.obs, byrow=T)
	obj.phe$log       <- FALSE;

	colnames(obj.phe$pheT) <- paste("T", obj.par$simu.times, sep="_");
	rownames(obj.phe$pheT) <- rownames( pheY );

	dat <- list(obj.gen=obj.gen, obj.phe=obj.phe, obj.par=obj.par, error=F);
	class(dat) <- "FM2.dat";

	cat("Data simulation is done![QTL=", obj.par$simu.qtlpos, "]\n");

	return(dat);
}

#--------------------------------------------------------------
# public: dat.load
#
# load a real data set for backcross and F2 experiment.
#
# input:
#  file : pheno_csv, geno_csv, marker_csv
#  cross: cross type, BC or F2
#--------------------------------------------------------------
dat.load<-function( pheno_csv, time_csv, geno_csv, marker_csv, log=FALSE, head=TRUE )
{
	dat<-list(
		obj.gen = list(
			geno.csv     = geno_csv,
			marker.csv   = marker_csv,
			marker.obj   = NULL,
			marker.table = NULL,
			genos.matrix = NULL),
		obj.phe = list(
			pheno.csv    = pheno_csv,
			time.csv     = time_csv,
			sample.obs   = 0,
			log          = F,
			sample.times = NULL,
			pheY         = NULL,
			pheX         = NULL,
			pheT         = NULL) );

	tb.gen <- read.csv( file=geno_csv, sep=",", header=head, row.names=1);
	tb.phe <- read.csv( file=pheno_csv, sep=",", header=head,row.names=1);
	if(!is.null(time_csv))
		tb.time <- read.csv( file=time_csv, sep=",", header=head, row.names=1)
	else
		tb.time <- NULL;

	tb.gen <- tb.gen[ order(rownames(tb.gen)),]
	tb.phe <- tb.phe[ order(rownames(tb.phe)),]
	if(!is.null(tb.time) )
		tb.time <- tb.time[ order(rownames(tb.time)),]

   	#if the ids in two files are not consistent, it would be bad data.
   	if ( any(rownames(tb.gen) != rownames(tb.phe) ) )
   		stop("Error: the ids in phenotype file and genotype fils are not consistent!");

	if(!is.null(tb.time))
   		if ( any(rownames(tb.time) != rownames(tb.phe) ) )
	   		stop("Error: the ids in phenotype file and genotype fils are not consistent!");

	rowCheck<-function(vec)
	{
		return( all(is.na(vec))  );
	}

	tb.phe.missing <- apply(tb.phe, 1, rowCheck);
	tb.gen.missing <- apply(tb.gen, 1, rowCheck);
	missing <- unique( which(tb.phe.missing==TRUE), which(tb.gen.missing==TRUE) );

	if ( length(missing) > 0)
	{
		cat("Removing missing individuals", length(missing), ".\n");
		tb.gen <- tb.gen[ -(missing),];
		tb.phe <- tb.phe[ -(missing),];
		if(!is.null(tb.time))
			tb.time <- tb.time[ -(missing),];
	}

   	dat$obj.phe$pheY <- as.matrix(tb.phe);
   	dat$obj.gen$genos.matrix  <- as.matrix(tb.gen);

	if(is.null(tb.time))
	{
		time.str <- colnames(tb.phe);
		if ( substr(time.str[1],1,1)=="X" )
		{
			time.str<-substring(time.str, 2 );
		}
		time.std <- as.numeric( time.str );
		if ( any( is.na(time.std) ) )
			time.std <- c(1:length(time.std))
		if ( max(time.std) <=1)
			time.std <- time.std* length(time.std);

		colnames( dat$obj.phe$pheY ) <- time.std;

		dat$obj.phe$pheT <- matrix(time.std, nrow=NROW(dat$obj.phe$pheY), ncol=NCOL(dat$obj.phe$pheY), byrow=T)
	}
	else
		dat$obj.phe$pheT <- as.matrix(tb.time);

	if (log )
	{
		dat$obj.phe$pheY <- log(dat$obj.phe$pheY);
		dat$obj.phe$log <- TRUE;
	}

   	dat$obj.gen$genos.matrix <- as.matrix(tb.gen);
   	dat$obj.gen$marker.table <- read.csv(file=marker_csv,sep=",", header=TRUE, row.names=1);
   	colnames(dat$obj.gen$marker.table) <- c("Marker", "Dist", "grp_idx", "Group");

   	#dat$obj.gen$marker.obj   <- fin.get_markerobj(dat$obj.gen$marker.table);

   	dat$obj.phe$sample.obs     <- NROW( dat$obj.phe$pheY )
   	dat$obj.phe$sample.times <- c( 1:NCOL( dat$obj.phe$pheY ) )
   	dat$error  <- FALSE;
	class(dat) <- "FM2.dat";

  	return( dat );
}

#--------------------------------------------------------------
# private: dat.summar_dat
#
# summarize the important information for the simulate data or
# real data.
#
# input : data object.
#--------------------------------------------------------------
dat.summary<-function( dat )
{
	str0 <- str1 <- str2 <- str3 <- str4 <- str5 <- str6 <- str7 <- str8 <- str9  <- str10 <- str11 <- str12 <- str13 <- "";

	strt <- sprintf("The data set for FunMap model:\n");
	stru <- sprintf("-----------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	 "Date", 	   Sys.time() );

	str1 <- sprintf("%15s: %s\n", 	 "Pheno. file",	   dat$obj.phe$pheno.csv );
	str2 <- sprintf("%15s: %s\n", 	 "Time  file",	   dat$obj.phe$time.csv );
	str3 <- sprintf("%15s: %s\n", 	 "Geno. file", 	   dat$obj.gen$geno.csv );
	str4 <- sprintf("%15s: %s\n", 	 "Maker file", 	   dat$obj.gen$marker.csv );
	str5 <- sprintf("%15s: %-10.0f\n", "Sample size",  dat$obj.phe$sample.obs );
	str6 <- sprintf("%15s: %-10.0f\n", "Sample times", NCOL(dat$obj.phe$pheY) );
	str7 <- sprintf("%15s: %-10.0f\n", "Marker count", NROW(dat$obj.gen$marker.table) );
	str8 <- sprintf("%15s: %-10.0f\n", "Chr/Group count", length(unique(dat$obj.gen$marker.table$grp_idx)) );
	str9 <- sprintf("%15s: %-10.0f\n", "Marker missing", sum(is.na(dat$obj.gen$genos.matrix) || dat$obj.gen$genos.matrix==-1)/length(dat$obj.gen$genos.mattix) );
	str10 <- sprintf("%15s: %-10.0f\n", "Phenotype missing", sum(is.na(dat$obj.phe$pheY))/length(dat$obj.phe$pheY));

	str11 <- sprintf("%15s: %s\n", 	 "Cross", 	        dat$obj.cross$type);
	if(!is.null(dat$obj.curve))
		str12 <- sprintf("%15s: %s\n", 	 "Curve",       dat$obj.curve@type );
	if(!is.null(dat$obj.covar))
		str13 <- sprintf("%15s: %s\n", 	 "Covariance", 	dat$obj.covar@type);

	strd <- sprintf("------------------------------------\n\n");

	str <- paste(strt, stru, str0, str1, str2, str3, str4, str5, 
	             str6, str7, str8, str9, str10, str11, str12, str13, strd, sep="" );

	return (str);
}

#--------------------------------------------------------------
# dat.plot
#
# Plot the Pharmacology curve for data object.
# Used by plot( EXP.dat object );
#--------------------------------------------------------------
dat.plot<- function( dat, plot_type=NULL, pdf_file=NULL  )
{
	if (is.null(pdf_file))  X11() else pdf(pdf_file);

	# Figure 1
	if( is.null(plot_type) || plot_type==2 )
	{
		err.fig2 <- try( fpt.plot_overlapping_curves( dat$obj.phe$pheY, dat$obj.phe$pheT, dat$obj.phe$pheX ) );
		if ( class(err.fig2) != "try-error" )
			title(paste("The ", dat$obj.curve@type," for all individuals.", sep=""));
	}

	# Figure 2
	if( is.null(plot_type) || plot_type==1)
	{
		err.fig1 <- try( fpt.plot_tiled_curves( dat$obj.phe$pheY, dat$obj.phe$pheT, dat$obj.phe$pheX, max_curves=8*8 ) );
		if ( class(err.fig1) != "try-error" )
			title(paste("The ",dat$obj.curve@type," for all individuals.", sep=""));
	}

	if (!is.null(pdf_file))
	{
		dev.off()
		cat( paste( "* The figures for all individuals are saved to ", pdf_file, ".\n", sep=""));
	}	
}

#--------------------------------------------------------------
# public: fin.get_traits_mu
#
#--------------------------------------------------------------
fin.get_traits_mu<-function( par.curve, obj.curve )
{
	mu    <- c();
	if (!is.null(par.curve$qq0) )
		mu    <- rbind( mu, obj.curve$get_mu( par.curve$qq0, par.curve$simu.times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par.curve$simu.times) ) );

	if (!is.null(par.curve$Qq1) )
		mu    <- rbind( mu, obj.curve$get_mu( par.curve$Qq1, par.curve$simu.times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par.curve$simu.times) ) );

	if (!is.null(par.curve$QQ2) )
		mu    <- rbind( mu, obj.curve$get_mu( par.curve$QQ2, par.curve$simu.times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par.curve$simu.times) ) );

	return(mu);
}

#--------------------------------------------------------------
# private: fin.recalc_marker
#
# make a marker object according to maker table of the data set
#
# input : marker table
# output: marker object
#--------------------------------------------------------------
fin.get_markerobj<-function( marker.table )
{
	marker_obj<-list(
			count = NULL,
			grps  = list() );

	marker_grp<-list(
			index     = -1,
			id        = "",
			count     = NULL,
			start_idx = NULL,
			dists     = list(),
			names     = list() );

	grp_idx <- -1;
	index   <- 1;
	index2  <- 1;
	dists   <- list();
	names   <- list();

	for (i in 1:length(marker.table[,3]))
	{
		if (marker.table[i,3]>grp_idx)
		{
			grp_idx <- marker.table[i,3];

			marker_obj$grps[[index]]       <- list();
			marker_obj$grps[[index]]$index <- marker.table[i,3];
			marker_obj$grps[[index]]$id    <- marker.table[i,4];
			marker_obj$grps[[index]]$start_idx <- i;

			if (index>1)
			{
				marker_obj$grps[[index-1]]$dists <- dists;
				marker_obj$grps[[index-1]]$names <- names;
				marker_obj$grps[[index-1]]$count <- index2-1;
				dists  <- list();
				names  <- list();
				index2 <- 1;
			}
			index <- index+1;
		}

		dists[[index2]] <- marker.table[i,2];
		names[[index2]] <- marker.table[i,1];
		index2 <- index2+1;
	}

	marker_obj$grps[[index-1]]$dists <- dists;
	marker_obj$grps[[index-1]]$names <- names;
	marker_obj$grps[[index-1]]$count <- index2-1;
	marker_obj$count <- index-1;

	return (marker_obj);
}


