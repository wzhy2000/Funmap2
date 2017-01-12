###############################################################
#
# Funmap System
#
# Routine:
#  1) FM2.simulate
#  2) FM2.load.data
#  3) print.FM2.dat
#  4) plot.FM2.dat
#  5) FM2.qtlscan
#  6) summary.FM2.qtl.mle
#  7) print.FM2.qtl.mle
#  8) plot.FM2.qtl.mle
#  9) FM2.permutation
# 10) print.FM2.qtl.mle.perm
# 11) summary.FM2.qtl.mle.perm
# 12) plot.FM2.qtl.mle.perm
# 13) FM2.estimate.data
# 14) FM2.report
# 15) FM2.select.qtl
# 16) FM2.simu.pipe
# 17) FM2.pipe
#
# History:
# 12/17/2016 Version 2.7
#
###############################################################

cur <- function(){
    return(FG_ENV);
}

#--------------------------------------------------------------
# FM2.get.curve
#
# Search Curve object
#--------------------------------------------------------------
FM2.get.curve <- function( curve.type )
{
	obj.curve <- fg.getCurve(curve.type);
	return(obj.curve);
}

#--------------------------------------------------------------
# FM2.get.covariance
#
# Search Curve object
#--------------------------------------------------------------
FM2.get.covariance <- function( covar.type )
{
	obj.covar <- fg.getCovariance(covar.type);
	return(obj.covar);
}


#--------------------------------------------------------------
# FM2.simulate
#
# Create a simulation data set by parameter object
#-------------------------------------------------------------- 
FM2.simulate <- function( cross.type = "BC", curve.type="Logistic", covar.type="AR1", simu.mrkdist=rep(20,10), simu.qtlpos=95, simu.obs=800, simu.times=8,
			par.X=NULL, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, phe.missing=0.01, marker.missing=0.01, pdf.file =NULL )
{
	obj.curve <- fg.getCurve(curve.type);
	if (is.null(obj.curve))
		stop(paste("Invalid curve type(", curve.type, ")", sep=""));

	obj.cross <- FM2.get_cross(cross.type);
	if (is.null(obj.cross))
		stop(paste("Invalid cross type(", cross.type, ")", sep=""));

	obj.covar <- fg.getCovariance(covar.type);
	if (is.null(obj.covar))
		stop(paste("Invalid covariance type(", covar.type, ")", sep=""));

	dat <- dat.simulate( obj.cross, obj.curve, obj.covar, simu.mrkdist, simu.qtlpos, simu.obs, simu.times,
			par.X=par.X, par0=par0, par1=par1, par2=par2, par.covar=par.covar, 
			phe.missing=phe.missing, marker.missing=marker.missing );

	if( is.null(dat) || dat$error )
		stop("Failed to simulate data.");

	class(dat) <- "FM2.dat";
	if(!dat$error)
		dat <- FM2.estimate.data( dat, curve.type = curve.type, covar.type = covar.type, pdf.file=pdf.file );

	dat$obj.cross <- obj.cross;
	return(dat);
}

#--------------------------------------------------------------
# FM2.load.data
#
# Load a real data, phenotype file, genotype file, marker file
# is necessary.
#--------------------------------------------------------------
FM2.load.data <- function( pheno.csv, time.csv, geno.csv, marker.csv, cross.type, curve.type=NULL, covar.type=NULL, pdf.file=NULL, log=FALSE )
{
	obj.cross <- FM2.get_cross(cross.type);
	if (is.null(obj.cross))
		stop(paste("Invalid cross type(", cross.type, ")", sep=""));

	dat <- dat.load( pheno.csv, time.csv, geno.csv, marker.csv, log);
	if( is.null(dat) || dat$error )
		stop("Failed to load data.");

	class(dat)<-"FM2.dat";
	dat <- FM2.estimate.data( dat, curve.type=curve.type, covar.type=covar.type, pdf.file=pdf.file );

	dat$obj.cross <- obj.cross;
	return(dat);
}

#--------------------------------------------------------------
# FM2.estimate.data
#
# Estimate the parameters of curve and covariance for given
#	data det, in the case of the cross.type being known
#-------------------------------------------------------------
FM2.estimate.data <- function(dat, curve.type=NULL, covar.type=NULL, pdf.file=NULL )
{
	if( class(dat)!="FM2.dat")
		stop("Invalid data object.");

	if(is.null(curve.type)) curve.type<-"auto";
	if(is.null(covar.type)) covar.type<-"auto";

	est <- fg_dat_est( dat$obj.phe, curve.type=curve.type, covariance.type=covar.type, file.plot.pdf=pdf.file, options=list() )
	if ( est$error )
		stop("Failed to estimate the data.\n")
	else
	{
		#copy the estimate results into data object
		dat$obj.phe$est.covar <- est$est.covar;
		dat$obj.phe$est.curve <- est$est.curve;

		dat$obj.phe$summary.curve <- est$summary.curve
		dat$obj.phe$summary.covar <- est$summary.covar

		dat$obj.curve <- fg.getCurve( est$est.curve$type );
		dat$obj.covar <- fg.getCovariance( est$est.covar$type );
	}

	return(dat);
}

#--------------------------------------------------------------
# print.FM2.dat
#
# used by print( FM2.dat ) or show( FM2.dat )
#--------------------------------------------------------------
print.FM2.dat<- function( x, ... )
{
	cat( dat.summary( x ) ) ;
	invisible();
}

#--------------------------------------------------------------
# summary.FM2.dat
#
# used by summary( FM2.dat )
#--------------------------------------------------------------
summary.FM2.dat<- function( object, ... )
{
	cat( dat.summary( object ) ) ;
	invisible();
}

#--------------------------------------------------------------
# plot.FM2.dat
#
# used by plot( FM2.dat.obj )
#--------------------------------------------------------------
plot.FM2.dat<- function( x, plot.type=NULL, pdf.file=NULL, ... )
{
	dat <- x;
	if( class(dat)!="FM2.dat" )
		stop("Invalid data object.");

	dat.plot( dat, plot.type, pdf.file );
	invisible();
}

#--------------------------------------------------------------
# FM2.qtlscan
#
# Hypothesis test for data object, the test methods should be
# (10,11,12,13,14,15)
#
# 1) scan.step, default=2, an interval distance used to scan flanking
#    marker, default is 2cm.
# 2) peak.count, default=5, a number shows how many significant QTL will
#    be selected.
# 3) plot_doctype, default=pdf, the figure output type for summary command.
#
#--------------------------------------------------------------
FM2.qtlscan<-function( dat, model="MLE", grp.idx=NULL, options=list() )
{
	if(class(dat)!="FM2.dat")
		stop("Invalid data object.");

	options <- merge.default.options(options);

	ret <- NULL;
	if( model == "MLE" )
	{
		r.time <- system.time( ret<- Qtlmle.qtlscan( dat, NULL, grp.idx, options) );
		if(is.null(ret) || ret$error )
			stop("Failed to scan QTL positions.");
		
		#copy the data into result object
		ret$obj.phe   <- dat$obj.phe;
		ret$obj.gen   <- dat$obj.gen;
		ret$obj.curve <- dat$obj.curve;
		ret$obj.covar <- dat$obj.covar;
		ret$obj.cross <- dat$obj.cross;
		ret$time      <- r.time;
		
		class( ret ) <- "FM2.qtl.mle";
	}

	return (ret);
}

#--------------------------------------------------------------
# print.FM2.qtl.mle
#
# used by show( FM2.qtl.mle )
#--------------------------------------------------------------
print.FM2.qtl.mle<-function( x, ... )
{
	if(class(x)!="FM2.qtl.mle")
		stop("Invalid QTL result object.");

	cat( Qtlmle.summary( x ) );
	invisible();
}

#--------------------------------------------------------------
# summary.FM2.qtl.mle
#
# used by summary( FM2.qtl.mle )
#--------------------------------------------------------------
summary.FM2.qtl.mle<-function( object, ... )
{
	if(class(object)!="FM2.qtl.mle") 
		stop("Invalid QTL result object.");
	
	cat( Qtlmle.summary( object ) );
	invisible();
}

#--------------------------------------------------------------
# plot.FM2.qtl.mle
#
# used by plot( FM2.qtl.mle )
#--------------------------------------------------------------
plot.FM2.qtl.mle<-function( x, plot.type=NULL, pdf.file=NULL, ... )
{
	res <- x;
	if(class(res)!="FM2.qtl.mle") 
		stop("Invalid QTL result object.");
		
	Qtlmle.plot( res, plot.type, pdf.file );
	invisible();
}

#--------------------------------------------------------------
# FM2.permutation
#
# Permutation tests.
#--------------------------------------------------------------
FM2.permutation<-function( dat, res, grp.idx=NULL, options=list() )
{
	if(class(dat)!="FM2.dat") stop("Invalid data object.");
	if(class(res)!="FM2.qtl.mle") stop("Invalid scan result object.");

	options <- merge.default.options(options);
	if (options$permu.loop < 100)
	{
		#stop("Permutation loop is too few(<100).\n");
		warning("Permutation loop is too few(<100).\n");
	}

	#update result object with permutation results.
	r.time <- system.time( res$obj.permu<- permu.execute( dat, grp.idx, options$permu.loop, options$permu.filter.ratio, options$scan.step, options$n.cores ) );
	if( !is.null( res$obj.permu) )
		res$obj.permu$time <- r.time;

	class( res$obj.permu ) <- "FM2.qtl.mle.perm";

	return( res );
}

#--------------------------------------------------------------
# print.FM2.qtl.mle.perm
#
# used by show( FM2.qtl.mle.perm )
#--------------------------------------------------------------
print.FM2.qtl.mle.perm<-function( x, ... )
{
	if(class(x)!="FM2.qtl.mle.perm")
		stop("Invalid permutation result.");

	cat( permu.summary( x ) );
	invisible();
}

#--------------------------------------------------------------
# summary.FM2.qtl.mle.perm
#
# used by summary( FM.ret.perm.obj )
#--------------------------------------------------------------
summary.FM2.qtl.mle.perm<-function( object, ... )
{
	if(class(object)!="FM2.qtl.mle.perm")
		stop("Invalid permutation result.");

	cat( permu.summary( object ) );
	invisible();
}

#--------------------------------------------------------------
# plot.FM2.qtl.mle.perm
#
# used by plot( FM.ret.perm.obj )
#--------------------------------------------------------------
plot.FM2.qtl.mle.perm<-function( x, pdf.file=NULL, ... )
{
	perm <- x;
	if(class(perm)!="FM2.qtl.mle.perm")
		stop("Invalid permutation result.");

	permu.plot( perm, pdf.file );
	invisible();
}

#--------------------------------------------------------------
# public: 
#
#--------------------------------------------------------------
FM2.select.qtl<-function( res,  threshold=0.05, threshold.type="pvalue")
{
	if(class(res)!="FM2.qtl.mle") 
		stop("Invalid scan result object.");
	
	if(is.null(threshold.type) || missing(threshold.type))
		threshold.type <- "pvalue";
		
	threshold.type <- tolower(threshold.type);
	if (!(threshold.type %in% c("pvalue","count","lr")))
		stop("'threshold.type' has 3 optional values: 'pvalue', 'LR', 'count'.");
	
	if( threshold.type=="pvalue" && is.null(res$obj.permu) )	
		stop("No permutation results.");

	if(is.null(threshold) && missing(threshold) )	
		stop("No cutoff or pvlaue are used for the criterion.");
	
	count <- cutoff<- pvalue <- NULL; 
	if(threshold.type=="pvalue") pvalue <- threshold;
	if(threshold.type=="lr") count <-  threshold;
	if(threshold.type=="count") count <- threshold;
	
	res <- fin.qtl_locate( res, pvalue, cutoff, count );
	return(res);
}
#--------------------------------------------------------------
# public: FM2.report()
#
#--------------------------------------------------------------

FM2.report<-function( file.report.pdf, dat, res=NULL, options=list( debug=F ) )
{
	if(class(dat)!="FM2.dat") 
		stop("Invalid data object.");
	
	if(class(res)!="FM2.qtl.mle") 
		stop("Invalid scan result object.");

	Report.new( file.report.pdf, options );
	Report.title( "Functional Mappping Report", "FunMap Report", "http://statgen.psu.edu" );
	Report.par( "dat.file", dat$pheno_file);

	Report.AddHeadline( "Data Summary", level=1 );
	fre.report_dat(dat);

	Report.AddHeadline( "QTL Profile", level=1 );

	if ( !is.null(res) )
	{
		#output the QTL profile
		fre.report_res( dat, res) ;

		Report.AddHeadline( "QTL Position", level=1 );
		#output the QTL curve(2);
		fre.report_res2( dat, res) ;

		#output the permutation;
		if(!is.null(res$obj.permu))
			fre.report_perm( res$obj.permu )		
	}
	
	Report.Output( file.report.pdf );
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function: FM2.pipe
#
# The result is saved into [pheno_csv.rdata]
#
# Options:
#1) n.cores, default=1, the cluster count for parallel permutation.
#2) scan.step, default=2, an interval distance used to scan flanking
#   marker, default is 2cm.
#3) peak.count, default=5, a number shows how many significant QTL will
#   be selected.
#4) permu.loop, default=1000, the count of permutation loop.
#5) permu.filter.ratio, default=5, a number shows how many significant QTL will
#   be selected.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FM2.pipe<-function(pheno.csv, time.csv, geno.csv, marker.csv, cross.type, curve.type=NULL, covar.type=NULL, model="MLE", 
			grp.idx=NULL, pdf.prefix=NULL, threshold=0.05, threshold.type="pvalue", options=list() )
{
	options <- merge.default.options(options);

	dat<- FM2.load.data( pheno.csv, time.csv, geno.csv, marker.csv, cross.type, curve.type, covar.type );
	if ( is.null(dat)  || dat$error )
		stop("Failed to load the phenotype or estimate the parameters.");

	p <- FM2.pipe.rest(dat, model, grp.idx, pdf.prefix, threshold, threshold.type, options); 
	return(list(dat=p$dat, ret=p$ret));
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FM2.simu.pipe
#
# Abstract: any model, any cross by simulation data
#
# The result is saved into [LC_simu_test_XX_XX.rdata]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FM2.simu.pipe<-function( 
				cross.type = "BC", curve.type="Logistic", covar.type="AR1", simu.mrkdist=rep(20,10), simu.qtlpos=95, 
				simu.obs=800, simu.times=8, par.X=NULL, par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, 
				phe.missing=0.01, marker.missing=0.01,  threshold=0.05, threshold.type="pvalue", 
				model="MLE", pdf.prefix=NULL, options=list() )
{
	options <- merge.default.options(options);
	
	dat <- FM2.simulate( cross.type, curve.type, covar.type, simu.mrkdist, simu.qtlpos, simu.obs, simu.times, 
			par.X, par0, par1, par2, par.covar, phe.missing, marker.missing );
	if (is.null(dat) || dat$error)
		stop("Failed to generate simulation data.");
	
	grp.idx <- NULL;
	p <- FM2.pipe.rest(dat, model, grp.idx, pdf.prefix, threshold, threshold.type, options); 
	return(list(dat=p$dat, ret=p$ret));
}

merge.default.options <- function( options )
{
	default <- list( scan.step=1, peak.count=5, debug=F, n.cores=1, permu.loop=100, permu.filter.ratio=1 );
	default[names(options)] <- options;
	options <- default;	

	# debug
	if ( is.null(options$debug) || is.na(options$debug) ) options$debug <- FALSE;
	#scan.step
	if (is.null(options$scan.step) || is.na(options$scan.step) ) options$scan.step <- 1
	#peak.count
	if (is.null(options$peak.count) || is.na(options$peak.count) ) options$peak.count <- 5
	# n.cores
	if ( is.null(options$n.cores) || is.na(options$n.cores) ) options$n.cores <- 1;
	# permu.loop
	if ( is.null(options$permu.loop) || is.na(options$permu.loop) ) options$permu.loop <- 100;
	# permu.filter.ratio
	if ( is.null(options$permu.filter.ratio) || is.na(options$permu.filter.ratio) ) options$permu.filter.ratio <- 1;

	.RW("debug", options$debug);
	
	return(options);
}

FM2.pipe.rest <- function(dat, model, grp.idx, pdf.prefix, threshold, threshold.type, options)
{
	show(dat);
	if(!is.null(pdf.prefix))
		plot(dat, pdf.file=paste(pdf.prefix, ".dat.pdf", sep=""))

	ret<- FM2.qtlscan( dat, model, grp.idx, options=options );

 	show( ret );
	if(!is.null(pdf.prefix))
		try( plot( ret, pdf.file=paste(pdf.prefix, ".mle.pdf", sep="") ) )

	try( FM2.report( paste(pdf.prefix, ".rpt.pdf", sep=""), dat, ret ) )

	if ( options$permu.loop > 1)
	{
		ret <- FM2.permutation( dat, ret, grp.idx=grp.idx, options );
		show( ret$obj.permu );
		show( ret )

		ret <- FM2.select.qtl( ret, threshold, threshold.type );
		if(!is.null(pdf.prefix))
			try( plot( ret, pdf.file=paste( pdf.prefix, ".mle.pdf", sep="") ) )
	}

	if(!is.null(pdf.prefix))
		try( FM2.report( paste(pdf.prefix, ".rpt.pdf", sep=""), dat, ret ) )
		
	return(list(dat=dat, ret=ret));	
}