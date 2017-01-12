###############################################################
#
# New Systems Mapping Application(SysMap1)
#
# report utility
#    1). fre.report_dat
#    2). fre.report_res
#    3). fre.report_res2
#
# History:
# 12/15/2011 Version 1.1
#
###############################################################

fre.report_dat<-function( dat )
{
	if ( class(dat) != "FM2.dat" )
		stop( "Error: Not a data set for Functional Mapping." );

	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Cross:", "", dat$obj.cross$type ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Curve:", "", dat$obj.curve@type) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Covariance:", "",  dat$obj.covar@type ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Pheno. file:", "", dat$obj.phe$pheno.csv ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Time. file:", "", dat$obj.phe$time.csv ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Geno. file:",  "", dat$obj.gen$geno.csv ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Marker file:", "", dat$obj.gen$marker.csv) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample size:", "", dat$obj.phe$sample.obs) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample times:","", length(dat$obj.phe$sample.times) ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Marker count:","", NROW(dat$obj.gen$marker.table) ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Chr/Group count:","", length(unique(dat$obj.gen$marker.table$grp_idx) ) ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Phenotype missing:","", sum(is.na(dat$obj.phe$pheY))/length(dat$obj.phe$pheY) ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Marker missing:","", sum(is.na(dat$obj.gen$genos.matrix) || dat$obj.gen$genos.matrix==-1)/length(dat$obj.gen$genos.matrix) ) );

	Report.AddParaBreak();

	c1 <- call("fpt.plot_tiled_curves", dat$obj.phe$pheY, dat$obj.phe$pheT, dat$obj.phe$pheX, max_curves=8*8 );
	Report.AddFigure( c1, " ", c(5, 5)*254, left.margin=0.1*254);

	c2 <- call("fpt.plot_overlapping_curves", dat$obj.phe$pheY, dat$obj.phe$pheT, dat$obj.phe$pheX );
	Report.AddFigure( c2, " ", c(5, 5)*254, left.margin=0.1*254);

	Report.AddParaBreak();

	return();
}


fre.report_res<-function( dat, res, perm=NULL )
{
	if ( is.null(res) )
		stop( "Error: Not a result for FunMap." );
	
	pheT <- dat$obj.phe$pheT;
	if (!is.null(res$res$qtl.peaks))
	{			
		QQ_par <- Qq_par <- qq_par <- cov_par <- NULL
		LR_peaks <- res$full.res[res$qtl.peaks,,drop=FALSE];
		if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 0, sep="_") %in% colnames(LR_peaks) ) )
			QQ_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 0, sep="_"), drop=F];
		if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 1, sep="_") %in% colnames(LR_peaks) ) )
			Qq_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 1, sep="_"), drop=F];
		if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 2, sep="_") %in% colnames(LR_peaks) ) )
			qq_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 2, sep="_"), drop=F];
		if (all(paste( "H1", get_param_info(res$obj.covar,, pheT )$names, sep="_") %in% colnames(LR_peaks) ) )
			cov_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.covar, pheT )$names, sep="_"), drop=F];

		get_simple_string<-function(par.table)
		{
			if(is.null(par.table)) return("");
			unlist( apply(par.table, 1, function(par.vec)
				{
					par.str0 <- "";
					for( k in 1:length(par.vec))
						par.str0 <- sprintf("%s%.2f,", par.str0, par.vec[k] );

					if (nchar(par.str0)>20)
						par.str0 = paste( strtrim(par.str0, 20), "...", sep="");

					return( par.str0 );
				}) );
		}

		sigs_list <- data.frame(
					Grp = c(LR_peaks[,1, drop=F]),
					Pos = sprintf("%.2f", LR_peaks[,2, drop=F] ),
					LR  = sprintf("%.2f", LR_peaks[,3, drop=F] ),
					Covariance = get_simple_string(cov_par),
					Curve_QQ  = get_simple_string(QQ_par),
					Curve_Qq  = get_simple_string(Qq_par),
					Curve_qq  = get_simple_string(qq_par) );

		Report.AddTable( sigs_list,
					title = "LR profile",
					frame = "wang",
					col.names = c("Grp.", "Pos.", "LR", "Covariance", "QQ", "Qq", "qq" ),
					col.width = c(100,    100,    200,   320,    400, 400, 400  ),
					col.align = c("L",    "L",    "R",  "L",  "L", "L", "L"),
					offset.x = 50,
					max.show = 20)
		Report.AddParaBreak(50);
	}
	
	perm.cutoff.05 <- fin.find_cutoff(perm, 0.05);
	perm.cutoff.01 <- fin.find_cutoff(perm, 0.01);

	grp.idx <- unique(res$full.res[,1]);
	marker.table <- res$obj.gen$marker.table[ res$obj.gen$marker.table$grp_idx %in% grp.idx,,drop=F]

	c1 <- call("fpt.plot_qtl_map", res$full.res, marker.table, perm.cutoff.05, perm.cutoff.01 );
	Report.AddFigure( c1, "QTL profile map", c(6, 6)*254, left.margin=0.5*254)

	Report.AddParaBreak();

	return();
}

fre.report_res2<-function( dat, res )
{
	if (is.null(dat$obj.curve) )
		stop("Invalid curve type");

	if(is.null(res$qtl.peaks))
		return();
	
	perm.cutoff.05 <- fin.find_cutoff(res$obj.permu, 0.05);
	perm.cutoff.01 <- fin.find_cutoff(res$obj.permu, 0.01);

	simu_QQ <- simu_Qq <- simu_qq <- NULL;
	QQ_par <- Qq_par <- qq_par <-  NULL
	
	LR_peaks <- res$full.res[res$qtl.peaks,,drop=F];
	pheT <- dat$obj.phe$pheT;

	if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 0, sep="_") %in% colnames(LR_peaks) ) )
		QQ_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 0, sep="_"), drop=F];
	if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 1, sep="_") %in% colnames(LR_peaks) ) )
		Qq_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 1, sep="_"), drop=F];
	if (all(paste( "H1", get_param_info(res$obj.curve, pheT )$names, 2, sep="_") %in% colnames(LR_peaks) ) )
		qq_par <- LR_peaks[,paste( "H1", get_param_info(res$obj.curve, pheT )$names, 2, sep="_"), drop=F];
	
	if(NROW(LR_peaks)>0)
		for(i in 1:NROW(LR_peaks) )
		{
			str <- sprintf("Groupp=%d, Postion=%.2f", LR_peaks[i,1], LR_peaks[i,2] );
			Report.AddHeadline( str, level=2 );
	
			grp.idx <- unique(res$full.res[,1]);
			marker.table <- res$obj.gen$marker.table[ res$obj.gen$marker.table$grp_idx %in% grp.idx,,drop=F]
			c1 <- call("fpt.plot_qtl_pos", res$full.res, marker.table,  LR_peaks[i,1], LR_peaks[i,2], cutoff.05=perm.cutoff.05, cutoff.01=perm.cutoff.01 );
			Report.AddFigure( c1, "", c(3.5, 3.5)*254, left.margin=0*254);
	
			c2 <- call("fpt.plot_com_curve", max(res$obj.phe$pheT, na.rm=T),
									round(max(res$obj.phe$pheT, na.rm=T) *1.1),
									res$obj.curve,
									pheY   = res$obj.phe$pheY,
									pheT   = res$obj.phe$pheT,
									QQ_par = QQ_par[i,],
									Qq_par = Qq_par[i,],
									qq_par = qq_par[i,],
									xlab="Time", ylab="Model");
	
			Report.AddFigure( c2, "", c(3.5, 3.5)*254, left.margin=0*254);
	
			Report.AddParaBreak();
		}
	
	return();
}

fre.report_perm<-function(perm)
{
	if ( is.null(perm) )
		stop( "Error: Not a permutation result for FunMap." );

	p.pv <- c( 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001, 0.00001);
	p.cut <- round( perm$pv.table[ match(p.pv, perm$pv.table[,1]), 2], digits = 2 );

	p.df <- data.frame(N0="Cutoff", N1=p.cut[1], N2=p.cut[2], N3=p.cut[3], N4=p.cut[4], N5=p.cut[5], N6=p.cut[6], N7=p.cut[7])

	Report.AddTable( p.df,
		 	title = "Permutation Cutoff",
			frame = "wang",
			col.names = c("pv", "0.1", "0.05", "0.01", "0.005", "0.001", "0.0001", "0.00001" ),
			col.width = c(100,   150,   150,    150,    150,      150,      150,      150  ),
			col.align = c("L",   "R",   "R",    "R",    "R",      "R",       "R",     "R"),
			offset.x = 50,
			max.show = 20)

	Report.AddParaBreak();

	c2 <- call("fpt.plot_permutation", perm$pv.table);

	Report.AddParaBreak();

	Report.AddFigure( c2, "", c(3.5, 3.5)*254, left.margin=0*254);

	Report.AddParaBreak();

	return();
}