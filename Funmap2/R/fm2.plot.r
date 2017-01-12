#################################################################
#
# New Systems Mapping Application(SysMap1)
#
# plot utility
#
#    1) plot_qtl_map
#    2) plot_qtl_pos
#    3) plot_tiled_curves
#    4) plot_overlapping_curves
#	 5) plot_com_curve
#    6) plot_permutation
#    7) plot_correlation
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# fpt.plot_to_doc
#   5 file formats are supported.
#
#   1) PDF
#   2) PS
#   3) WMF
#   4) JPG
#   5) PNG (default)
#--------------------------------------------------------------
fpt.plot_to_doc<-function( filename, doctype, h=800, w=800)
{
	strFile <- "";

	if (doctype=="pdf")
	{
		strFile <- paste(filename, FM2.seqid(), "pdf", sep=".");
		pdf( strFile );
	}
	else if (doctype=="jpg")
	{
		strFile <- paste(filename, FM2.seqid(), "jpg", sep=".");
		jpeg( strFile );
	}
	else if (doctype=="wmf")
	{
		strFile <- paste(filename, FM2.seqid(), "wmf", sep=".");
		win.metafile( strFile );
	}
	else if (doctype=="ps")
	{
		strFile <- paste(filename, FM2.seqid(), "ps", sep=".");
		postscript( strFile );
	}
	else if (doctype=="png")
	{
		strFile <- paste(filename, FM2.seqid(), "png", sep=".");
		png( strFile, width=300*9, height=300*9, res=300 );
	}
	else
	{
		strFile <- paste(filename, FM2.seqid(), "png", sep=".");
		png( strFile, height=h, width=w);
	}
	return (strFile);
}

#--------------------------------------------------------------
# plot_qtl_pos
#
# draw the profile of the likelihood ratio(LR) for for one
# chromosome, not for all chromosoms.
#
# In this plot, the marker information will be given at the
# bottom of curve.
#
# Input:
#     index: chromosome index
#       dat: data object, the marker information will be used.
#       res: LR table,like as [,1]=index, [,2]=QTL pos. [,3]=LR
#       cutoff05: the significance levels P=0.05
#       cutoff01: the significance levels P=0.01
#--------------------------------------------------------------
fpt.plot_qtl_pos<-function( qtl.map.table, marker.table, index, qtl.pos=NULL, cutoff.05=NULL, cutoff.01=NULL )
{
	chr_logs <- qtl.map.table[, c(1,2,3) ];
	chr_logs <- chr_logs[ chr_logs[,1]==index,]
	marker_list <- marker.table[which(marker.table$grp_idx==index), ];

	fpt.internal_plot_qtl_pos( chr_logs, marker_list, cutoff.05=cutoff.05, cutoff.01=cutoff.01, qtl_ps=qtl.pos );
}

#--------------------------------------------------------------
# internal_plot_qtl_pos
#
#--------------------------------------------------------------
fpt.internal_plot_qtl_pos<-function(chr_logs, marker_list=NA, threshold=NA, cutoff.05=NULL, cutoff.01=NULL, qtl_ps=NULL)
{
	op <- par(mar=c( 0, 0, 1, 0.2), bty="n" );

	plot( chr_logs[,2], chr_logs[,3], xlim=c( -10, 90), ylim=c(-25, 75),
			type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");

	#========Draw LR Curve in a range(-10, 0, 90, 75)========
	# Border(0,0,85,75);
	rect(0, 0, 85, 75, col = c(NA,"midnightblue"));

	max_xlim <- max( chr_logs[,2] );
	min_xlim <- min( chr_logs[,2] );
	max_ylog <- max( chr_logs[,3] )*1.2;
	#min_ylog <- min( chr_logs[,3] );
	min_ylog <- 0;

	if (!any(is.na(marker_list))  )
	{
		max_xlim <- max( marker_list[,2] );
		min_xlim <- min( marker_list[,2] );
	}

	#Draw Curve in (0,0,85,75)
	lines( chr_logs[,2]/(max_xlim - min_xlim)*85, (chr_logs[,3]- min_ylog)/(max_ylog - min_ylog)*75, lwd=2, pch=19);
	points( chr_logs[,2]/(max_xlim - min_xlim)*85, (chr_logs[,3]- min_ylog)/(max_ylog - min_ylog)*75,type = "p", pch=21);

	#Draw QTL position
	if ( !is.null(qtl_ps) )
		for(i in 1:length(qtl_ps) )
		{
			idxs <- which(chr_logs[,2] == qtl_ps[i] );
			if (length(idxs)<1)
				next;

			qtl_pv <- chr_logs[ idxs[1], 3 ] ;
			segments(
				qtl_ps[i]/(max_xlim - min_xlim)*85, 0,
				qtl_ps[i]/(max_xlim - min_xlim)*85,
				qtl_pv/(max_ylog - min_ylog)*75, lwd=2,lty="dotted");

			points( qtl_ps[i]/(max_xlim - min_xlim)*85,
					chr_logs[ idxs[1], 3 ]/(max_ylog - min_ylog)*75,
					type = "p", pch=19);
		}

	if ( !any(is.na(marker_list)) )
		for (j in 2:length(marker_list[,1])-1 )
		{
			x0 <- marker_list[j,2]/(max_xlim - min_xlim)*85;
			segments( x0, 0, x0, 75, lty="dotted");
		}


	if (!is.null(cutoff.05))
		segments( 0, cutoff.05/(max_ylog - min_ylog)*75, 85, cutoff.05/(max_ylog - min_ylog)*75, lwd=1, lty="22", col="red");

	if (!is.null(cutoff.01))
		segments( 0, cutoff.01/(max_ylog - min_ylog)*75, 85, cutoff.01/(max_ylog - min_ylog)*75, lwd=1, lty="22", col="orange");

	x_temp <- (max_xlim - min_xlim)/6;
	x_temp <- 10*ceiling(x_temp/10);

	for (j in 1:(max_xlim%/%x_temp) )
	{
		if( j*x_temp>=min_xlim && j*x_temp<=max_xlim)
		{
			segments( j*x_temp/(max_xlim - min_xlim)*85, 0, j*x_temp/(max_xlim - min_xlim)*85, 1 );
			segments( j*x_temp/(max_xlim - min_xlim)*85, 75-1, j*x_temp/(max_xlim - min_xlim)*85, 75 );
			text( j*x_temp/(max_xlim - min_xlim)*85, -0.6, (j*x_temp), adj=c(0.5, 1), cex=0.8 );
		}

	}

	text( min_xlim/(max_xlim - min_xlim)*85, -0.3, "0", adj=c(1, 0.5) );


	y_temp <- (max_ylog - min_ylog)/12;
	y_temp <- 2*ceiling(y_temp/2);
	for (j in 1:((max_ylog-min_ylog)%/%y_temp) )
	{
		if( j*y_temp >= min_ylog && j*y_temp <= max_ylog )
		{
			segments( 0, j*y_temp/(max_ylog - min_ylog)*75, 1, j*y_temp/(max_ylog - min_ylog)*75, 1 );
			segments( 85, j*y_temp/(max_ylog - min_ylog)*75, 84, j*y_temp/(max_ylog - min_ylog)*75, 1 );
			text( -1, j*y_temp/(max_ylog - min_ylog)*75, (j*y_temp), adj=c(1, 0.5), cex=0.8 );
		}

	}
	
	sticker_h <- 0;
	if ( any(is.na(marker_list)))
	{
		for (j in 2:length(chr_logs[,2]) )
		{
			segments(chr_logs[j,2], min_ylog, chr_logs[j,2], min_ylog- sticker_h*3);
		}
   	}

	#======= Draw marker ruler in a range(-10,-25,90,25)=======
	segments( 0, -17, 85, -17, lwd=2, col="black");
	draw_last_x <- -100;
	for (j in 1:length(marker_list[,1]) )
	{
		x0 <- marker_list[j,2]/(max_xlim - min_xlim)*85;
		segments( x0, -17, x0, -16);
		if ( x0 - draw_last_x > strheight(marker_list[j,3], srt=90, cex=0.8) )
		{
			text( x0, -16, marker_list[j,3], adj=c(0,0.5), srt=90, cex=0.8);
			draw_last_x <- x0;
		}
	}

	draw_last_x <- -100;
	for (j in 1:( length(marker_list[,1])-1) )
	{
		x0 <- 0.5*(marker_list[j+1,2] + marker_list[j,2])/(max_xlim - min_xlim)*85;
		segments( x0, -17, x0, -18);

		str<- sprintf("%6.2f", (marker_list[j+1,2]-marker_list[j,2]));
		if ( x0 - draw_last_x > strheight(str, srt=90, cex=0.8) )
		{
			text( x0, -18, str, adj=c(1,0.5), srt=90,cex=0.8 );
			draw_last_x <- x0;
		}
	}

	text( -8, 35, "LR2", adj=c(1,0.5), srt=90,cex=1.2);

	par(op);
}

#--------------------------------------------------------------
# fpt.plot_qtl_map
#
# This plot will draw a map of the LR profilefor full chromosomes.
# It can layout the chromosomes automatically. If the layout need to
# be customized,  plot_chromosom_map is better than this one
#
# Input:
#      dat:  data object, the marker information will be used.
#      res:  LR table.
#--------------------------------------------------------------
fpt.plot_qtl_map<-function( qtl.map.table, marker.table, cutoff.05=NULL, cutoff.01=NULL )
{
	chr_nums <- max( marker.table[, c("grp_idx")] );
	x <- (chr_nums)^0.5;
	level_cnt <- floor(x);
	if(chr_nums<4) level_cnt <-1;

	chr_logs <- qtl.map.table[, c(1,2,3) ];
	marker_list<- marker.table[, c("grp_idx", "Dist", "Marker")];

	fpt.internal_plot_qtl_map( chr_nums, level_cnt, chr_logs, marker_list=marker_list, cutoff.05=cutoff.05, cutoff.01=cutoff.01 );
}

#--------------------------------------------------------------
# internal_plot_qtl_map
#
#--------------------------------------------------------------
fpt.internal_plot_qtl_map<-function(chr_nums, level_cnt, chr_logs, marker_list=NULL, cutoff.05=NULL, cutoff.01=NULL )
{
	max_log <- 0;
	min_log <- 0;
	#Chromosome Length;max LR2; markernumber
	ch_ev <- matrix(0, nrow=chr_nums, ncol=3);

	for ( i in  1:length(chr_logs[,1]) )
	{
		n <- chr_logs[i,1];
	    if ( chr_logs[i,2]>ch_ev[n,1] )
			ch_ev[n,1] <- chr_logs[i,2];

	    if ( chr_logs[i,3]>ch_ev[n,2] )
			ch_ev[n,2] <- chr_logs[i,3];

	    if (chr_logs[i,3]>max_log)
			max_log <- chr_logs[i,3];

	    if (chr_logs[i,3]<min_log)
			min_log <- chr_logs[i,3];
	    ch_ev[n,3] <- ch_ev[n,3] + 1;
	}

	nRowChrs <- array(1, level_cnt);
	nRowChrs[1]<- chr_nums - (level_cnt-1)
	delt_sd <- Inf;
	delt_old_sd<-0;
	if (level_cnt>1)
	{
		while(delt_sd>10)
		{
			for(i in 1:(level_cnt-1))
			{
				dm <- fpt.get_SmallestVar2( ch_ev, nRowChrs, i );
				nRowChrs[i]<-dm[2];
				nRowChrs[i+1]<-dm[3];
			}

			dm <- fpt.get_deltvar( ch_ev, nRowChrs, 1 , level_cnt );
			delt_sd <- abs( dm[1]-delt_old_sd );
			delt_old_sd<-dm[1];
		}
	}
	else
	{
		dm<-c(0, sum(ch_ev[,1]), 0 )
	}

	pos <- matrix(0, nrow=chr_nums, ncol=4);
	x0 <- 0;
	h  <- 1/level_cnt;
	y0 <- 1-h;

	nChr <- 1;
	for (i in 1:level_cnt)
	{
   		for ( j in 1: nRowChrs[i] )
    		{
    			pos[nChr,1] <- x0;
    			pos[nChr,2] <- y0;
    			pos[nChr,3] <- ch_ev[nChr,1]/dm[i+1];
    			pos[nChr,4] <- h;
    			if (pos[nChr,3] < 0.001) pos[nChr,3] <- 0.001;

    			x0 <- x0 +pos[nChr,3];
    			nChr <- nChr+1;
   		}
   		x0 <- 0;
   		y0 <- y0-h;
	}

	op <- par(mar=c(5, 4, 3, 2), mgp=c(1.8, 0.5, 0.00 )  );

	plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=T, xlab="",ylab="", xlim=c(0, 1000), ylim=c(0, 1000), bty="o", xaxs="i", yaxs="i" );
	title(xlab="", ylab=" -2 Log Likelihood Ratio", cex.lab=1.2);

	for (i in 1:chr_nums )
	{
		sub_plot <- c(pos[i,1], pos[i,1]+pos[i,3], pos[i,2], pos[i,2]+pos[i,4]);
		sub_rc <- c(sub_plot[1], sub_plot[3], sub_plot[2], sub_plot[4])*1000;

		rect( sub_rc[1], sub_rc[2], sub_rc[3], sub_rc[4], border="black");

		ll <- matrix(0, nrow=ch_ev[i,3], ncol=2 );
		n <-1;
		for (j in 1:length(chr_logs[,1]) )
		{
			if ( chr_logs[j,1] == i )
			{
				ll[n,1] <- chr_logs[j,2];
				ll[n,2] <- chr_logs[j,3];
				n <- n+1;
			}
		}

		xlim.cur <- c(0, ch_ev[i,1] );
		ylim.cur <- c(min_log*0.9, max_log*1.3);

		x0 <- (ll[,1])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
		y0 <- (ll[,2])*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
		lines( x0, y0, lwd=1);

		x0 <- min(ll[,1])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
		x1 <- max(ll[,1])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];

		if(!is.null(cutoff.05))
		{
			y0 <- cutoff.05*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
			if( y0 <= sub_rc[4] && y0 >= sub_rc[2] )
				segments( x0, y0, x1, y0, col="red", lwd=1, lty="22");
		}

		if(!is.null(cutoff.01))
		{
			y0 <- cutoff.01*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
			if( y0 <= sub_rc[4] && y0 >= sub_rc[2] )
				segments( x0, y0, x1, y0, col="orange", lwd=1, lty="22");
		}

		#marker
		sticker_h<- ( max_log*1.3- min_log*0.9)/20;
		if ( !is.null(marker_list))
		{
			for (j in 2:length(ll[,1]) )
			{
				#segments( ll[j,1], min_log*0.9, ll[j,1], min_log*0.9+ sticker_h);
				x0 <- (ll[j,1])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
				y0 <- (min_log*0.9)*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
				x1 <- (ll[j,1])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
				y1 <- (min_log*0.9+ sticker_h)*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
				segments( x0, y0, x1, y1 );
			}
    		}
		else
		{
			for (j in 1:length(marker_list[,1]) )
			{
				if (marker_list[j,1]==i)
				{
					x0 <- (marker_list[j,2])*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
					y0 <- (min_log*0.9)*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
					y1 <- (min_log*0.9+ sticker_h)*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];
					segments( x0, y0, x0, y1 );

				}
			}
		}

		#chrom no
		x0 <- (5)*(sub_rc[3]-sub_rc[1])/(xlim.cur[2]-xlim.cur[1]) + sub_rc[1];
		y0 <- (max_log*1.2)*(sub_rc[4]-sub_rc[2])/(ylim.cur[2]-ylim.cur[1]) + sub_rc[2];

		text(x0, y0, paste("",i) , font=4);
	}


	ticker <-  round(seq.int(ylim.cur[1], ylim.cur[2], length.out=5)[1:2]);
	ticker5 <- seq.int(ticker[1], ylim.cur[2], ticker[2]-ticker[1])

	for(i in 1:level_cnt)
	{
		if(i==level_cnt)
		{
			y0 <- ticker5[1:5]/(ylim.cur[2]-ylim.cur[1])*(1/level_cnt) + (i-1)*(1/level_cnt)
			axis(2, at=y0*1000, labels=ticker5[1:5], cex.axis=0.9);
		}
		else
		{
			y0 <- ticker5[1:4]/(ylim.cur[2]-ylim.cur[1])*(1/level_cnt) + (i-1)*(1/level_cnt)
			axis(2, at=y0*1000, labels=ticker5[1:4], cex.axis=0.9);
		}
	}

	par(op);
	return(pos);
}

fpt.get_SmallestVar2<-function( ch_ev, nRowChrs, iStart)
{
	nSum<-sum(nRowChrs[c(iStart,iStart+1)])
	nRowBack <- nRowChrs;
	x1<-1;
	x2<-nSum - 1
	nRowBack[iStart]<-x1;
	nRowBack[iStart+1]<- x2;

	delt_var <- fpt.get_deltvar(ch_ev, nRowBack, iStart , 2);
	delt.set <- c( 2:(nSum-1) )
	if (length(delt.set)==0) delt.set<-c(2);

	for ( i in delt.set)
	{
		nRowBack[iStart]   <- i;
		nRowBack[iStart+1] <- nSum-i;
		var2 <- fpt.get_deltvar(ch_ev, nRowBack, iStart , 2);
		if (var2[1]<delt_var)
		{
			delt_var<-var2[1];
			x1<-i;
			x2<-nSum-i;
		}
	}

	return( c(delt_var[1], x1, x2))
}

fpt.get_deltvar<-function( ch_ev, nRowChrs, iStart , iCount )
{
	nWidths <- array(0, iCount );

	index<-1;
	ns<-1;

	for ( i in 1:length(nRowChrs) )
	{
		ne <- nRowChrs[i]+ns-1;
		if ( (i>=iStart) && i<(iStart+iCount) )
		{
			nWidths[index] <- sum( ch_ev[ c(ns:ne), 1] );
			index<-index+1;
		}
		ns <- ns + nRowChrs[i];
	}

	return( c( var(nWidths), nWidths ) );
}


#--------------------------------------------------------------
# ftp.plot_tiled_curves
#
# Draw the growth curve of the individuals in some tiled plots,
# the plot's layout can be automatically calculated.
#
# Input:
#       pheY : data object( pheY );
#--------------------------------------------------------------
fpt.plot_tiled_curves<-function( pheY, pheT, pheX=NULL, rows=NA, cols=NA, max_curves = NULL, selected=NULL )
{
	max_log <- 0;
	min_log <- 0;

	if (is.null(max_curves))
		max_curves <- min(100, NROW(pheY));

	if ( is.null(selected) )
		selected <- c(1:max_curves)

	selected.greater <- which( selected > NROW(pheY) );
	if (length (selected.greater) )
		selected <- selected[- selected.greater ];

	if ( is.na (rows) )
	{
		cols <- ceiling(  (length(selected))^0.5 );
		rows <- ceiling( (length(selected))/cols);
	}

	minv <- min(as.matrix(pheY[selected,]), na.rm=TRUE)*0.9;
	maxv <- max(as.matrix(pheY[selected,]), na.rm=TRUE)*1.1;

	xmax <- max(pheT , na.rm=TRUE);
	xmin <- min(pheT , na.rm=TRUE);

	py <- length ( pheY[1, ] );
	rx <- ceiling( py/8 );
	xunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
	rx <- xunit[min(which(xunit>rx))];
	if (rx>2000) rx<-as.integer(py/8);

	ry<- ceiling((maxv - minv)/8)
	yunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
	ry <- yunit[min(which(yunit>ry))];
	if (ry>2000) ry<-as.integer((maxv - minv)/8);

	p.width <- (xmax-xmin)*cols/9*10;
	p.height <- (maxv-minv)*rows/9*10;

	plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=T, xlab="",ylab="", xlim=c(0, p.width), ylim=c(0, p.height) );
	for (i in 0:(rows-1) )
	{
		for (j in 0:(cols-1) )
		{
			sub_fig<- c( 1/cols*j, 1/cols*(j+1), 1/rows*i,  1/rows*(i+1) )*0.98+0.02;
			sub_rc <- c( p.width * sub_fig[1], p.height* sub_fig[3],
			             p.width * sub_fig[2], p.height* sub_fig[4] );
			rect(sub_rc[1], sub_rc[2], sub_rc[3], sub_rc[4], border="black", lwd=1);
			sub_fc <- c( sub_rc[1]-xmin, sub_rc[2]-minv );

			if ( (i*cols+j+1) <= length(selected) )
			{
				idx <- i*cols+j+1;
				px<-as.numeric( pheT[selected[idx], ] );
				py<-as.numeric( pheY[selected[idx], ] );
				if (length(which(py==-1))>0)
				{
					px<-px[-(which(py==-1))];
					py<-py[-(which(py==-1))];
				}

				if (length(which(is.na(py)))>0)
				{
					px<-px[-(which(is.na(py)))];
					py<-py[-(which(is.na(py)))];
				}

				s01 <- smooth.spline(px, y=py);
				xx1 <- seq(xmin,xmax,by=0.2)

				#plot(c(1:2), c(1:2),xlim=c(xmin, xmax), ylim=c(minv, maxv),
				#	type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i", );

				if (.RR("curve_smooth", def=0) == 1)
					lines( predict(s01,xx1) + sub_fc, lwd=1)
				else
					lines(px+sub_fc[1], py+sub_fc[2], lwd=1);

				for (k in 1:length(px))
					points(px[k]+sub_fc[1], py[k]+sub_fc[2] ,type="o", pch=19, cex=0.5)

				h<-( maxv - minv )
				if (i==0)
					for (k in as.integer(xmin/rx):as.integer(xmax/rx))
					{
						if ( k*rx+sub_fc[1]+xmin <= sub_rc[1] || k*rx+sub_fc[1]+xmin >= sub_rc[3]) next;
						segments(k*rx+sub_fc[1]+xmin, minv+sub_fc[2], k*rx+sub_fc[1]+xmin, minv+h/40+sub_fc[2]);
					}

				if (j==0)
			  		for (k in as.integer(minv/ry):as.integer(maxv/ry) )
					{
						if ( k*ry+sub_fc[2] <= sub_rc[2] || k*ry+sub_fc[2] >= sub_rc[4]) next;
						segments(xmin+sub_fc[1], k*ry+sub_fc[2], xmin+(px-1)/40+sub_fc[1], k*ry+sub_fc[2]);
					}
			}
			else
			{
				#plot(c(1:2), c(1:2),xlim=c(1, 10), ylim=c(minv, maxv),
		    	#		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");
			}

			if (i*cols+j+1 <= NROW(pheY))
				text( 1+(xmax-1)/5 + sub_fc[1], h*8.5/10+minv + sub_fc[2], paste(i*cols+j+1, "", sep=""), font=4);

			if (j==0)
			{
		  		for (k in as.integer(minv/ry):(as.integer(maxv/ry)-1) )
				{
					if ((k*ry)<maxv && (k*ry)>minv && (k*ry-minv)/(maxv-minv)<0.95 && (k*ry-minv)/(maxv-minv)>0.01)
						text(0.5+sub_fc[1], k*ry+ sub_fc[2], k*ry, adj=c(1, 0.5),cex=0.5);
				}
			}

			if (i==0)
			{
				for (k in as.integer(xmin/rx):as.integer(xmax/rx))
				{
					if ((k*rx)/xmax<0.9 && (k*rx)/xmax>0.05)
						text(k*rx+ sub_fc[1]+xmin, 0.5+ sub_fc[2]+minv, (k*rx), adj=c(0.5,1),cex=0.5)
				}
			}
		}
	}
}

#--------------------------------------------------------------
# ftp.plot_overlapping_curves
#
# Draw the growth curve of the individuals in one plots, the curves
# will be overlapped.
#
# Input:
#       pheY : data object( pheY );
#--------------------------------------------------------------
fpt.plot_overlapping_curves<-function( pheY, pheT, pheX=NULL, options=list() )
{
	max_log <- 0;
	min_log <- 0;

	xmax <- max(pheT, na.rm=TRUE);
	xmin <- min(pheT, na.rm=TRUE);
	minv <- floor( min(as.matrix(pheY), na.rm=TRUE) )*0.9;
	maxv <- ceiling( max(as.matrix(pheY), na.rm=TRUE) )*1.1;

	rx <- xmax/16 ;
	xunit <- c(1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000);
	rx <- xunit[min(which(xunit>rx))];
	if (rx>2000) rx<-as.integer(xmax/16);

	ry<- (maxv - minv)/25 ;
	yunit <- c(1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000);
	ry <- yunit[min(which(yunit>ry))];
	if (ry>2000) ry<-as.integer((maxv - minv)/25);

	if (!is.null(options$fixed_yticks))
	{
		minv <-0;
		ry <- ceiling( (maxv - minv)/options$fixed_yticks );
		maxv <- minv + ry * options$fixed_yticks*1.1
	}

	#par(bty="n")
	#if (is.null(options$par_mai))
	#	par(mai=c(0.5,0.5,0.5,0.5))
	#else
	#	par(mai=options$par_mai);

	op <- par(mar=c( 3, 2, 1, 0), mgp=c(1,1,0), bty="n" );

	plot(c(1:2), c(1:2),xlim=c(xmin-abs(xmax-xmin)*1/100, xmax+abs(xmax-xmin)*1/100), ylim=c(minv, maxv ),
		type="n", main="", xlab="Time", ylab="Phenotypic Values", xaxt="n", yaxt="n", frame=F); #xaxs="i", yaxs="i",

	rect(xmin-abs(xmax-xmin)*1/100, minv, xmax+abs(xmax-xmin)*1/100, (maxv-minv)*0.99+minv, col = c(NA,"black"), lwd=1);

	for (i in 1:length(pheY[,1]) )
	{
		py<-as.numeric( pheY[i, ] );
		px<-as.numeric( pheT[i, ] );
		if (length(which(is.na(py)))>0)
		{
			px<-px[-(which(is.na(py)))];
			py<-py[-(which(is.na(py)))];
		}

		s01<-smooth.spline(px, y=py );
		xx1<-seq(xmin,xmax,by =0.2 )

		if (.RR("curve_smooth",def=0) == 1)
		{
			pt<-predict(s01,xx1);
			lines(pt$x, pt$y, lwd=1, col="gray")
		}
		else
			lines(px, py, lwd=1, col="gray")


		for (k in 1:length(px))
			points(px[k], py[k], type="o", pch=19, cex=0.5,col="darkgray")
    }

	k <- c( as.integer(minv/ry):as.integer(maxv/ry) )
	axis(side=2, at=k*ry, col="black", col.axis="black", pos=xmin-abs(xmax-xmin)*1/100, col.lab="black", font.axis=1, font= 1, cex.axis = ifelse(is.null(options$cex_axis), 1, options$cex_axis)  )

	k <- c( as.integer(xmin/rx):as.integer(xmax/rx) )
	axis(side=1, at=k*rx, col="black", col.axis="black", pos=minv, col.lab="black", font.axis=1, font= 1, cex.axis = ifelse(is.null(options$cex_axis), 1, options$cex_axis), padj=0.5 )


	#if ( dat.phe$curve.type == "CM" )
	#{
	#	par <- CM.get_est_value(dat);
	#	x10 <- seq(0,xmax,0.1);
	#	y10  <- CM.get_model_mu(x10, par[3:length(par)], tmin=1, tmax=xmax);
	#	lines(x10, y10, lwd=2, col="black" , lty="solid");
	#
	#	points(xmin:xmax, colMeans(dat.phe$pheY), type="o", pch=19, cex=0.5, col="red" );
	#}

	par(op);

}

#--------------------------------------------------------------
# ftp.plot_com_curve
#
# Draw the np curve for QQ, Qq, qq genotypes. the paramter of
# the pharmacology curve is E0, E50, Emax.
#
# Input:
#     nMesa : Mesuared times
#     nLong : total times
#     data  : data
#     QQ_Par: E0, E50 and Emax for QQ type
#     Qq_Par: E0, E50 and Emax for Qq type
#     qq_Par: E0, E50 and Emax for qq type
#--------------------------------------------------------------

fpt.plot_com_curve<-function( nMesa, nLong, obj.curve, pheY = NULL, pheT=NULL, QQ_par=NULL, Qq_par=NULL, qq_par=NULL, simu_QQ=NULL, simu_Qq=NULL, simu_qq=NULL, xlab="Time", ylab="Model" )
{
	if ( nMesa > nLong)
	 	nLong <- nMesa + 4;

	limit1 <- 10;
	if (!is.null(pheY)) limit1 <- max( pheY, na.rm=TRUE );

	if (!is.null(QQ_par))
		limit1 <- max( limit1, get_curve( obj.curve, QQ_par, 1:nMesa, options=list(tmin=QQ_par[1], tmax=nMesa ) ) );
	if (!is.null(Qq_par))
		limit1 <- max( limit1, get_curve( obj.curve, Qq_par, 1:nMesa, options=list(tmin=Qq_par[1], tmax=nMesa ) ) );
	if (!is.null(qq_par))
		limit1 <- max( limit1, get_curve( obj.curve, qq_par, 1:nMesa, options=list(tmin=qq_par[1], tmax=nMesa ) ) );

	op <- par(mar=c(2,3,1,1), bty="o");


	plot(c(1:2), c(1:2),xlim=c(1, nLong+1), ylim=c(0, limit1*1.1),
		type="n", main="",xlab=xlab, ylab=ylab, xaxt="s", yaxt="s", xaxs="i", yaxs="i", frame=T);

	if (!is.null(pheY))
	{
		tLen <- length( pheY[1,] );
		for (i in 1:length(pheY[,1]) )
		{
			#lines(c(1:tLen), pheY[i,], col=rgb(0.75,0.75,0.75), lwd=1);

			py <- as.numeric( pheY[i,] );
			px <- as.numeric( pheT[i, ]  ) ;

			if (length(which(py==-1))>0)
			{
				px<-px[-(which(py==-1))];
				py<-py[-(which(py==-1))];
			}

			lines(px, py, col=rgb(0.75,0.75,0.75), lwd=1);
		}
	}

	x10 <- seq(0,nLong,0.1)	;
	p0  <- which(x10>=nMesa)[1];
	p1  <- length(x10);

	if (!is.null(QQ_par))
	{
		y20 <- get_curve( obj.curve, QQ_par, x10, options=list(tmin=QQ_par[1], tmax=nMesa) );
		lines(x10[1:p0],  y20[1:p0], lwd=2, col="darkgreen", lty="solid");
		lines(x10[p0:p1], y20[p0:p1], lwd=2, col="darkgreen", lty="dotted");

		#t1<-as.integer(QQ_par[1]*nMesa);
		#segments(t1, 0, t1, limit1*1.1, lty="solid");
	}

	if (!is.null(simu_QQ))
	{
		y20 <- get_curve( obj.curve, simu_QQ, x10, options=list(tmin=simu_QQ[1], tmax=nMesa)) ;
		lines(x10[1:p0],  y20[1:p0], lwd=2, col="black", lty="dotted");
		lines(x10[p0:p1], y20[p0:p1], lwd=2, col="black", lty="dotted");
	}


	if (!is.null(Qq_par))
	{
		y10 <- get_curve( obj.curve, Qq_par, x10, options=list(tmin=Qq_par[1], tmax=nMesa) );
		lines(x10[1:p0],  y10[1:p0], lwd=2, col="blue" , lty="solid");
		lines(x10[p0:p1], y10[p0:p1], lwd=2, col="blue" , lty="dotted");

		#t1<-as.integer(Qq_par[1]*nMesa);
		#segments(t1, 0, t1, limit1*1.1, lty="solid");
	}

	if (!is.null(simu_Qq))
	{
		y10 <- get_curve( obj.curve, simu_Qq, x10, options=list(tmin=simu_Qq[1], tmax=nMesa) );
		lines(x10[1:p0],  y10[1:p0], lwd=2, col="black", lty="dotted");
		lines(x10[p0:p1], y10[p0:p1], lwd=2, col="black", lty="dotted");
	}

	if (!is.null(qq_par))
	{
		y00 <- get_curve( obj.curve, qq_par, x10, options=list(tmin=qq_par[1], tmax=nMesa)) ;
		lines(x10[1:p0],  y00[1:p0], lwd=2, col="red" , lty="solid");
		lines(x10[p0:p1], y00[p0:p1], lwd=2, col="red" , lty="dotted");

		#t2<-as.integer(qq_par[1]*nMesa);
		#segments(t2, 0, t2, limit1*1.1, lty="solid");
	}


	if (!is.null(simu_qq))
	{
		y00 <- get_curve( obj.curve, simu_qq, x10, options=list(tmin=simu_qq[1], tmax=nMesa));
		lines(x10[1:p0],  y00[1:p0], lwd=2, col="black" , lty="dotted");
		lines(x10[p0:p1], y00[p0:p1], lwd=2, col="black" , lty="dotted");
	}

	segments( nMesa, 0, nMesa, limit1*1.1, lty="dotted", col="black");

	plotchar <-c(20,21);
	sLegend<- c();
	colors <- c();

	if (!is.null(QQ_par))
	{
		sLegend<- c( sLegend, sprintf("(QQ%i):%3.2f,%3.2f,...",2, QQ_par[1], QQ_par[2]) );
		colors <- c( colors, "darkgreen");
	}
	if (!is.null(Qq_par))
	{
		sLegend<- c( sLegend, sprintf("(Qq%i):%3.2f,%3.2f,...",1, Qq_par[1], Qq_par[2]) );
		colors <- c( colors, "blue");
	}
	if (!is.null(qq_par))
	{
		sLegend<- c( sLegend, sprintf("(qq%i):%3.2f,%3.2f,...",0, qq_par[1], qq_par[2]) );
		colors <- c( colors, "red");
	}

	if (length(sLegend)>0)
		legend(x="topright", y=NULL, sLegend, cex=0.6, col=colors,  pch=plotchar, lty=1, title="Legend")

	par(op);

}

#--------------------------------------------------------------
# ftp.plot_permutation
#
# Draw the permutation graph.
#
# Input:
#     pv.table :
#--------------------------------------------------------------
fpt.plot_permutation<-function( pv.table )
{
	op <- par(mar=c(2,3,1,1), bty="o");

	plot(c(1:2), c(1:2),xlim=c(-6,1), ylim=c(0, max(pv.table[,2]) ),
		type="n", main="",xlab="p-value", ylab="cutoff", xaxt="n", yaxt="s", xaxs="i", yaxs="i");

	title("Permutation result");
	lines( log10(pv.table[,1]), pv.table[,2], lty=1, col="blue")

	legend_left <- c();
	idx05 <- which(pv.table[,1]==0.05);
	str_legend <- c();
	legend_left <- c();

	if (length(idx05)>0)
	{
		segments(log10(0.05), 0, log10(0.05), pv.table[idx05[1], 2], lty=2 );
		str_legend <- sprintf("%-7.6f: %6.2f", 0.05, pv.table[idx05[1], 2]);
		legend_left <- c( str_legend );
	}

	sig_lines<-c( 0.01, 0.001, 0.0001, 0.00001, 0.000001);
	for (i in 1:length(sig_lines) )
	{
		idx01 <- which(pv.table[,1]==sig_lines[i]);
		if (length(idx01)>0)
		{
			segments(log10(sig_lines[i]), 0, log10(sig_lines[i]), pv.table[idx01[1], 2] , lty=2);
			str_legend <- sprintf("%-7.6f: %6.2f", sig_lines[i], pv.table[idx01[1], 2]);
			legend_left <- c( legend_left, str_legend );
		}
	}

	sig_lables<-c( 1, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001);
	axis(1, at=log10(sig_lables),labels=sig_lables, las=2);

	if (length(legend_left)>0)
	{
		str_max <- max(strwidth(str_legend));
		legend("topright", legend = legend_left,text.width =str_max , title = "Cutoffs");
	}

	par(op);
}


#--------------------------------------------------------------
# ftp.plot_correlation
#
# Draw the corelation graph.
#
# Input:
#     dat :
#     sTitle :
#     labels:
#     log10:
#     significance:
#--------------------------------------------------------------

fpt.plot_correlation<-function( dat, sTitle, labels=NULL, log10=TRUE, significance=0.05)
{
	nLenArm <- max( max(dat[,1]), max(dat[,2]) ) ;

	len <- length(dat[,1]);
	rm <- dat[,3 ];
	if (log10)
	{
		rm[ rm>1 ] <- 1;
		rm <- -log10 ( dat[,3 ] );
		rm[ rm<0 ] <- 0;
	}

	maxv <- max( rm );
	if (maxv <= -log10(0.05))
		maxv <- -log10(0.04);
	minv <- min( rm );
	minv <- 0;

	ramp <- colorRamp(c( "green", "blue" ));
	cols <- rgb ( ramp(seq(0, 1, length = 1000) ), max=255 );
	mc <- array("#FFFFFF", dim = c( nLenArm, nLenArm));
	cls <- round( (rm-minv)/(maxv-minv)*1000 ) +1;
	cls[which(cls>1000)]<-1000;

	for (n in 1:len)
		mc[ dat[n,1], dat[n,2] ] <- cols[  cls[n] ];
	mc[is.na(mc)]<-"#FFFFFF";

	par(mar=c(2,2,2,2)+0.1);
	plot(c(0, nLenArm*1.45), c(0, nLenArm*1.45), type= "n", xlab="", ylab="", xaxt="n", yaxt="n");
	for (x in 1:nLenArm)
	for (y in 1:nLenArm)
	{
		ox <- (nLenArm-x+1);
		oy <- y;
		if (oy+ox<(nLenArm+1)) next;

		x0 <- -sqrt(2)/2*(ox-1) + sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
		y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4;
		xs <- c(x0, x0+sqrt(2)/2, x0, x0-sqrt(2)/2, x0);
		ys <- c(y0, y0+sqrt(2)/2, y0+sqrt(2), y0+sqrt(2)/2, y0);

		if ( mc[x,y] != "#FFFFFF" )
			polygon(xs, ys, col=mc[x,y], border="gray", angle=-45)
		else
			polygon(xs, ys, col=mc[x,y], border="white", angle=-45);
	}


	l <- nLenArm * sqrt(2)/2;
	x0 <- nLenArm*sqrt(2)/4;
	for (i in 0:100)
	{
		rect(x0+i*(l/100), nLenArm*sqrt(2)/2*0.1,
		     x0+(i+1)*l/100, nLenArm*sqrt(2)/2*0.2,
		     col=cols[i*10+1], border=cols[i*10+1])
	}

	text( x0+0* (l/100), nLenArm*sqrt(2)/2*0.1-strheight("1"), "1",cex=0.5 );
	x50 <- round( (-log10(0.5) - minv)/(maxv-minv)*100 ) +1;
	if (x50<100)
	{
		segments(x0+x50*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x50*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x50*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".5",cex=0.5 );
	}
	x005 <- round( (-log10(0.05) - minv)/(maxv-minv)*100 ) +1;
	if (x005<100)
	{
		segments(x0+x005*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x005*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x005*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".05",cex=0.5 );
	}
	x001 <- round( (-log10(0.01) - minv)/(maxv-minv)*100 ) +1;
	if (x001)
	{
		segments(x0+x001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".01",cex=0.5 );
	}
	x0001 <- round( (-log10(0.001) - minv)/(maxv-minv)*100 ) +1;
	if (x0001<100)
	{
		segments(x0+x0001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x0001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x0001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".001",cex=0.5 );
	}
	x00001 <- round( (-log10(0.0001) - minv)/(maxv-minv)*100 ) +1;
	if (x00001<100)
	{
		segments(x0+x00001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x00001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x00001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".0001",cex=0.5 );
	}

	for (n in 1:len )
	{
		if ( dat[n,3] <= significance )
		{
			ox <- (nLenArm - dat[n,1] + 1);
			oy <- dat[n,2];
			if (oy+ox<(nLenArm+1)) next;

			x0 <- -sqrt(2)/2*(ox-1) + sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
			y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4;
			if (nLenArm<10)
			{
				str <- sprintf("%5.3f", dat[n,3]);
				text( x0, y0+sqrt(2)/2, str, col="black", srt=-45);
			}
			else
			{
				points( x0, y0+sqrt(2)/2, pch=3,cex=0.8, col="yellow");
			}
		}
	}

	if (!is.null( labels ))
	{
		for (n in 1:length(labels) )
		{
			x <- (n-1)*sqrt(2)+sqrt(2)/2;
			y <- sqrt(2)/2*nLenArm * (1.4) ;
			text(x, y-0.5, labels[n], adj=c(1, 0), srt=-90);
		}
	}


	title( sTitle );
}

