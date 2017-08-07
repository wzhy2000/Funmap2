#fpt.plot_genotype( dat )

fpt.plot_genotype<-function( dat )
{
	marker.table <- dat$obj.gen$marker.table;

	chr_nums <- max( marker.table[, c("grp_idx")] );
	x <- (chr_nums)^0.5;
	level_cnt <- floor(x):(floor(x)*2);
	if(chr_nums<4) level_cnt <-1;

	group <- unique(marker.table[, c("grp_idx")])
	group.id <- unique(marker.table[, c("Group")])

	ch_ev <- as.data.frame(matrix(0, nrow=chr_nums, ncol=3));
	colnames(ch_ev) <- c("Dist", "Group.idx", "Group")
	ch_ev$Group.idx <- group
	ch_ev$Group <- group.id

	for(i in 1:length(group))
	{
		group.pos <- which(marker.table[, c("grp_idx")]==group[i]);
		ch_ev[i, 1] <- c(max(marker.table[group.pos, c("Dist")]));
	}


	fpt.internal_plot_genotype( chr_nums, level_cnt, ch_ev, marker.table );
}


fpt.internal_plot_genotype<-function(chr_nums, 	Level_cnt, ch_ev, marker.table)
{

	max_log <- 0;
	min_log <- 0;

	dm.init <- c(Inf)
	nColChrs.init <- c()

	for(r in 1:length(Level_cnt))
	{
		level_cnt <- Level_cnt[r]
		nColChrs <- array(1, level_cnt);
		nColChrs[1]<- chr_nums - (level_cnt-1);
		delt_sd <- Inf;
		delt_old_sd<-0;

		if (level_cnt>1)
		{
			while(delt_sd>10)
			{
				for(i in 1:(level_cnt-1))
				{
					dm <- fg.get_SmallestVar2( ch_ev, nColChrs, i );
					nColChrs[i]<-dm[2];
					nColChrs[i+1]<-dm[3];

					## cat("nColChrs", nColChrs, "\n")
				}

				dm <- fg.get_deltvar( ch_ev, nColChrs, 1 , level_cnt );
				delt_sd <- abs( dm[1]-delt_old_sd );
				delt_old_sd<-dm[1];
			}
		}

		else
		{
			dm<-c(0, sum(ch_ev[,1]), 0 )
		}

		## cat("dm", dm, "\n");
		## cat("nColChrs.init", nColChrs.init, "\n")
		if(dm[1] < dm.init[1])
		{
			dm.init <- dm;
			nColChrs.init <- nColChrs;
		}

	}

	dm <- dm.init;
	nColChrs <- nColChrs.init;

	width.bar <- 0.3;
	y.lim <- max(dm[-1])+max(dm[-1])*0.1*(max(nColChrs)-1)
	x.lim <- length(nColChrs);

	par(mar=c( 1, 1, 1, 1) );
	plot(c(-0.05*x.lim, x.lim*0.95), c(-y.lim, y.lim*0.05),type="n",cex=0.6,axes=FALSE,xlab="",ylab="")
	rect(-0.05*x.lim, -y.lim*1.02, x.lim*0.95, y.lim*0.05);

	N.CHR <- 0;
	left.start <- 0;
	for(n.chrom in nColChrs)
	{
		top.start <- 0;
		for(i in 1:n.chrom)
		{
			sub.mar_table <- marker.table[which(marker.table[, "grp_idx"]==N.CHR+i),];
			fpt.chrom.genotype(sub.mar_table, left.start, top.start, ch_ev[N.CHR+i, , drop=F], x.lim)

			top.start <- top.start + max(dm[-1])*0.1 + ch_ev[N.CHR+i, 1]
		}
		N.CHR <- N.CHR+n.chrom;

		#left.start <- left.start+ (strwidth("G")*max(nchar(as.character(marker.table[, c("Marker")]))))*1.5
		left.start <- left.start+ 1
	}
}

fpt.chrom.genotype <- function(marker.table, left.start, top.start, ch.ev, x.lim)
{
	bar.width <- 0.015*x.lim;
	rect( left.start, -(top.start+ch.ev[1, 1]), left.start+bar.width, -top.start, col="Blue", border = "black")
	segments(left.start-0.005*x.lim, -top.start, left.start+bar.width+0.005*x.lim,  -top.start, lwd=0.8);

	text( left.start+bar.width+0.01*x.lim, -top.start, labels=sprintf( "%4s",as.character(marker.table[1, "Marker"])), adj=0, cex=0.4, font=2 );
	text( left.start-0.01*x.lim, -(top.start+ch.ev[1, 1]), pos=1, labels=ch.ev[1, c("Group")], adj=c(0.5,0.5),cex=0.7,col="black", font=3);

	dis2 <- marker.table[2, "Dist"]
	dis1 <- marker.table[1, "Dist"]

	for(i in 2:length(marker.table[,1]))
	{
		segments(left.start-0.005*x.lim, -(top.start+marker.table[i, "Dist"]), left.start+bar.width+0.005*x.lim,  -(top.start+marker.table[i, "Dist"]), lwd=0.8);

		if(strheight("G")/2<= dis2-dis1)
		{
			text( left.start+bar.width+0.01*x.lim, -(top.start+marker.table[i, "Dist"]), labels=sprintf( "%4s",as.character(marker.table[i, "Marker"])), adj=0, cex=0.4, font=2 );
			text( left.start-0.01*x.lim,  -(top.start+ (marker.table[i, "Dist"]+marker.table[i-1, "Dist"])/2), labels=sprintf("%5.1f", marker.table[i, "Dist"]- marker.table[i-1, "Dist"]), adj=1, cex=0.4 );
			dis2 <- marker.table[i+1, "Dist"];
			dis1 <- marker.table[i, "Dist"];

		}else{
			dis2 <- marker.table[i+1, "Dist"];
			dis1 <- marker.table[i-1, "Dist"]
		}

		if(i%%2 ==0)
		{
			rect(left.start , -(top.start+marker.table[i, "Dist"]), left.start+bar.width,  -(top.start+marker.table[i-1, "Dist"]), col="Gray20" );
		}
	}
}


fg.get_SmallestVar2<-function( ch_ev, nColChrs, iStart)
{
	nSum<-sum(nColChrs[c(iStart,iStart+1)])
	nColBack <- nColChrs;
	x1<-1;
	x2<-nSum - 1
	nColBack[iStart]<-x1;
	nColBack[iStart+1]<- x2;

	delt_var <- fg.get_deltvar(ch_ev, nColBack, iStart , 2);
	delt.set <- c( 2:(nSum-1) )
	if (length(delt.set)==0) delt.set<-c(2);

	for ( i in delt.set)
	{
		nColBack[iStart]   <- i;
		nColBack[iStart+1] <- nSum-i;
		var2 <- fg.get_deltvar(ch_ev, nColBack, iStart , 2);
		if (var2[1]<delt_var[1])
		{
			delt_var[1]<-var2[1];
			x1<-i;
			x2<-nSum-i;
		}
	}

	return( c(delt_var[1], x1, x2))
}


fg.get_deltvar<-function( ch_ev, nColChrs, iStart , iCount )
{
	nWidths <- array(0, iCount );

	index<-1;
	ns<-1;

	for ( i in 1:length(nColChrs) )
	{
		ne <- nColChrs[i]+ns-1;
		if ( (i>=iStart) && i<(iStart+iCount) )
		{
			nWidths[index] <- sum( ch_ev[ c(ns:ne), 1] );
			index<-index+1;
		}
		ns <- ns + nColChrs[i];
	}

	return( c( var(nWidths), nWidths ) );
}



#source("plot.marker.r")
#fpt.gene.map(dat)

fpt.gene.map <- function (dat)
{
	genos_table <- dat$obj.gen$genos.matrix;
	marker_table <- dat$obj.gen$marker.table;
	len.marker <- dim( genos_table)[2];
	len.sample <- dim( genos_table)[1];

	Group <- unique(marker_table$grp_idx);
	xlim <- len.marker + length(Group)*3;
	par(mai=c(mai=c(0.5,0.5,0.5,0.5)));
	plot(c(0,xlim),c(-len.sample*0.5,len.sample*1.5),type="n",cex=0,axes=FALSE,xlab="",ylab="");


	for(i in 1:length(Group))
	{
		pos.mar <- which( marker_table$grp_idx==Group[i]);

		geno <- genos_table[, pos.mar, drop=F];

		for(j in 1:len.sample)
		{
			col <- rep("Gray25",length(pos.mar));
			col[which(geno[j,]==2)] <- "yellow";
			col[which(geno[j,]==1)] <- "red";
			col[which(geno[j,]==0)] <- "blue";
			m <- pos.mar + (i-1)*3;
			rect(m,j,m+1,j+1,col = col,border=NA);
		}
		text( (pos.mar[1]+pos.mar[length(pos.mar)])/2+(i-1)*3, 0, labels=marker_table$Group[pos.mar[1]] ,adj=c(1,1), cex=0.7,srt=90,font=2,lwd=0.001);
	}

	m <- 0.55
	rect(xlim*(m-0.01), len.sample*1.05, xlim, len.sample*1.15);


	rect(xlim*m,len.sample*1.08, xlim*(m+0.04), len.sample*1.12,col="yellow");
	text(xlim*(m+0.07),len.sample*1.1,"AA", cex=0.7);
	rect(xlim*(m+0.1),len.sample*1.08, xlim*(m+0.14), len.sample*1.12,,col="red");
	text(xlim*(m+0.17),len.sample*1.1,"Aa", cex=0.7);
	rect(xlim*(m+0.2),len.sample*1.08, xlim*(m+0.24), len.sample*1.12,,col="blue");
	text(xlim*(m+0.27),len.sample*1.1,"aa", cex=0.7);
	rect(xlim*(m+0.3),len.sample*1.08, xlim*(m+0.34), len.sample*1.12,col="black");
	text(xlim*(m+0.395),len.sample*1.1,"missing", cex=0.7);

}



