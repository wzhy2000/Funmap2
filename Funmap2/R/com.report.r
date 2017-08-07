###############################################################
#
# New Systems Mapping Application(SysMap1)
#
# Basic Report Functions(Report Library)
#    1) Report.new
#    2) Report.par
#    3) Report.title
#    4) Report.AddHeadline
#    5) Report.AddParaLine
#    6) Report.AddTable
#    7) Report.AddFigure
#    8) Report.AddPageBreak
#    9) Report.AddParaBreak
#   10) Report.Output
#
# History:
# 12/15/2011 Version 1.1
#
###############################################################

Report.new<-function(report.file, options)
{
	FG_ENV$g_rpt <- list(
		doc.title="title",
		page.header="",
		page.footer="http://statgen.psu.edu/",
		page.offset = -50,
		dat.file = "",
		page.width = 0,
		page.height = 0,
		para.width = 0,
		line.height = 33,
		matr=list(),
		matr.num=0,
		layout=list(),
		pages=list(),
		tmpfile="temp.reprt.9988.pdf",
		options=options,
		report.file= report.file,
		heading.no = c(0) );
}

Report.par<-function( par_name, par_val )
{
	if (par_name=="dat.file")
		FG_ENV$g_rpt$dat.file <- par_val;

}

Report.title<-function( doc.title, page.header, page.footer )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	FG_ENV$g_rpt$doc.title   <- doc.title;
	FG_ENV$g_rpt$page.header <- page.header;
	FG_ENV$g_rpt$page.footer <- page.footer;
}

#level:c(1,2,3,4)
Report.AddHeadline<-function( str, level=1, options=list() )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");


	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="H", str=str, level=level, heading.no=get_headingno(level),options=options);
}


get_headingno<-function(level)
{
	if ( level == length(FG_ENV$g_rpt$heading.no)+1 )
	{
		FG_ENV$g_rpt$heading.no <- c( FG_ENV$g_rpt$heading.no, 1 )
		return(paste(FG_ENV$g_rpt$heading.no, sep="",collapse="."));
	}

	if ( level <= length(FG_ENV$g_rpt$heading.no) )
	{
		FG_ENV$g_rpt$heading.no[level] <- FG_ENV$g_rpt$heading.no[level]+1;
		if (level<length(FG_ENV$g_rpt$heading.no))
			for(i in length(FG_ENV$g_rpt$heading.no):(level+1))
				FG_ENV$g_rpt$heading.no <- FG_ENV$g_rpt$heading.no[-i];

		return(paste(FG_ENV$g_rpt$heading.no, sep="",collapse=".") );
	}

	return (" ");
}

Report.AddParaLine<-function( str, tab.pos=c(), tab.align=c(), options=list(), offset=0 )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="P", str=str, tab.pos=tab.pos, tab.align=tab.align, offset=offset);
}

Report.AddTable<-function( df, title=NULL, frame=c("full","wang", "san"), col.names=NULL, col.width=NULL, col.align=NULL, offset.x =NULL, max.show=0, options=list() )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="TB", data=df,
				title = title,
				frame = frame,
				col.names=col.names,
				col.width=col.width,
				col.align=col.align,
				options=options,
				offset.x = offset.x,
				max.show = max.show);
}

Report.AddFigure<-function( fCallDraw, title, size, left.margin=0, options=list(border=0) )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="F", title=title, size=size, fCallDraw=fCallDraw, fCallDraw2=NULL, left.margin=left.margin);
}

Report.AddPageBreak<-function()
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="NEWPAGE");
}


Report.AddParaBreak<-function(height=48)
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	if (height<48) height = 48;

	FG_ENV$g_rpt$matr.num <- FG_ENV$g_rpt$matr.num+1;
	FG_ENV$g_rpt$matr[[ FG_ENV$g_rpt$matr.num ]] <- list(id="NEWLINE", height=height);
}


Report.Output<-function( file )
{
	if (is.null(FG_ENV$g_rpt))
		stop("No report object in current context.");

	ret <- FALSE;
	pdf( FG_ENV$g_rpt$tmpfile, paper="letter", width=7, height=9.5);

	par(lwd=1) ;

	#try(ret <- Report.do_layout());
	ret <- Report.do_layout();
	dev.off();
	file.remove( FG_ENV$g_rpt$tmpfile)

	if (!ret)
		stop("Failed to layout this report.");

	pdf(file, paper="letter", width=7, height=9.5);

	#try( Report.do_pagecover() );
	#for (i in 1:length(FG_ENV$g_rpt$pages))
	#{
	#	try( Report.do_pagebody(i) );
	#
	#}

	Report.do_pagecover();
	for (i in 1:length(FG_ENV$g_rpt$pages))
	{
		Report.do_pagebody(i);

	}

	dev.off();
}


Report.get_cex<-function(level)
{
	cex <-c( 2,1.5,1.2,1 );
	return(cex[level]);
}

Report.do_pagebody<-function( page_no )
{
	if (.RR("debug", 0))
		cat("PageNo:", page_no, "Width,Height",
			length(FG_ENV$g_rpt$pages[[page_no]]$LayoutWidth),
			length(FG_ENV$g_rpt$pages[[page_no]]$LayoutHeight), "\n" );

	cat("The page", page_no, " is being plotted.............\n");

	nLayRange <- FG_ENV$g_rpt$pages[[page_no]]$LayoutRange;
	nf <- layout( FG_ENV$g_rpt$pages[[page_no]]$LayoutMatrix,
				 widths =FG_ENV$g_rpt$pages[[page_no]]$LayoutWidth,
				 heights=FG_ENV$g_rpt$pages[[page_no]]$LayoutHeight, TRUE);
	Report.do_pageheader(page_no);

	for(i in nLayRange[1]:nLayRange[2])
	{
		if ( FG_ENV$g_rpt$layout[[i]]$id=="H" )
			Report.do_heading( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="P" )
			Report.do_paragraph( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="TH" )
			Report.do_th( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="TR" )
			Report.do_tr( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="F")
			Report.do_figures( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="BR" )
			Report.do_br( FG_ENV$g_rpt$layout[[i]] );

		if ( FG_ENV$g_rpt$layout[[i]]$id=="BRECT" )
			Report.do_brect( FG_ENV$g_rpt$layout[[i]] );
	}

	Report.do_pagefooter();
}

Report.do_pagecover<-function()
{
	#cover
	par(mar=c(0,0,0,0))
	plot(1:20, 1:20, type="n", xlim=c(0, 7*254), ylim=c(0, 9.5*254), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n");

	text(7*254/2, 1600, FG_ENV$g_rpt$doc.title, cex=2, adj=0.5);
	text(700, 1200, paste("Data File: ",  FG_ENV$g_rpt$dat.file,  sep=""), cex=0.8, adj=0);
	text(700, 1130, paste("Report Date: ", Sys.time(), sep=""), cex=0.8, adj=0);
}

Report.do_pageheader<-function( pageno )
{
	par(mai=c(0,0,0,0));
	plot(1:20, 1:20, type="n", xlim=c(0, 7*254), ylim=c(0, 0.25*254), xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n");
	segments(0, 2, 7*254, 2, lwd=1, col="black");
	text(20, 20, FG_ENV$g_rpt$page.header, cex=0.8, adj=0);
	text(7*254-50, 20, paste("Page: ", pageno, sep=""), cex=0.8, adj=1);
}

Report.do_pagefooter<-function()
{
	plot(1:20, 1:20, type="n", xlim=c(0, 7*254), ylim=c(0, 0.25*254), xlab="", ylab="", xaxt="n", yaxt="n", bty="n",xaxs="i", yaxs="i");
	segments(0, 60, 7*254, 60, lwd=1, col="black");
	text(7*254-50, 30, FG_ENV$g_rpt$page.footer, cex=0.8, adj=1);
}

#LEVEL:1 cex=1.3
#LEVEL:2 cex=1
#LEVEL:3 cex=0.8
#LEVEL:4 cex=0.5
Report.do_heading<-function( lay )
{
	par(mar=c(0,0,0,0))
	plot(NA, NA, type="l", xlim=c(0, FG_ENV$g_rpt$page.width), ylim=c(0, lay$size[2]), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");

	text(0 + FG_ENV$g_rpt$page.offset, lay$size[2]/4, lay$str, cex=lay$cex, adj = c(0, 0) );
}

Report.do_paragraph<-function( lay )
{
	par(mar=c(0,0,0,0))

	plot(1:20, 1:20, type="n", xlim=c(0, FG_ENV$g_rpt$page.width), ylim=c(0, ifelse(lay$size[2]>0, lay$size[2],100)), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");

	for (i in 1:length(lay$sublay))
	{
		sub <- lay$sublay[[i]];
		adj <- c(0,0);
		pos.x<-	sub$tab.pos[1];

		if (sub$tab.align=="C")
			adj <-c(0.5,0);
		if (sub$tab.align=="R")
		{
			adj <-c(1,0);
			sub.next <- lay$sublay[[i+1]];
			if (is.null(sub.next))
				pos.x <- FG_ENV$g_rpt$page.width
			else
				pos.x <- sub.next$tab.pos[1];
		}

		text(pos.x+ FG_ENV$g_rpt$page.offset+lay$offset, lay$size[2]/4, sub$str, cex=1, adj = adj, vfont=c("serif", "plain")  );
	}
}

Report.do_th<-function( lay )
{
	par(mar=c(0,0,0,0))
	plot(NA, NA, type="n", xlim=c(0, FG_ENV$g_rpt$page.width), ylim=c(0, lay$size[2]), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");

	tab.pos <- cumsum(c(0, lay$col.width));

	for (i in 1:length(lay$col.names))
	{
		adj <-c(0,0);
		pos.x <- tab.pos[i];

		if (lay$col.align[i]=="C")
		{
			adj <-c(0.5,0);
			pos.x <- (tab.pos[i]+tab.pos[i+1])/2;
		}
		if (lay$col.align[i]=="R")
		{
			adj <-c(1,0);
			pos.x <- tab.pos[i+1];
		}
		if (lay$col.align[i]=="L")
		{
			pos.x <- pos.x + 8;
		}

		layoff <- 0;
		if (!is.null(lay$offset.x)) layoff = lay$offset.x;
		if (lay$frame=="full")
		{
			rect(tab.pos[i]  +FG_ENV$g_rpt$page.offset+layoff, -1,
				 tab.pos[i+1]+FG_ENV$g_rpt$page.offset+layoff, lay$size[2]+1, col=NA, border="black");
		}

		if (lay$frame=="wang")
		{
			x0 <- tab.pos[i]  +FG_ENV$g_rpt$page.offset+layoff;
			x1 <- tab.pos[i+1]+FG_ENV$g_rpt$page.offset+layoff;
			y0 <- -1;
			y1 <- lay$size[2]+1;

			segments(x0, y0, x1, y0 )
			segments(x0, y1, x1, y1 )
			if (i==1)
				segments(x1, y0, x1, y1 );

		}

		text(pos.x+ FG_ENV$g_rpt$page.offset + layoff, lay$size[2]/4, labels=lay$col.names[i], cex=1, adj = adj, vfont=c("serif", "plain")  );
	}
}

Report.do_tr<-function( lay )
{
	par(mar=c(0,0,0,0))
	plot(NA, NA, type="n", xlim=c(0, FG_ENV$g_rpt$page.width), ylim=c(0, lay$size[2]), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");

	tab.pos <- cumsum(c(0, lay$col.width));
	for (i in 1:length(lay$data))
	{
		adj <- c(0,0);
		pos.x <- tab.pos[i];

		if (lay$col.align[i]=="C")
		{
			adj <- c(0.5,0);
			pos.x <- (tab.pos[i]+tab.pos[i+1])/2;
		}
		if (lay$col.align[i]=="R")
		{
			adj <- c(1,0);
			pos.x <- tab.pos[i+1];
		}

		if (lay$col.align[i]=="L")
		{
			pos.x <- pos.x + 8;
		}


		layoff <- 0;
		if (!is.null(lay$offset.x)) layoff = lay$offset.x;

		if (lay$frame=="full")
		{
			rect(tab.pos[i]  +FG_ENV$g_rpt$page.offset+layoff, 0-1,
				 tab.pos[i+1]+FG_ENV$g_rpt$page.offset+layoff, lay$size[2]+1, col=NA, border="black");
		}

		if (lay$frame=="wang")
		{
			x0 <- tab.pos[i]  +FG_ENV$g_rpt$page.offset+layoff;
			x1 <- tab.pos[i+1]+FG_ENV$g_rpt$page.offset+layoff;
			y0 <- -1;
			y1 <- lay$size[2]+1;

			if(i==1)
				segments(x1, y0, x1, y1 );
			if (lay$bLastRow)
				segments(x0, y0, x1, y0 );
		}


		text(pos.x+FG_ENV$g_rpt$page.offset+layoff, lay$size[2]/4, labels=as.character(unlist(lay$data[i])), cex=0.9, adj = adj, vfont=c("serif", "plain")  );
	}
}

Report.do_brect<-function( lay )
{
	par(mar=c(0,0,0,0))
	plot(1:20, 1:20, type="n", xlim=c(0, lay$size[1]), ylim=c(0, lay$size[2]), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");
}

Report.do_br<-function( lay )
{
	par(mar=c(0,0,0,0))
	plot(1:20, 1:20, type="n", xlim=c(0, FG_ENV$g_rpt$page.width), ylim=c(0, lay$size[2]), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n", bty="n");
}

Report.do_figures<-function( lay )
{
	par(mar=c(0,0,0,0))
	eval(lay$fCallDraw)
}

Report.get_textsize<-function( str, text.cex, text.widthlim, vfont=c("serif", "plain") )
{
	w <- strwidth(  str, cex=text.cex, adj=0, vfont=vfont );
	h <- strheight( str, cex=text.cex, adj=0, vfont=vfont )*1.3;
	if (w <= text.widthlim )
		return(c(w,h));

	str.words <- strsplit(str, " ");
	str.newwords <- c();
	str.testline <- c(str.words[1]);

	if(length(str.words ) > 1)
	{
		for (i in 2:length(str.words ) )
		{
			str.testwidth <- strwidth( paste( c( str.testline, str.words[i] ) , sep=" " ),
							cex=text.cex, adj=0, vfont=c("serif", "plain") );
			if ( str.testwidth > text.widthlim )
			{
				str.newwords <- c(str.newwords, str.testline, "\n");
				str.testline <- c(str.words[i]);
			}
			else
				str.testline <- paste( str.testline, str.words[i], " ");

		}
		w <- strwidth(  str.newwords, cex=text.cex, adj=0, vfont=vfont );
		h <- strheight( str.newwords, cex=text.cex, adj=0, vfont=vfont )*1.3;
	}

	return(c( min(w, text.widthlim), sum(h) ) );
}

Report.do_layout<-function()
{
	FG_ENV$g_rpt$page.width  <- 7*254;
	FG_ENV$g_rpt$page.height <- 9.5*254;

	page.height = 9.5*254-0.5*254;
	page.width  = 7*254;
	page.pos    = 0;
	page.figx   = 0;
	page.figy   = 0;
	page.no     = 1;
	last.lay    = 1;


	plot(0, 0, type="n", xlim=c(0, page.width), ylim=c(0, page.height), xlab="", ylab="", lwd=2, xaxt="n", yaxt="n");

	nlay <-1;
	for(i in 1:length(FG_ENV$g_rpt$matr))
	{
		if (.RR("debug", 0))
			cat(FG_ENV$g_rpt$matr[[i]]$id, "\n");

		if ( FG_ENV$g_rpt$matr[[i]]$id=="H" )
		{
			elm_h = FG_ENV$g_rpt$matr[[i]];
			size=Report.get_textsize( elm_h$str, Report.get_cex(elm_h$level), page.width );

			if ( page.pos + size[2] > page.height)
			{
				#append a blank area.
				page.break<- page.height - page.pos;
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, page.break), pos=page.pos, page=page.no );
				FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
				page.pos <- 0;
				nlay     <- nlay+1;
				last.lay <- nlay;
				page.no  <- page.no + 1;

				#<--Add a space between the headline and the page body
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
				page.pos <- page.pos +  size[2];
				page.figy <- page.pos;
				nlay <- nlay+1;
				#-->End of space adding
			}


			FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="H", size=size, pos=page.pos, page=page.no, str=paste(elm_h$heading.no, elm_h$str, sep=" "), cex=Report.get_cex(elm_h$level) );
			page.pos <- page.pos +  size[2];
			page.figy <- page.pos;
			nlay <- nlay+1;
		}

		if ( FG_ENV$g_rpt$matr[[i]]$id=="P" )
		{
			elm_p = FG_ENV$g_rpt$matr[[i]];
			offset = FG_ENV$g_rpt$matr[[i]]$offset;

			p_str <- c();
			p_pos <- c();
			p_alg <- c();

			for(i in 1:length(elm_p$str))
			{
				p_str <- c(p_str, elm_p$str[i]);
				if (i <= length(elm_p$tab.pos))
					p_pos <- c( p_pos, elm_p$tab.pos[i] )
				else
					p_pos <- c( p_pos, -1 );
				if (i <= length(elm_p$tab.align))
					p_alg <- c( p_alg, elm_p$tab.align[i] )
				else
					p_alg <- c( p_alg, "L");
			}

			if ( page.pos + FG_ENV$g_rpt$line.height > page.height)
			{
				page.break<- page.height - page.pos;
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
				FG_ENV$g_rpt$pages[[page.no]] <-  list(LayoutRange = c( last.lay, nlay ) );
				last.lay <- nlay;
				nlay     <- nlay+1;
				page.pos <- 0;
				page.no  <- page.no + 1;

				#<--Add a space between the headline and the page body
				FG_ENV$g_rpt$layout[[ nlay ]]  <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
				page.pos <- page.pos +  size[2];
				page.figy <- page.pos;
				nlay <- nlay+1;
				#-->End of space adding

			}

			for(i in 1:length(p_str))
			{
				if (p_pos[i] == -1)
				{
					size=Report.get_textsize( p_str[i], 1, page.width - 5 );
					p_pos[i]<- size[1];
				}
			}

			FG_ENV$g_rpt$layout[[nlay]] <- list(id="P", size=size, pos=page.pos, page=page.no, offset=offset, sublay=list())
			sublay <- list();
			for (i in 1:length(p_str))
			{
				if (i==1)
					sublay[[i]] <- list(id="T", str=p_str[i], tab.pos=c(0, p_pos[i]), tab.align = p_alg[i])
				else
					sublay[[i]] <- list(id="T", str=p_str[i], tab.pos=c(p_pos[i-1], p_pos[i]), tab.align = p_alg[i]);
			}

			FG_ENV$g_rpt$layout[[nlay]]$sublay <- sublay;
			page.pos <- page.pos +  size[2];
			page.figy <- page.pos;
			nlay <- nlay+1;
		}

		if ( FG_ENV$g_rpt$matr[[i]]$id=="TB" )
		{
			elm_tb = FG_ENV$g_rpt$matr[[i]];

			size = Report.get_textsize( elm_tb$title, 1, page.width- 5 );
			if ( page.pos + size[2] > page.height)
			{
				page.break<- page.height - page.pos;
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
				FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
				page.pos <- 0;
				nlay <- nlay + 1;
				last.lay <- nlay;
				page.no  <- page.no + 1;

				#<--Add a space between the headline and the page body
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
				page.pos <- page.pos +  size[2];
				page.figy <- page.pos;
				nlay <- nlay+1;
				#-->End of space adding

			}

			FG_ENV$g_rpt$layout[[nlay]] <- list(id="P", size=size, pos=page.pos, page=page.no, sublay=list(), offset=0 );
			sublay <- list(str=elm_tb$title, tab.pos=c(0.25*254), tab.align="L" );
			FG_ENV$g_rpt$layout[[nlay]]$sublay[[1]] <- sublay;
			page.pos <- page.pos +  size[2];
			nlay<-nlay+1;

			if ( page.pos + 3.3 > page.height)
			{
				page.break<- page.height - page.pos;
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
				FG_ENV$g_rpt$pages[[page.no]] <-  list(LayoutRange = c( last.lay, nlay) );
				page.pos <- 0;
				page.no  <- page.no + 1;
				nlay <- nlay + 1;
				last.lay <- nlay;

				#<--Add a space between the headline and the page body
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
				page.pos <- page.pos +  size[2];
				page.figy <- page.pos;
				nlay <- nlay+1;
				#-->End of space adding

			}

			FG_ENV$g_rpt$layout[[nlay]] <- list(id="TH", size=size, pos=page.pos, page=page.no,
			            col.names=elm_tb$col.names,
			            col.width = elm_tb$col.width,
			            col.align = elm_tb$col.align,
			            offset.x=elm_tb$offset.x,
			            frame=elm_tb$frame );
			page.pos <- page.pos +  size[2];
			nlay<-nlay+1;
			page.figx <- 0;

			dat.nrow <- dim(FG_ENV$g_rpt$matr[[i]]$data)[1];
			if (FG_ENV$g_rpt$matr[[i]]$max.show!=0 && dat.nrow > FG_ENV$g_rpt$matr[[i]]$max.show)
				dat.nrow <- FG_ENV$g_rpt$matr[[i]]$max.show;

			for (d in 1:dat.nrow)
			{
				if ( page.pos + 48 > page.height)
				{
					page.break <- page.height - page.pos;
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
					FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
					page.pos <- 0;
					page.no  <- page.no + 1;
					nlay <- nlay + 1;
					last.lay <- nlay;

					#<--Add a space between the headline and the page body
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
					page.pos <- page.pos +  size[2];
					page.figy <- page.pos;
					nlay <- nlay+1;
					#-->End of space adding

				}

				dat.vec <- elm_tb$data[d,];
				if (d==dat.nrow && dat.nrow<dim(FG_ENV$g_rpt$matr[[i]]$data)[1])
					for(dr in 1:length(dat.vec) )
						dat.vec[dr] <- "...";

				FG_ENV$g_rpt$layout[[nlay]] <- list(id="TR", size=c(page.width, 48), pos=page.pos, page=page.no,
						data      = dat.vec,
						col.width = elm_tb$col.width,
						col.align = elm_tb$col.align,
						offset.x  = elm_tb$offset.x,
			            frame     = elm_tb$frame,
			            bLastRow  = (d==dat.nrow) );
				page.pos <- page.pos +  48;
				nlay <- nlay+1;
			}


			page.figy <- page.pos;

		}
		if ( FG_ENV$g_rpt$matr[[i]]$id=="F")
		{
			elm_f = FG_ENV$g_rpt$matr[[i]];

			if (page.figx >0)
			{
				if (page.figx + elm_f$left.margin + elm_f$size[1] <= page.width )
				{
					if (elm_f$left.margin>0)
					{
						FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BRECT",
								page=page.no,
								page.figx = page.figx,
								page.figy = page.figy,
								size = c(elm_f$left.margin, elm_f$size[2]) );
						nlay <- nlay + 1;
					}

					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="F", page=page.no,
								page.figx = page.figx+elm_f$left.margin,
								page.figy = page.figy,
								size = elm_f$size,
								title=elm_f$title,
								fCallDraw=elm_f$fCallDraw);

					nlay <- nlay + 1;
					page.figx <- page.figx + elm_f$left.margin + elm_f$size[1];
					pagep.pos <- page.figy + elm_f$size[2];
				}
				else
				{
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BRECT",
								page=page.no,
								page.figx = page.figx,
								page.figy = page.figy,
								size = c(page.width - page.figx, page.pos - page.figy) );
					nlay <- nlay + 1;
					page.figx <- 0;
					page.figy <- page.pos;
				}
			}

			if (page.figx ==0)
			{
				if ( page.pos + elm_f$size[2] > page.height)
				{
					page.break<- page.height - page.pos;
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
					FG_ENV$g_rpt$pages[[page.no]] <-  list(LayoutRange = c( last.lay, nlay ) );

					page.no  <- page.no + 1;
					nlay <- nlay + 1;
					last.lay <- nlay;
					page.pos <- 0;
					page.figy <- 0;
					page.figx <- 0;

					#<--Add a space between the headline and the page body
					FG_ENV$g_rpt$layout[[ nlay ]]  <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
					page.pos <- page.pos +  size[2];
					page.figy <- page.pos;
					nlay <- nlay+1;
					#-->End of space adding

				}


				if (elm_f$left.margin>0)
				{
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BRECT",
								page=page.no,
								page.figx = page.figx,
								page.figy = page.figy,
								size = c(elm_f$left.margin, elm_f$size[2]) );

					nlay <- nlay + 1;
					page.figx <-  page.figx + elm_f$left.margin;
				}

				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="F" , page=page.no,
								page.figx = page.figx,
								page.figy = page.figy,
								size = elm_f$size,
								title=elm_f$title,
								fCallDraw=elm_f$fCallDraw,
								fCallDraw2=elm_f$fCallDraw2);

				nlay <- nlay + 1;
				page.figx <- page.figx + elm_f$size[1];
				page.figy <- page.pos;
				page.pos <- page.pos + elm_f$size[2];
			}
		}


		if ( FG_ENV$g_rpt$matr[[i]]$id=="NEWPAGE" )
		{
			page.break<- page.height - page.pos;
			FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
			FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
			page.no  <- page.no + 1;
			page.pos <- 0;
			nlay <- nlay + 1;
			last.lay <- nlay;
			page.figx <- 0;
			page.figy <- page.pos;

			#<--Add a space between the headline and the page body
			FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
			page.pos <- page.pos +  size[2];
			page.figy <- page.pos;
			nlay <- nlay+1;
			#-->End of space adding
		}

		if ( FG_ENV$g_rpt$matr[[i]]$id=="NEWLINE" )
		{
 			if (page.figx>0)
			{
				if (.RR("debug", 0))
					cat("page.figx:", page.figx, "\n");

				if (page.width - page.figx>0)
				{
					FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BRECT",
						page=page.no,
						page.figx = page.figx,
						page.figy = page.figy,
						size = c(page.width - page.figx, page.pos - page.figy) );
					nlay <- nlay+1;
				}

				page.figx <- 0;
				page.figy <- page.pos;
			}

			if ( page.pos + FG_ENV$g_rpt$matr[[i]]$height > page.height)
			{
				page.break<- page.height - page.pos;
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
				FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
				page.no  <- page.no + 1;
				page.pos <- 0;
				nlay <- nlay + 1;
				last.lay <- nlay;


				#<--Add a space between the headline and the page body
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width, 30), pos=page.pos, page=page.no );
				page.pos <- page.pos +  size[2];
				page.figy <- page.pos;
				nlay <- nlay+1;
				#-->End of space adding

			}
			else
			{
				FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,FG_ENV$g_rpt$matr[[i]]$height), pos=page.pos, page=page.no );
				page.pos <- page.pos +  FG_ENV$g_rpt$matr[[i]]$height;
				nlay <- nlay+1;
				page.figx <- 0;
				page.figy <- page.pos;
			}
		}

	}

	if (last.lay < nlay )
	{
		page.break <- page.height - page.pos;
		FG_ENV$g_rpt$layout[[ nlay ]] <- list(id="BR", size=c(page.width,page.break), pos=page.pos, page=page.no );
		FG_ENV$g_rpt$pages[[page.no]] <- list(LayoutRange = c( last.lay, nlay ) );
	}

	for (i in 1:length(FG_ENV$g_rpt$pages))
	{
		layMat <- c();
		layHeight <- c();
		layWidth  <- c();
		lay.range <- FG_ENV$g_rpt$pages[[i]]$LayoutRange;

		cat("The layout for page", i, " is being generated.............\n");

		laymat <- list( mat=c(1), width=c(FG_ENV$g_rpt$page.width), height=c(0.25*254), col=1, row=1, elm_no=1, rows=list(c(1)));

		for (j in lay.range[1]:lay.range[2])
		{
			if (FG_ENV$g_rpt$layout[[j]]$id != "F" && FG_ENV$g_rpt$layout[[j]]$id != "BRECT")
				laymat <- laymat.add( laymat, FG_ENV$g_rpt$layout[[j]]$id,  FG_ENV$g_rpt$layout[[j]]$size[1], FG_ENV$g_rpt$layout[[j]]$size[2] )
			else
				laymat <- laymat.add( laymat, FG_ENV$g_rpt$layout[[j]]$id,
									  FG_ENV$g_rpt$layout[[j]]$size[1],
									  FG_ENV$g_rpt$layout[[j]]$size[2],
									  page.figx = FG_ENV$g_rpt$layout[[j]]$page.figx,
									  page.figy = FG_ENV$g_rpt$layout[[j]]$page.figy );
		}

		laymat <- laymat.add( laymat, "BR", FG_ENV$g_rpt$page.width, 0.25*254, done=TRUE );

		FG_ENV$g_rpt$pages[[i]]$LayoutMatrix <- laymat$mat;
		FG_ENV$g_rpt$pages[[i]]$LayoutWidth  <- laymat$mat.width;
		FG_ENV$g_rpt$pages[[i]]$LayoutHeight <- laymat$mat.height;
	}

	return(TRUE)
}

laymat.add<-function(laymat, id, width, height, page.figx=-1, page.figy =-1,done=FALSE )
{
	if (.RR("debug", 0))
		cat(id, width, height, page.figx, page.figy, "\n");

	laymat$elm_no <- laymat$elm_no+1;
	laymat$width  <- c(laymat$width, width);
	laymat$height <- c(laymat$height, height);
	addno <- laymat$elm_no;

	if (page.figy == -1)
	{
		laymat$row <- laymat$row+1;
		laymat$rows[[ laymat$row ]] <- c( addno );
	}
	else
	{
		if (page.figx==0)
		{
			laymat$row <- laymat$row+1;
			laymat$rows[[ laymat$row ]] <- c( addno );
		}
		else
			laymat$rows[[ laymat$row ]] <- c( laymat$rows[[ laymat$row ]], addno );

	}

	if (done)
	{
		nmax <- 1;
		nmax.i <- 0;
		nmax.m <- 1;
		nmax.h <- c();
		nmax.w <- c(FG_ENV$g_rpt$page.width);

		for (i in 1:length(laymat$rows))
		{
			if ( length(laymat$rows[[i]])> nmax )
			{
				nmax <- length(laymat$rows[[i]]);
				nmax.m <- i;
				nmax.h <- c( nmax.h, laymat$height[nmax.i+1] );
				nmax.w <- laymat$width[ (nmax.i+1): (nmax.i+length(laymat$rows[[i]])) ];
			}
			else
			{
				nmax.h <- c( nmax.h, laymat$height[nmax.i+1] );
			}

			nmax.i <- nmax.i + length(laymat$rows[[i]]);
		}

		laymat$mat <- array(0, dim=c(laymat$row, nmax));
		for (i in 1:length(laymat$rows))
		{
			if ( length(laymat$rows[[i]])>1)
				laymat$mat[i,] <- laymat$rows[[i]]
			else
				laymat$mat[i,] <- rep(laymat$rows[[i]], nmax);
		}

		laymat$mat.width <- nmax.w
		laymat$mat.height <- nmax.h
	}

	return( laymat );
}
