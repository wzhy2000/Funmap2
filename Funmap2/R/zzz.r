FG_ENV <- new.env()

msg <- function(...)
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

    #packageStartupMessage("## Funmap2 Package v.2.5")
    #packageStartupMessage("## Build date: ", date(), "")
    #packageStartupMessage("## Copyright (C) 2011-", yr, ", http://ccb.bjfu.edu.cn", sep="")
    #packageStartupMessage("## Written by Zhong Wang(zhongwang@bjfu.edu.cn)\n")
}

.onAttach <- function(...)
{
	#msg();

	FG_ENV$i_hash <- new.env( hash=TRUE, parent=emptyenv(), size=100L);

	.RW("n.seed", 1);
	.RW("debug", FALSE);
	.RW("try.silent", TRUE);
	.RW("loop.est", 10);
	.RW("loop.mle", 4);
	.RW("mu.range.loop", 10);
}

.onLoad <- function(libname, pkgname)
{
	#Do Nothing
}

