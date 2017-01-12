###############################################################
# 
# System management utility
#
# History:
# 03/17/2010 Version 0.1
#
###############################################################

task_start<-function(...)
{
	p <- list(start_timer = proc.time());
		
	msgs <- list(...);
	if (length(msgs)==0)
		cat("The task is started...\r\n")
	else	
	{
		cat(..., sep="");
	}
		
	flush.console();
	return(p);
}

task_elapsed<-function(p, finished = NA, ...)
{
	nt <- proc.time() - p$start_timer;
		
	if(!is.na(finished))
		sys_t <-  round(c(nt[3], nt[3]*( 1 - finished)/finished ));
		
	sTime1<-"";
	sTime2<-"";
	if (sys_t[1]>60)
		sTime1 <- sprintf( "%02g:%02g:%02g", sys_t[1]%/%3600, (sys_t[1]- sys_t[1]%/%3600*3600)%/%60, sys_t[1]%%60 )
	else
		sTime1 <- sprintf( "%02g seconds", sys_t[1] );

	if (sys_t[2]>60)
		sTime2 <- sprintf( "%02g:%02g:%02g", sys_t[2]%/%3600, (sys_t[2]- sys_t[2]%/%3600*3600)%/%60, sys_t[2]%%60 )
	else
		sTime2 <- sprintf( "%02g seconds", sys_t[2] );
	sTime<-paste( sTime1, " has elapsed, left time: ", sTime2, ". ");
		
	msgs <- list(...);
	if ( length(msgs) != 0 )
	{
		for (i in 1:length(msgs))
			if (msgs[[i]]=="$SYS_PROMPT$")
				msgs[[i]] <- sTime;
	}
	else
		msgs[[1]]<-paste( sTime, "\r\n");

	cat( unlist(msgs), sep="" );
	flush.console();
}

task_stop<-function(...)
{
	msgs <- list(...);
	if (length(msgs)==0)
		cat("The task is stopped...\r\n")
	else	
		cat( ... , sep="" );
	flush.console();
}