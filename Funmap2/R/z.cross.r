
#--------------------------------------------------------------
# private: fin.generate_bc_marker;
#
# genarate N Backcross Markers from marker disttance (cM): dist.
#
# input
#     dist : the vector for the marker distances
#   samp_N : the sample size
#--------------------------------------------------------------
fin.generate_bc_marker<-function(samp_N, dis)
{
	if (dis[1] != 0)
		Cm=c(0, dis)/100
	else
		Cm=dis/100;

	n <- length(Cm);
	rs <- 1/2*( exp(2*Cm) - exp(-2*Cm) ) / (exp(2*Cm)+exp(-2*Cm));
	mk<-array( 0, dim=c(samp_N,n) );

	for (j in 1:samp_N)
		mk[j,1] <- ( runif(1)>0.5 );

	for (i in 2:n)
		for (j in 1:samp_N)
		{
			if (mk[j,i-1]==1)
				mk[j,i] <- ( runif(1)>rs[i] )
			else
				mk[j,i] <- ( runif(1)<rs[i] );
    	}
	return(mk);
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Backcross
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

BC.sum_par <- function(par, fmt)
{
	return("");
}

BC.get_qtl_prob <- function( lmarker, rmarker, dist, d.qtl)
{
	prob     <- BC.get_cond_prob( dist, d.qtl );
	nInx     <- 4 - (lmarker*2 + rmarker);
	allprob  <- prob[ nInx ];

	return ( cbind(allprob, 1-allprob) );
}

BC.get_cond_prob <- function(dist, d.qtl, unit="cm")
{
	if (unit=="cm")
	{
		d.qtl <- d.qtl/100;
		dist  <- dist/100;
	}

	theta  <- d.qtl/dist;
	prob   <- c(1,1-theta,theta,0);

	return(prob);
}

BC.get_simu_marker<-function(simu_N, simu_mrkdist, cross_par)
{
	marker <- fin.generate_bc_marker( simu_N, simu_mrkdist );
	return(marker);
}


BC.get_simu_qtl<-function( simu_N, lmarker, rmarker, qtl.pos, lmarker.pos, rmarker.pos, option)
{
	if (is.na(rmarker.pos))
		return(rmarker);

	qtl  <- rep(0, simu_N );
	prob <- BC.get_cond_prob(rmarker.pos - lmarker.pos, qtl.pos - lmarker.pos );

	allprob <- prob[ 4 - ( lmarker*2 + rmarker ) ];
	qtl<-array(0, dim=c( simu_N ))
	for (i in 1:simu_N)
	{
		#return 0 or 1
		qtl[i] <- ( runif(1)>allprob[i] );
	}

	return(qtl);
}

fin.get_gen_par<-function(par, fm2_cross_gen_num, fm2_cross_gen_QQ, fm2_cross_gen_Qq, fm2_cross_gen_qq)
{
	mu_QQ <- c();
	mu_Qq <- c();
	mu_qq <- c();

	nstart <- 0;
	wlen <- length(par)/fm2_cross_gen_num;
	if(fm2_cross_gen_QQ)
	{
		mu_QQ <- par[(nstart+1):(wlen+nstart)];
		nstart <- nstart+wlen;
	}

	if(fm2_cross_gen_Qq)
	{
		mu_Qq <- par[(nstart+1):(wlen+nstart)];
		nstart <- nstart+wlen;
	}

	if(fm2_cross_gen_qq)
	{
		mu_qq <- par[(nstart+1):(wlen+nstart)];
		nstart <- nstart+wlen;
	}

	return(list(QQ=mu_QQ, Qq=mu_Qq, qq=mu_qq));
}

cross_BC<-list(
	type 	= "BC",
	desc 	= "Backcross",
	hp_desc = "mu_Qq == mu_qq",
	gen_num = 2,
	par_num = 2,
	gen_QQ  = F,
	gen_Qq  = T,
	gen_qq  = T,
	single_marker	= F,
	get_gen_par  = function(par){ fin.get_gen_par(par, 2, F, T, T)},
	get_qtl_prob = BC.get_qtl_prob,
	get_sum_par  = BC.sum_par,
	get_simu_marker = BC.get_simu_marker,
	get_simu_qtl = BC.get_simu_qtl );


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  F2
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

F2.get_qtl_prob <- function( lmarker, rmarker, dist, d.qtl)
{
	prob.mat <- F2.get_cond_prob( dist, d.qtl );
	nInx     <- lmarker*3 + rmarker +1;
	allprob  <- prob.mat[ nInx, ];

	return (allprob);
}
#--------------------------------------------------------------
# private: F2.get_cond_prob
#
# Create a matrix for conditional probility of F2 cross.
# See also Statistical Genetics of Quantitative Traits.
#          Page 255, Table11.7
#
# Input : d: marker distance,
#        d1: QTL pos. between two markers.
# output: the matrix[9,3] of conditional probility.
#--------------------------------------------------------------
F2.get_cond_prob <- function(dist, d.qtl, unit="cm")
{
	if (unit=="cm")
	{
		dist  <- dist/100;
		d.qtl <- d.qtl/100;
	}

	r  = (1-exp(-2*dist))/2;
	r1 = (1-exp(-2*d.qtl))/2;

	nita <- r^2/( (1-r)^2 + r^2 );
	seta <- r1/r;

	M<-array(0, dim=c(9,3));
	#mm nn (n00)
	M[1,] <- c( 0, 		0, 		1);
	#mm Nn (n01)
	M[2,] <- c( 0, 		seta, 	1-seta);
	#mm NN (n02)
	M[3,] <- c( seta^2, 2*seta*(1-seta), (1-seta)^2);

	#Mm nn (n10)
	M[4,] <- c( 0, 		1-seta, seta);
	#Mm Nn (n11)
	M[5,] <- c( nita*seta*(1-seta), 1-2*nita*seta*(1-seta), nita*seta*(1-seta));
	#Mm NN (n12)
	M[6,] <- c( seta, 	(1-seta), 0);

	#MM nn (n20)
	M[7,] <- c( (1-seta)^2, 2*seta*(1-seta), seta^2);
	#MM Nn (n21)
	M[8,] <- c( (1-seta), seta, 0);
	#MM NN (n22)
	M[9,] <- c( 1, 		0, 		0);

	return (M);
}

F2.get_simu_marker<-function(simu_N, simu_mrkdist, cross_par)
{
	marker1 <- fin.generate_bc_marker( simu_N, simu_mrkdist );
	marker2 <- fin.generate_bc_marker( simu_N, simu_mrkdist );
	return( marker1 + marker2 );
}


F2.get_simu_qtl<-function( simu_N, lmarker, rmarker, qtl.pos, lmarker.pos, rmarker.pos, option)
{
	if (is.na(rmarker.pos))
		return(rmarker);

	prob <- F2.get_cond_prob(rmarker.pos - lmarker.pos, qtl.pos - lmarker.pos );
	nInx <- lmarker*3 + rmarker +1;
	allprob<- prob[ nInx, ];

	qtl<-array(0, dim=c( simu_N ))
	for (i in 1:simu_N)
	{
		qtl_prob <- cumsum( allprob[i,] );
		#return 0/1/2
		qtl[i] <- min( which(runif(1) < qtl_prob )  ) - 1  ;
	}

	return(qtl);
}

F2.sum_par <- function(par, fmt)
{
	return("");
}

cross_F2<-list(
	type 	= "F2",
	desc 	= "F2",
	hp_desc = "mu_QQ == mu_Qq == mu_qq",
	gen_num = 3,
	gen_QQ  = T,
	gen_Qq  = T,
	gen_qq  = T,
	single_marker	= F,
	get_gen_par  = function(par){ fin.get_gen_par(par, 3, T, T, T)},
	get_simu_marker = F2.get_simu_marker,
	get_simu_qtl = F2.get_simu_qtl,
	get_qtl_prob = F2.get_qtl_prob,
	get_cond_mat = F2.get_cond_prob,
	get_sum_par  = F2.sum_par);

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  RIL
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RIL.get_qtl_prob <- function( lmarker, rmarker, dist, d.qtl)
{
	prob   <- RIL.get_cond_prob( dist, d.qtl );
	lmarker[ which(lmarker == 2) ] <- 1;
	rmarker[ which(rmarker == 2) ] <- 1;
	allprob<- prob[lmarker*2 + rmarker+1 ];
	return ( cbind(allprob, 1-allprob) );
}

#--------------------------------------------------------------
# private: get_cond_prob
#
# Create a matrix for conditional probility of RIL cross.
# See also Statistical Genetics of Quantitative Traits.
#          Page 235, Table10.7
#
# Input : d: marker distance,
#         d1: QTL pos. between two markers.
# output: the matrix[4,1] of conditional probility.
#--------------------------------------------------------------
RIL.get_cond_prob <- function(dist, d.qtl, unit="cm")
{
	if (unit=="cm")
	{
		d.qtl <- d.qtl/100;
		dist  <- dist/100;
	}

	r1 = (1-exp(-2*d.qtl))/2;
	r  = (1-exp(-2*dist))/2;
	#r2 = r - r1;
	r2 = (r-r1)/(1-2*r1);

	R1 <- 2*r1/( 1+2*r1);
	R2 <- 2*r2/( 1+2*r2);

	M<-array(0, dim=c(4));

	#mm nn (n00)
	x <- R1*R2*(3-2*R1-2*R2)/(2*(1-R1)*(1-R2));
	#genotype=1(QQ)
	M[1] <- x;

	#mm NN (n01)
	x <- (2*R1 - R1*R2*(3+2*R1-2*R2))/( 2*R2 + R1*(2-6*R2) );
	#genotype=1(QQ)
	M[2] <- x;

	#MM nn (n10)
	x <- (2*R1 - R1*R2*(3+2*R1-2*R2))/( 2*R2 + R1*(2-6*R2) );
	#genotype=1(QQ)
	M[3] <-  1-x;

	#MM NN (n11)
	x<-R1*R2*(3-2*R1-2*R2)/(2*(1-R1)*(1-R2));
	#genotype=1(QQ)
	M[4] <-  1-x;


	return (M);
}

RIL.get_simu_marker<-function(simu_N, simu_mrkdist, cross_par)
{
	marker <- fin.generate_bc_marker( simu_N, simu_mrkdist );
	marker[which(marker==1)]<-2;
	return(marker);
}

RIL.get_simu_qtl<-function( simu_N, lmarker, rmarker, qtl.pos, lmarker.pos, rmarker.pos, option)
{
	if (is.na(rmarker.pos))
		return(rmarker);

	qtl  <- rep(0, simu_N );
	prob <- RIL.get_cond_prob(rmarker.pos - lmarker.pos, qtl.pos - lmarker.pos );

	lmarker[ which(lmarker == 2) ] <- 1;
	rmarker[ which(rmarker == 2) ] <- 1;
	allprob <- prob[ lmarker*2 + rmarker+1];
	qtl<-array(0, dim=c( simu_N ))

	for (i in 1:simu_N)
	{
		# 0: homozygote(qq)
		# 2: homozygote+additive(QQ)
		qtl[i] <- ( runif(1)>allprob[i] ) *2  ;
	}

	return(qtl);
}

RIL.sum_par <- function(par, fmt)
{
	return("");
}

cross_RIL<-list(
	type 	= "RIL",
	desc 	= "RIL",
	hp_desc = "mu_QQ == mu_qq",
	gen_num = 2,
	gen_QQ  = T,
	gen_Qq  = F,
	gen_qq  = T,
	single_marker	= F,
	get_gen_par  = function(par){ fin.get_gen_par(par, 2, T, F, T)},
	get_simu_marker = RIL.get_simu_marker,
	get_qtl_prob = RIL.get_qtl_prob,
	get_cond_mat = RIL.get_cond_prob,
	get_sum_par  = RIL.sum_par,
	get_simu_qtl = RIL.get_simu_qtl );


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Natural Population
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NP1.get_qtl_prob <- function( lmarker, rmarker, dist, d.qtl, option)
{
	if (class(option)=="list")
	{
		p <- option$p;
		q <- option$q;
		D <- option$D;
	}
	else
	{
		p <- option[1];
		q <- option[2];
		D <- option[3];
	}

	prob   <- NP1.get_cond_prob( p, q, D );
	allprob<- prob[ lmarker + 1, ];

	return (allprob);
}

#--------------------------------------------------------------
# private: get_cond_prob
#
# Create a matrix for conditional probility of RIL cross.
# See also Statistical Genetics of Quantitative Traits.
#          Page 235, Table10.7
#
# Input : d: marker distance,
#         d1: QTL pos. between two markers.
# output: the matrix[4,1] of conditional probility.
#--------------------------------------------------------------
NP1.get_cond_prob <- function(p, q, D)
{
  	p11 <- p*q+D
  	p10 <- p*(1-q)-D
  	p01 <- (1-p)*q-D
  	p00 <- (1-p)*(1-q)+D

	M  = matrix(0,3,3);
	MM = 3
	Mm = 2;
	mm = 1;
	aa = 1
	Aa = 2
	AA = 3

	M[mm, aa] <- p00^2
	M[mm, Aa] <- 2*p01*p00
	M[mm, AA] <- p01^2

	M[Mm, aa] <- 2*p10*p00
	M[Mm, Aa] <- 2*(p11*p00+p10*p01)
	M[Mm, AA] <- 2*p11*p01

	M[MM, aa] <- p10^2
	M[MM, Aa] <- 2*p11*p10
	M[MM, AA] <- p11^2

	M<- M/rowSums(M);

	return (M);
}

NP1.sum_par <- function(par, fmt)
{
	return("");
}

NP1.get_simu_marker<-function(simu_N, simu_mrkdist, option)
{
	marker <- array(0, dim=c(simu_N, 2));
	prob <- NP1.get_cond_prob( option$p, option$q, option$D );

	#MM, Mm, mm
	prob.1 <- cumsum( c(option$p^2, 2*option$p*(1-option$p), (1-option$p)^2));
	for (i in 1:simu_N)
	{
		x<-runif(1)
		x.ord <- min(which( x < prob.1 ) );
		marker[i,1]<-c(2,1,0)[x.ord];
	}

	#aa, Aa, AA
	for (i in 1:simu_N)
	{
		x<-runif(1)
		x.p <- cumsum( prob[ marker[i,1] + 1,  ] );
		x.ord <- min(which( x < x.p) );
		marker[i,2]<-c(0, 1, 2)[x.ord];
	}

	return(marker);
}

NP1.get_simu_qtl<-function( simu_N, lmarker, rmarker, qtl.pos, lmaker.pos, rmaker.pos, option)
{
	if (is.na(rmaker.pos))
		return(rmarker);

	qtl  <- rep(0, simu_N );
	prob <- NP1.get_cond_prob( option$p, option$q, option$D );

	for (i in 1:simu_N)
	{
		x<-runif(1)
		x.p <- cumsum( prob[ lmarker[i] + 1,  ] );
		x.ord <- min(which( x < x.p) );
		qtl[i]<-c(0, 1, 2)[x.ord];
	}

	return(qtl);

}

NP1.get_init_rand<-function(dat)
{
	p=runif(1, min=0.51, max=1);
	q=runif(1, min=0.51, max=1);
	D=runif(1, min=0, max=0.10);
	#p<-0.6;
	#q<-0.7;
	#D<-0.05;

	return(list(p=p,q=q,D=D));
}

NP1.get_est_param<-function( dat, par.cross, par.covar, par.curve, obj.curve)
{
	if (!SM.covar$is_valid(par.covar))
		return(NaN);

	qtl.prob <- NP1.get_qtl_prob( dat$obj.gen$genos.matrix[,1], NULL, 0, 0, par.cross);
	time.std <- dat$obj.phe$sample_times;
	y <- as.matrix( dat$obj.phe$pheY );
	sig.inv <- FM.covar$get_inv_mat(par.covar, dat$obj.phe$sample_times, FM.curve$trait_num);
	sig.det <- FM.covar$get_mat_det(par.covar, dat$obj.phe$sample_times, FM.curve$trait_num);

	marker <- dat$obj.gen$genos.matrix[,1];
	nn <- nrow(y);

	p  <- par.cross$p;
	q  <- par.cross$q;
	D  <- par.cross$D;

	cat("NP1:", p, q, D, "\n");

	op <- 0;
	oq <- 0;
	oD <- 0;

    p11 <- p*q+D
    p10 <- p*(1-q)-D
    p01 <- (1-p)*q-D
    p00 <- (1-p)*(1-q)+D
	loop_k<-1

	while( (abs(p-op)>1e-5 || abs(q-oq)>1e-5 || abs(D-oD)>1e-5 )&& loop_k < 100)
	{
		op <- p;
		oq <- q;
		oD <- D;

		qtl.prob <- NP1.get_qtl_prob( marker, NULL, 0, 0, list(p=p, q=q, D=D) );
		A <- obj.curve$mlefunc(y, par.curve, sig.inv, sig.det, qtl.prob, time.std);
		if (is.na(A) || is.nan(A) || A$val==Inf || A$val==-Inf)
			return(NaN)

		QQ <- A$pf/rowSums(A$pf);

		#MM
        ii2 <- seq(nn)[marker==2]
        mm2 <- length(ii2)

		#Mm
        ii1 <- seq(nn)[marker==1]
        mm1 <- length(ii1)

		#mm
        ii0 <- seq(nn)[marker==0]
        mm0 <- length(ii0)

        theta<- p11*p00/(p11*p00+p10*p01)
        #p11 <- 1/(2*nn)*(sum( 2*QQ[ii2,1] +   QQ[ii1,2]) + sum(QQ[ii2,1] + theta * QQ[ii2,2]))
        #p10 <- 1/(2*nn)*(sum(   QQ[ii2,2] + 2*QQ[ii1,3]) + sum(QQ[ii2,3] + (1-theta) * QQ[ii2,2]))
        #p01 <- 1/(2*nn)*(sum( 2*QQ[ii1,1] +   QQ[ii3,2]) + sum(QQ[ii2,1] + (1-theta) * QQ[ii2,2]))
        #p00 <- 1/(2*nn)*(sum(   QQ[ii3,2] + 2*QQ[ii3,1]) + sum(QQ[ii2,3] + theta * QQ[ii2,2]))

        p11 <- 1/(2*nn)*(sum(  2*QQ[ii2,3] +   QQ[ii2,2]) + sum(QQ[ii1,3] + theta * QQ[ii1,2]))
        p10 <- 1/(2*nn)*(sum(    QQ[ii2,2] + 2*QQ[ii2,1]) + sum(QQ[ii1,1] + (1-theta) * QQ[ii1,2]))
        p01 <- 1/(2*nn)*(sum(  2*QQ[ii0,3] +   QQ[ii0,2]) + sum(QQ[ii1,3] + (1-theta) * QQ[ii1,2]))
        p00 <- 1/(2*nn)*(sum(    QQ[ii0,2] + 2*QQ[ii0,1]) + sum(QQ[ii1,1] + theta * QQ[ii1,2]))

        p   <- p11+p10
        q   <- p11+p01
        D   <- p11*p00-p01*p10

		loop_k <- loop_k + 1
		if (p<.Machine$double.eps) p<-.Machine$double.eps;
		if (q<.Machine$double.eps) q<-.Machine$double.eps;
		if (D<.Machine$double.eps) D<-.Machine$double.eps;

		if (p>=1) p<-0.99;
		if (q>=1) q<-0.99;
		if (D>=1) D<-0.99;

		cat("NP1(", loop_k, ")", p, q, D, p11,p10,p01,p00, "\n");
	}

	return(list(p=p, q=q, D=D));
}


cross_NP1<-list(
	type 	= "NP1",
	desc 	= "NP1",
	hp_desc = "mu_QQ == mu_Qq == mu_qq",
	gen_QQ  = T,
	gen_Qq  = T,
	gen_qq  = T,
	gen_num = 3,
	par_num = 3,
	single_marker	= T,
	get_gen_par  = function(par){ fin.get_gen_par(par, 3, T, T, T)},
	get_init_rand= NP1.get_init_rand,
	get_simu_marker = NP1.get_simu_marker,
	get_simu_qtl = NP1.get_simu_qtl,
	get_qtl_prob = NP1.get_qtl_prob,
	get_cond_prob= NP1.get_cond_prob,
	get_sum_par  = NP1.sum_par,
	get_est_param= NP1.get_est_param);


FM2.get_cross <- function(cross.type)
{
	FM2.crosss <- list();

	FM2.crosss[[length(FM2.crosss)+1]] <- cross_BC;
	FM2.crosss[[length(FM2.crosss)+1]] <- cross_F2;
	FM2.crosss[[length(FM2.crosss)+1]] <- cross_RIL;
	#FM2.crosss[[length(FM2.crosss)+1]] <- cross_NP1;

	for (i in 1:length(FM2.crosss))
		if (FM2.crosss[[i]]$type==cross.type)
			return(FM2.crosss[[i]])

	return(NULL);
}