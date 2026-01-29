library(mfnode)

gpfre = sprintf("%s|DIV\\w*|VOR\\w*|E(TA)?DOT",gpfre)
gpfre = sprintf("%s|PRESS\\. DEPARTURE|VERT\\. DIVERGENCE|TERM X",gpfre)
gpfre = sprintf("%s|\\w+_?NL",gpfre)

gpnormerr = function(nd)
{
	ind = grep(sprintf("^ *\\w+:( +%s){3} *$",Gnum),nd)
	if (length(ind) == 0) return(NULL)

	indth = grep("OpenMP threads",nd)
	nth = 0
	if (length(indth) > 0) {
		nth = as.integer(sub(".+ OpenMP threads *= *(\\d+).*","\\1",nd[indth[1]]))
	}

	gpn = numlines(nd[ind])
	noms = unique(sub("^ *(\\w+): +.+","\\1",nd[ind]))
	indp = grep("\\| +stderr\\.\\d+ +\\|",nd)
	proc = as.integer(sub("\\| +stderr\\.(\\d+) +\\| *","\\1",nd[indp]))+1
	nt = length(gpn)/(3*length(noms)*length(proc))
	stopifnot(nt == as.integer(nt))
	if (length(proc) == 1) {
		gpn = array(gpn,c(3,nt,1,length(noms)),
			dimnames=list(c("avg","min","max"),NULL,0,noms))
		gpn = aperm(gpn,c(2,3,1,4))
	} else {
		dim(gpn) = c(3,nt,length(noms),length(proc))
		m = apply(gpn[1,,,],1:2,mean)
		n = apply(gpn[2,,,],1:2,min)
		x = apply(gpn[3,,,],1:2,max)
		gpn = array(c(m,n,x),c(nt,1,length(noms),3),
			dimnames=list(NULL,0,noms,c("avg","min","max")))
		gpn = aperm(gpn,c(1,2,4,3))
	}

	indt = grep("^ *\\d+:\\d+:\\d+ +STEP +\\d+ +H *=",nd)
	step = as.integer(sub("^ *\\d.+ +STEP +(\\d+) +H *=.*","\\1",nd[indt]))
	if (length(step) == 1) {
		if (nt > 1) step = rep(step,nt)
	} else {
		# account for stepX or DFI before step0
		ii = which(ind < indt[1])
		nevent0 = length(ii)/length(noms)

		# avoid step0 because of potential pre-CNT4 printouts (eg DFI)
		ii = which(indt[1] < ind & ind < indt[2])
		nevent1 = length(ii)/length(noms)

		if (nevent0 %% nevent1 == 0) {
			nx = nevent0/nevent1-1
			step = c(rep(-1,nx),step)
		}

		step = rep(step,each=nevent1)

		if (dim(gpn)[1] > length(step)) {
			nlast = dim(gpn)[1]-length(step)
			if (nlast == nevent1) step = c(step,rep(max(step)+1,nevent1))
		} else {
			# tailor last step because it often has no printouts
			ii = which(ind < max(indt))
			nevent = length(ii)/length(noms)
			if (length(step) > nevent) length(step) = nevent
		}
	}

	dimnames(gpn)[[1]] = step

	gpn
}

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])
cargs[sapply(cargs,is.null)] = ""

if (interactive()) browser()

nd = readLines(cargs$fic1,skipNul=TRUE)
ts1 = getvar("TSTEP",nd)
has.fc = any(regexpr("^ *START CNT4",nd) > 0)
if (! has.fc) cat("--> no forecast conf (cnt4) in 1st file\n")

if (is.null(cargs$re)) {
	nd = grep("^ *gpnorm gflt0",nd,invert=TRUE,ignore.case=TRUE,value=TRUE)
	gp1 = gpnorm(nd,lev=0,gpout=gpfre)
} else if (! nzchar(cargs$re) || regexpr("gpnorm gfl",cargs$re,ignore.case=TRUE) > 0) {
	gp1 = gpnorm(nd,lev=0,cargs$re,gpout=gpfre)
} else {
	gp1 = gpnorm(nd,lev=0,cargs$re,gpfre)
}

if (is.null(gp1)) {
	cat("--> no GP norms in 1st file\n")
	q("no")
}

nfrgdi = getvar(".+ NFRGDI",nd)
cat("nb of steps, file 1:",dim(gp1)[1],"- norms frequency:",nfrgdi,"\n")
step1 = dimnames(gp1)[[1]]

nd = readLines(cargs$fic2,skipNul=TRUE)
ts2 = getvar("TSTEP",nd)

if (ts1 != ts2) stop("different TSTEP")

if (is.null(cargs$re)) {
	nd = grep("^ *gpnorm gflt0",nd,invert=TRUE,ignore.case=TRUE,value=TRUE)
	gp2 = gpnorm(nd,lev=0,gpout=gpfre)
} else if (! nzchar(cargs$re) || regexpr("gpnorm gfl",cargs$re,ignore.case=TRUE) > 0) {
	gp2 = gpnorm(nd,lev=0,cargs$re,gpout=gpfre)
} else {
	gp2 = gpnorm(nd,lev=0,cargs$re,gpfre)
}

if (is.null(gp2)) {
	cat("--> no GP norms in 2nd file\n")
	q("no")
}

nfrgdi = getvar(".+ NFRGDI",nd)
cat("nb of steps, file 2:",dim(gp2)[1],"- norms frequency:",nfrgdi,"\n")
step2 = dimnames(gp2)[[1]]

noms1 = dimnames(gp1)[[4]]
noms2 = dimnames(gp2)[[4]]

indv = match(noms1,noms2)
if (any(is.na(indv))) cat("missing variables in 2nd file :",noms1[is.na(indv)],"\n")

indv = match(noms2,noms1)
iv = which(noms2 %in% noms1)
if (any(is.na(indv))) cat("new variables :",noms2[is.na(indv)],"\n")
if (length(iv) == 0) {
	cat("variables (1):",noms1,"\n")
	cat("variables (2):",noms2,"\n")
	stop("no variables in common to compare\n")
}

indt = match(step2,step1)
it = which(step2 %in% step1)
if (length(it) == 0) {
	cat("steps:",length(step1),length(step2),"\n")
	stop("no steps in common to compare\n")
}

if (length(step1) != length(step2)) {
	nt = min(length(step1),length(step2))
	cat("--> different number of steps in files, limiting norms to",nt,"1st ones\n")
}

it = which(! is.na(indt))
gp1 = gp1[na.omit(indt),,,na.omit(indv),drop=FALSE]
gp2 = gp2[it,,,iv,drop=FALSE]
step1 = step1[na.omit(indt)]

env = Sys.getenv("DEBUG_EGPNORM")
if (nzchar(env) && env != "0") {
	ndglg = getvar("NDGLG",nd)
	ndlon = getvar("NDLON",nd)
	cat("apply avg norm ratio on ref:",ndlon/ndglg,"\n")
	gp1[,1,1,] = gp1[,1,1,]*ndlon/ndglg
}

fpm = function(x) {
	xx = max(abs(x))
	x[1] = (x[1]+xx)-xx
	x
}

if (! is.null(cargs$ftz) && as.logical(cargs$ftz)) {
	gp1 = apply(gp1,c(1,2,4),fpm)
	gp1 = aperm(gp1,c(2,3,1,4))
	gp2 = apply(gp2,c(1,2,4),fpm)
	gp2 = aperm(gp2,c(2,3,1,4))
}

ndiff = array(round(digitsdiff(gp1,gp2)),dim=dim(gp1))

mnx = "mnx" %in% names(cargs) && as.logical(cargs$mnx)

noms1 = noms1[na.omit(indv)]
nt = dim(gp1)[1]
nvar = dim(ndiff)[4]

for (j in seq(1,nvar,by=10)) {
   indv = seq(j,min(j+9,nvar))

   noms = noms1[indv]

	if (max(nchar(noms))*length(noms) > 80) noms = abbreviate(noms,8)
	if (max(nchar(noms)) > 10) {
		fmt = "%10s"
	} else if (max(nchar(noms)) > 7 || mnx) {
		fmt = "%8s"
	} else if (max(nchar(noms)) > 5) {
		fmt = "%6s"
	} else {
		fmt = "%5s"
	}

	cat(" step",sprintf(fmt,noms),"\n")
	ndf = ndiff[,,,indv,drop=FALSE]
	if (all(ndf == 0,na.rm=TRUE)) {
		ind = seq(min(nt,5))
	} else {
		ind = seq(min(nt,15))
	}

	if (mnx) {
		sdiff = apply(ndf,c(1,4),function(x) paste(sprintf("%g",x[1,]),collapse="/"))
		for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,sdiff[i,]),"\n")
	} else {
		for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,ndf[i,1,1,]),"\n")
	}

	if (all(ndf == 0,na.rm=TRUE)) {
		if (nt > length(ind)) cat("...",nt-length(ind),"more 0 lines\n")
	} else if (length(ind) > nt) {
		if (nt > 30) {
			cat("... (every",nt%/%30,"printed time-step)\n")
			ind = seq(length(ind),nt,by=nt%/%30)[-1]
		} else {
			ind = seq(length(ind),nt)[-1]
		}

		if (mnx) {
			for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,sdiff[i,]),"\n")
		} else {
			for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,ndf[i,1,1,]),"\n")
		}
	}

	if (any(is.na(ndf))) {
		ind = apply(ndf,4,function(x) any(is.na(x)))
		cat("Warning: NaN for variables",noms[ind],"\n")
	}
}

if (any(regexpr("TL|AD",step1) > 0) && ! is.null(cargs$re) &&
	regexpr("gpnorm g(mv|fl)t0 traj",cargs$re) > 0) {
	gpnl = gp1[grep("TL|AD",step1,invert=TRUE),,,,drop=FALSE]
	gptl = gp1[regexpr("TL",step1) > 0,,,,drop=FALSE]
	gpad = gp1[rev(regexpr("AD",step1) > 0),,,,drop=FALSE]

	cat("+ NL/TL comparison:\n")
	ndiff = array(round(diffnorm(gpnl,gptl)),dim=dim(gpnl))
	cat(" step",sprintf(fmt,noms1),"\n")
	nt = dim(gpnl)[1]
	if (all(ndiff == 0)) {
		for (i in seq(min(5,nt))) {
			cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
		}
	} else if (nt < 30) {
		for (i in seq(nt)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
	} else {
		for (i in seq(15)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
		cat("...\n")
		ind = seq(15,nt,by=nt%/%30+1)[-1]
		for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
	}

	cat("+ NL/AD comparison:\n")
	ndiff = array(round(diffnorm(gpnl,gptl)),dim=dim(gpnl))
	cat(" step",sprintf(fmt,noms1),"\n")
	nt = dim(gpnl)[1]
	if (all(ndiff == 0)) {
		for (i in seq(min(5,nt))) {
			cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
		}
	} else if (nt < 30) {
		for (i in seq(nt)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
	} else {
		for (i in seq(15)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
		cat("...\n")
		ind = seq(15,nt,by=nt%/%30+1)[-1]
		for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,1,]),"\n")
	}
}

