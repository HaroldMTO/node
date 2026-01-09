library(mfnode)

corrpc = function(step,nsiter)
{
	corr = sprintf("C%d",1:nsiter)
	# C[jstep] or C[jiter]: better look for predictor
	indp = which(regexpr("^C",step) < 0)
	# no corrector at last step
	indp = indp[indp < length(step)]
	for (i in seq(nsiter)) step[indp+i] = paste(step[indp],corr[i],sep="_")

	step
}

corrtl = function(step)
{
	lambda = getvar("LAMBDA",nd)
	if (length(lambda) > 1) {
		lambda = abs(lambda[lambda <= 0 & lambda == as.integer(lambda)])
	}

	if (length(lambda) == 0) return(step)

	id = which(duplicated(step))
	nt = length(id)/length(lambda)
	stopifnot(nt == as.integer(nt))
	step[id] = sprintf("%s_%s",paste("L",rep(lambda,each=nt),sep=""),step[id])

	step
}

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])
cargs[sapply(cargs,is.null)] = ""

nd = readLines(cargs$fic1,skipNul=TRUE)
has.fc = any(regexpr("START CNT4",nd) > 0)

if (is.null(cargs$re)) {
	sp1 = spnorm(nd,lev=0)
} else {
	sp1 = spnorm(nd,lev=0,cargs$re)
}

if (is.null(sp1)) {
	cat("--> no SP norm in 1st file, pattern",cargs$re,"\n")
	q("no")
}

nfrsdi = getvar(".+ NFRSDI",nd)
step1 = dimnames(sp1)[[1]]
ts1 = getvar("TSTEP",nd)
cat(". file1 - nb of steps:",dim(sp1)[1],"- norms frequency:",nfrsdi,"\n")

if (any(duplicated(step1))) {
	nsiter = getvar("NSITER",nd)
	if (nsiter > 0) {
		step1 = corrpc(step1,nsiter)
	} else {
		step1 = corrtl(step1)
	}
}

nd = readLines(cargs$fic2,skipNul=TRUE)

if (is.null(cargs$re)) {
	sp2 = spnorm(nd,lev=0)
} else {
	sp2 = spnorm(nd,lev=0,cargs$re)
}

if (is.null(sp2)) {
	cat("--> no SP norm in 2nd file, pattern",cargs$re,"\n")
	q("no")
}

nfrsdi = getvar(".+ NFRSDI",nd)
step2 = dimnames(sp2)[[1]]
ts2 = getvar("TSTEP",nd)
cat(". file 2 - nb of steps:",dim(sp2)[1],"- norms frequency:",nfrsdi,"\n")

if (any(duplicated(step2))) {
	nsiter = getvar("NSITER",nd)
	if (nsiter > 0) {
		step2 = corrpc(step2,nsiter)
	} else {
		step2 = corrtl(step2)
	}
}

if (ts2 != ts1) stop("different TSTEP values")

if (! has.fc) stopifnot(dim(sp1)[1] == dim(sp2)[1])

noms1 = dimnames(sp1)[[3]]
noms2 = dimnames(sp2)[[3]]
stopifnot(all(! duplicated(noms2)))

indv = match(noms1,noms2)
if (any(is.na(indv))) cat("missing variables :",noms1[is.na(indv)],"\n")

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

sp1 = sp1[na.omit(indt),,na.omit(indv),drop=FALSE]
sp2 = sp2[it,,iv,drop=FALSE]
ndiff = array(round(digitsdiff(sp1,sp2)),dim=dim(sp1))

noms = noms1[na.omit(indv)]
if (max(nchar(noms))*length(noms) > 80) noms = abbreviate(noms,6)
if (max(nchar(noms))*length(noms) > 80) noms = abbreviate(noms,5)
if (max(nchar(noms)) > 9) {
	fmt = "%9g"
	cat(" step",sprintf("%9s",noms),"\n")
} else if (max(nchar(noms)) > 7) {
	fmt = "%7g"
	cat(" step",sprintf("%7s",noms),"\n")
} else {
	fmt = "%6g"
	cat(" step",sprintf("%6s",noms),"\n")
}

nt = dim(sp1)[1]
if (all(! is.na(ndiff)) && all(ndiff == 0)) {
	for (i in seq(min(5,nt))) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,]),"\n")
	if (nt > 5) cat("...",nt-min(5,nt),"more 0 lines\n")
} else if (nt < 30) {
	for (i in seq(nt)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,]),"\n")
} else {
	for (i in seq(15)) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,]),"\n")
	cat("... (every",nt%/%30,"printed time-step)\n")
	ind = seq(15,nt,by=nt%/%30+1)[-1]
	for (i in ind) cat(format(step1[i],width=5),sprintf(fmt,ndiff[i,1,]),"\n")
}

