library(mfnode)

getarg = function(x,args)
{
	ind = grep(sprintf("\\<%s=",x),args)
	if (length(ind) == 0) return(NULL)

	strsplit(sub(sprintf("\\<%s=",x),"",args[ind]),split=":")[[1]]
}

abh = function(nd,nflevg)
{
	ih = grep("A and B (at half levels|on half layers)",nd)
	ind = ih+1+seq(nflevg+1)
	snum = "-?\\d+\\.\\d+"
	if (regexpr("JLEV +A +B +ETA +ALPHA",nd[ih+1]) > 0) {
		# take only 4 values (DELB is not present at level 0)
		re = sprintf(" *\\d+ +(%s) +(%s) +(%s) +(%s)",snum,snum,snum,snum)
		ire = regexec(re,nd[ind])
		alh = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))
		ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	} else {
		ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",snum,snum,snum),nd[ind])
		alh = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
		ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	}

	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	data.frame(Ah=ah,Bh=bh,alpha=alh)
}

abhfp = function(nd,nfplev)
{
	ih = grep("Set up F-post processing, vertical geometry",nd)
	ind = grep("FPVALH.+FPVBH",nd)
	ind = ind[ind > ih]
	stopifnot(length(ind) == nfplev+1)

	snum = "-?\\d+\\.\\d+"
	ire = regexec(sprintf(" *FPVALH\\(( *\\d+|\\*)+\\) *= *(%s) +FPVBH\\(( *\\d+|\\*+)\\) *= *(%s)",
		snum,snum),nd[ind])
	ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))
	data.frame(Ah=ah,Bh=bh)
}

args = commandArgs(trailingOnly=TRUE)

hasx11 = is.null(getarg("png",args)) && capabilities("X11")
ask = hasx11 && interactive()
if (hasx11) {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) {
		options(device.ask.default=TRUE)
	} else {
		cat("--> sending plots to Rplots.pdf\n")
	}
}

xaxis = data.frame(unit=c(1,60,3600,86400),label=sprintf("fc time (%s)",
	c("s","mn","h","days")),mindiff=c(0,240,14400,6*86400),freq=c(6,1,6,1))

pngd = getarg("png",args)
if (is.null(pngd)) pngd = "."
lev = as.integer(getarg("lev",args))
if (is.null(lev)) lev = 0L
hmin = as.numeric(getarg("hmin",args))
hmax = as.numeric(getarg("hmax",args))
ptype = getarg("type",args)
fpre = getarg("fpre",args)
fpref = getarg("fpref",args)
if (is.null(fpref)) fpref = "full-pos spnorms"

fnode = grep("=",args,invert=TRUE,value=TRUE)

cat("Read file",fnode[1],"\n")
nd = readLines(fnode[1])
nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
nflevg = getvar("NFLEVG",nd)
nfp3s = getvar("NFP3S",nd)
has.levels = getvar(".*NSPPR",nd) > 0

grouplev = NULL
ab = abh(nd,nflevg)
abfp = abhfp(nd,nfp3s)
vp00 = 101325
etah = ab$alpha+ab$Bh
eta = (etah[-1]+etah[-(nflevg+1)])/2
etahfp = abfp$Ah/vp00+abfp$Bh
etafp = (etahfp[-1]+etahfp[-(nfp3s+1)])/2

lev1 = lev
if (length(lev) > 1) {
	grouplev = lev
	lev1 = sapply(etafp,function(x) which.min(abs(x-eta)))
	lev = seq(nfp3s)
} else if (lev == -1) {
	if (! has.levels) {
		cat("--> no level norms\n")
	}

	lev0 = seq(nflevg)
	lev1 = sapply(etafp,function(x) which.min(abs(x-eta)))
	lev = seq(nfp3s)
}

tstep = getvar("TSTEP",nd)
nstop = getvar("NSTOP",nd)

if (interactive()) browser()

cat("Parse FP norms of type",ptype,"\n")
fp1 = fpspnorm(nd,lev)

if (is.null(fp1)) {
	cat("--> no FP norms, quit\n")
	quit("no")
}

i0 = apply(fp1,1,function(x) all(x==0))
if (any(i0) && ! all(i0)) {
	cat("--> FP norms all 0 for some steps, removed\n")
	fp1 = fp1[-which(i0),,,,drop=FALSE]
}

step = dimnames(fp1)[[1]]
cat(". steps:",head(step[-length(step)]),"...",step[length(step)],"\n")
if (nstop == 0 && dim(fp1)[1] > 1) cat("--> steps are events of the job\n")

ix = grep("^X",step)
if (length(ix) == length(step)) {
	cat("--> FP norms for STEPX only, quit\n")
	quit("no")
} else if (length(ix) > 0) {
	cat("--> FP norms present for STEPX, remove",length(ix),"steps\n")
	fp1 = fp1[-ix,,,,drop=FALSE]
	step = dimnames(fp1)[[1]]
}

istep = as.numeric(gsub("C(\\d+)","\\1.5",step))

times = tstep*istep

fpl = list(fp1)

ndim = length(dim(fp1))
fpnoms = dimnames(fp1)[[ndim]]

if (is.null(fpre) && length(fnode) == 1) {
	leg = "FPSP"
	cat("Select times and variables among:
times:",times,"
vars:",fpnoms,"\n")
	if (length(lev) > 1) {
		sp1 = spnorm(nd,lev0,abbrev=FALSE)
	} else {
		sp1 = spnorm(nd,lev1,abbrev=FALSE)
	}

	if (! is.null(sp1)) {
		indt = match(times,dimnames(sp1)[[1]])
		indv = match(fpnoms,dimnames(sp1)[[3]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"SP")
			sp1 = sp1[indt,,indv,drop=FALSE]
			if (length(lev) > 1) {
				cat("--> augment fp1 to",nflevg,"levels\n")
				fp0 = array(dim=c(dim(fp1)[1],nflevg,dim(fp1)[-(1:2)]))
				dimnames(fp0) = c(dimnames(fp1)[1],list(seq(nflevg)),dimnames(fp1)[-(1:2)])
				fp0[,lev1,] = fp1
				fp1 = fp0
				fpl = list(fp1)
				lev = lev0
			}

			fpl = c(fpl,list(sp1))
		} else {
			cat("--> SP times:",dimnames(sp1)[[1]],"\n")
			cat("--> SP vars:",dimnames(sp1)[[3]],"\n")
		}
	}

	if (length(leg) == 1) {
		cat("--> no mixed spectral norms, quit\n")
		quit("no")
	}
} else {
	if (length(fnode) > 1) {
		leg = sub("node\\.?","",fnode,ignore.case=TRUE)
		fpre = rep(fpref,length(fnode)-1)
	} else {
		leg = sub("spnorm +","",c(fpref,fpre))
	}

	for (i in seq(along=fpre)) {
		if (length(fnode) > 1) {
			cat("Read file",fnode[i+1],"\n")
			nd = readLines(fnode[i+1])
		}

		cat("Parse FP norms, pattern",fpre[i],"\n")
		fpi = fpspnorm(nd,lev,fpre[i])

		if (is.null(fpi)) {
			cat("--> no norms for pattern",fpre[i],"\n")
			next
		}

		i0 = apply(fpi,1,function(x) all(x==0))
		if (all(i0)) {
			cat("--> FP all 0 for all steps, continue\n")
			next
		}

		noms = dimnames(fpi)[[ndim]]
		indv = match(fpnoms,noms)
		nc = max(nchar(noms))
		if (max(nchar(fpnoms)) != nc) {
			fpnoms1 = substr(fpnoms,1,nc)
			indv1 = match(fpnoms1,noms)
			if (length(which(is.na(indv1))) < length(which(is.na(indv)))) {
				indv = indv1
				fpnoms = fpnoms1
			}
		}
		stopifnot(any(! is.na(indv)))
		iv = ! noms %in% fpnoms
		if (any(iv)) cat("--> dropped variables:",noms[iv],"\n")

		stepi = dimnames(fpi)[[1]]
		ix = grep("^X",stepi)
		if (length(ix) > 0) {
			fpi = fpi[-ix,,,,drop=FALSE]
			stepi = dimnames(fpi)[[1]]
		}

		inds = match(step,stepi)
		stopifnot(any(! is.na(inds)))
		fpl[[i+1]] = fpi[inds,,,indv,drop=FALSE]
	}
}

leg = leg[which(! sapply(fpl,is.null))]
fpl = fpl[! sapply(fpl,is.null)]

nf = length(fpnoms)
nt = dim(fp1)[1]

if (length(lev) > 1) {
	cat("Produce vertical profiles\n")
	tt = sprintf("FP norm of %s",fpnoms)
	indt = c(1,(nt+1)%/%2,nt)
	ts = sprintf("t%s",paste(istep[indt],collapse="/"))
	nc = min(nf,2)
	nr = min(length(fpl),2)
	nj = nc

	for (i in seq((nf-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/%snormv%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfcol=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),cex=.83)

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			cat(". FP field",fpnoms[j],"in",ficpng,"\n")
			y = t(fp1[indt,,j])
			plotvmean(y,lev,type="o",pch=".",cex=1.1,xlab=fpnoms[j],main=c(tt[j],ts))

			if (nr == 1) next

			y = sapply(fpl,function(x) x[indt[2],,j],simplify="array")
			il = which(apply(y,2,function(x) any(! is.na(x))))
			ts2 = sprintf("t%s",istep[indt[2]])
			plotvmean(y,lev,xlab=fpnoms[j],main=c(tt[j],ts2),lty=c(1,2,2),col=il,
				type="o",pch=".",cex=1.1,legend=leg[il])
		}

		dev.off()
	}
} else if (length(lev) == 1 && nt == 1) {
	cat("1 time-step only\n")
} else {
	cat("Produce time-series\n")
	tt = sprintf("FP norm of %s",fpnoms)
	nc = 2
	nr = min(1+(nf-1)%/%2,2)
	nj = nr*nc

	for (i in seq((nf-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/%snorm%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),cex=.83)

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			cat(". FP field",fpnoms[j],"\n")
			y = sapply(fpl,function(x) x[,j],simplify="array")
			plotmean(times,y,tt[j],ylab=fpnoms[j])
		}

		dev.off()
	}
}
