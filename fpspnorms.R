library(mfnode)

getarg = function(x,args)
{
	ind = grep(sprintf("\\<%s=",x),args)
	if (length(ind) == 0) return(NULL)

	strsplit(sub(sprintf("\\<%s=",x),"",args[ind]),split=":")[[1]]
}

augmentlev = function(fp,nlev,ind)
{
	stopifnot(all(ind %in% seq(nlev)))

	if (dim(fp)[2] == nlev) return(fp)

	fpn = array(dim=c(dim(fp)[1],nlev,dim(fp)[-(1:2)]))
	dimnames(fpn) = c(dimnames(fp)[1],list(seq(nlev)),dimnames(fp)[-(1:2)])
	fpn[,ind,] = fp

	fpn
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
vp00 = 101325
etah = ab$alpha+ab$Bh
eta = (etah[-1]+etah[-(nflevg+1)])/2
if (is.null(nfp3s)) {
	etafp = eta
	nfp3s = nflevg
} else {
	abfp = abhfp(nd,nfp3s)
	etahfp = abfp$Ah/vp00+abfp$Bh
	etafp = (etahfp[-1]+etahfp[-(nfp3s+1)])/2
}

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

	# some small fix for different names (but same variable)
	spnoms = dimnames(sp1)[[3]]
	ind = match("HUMIDITY",spnoms)
	if (! is.na(ind)) dimnames(sp1)[[3]][ind] = "HUMI.SPECIFI"

	if (! is.null(sp1)) {
		indt = match(times,dimnames(sp1)[[1]])
		indv = match(fpnoms,dimnames(sp1)[[3]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"SP")
			sp1 = sp1[indt,,indv,drop=FALSE]
			if (length(lev) > 1) {
				cat("--> augment fpl to",nflevg,"levels\n")
				fpl = lapply(fpl,augmentlev,nflevg,lev1)
				lev = lev0
			}

			fpl = c(fpl,list(sp1))
		} else {
			cat("--> SP times:",dimnames(sp1)[[1]],"\n")
			cat("--> SP vars:",dimnames(sp1)[[3]],"\n")
		}
	}

	if (length(lev) > 1) {
		gp1 = gpnorm(nd,lev0,gpout=gpfre,abbrev=FALSE)
	} else {
		gp1 = gpnorm(nd,lev1,gpout=gpfre,abbrev=FALSE)
	}

	if (! is.null(gp1)) {
		indt = match(times,dimnames(gp1)[[1]])
		indv = match(fpnoms,dimnames(gp1)[[4]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"GFL")
			# let's drop dim 3 only
			dd = dim(gp1)
			dn = dimnames(gp1)
			gp1 = gp1[indt,,1,indv]
			if (length(dim(gp1)) != 3) gp1 = array(gp1,dd[-3],dn[-3])
			if (length(lev) > 1) {
				cat("--> augment fpl to",nflevg,"levels\n")
				fpl = lapply(fpl,augmentlev,nflevg,lev1)
				lev = lev0
			}

			fpl = c(fpl,list(gp1))
		} else {
			cat("--> GP times:",dimnames(gp1)[[1]],"\n")
			cat("--> GP vars:",dimnames(gp1)[[4]],"\n")
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
nt = dim(fpl[[1]])[1]

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
			y = t(fpl[[1]][indt,,j])
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

	if (! is.null(grouplev)) {
		cat("--> groups of levels:",grouplev,"\n")
		indg = 0
		leg = character()
		fp1 = fpl[[1]]
		for (i in seq(along=grouplev)) {
			indi = which(lev <= grouplev[i])
			indg = indi[indi > max(indg)]
			cat(".. level indices:",range(indg),"\n")
			leg[i] = sprintf("L %d-%d",min(indg),max(indg))
			fpi = fp1[,1,,,drop=FALSE]
			fpi[,,1,] = apply(fp1[,indg,1,,drop=FALSE],c(1,3,4),mean)
			fpi[,,2,] = apply(fp1[,indg,2,,drop=FALSE],c(1,3,4),min)
			fpi[,,3,] = apply(fp1[,indg,3,,drop=FALSE],c(1,3,4),max)
			fpl[[i]] = fpi
		}

		lev = 0
	}
}

if (length(lev) == 1) {
	cat("Produce time-series, level",lev,"\n")
	if (nt == 1) {
		cat("1 time-step only\n")
		quit("no")
	}

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
			plotmean(times,y,main=tt[j],ylab=fpnoms[j])
		}

		dev.off()
	}
}
