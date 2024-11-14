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
	fpn[,ind,,] = fp

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
detail = regexpr("fp.*detail",ptype) > 0
fpre = getarg("fpre",args)
fpref = getarg("fpref",args)

fnode = grep("=",args,invert=TRUE,value=TRUE)

cat("Read file",fnode[1],"\n")
nd = readLines(fnode[1])
nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
nflevg = getvar("NFLEVG",nd)
nfplev = getvar("NFPLEV",nd)
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
	abfp = abhfp(nd,nfplev)
	etahfp = abfp$Ah/vp00+abfp$Bh
	etafp = (etahfp[-1]+etahfp[-(nfp3s+1)])/2
}

levf = lev
if (length(lev) > 1) {
	grouplev = lev
	levf = sapply(etafp,function(x) which.min(abs(x-eta)))
	lev = seq(nfp3s)
} else if (lev == -1) {
	if (! has.levels) {
		cat("--> no level norms\n")
		q("no")
	}

	lev0 = seq(nflevg)
	levf = sapply(etafp,function(x) which.min(abs(x-eta)))
	lev = seq(nfp3s)
}

if (! identical(lev,0L)) stopifnot(has.levels)

tstep = getvar("TSTEP",nd)
nstop = getvar("NSTOP",nd)

if (interactive()) browser()

cat("Parse FP norms of type",ptype,"\n")
if (is.null(fpref)) {
	# take all FP events (only possible if no 'gpnorm' present, ie mixed levels)
	fp1 = fpgpnorm(nd,lev,quiet=TRUE)
} else if (regexpr("dynfpos z",fpref) > 0) {
	fp1 = gpnorm(nd,levf,fpref,abbrev=FALSE)
} else if (regexpr("dynfpos",fpref) > 0) {
	fp1 = fpgpnorm(nd,levf,fpref,quiet=TRUE)
} else if (regexpr("allfpos",fpref) > 0) {
	# special case when 'gpnorm' present (mixed levels)
	fp1 = fpgpnorm(nd,lev,tag="gpnorm",invert=TRUE,quiet=TRUE)
} else {
	fp1 = fpgpnorm(nd,lev,fpref,quiet=TRUE)
}

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
fpnoms = dimnames(fp1)[[4]]

if (is.null(fpre) && length(fnode) == 1) {
	leg = "FPGP"
	cat("Select times and variables among:
times:",times,"
vars:",fpnoms,"\n")
	fpsp1 = fpspnorm(nd,lev)
	if (! is.null(fpsp1)) {
		indt = match(times,dimnames(fpsp1)[[1]])
		indv = match(fpnoms,dimnames(fpsp1)[[3]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"FPSP")
			fpsp1 = fpsp1[indt,,indv,drop=FALSE]
			fpl = c(fpl,list(fpsp1))
		} else {
			cat("--> FPSP times:",dimnames(fpsp1)[[1]],"\n")
			cat("--> FPSP vars:",dimnames(fpsp1)[[3]],"\n")
		}
	}

	if (length(lev) > 1) {
		sp1 = spnorm(nd,lev0,abbrev=FALSE)
	} else {
		sp1 = spnorm(nd,levf,abbrev=FALSE)
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
				fpl = lapply(fpl,augmentlev,nflevg,levf)
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
		gp1 = gpnorm(nd,levf,gpout=gpfre,abbrev=FALSE)
	}

	if (! is.null(gp1)) {
		indt = match(times,dimnames(gp1)[[1]])
		indv = match(fpnoms,dimnames(gp1)[[4]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"GFL")
			gp1 = gp1[indt,,,indv,drop=FALSE]
			if (length(lev) > 1) {
				cat("--> augment fpl to",nflevg,"levels\n")
				fpl = lapply(fpl,augmentlev,nflevg,levf)
				lev = lev0
			}

			fpl = c(fpl,list(gp1))
		} else {
			cat("--> GP times:",dimnames(gp1)[[1]],"\n")
			cat("--> GP vars:",dimnames(gp1)[[4]],"\n")
		}
	}

	if (length(leg) == 1) {
		cat("--> no mixed norms, quit\n")
		quit("no")
	}
} else {
	if (length(fnode) > 1) {
		leg = sub("node\\.?","",fnode,ignore.case=TRUE)
		fpre = rep(fpref,length(fnode)-1)
	} else {
		leg = sub("gpnorm +","",c(fpref,fpre))
	}

	for (i in seq(along=fpre)) {
		if (length(fnode) > 1) {
			cat("Read file",fnode[i+1],"\n")
			nd = readLines(fnode[i+1])
		}

		cat("Parse FP norms, pattern",fpre[i],"\n")
		if (regexpr("dynfpos z",fpre[i]) > 0) {
			fpi = gpnorm(nd,levf,fpre[i],abbrev=FALSE)
		} else if (regexpr("dynfpos",fpre[i]) > 0) {
			fpi = fpgpnorm(nd,levf,fpre[i],quiet=TRUE)
		} else if (regexpr("allfpos",fpre[i]) > 0) {
			# special case when 'gpnorm' present (mixed levels)
			fpi = fpgpnorm(nd,lev,tag="gpnorm",invert=TRUE,quiet=TRUE)
		} else {
			fpi = fpgpnorm(nd,lev,fpre[i],quiet=TRUE)
		}

		if (is.null(fpi)) {
			cat("--> no norms for pattern",fpre[i],"\n")
			next
		}

		i0 = apply(fpi,1,function(x) all(x==0))
		if (all(i0)) {
			cat("--> FP all 0 for all steps, continue\n")
			next
		}

		noms = dimnames(fpi)[[4]]
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
	nc = 3
	nr = min(length(fpl),2)
	nj = 1

	for (i in seq((nf-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/%snormv%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			cat(". FP field",fpnoms[j],"in",ficpng,"\n")
			y = aperm(fpl[[1]][indt,,,j],c(2,1,3))
			plotvmnx(y,lev,xlab=fpnoms[j],main=c(tt[j],ts))

			if (nr == 1) next

			yn = sapply(fpl,function(x) x[indt[2],,,j],simplify="array")
			yn = aperm(yn,c(1,3,2))
			il = which(apply(yn,2,function(x) any(! is.na(x))))
			ts2 = sprintf("t%d",indt[2])
			plotvmnx(yn[,il,,drop=FALSE],lev,xlab=fpnoms[j],main=c(tt[j],ts2),legend=leg[il],
				lty=2,col=il)
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

	if (length(fpl) > 1) {
		con = file(sprintf("%s/%s.txt",pngd,ptype),"wt")
		for (i in seq(along=fpl)[-1]) {
			cat("\nNorm diff. tendency for patterns:",sprintf("'%s'",leg[(i-1):i]),"\n",
				file=con)
			iv = which(! is.na(dimnames(fpl[[i]])[[4]]))
			gpdiff = fpl[[i]][,1,1,iv,drop=FALSE]-fpl[[i-1]][,1,1,iv,drop=FALSE]
			a = apply(gpdiff,4,function(x) coef(line(x))[2])
			s = apply(gpdiff,4,function(x) sd(residuals(line(x))))
			evo = rbind(tend=a,sd=s)
			evo = format(evo,digits=3,width=12)
			cat("field",sprintf("% 12s",dimnames(evo)[[2]]),"\n",file=con)
			write.table(evo,con,quote=FALSE,row.names=format(dimnames(evo)[[1]],width=5),
				col.names=FALSE)
		}

		close(con)
	}

	iu = findInterval(diff(range(times)),xaxis$mindiff)
	tunit = xaxis$unit[iu]
	xlab = xaxis$label[iu]
	tfreq = xaxis$freq[iu]

	ttime = times/tunit

	xlim = range(ttime)
	if (length(hmin) == 1) xlim[1] = max(xlim[1],hmin*3600/tunit)
	if (length(hmax) == 1) xlim[2] = min(xlim[2],hmax*3600/tunit)

	it = which(ttime >= xlim[1] & ttime <= xlim[2])
	ttime = ttime[it]
	x = pretty(ttime[it]/tfreq,7)*tfreq
	xaxp = c(range(x),length(x)-1)
	fpl = lapply(fpl,function(x) x[it,,,,drop=FALSE])

	cat("FP norms",fpref,"\n")
	fpnoms = dimnames(fpl[[1]])[[4]]
	titre = paste("FP norm of",fpnoms)
	if (lev > 0) titre = paste(titre,"- lev",lev)

	nf = length(fpnoms)
	nc = 2
	ltymnx = c(1,3,3)
	if (detail) {
		# 2 rows, 2 columns, 1 field per figure
		nr = 2
		nj = 1
	} else {
		# 3 rows max, 2 columns, nr fields per figure (1 row by field)
		nr = max(2,min(3,nf))
		nj = nr
	}

	for (i in seq((nf-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/%snorm%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			y = sapply(fpl,function(x) x[,1,,j],simplify="array")
			il = which(apply(y,3,function(x) any(! is.na(x))))
			scal = 1/scale10(y[,1,])
			if (.001 <= scal && scal < 1 || is.infinite(scal)) scal = 1
			plotmean(ttime,y[,1,il],main=titre[j],leg[il],tunit,xlim=xlim,
				xlab=xlab,ylab=fpnoms[j],xaxp=xaxp,scale=scal,col=il)

			plotmnx(ttime,y,titre[j],xlim=xlim,xlab=xlab,ylab=fpnoms[j],xaxp=xaxp)

			if (! detail) next

			plotmnx(ttime,y,titre[j],imnx=2,xlim=xlim,xlab=xlab,ylab=fpnoms[j],xaxp=xaxp)
			plotmnx(ttime,y,titre[j],imnx=3,xlim=xlim,xlab=xlab,ylab=fpnoms[j],xaxp=xaxp)
		}

		dev.off()
	}
}

