Gnum = "[-+]?[0-9][-+0-9.E]+[0-9]"
Gnum = "-?\\d*\\.\\d+([eE]?[-+]?\\d+)?\\>"
Gint = "-?\\d+\\>"

getarg = function(x,args)
{
	ind = grep(sprintf("\\<%s=",x),args)
	if (length(ind) == 0) return(NULL)

	strsplit(sub(sprintf("\\<%s=",x),"",args[ind]),split=":")[[1]]
}

getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(%s|%s).*",var,sep,Gint,Gnum)
	unique(as.numeric(gsub(re,"\\1",grep(re,nd,value=TRUE))))
}

timestep = function(nd)
{
	s = grep("^ *TSTEP += *[-+]?[0-9][-+0-9E.]+",unlist(nd),value=TRUE)
	stopifnot(length(s) == 1)

	tstep = as.numeric(sub("TSTEP += *([-+]?[0-9][-+0-9E.]+) .+","\\1",s))
}

spnorm = function(fnode,time,nmax=NA,lev)
{
	nd = readLines(fnode)
	nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
	ind = grep("SPECTRAL NORMS",nd)

	inds = grep("NORMS AT NSTEP CNT4",nd[ind-1])
	ind = ind[inds]
	istep = as.integer(sub("NORMS AT NSTEP CNT4( \\(PREDICTOR\\))? +(\\d+)","\\2",nd[ind-1]))

	#ioro = grep("OROGRAPHY",nd[ind])
	#if (length(ioro) > 0) ind = ind[-ioro]

	indl = grep("^ *NFLEVG += *[0-9]+",nd)
	nl = as.integer(sub("^ *NFLEVG += *([0-9]+) .+","\\1",nd[indl]))
	has.levels = regexpr("^ *1 +[-+]?[0-9][-+0-9.E]+[0-9]",nd[ind[1]+3]) > 0 &&
		regexpr("^ *\\d+ +[-+]?[0-9][-+0-9.E]+[0-9]",nd[ind[1]+nl+2]) > 0

	tstep = timestep(nd)
	times = tstep*istep

	if (is.null(time)) {
		n = length(ind)
		if (! is.na(nmax) && n > nmax) {
			inds = c(1,sort(sample.int(n-1,nmax))+1)
		} else {
			inds = seq(n)
		}
	} else {
		inds = which(times %in% time)
		if (length(inds) == 0) stop(sprintf("time %g not found in steps",time))
	}

	spsp = as.numeric(gsub("SPECTRAL NORMS.+ ([-0-9.E+]+|NaN)$","\\1",nd[ind[inds]]))

	if (is.null(lev)) {
		indi = ind[inds]+2
	} else {
		stopifnot(has.levels)
		if (length(lev) == 1 && lev == 0) lev = seq(nl)

		stopifnot(all(! duplicated(lev)) && all(0 < lev & lev <= nl))
		stopifnot(length(lev) == 1 || length(inds) == 1)

		indi = ind[inds]+lev+2
	}

	spre = regmatches(nd[indi],gregexpr("([-+]?[0-9][-+0-9.E]+[0-9]|NaN)",nd[indi]))
	spre = lapply(spre,function(x) gsub("(\\d+)(\\-\\d+)","\\1E\\2",x))
	sphyd = t(sapply(spre,as.numeric))
	stopifnot(is.matrix(sphyd))

	noms = strsplit(nd[ind[1]+1]," {2,}")[[1]][-1]

	spnh = NULL
	if (any(regexpr("\\<LNHDYN *= *T",nd) > 0)) {
		if (has.levels) {
			indi = indi+nl+2
			indn = ind[1]+nl+3
		} else {
			indi = indi+2
			indn = ind[1]+3
		}

		spre = regmatches(nd[indi],gregexpr("([-+]?[0-9][-+0-9.E]+[0-9]|NaN)",nd[indi]))
		spnh = t(sapply(spre,as.numeric))
		stopifnot(is.matrix(spnh))

		nomsnh = strsplit(nd[indn]," {2,}")[[1]][-1]
		nomsnh[nomsnh == "LOG(PRE/PREHYD)"] = "LOG(P/P_hyd)"
		nomsnh[nomsnh == "d4 = VERT DIV + X"] = "d4 (= vdiv+X)"
	}

	noms[noms == "KINETIC ENERGY"] = "TKE"

	ina = apply(sphyd,2,function(x) all(is.na(x)))
	if (any(ina)) {
		for (i in which(ina)) sphyd[,i] = 0
	}

	if (length(inds) > 1) {
		spvar = data.frame(sphyd,spsp)
		names(spvar) = c(noms,"SP")
	} else {
		spvar = data.frame(sphyd)
		names(spvar) = noms
	}

	if (! is.null(spnh)) {
		spvar = cbind(spvar,spnh)
		names(spvar) = c(names(spvar),nomsnh)
	}

	attr(spvar,"step") = istep[inds]
	attr(spvar,"time") = istep[inds]*tstep
	attr(spvar,"tstep") = tstep
	attr(spvar,"level") = lev

	spvar
}

plott = function(xaxis,col=par("col"),reg,...)
{
	plot(xaxt="n",col=col,...)
	axis(1,xaxis)

	if (! is.null(reg)) {
		mtext(sprintf("Reg. line tend.: %.4g",coef(reg)[2]),lty=2,col=col,cex=par("cex"))
		abline(reg,lty=2,col=col)
	}
}

args = commandArgs(trailingOnly=TRUE)
lev = NULL
ind = grep("lev=",args)
if (length(ind) > 0) {
	if (length(grep("lev=\\d+",args)) == 0) stop("level is not a positive integer value")

	lev = strsplit(sub("lev=","",args[ind[1]]),split=":")[[1]]
	lev = as.integer(lev)
	if (length(lev) > 1) {
		warning("several levels in 'lev=...', only 1st one will be used")
		lev = lev[1]
	}
}

time = NULL
ind = grep("time=\\d+",args)
if (length(ind) > 0) time = as.numeric(sub("time=","",args[ind[1]]))

nmax = NA
ind = grep("nmax=\\d+",args)
if (length(ind) > 0) nmax = as.integer(sub("nmax=","",args[ind[1]]))

files = grep("=",args,invert=TRUE,value=TRUE)

if (regexpr("\\.txt",files[1]) > 0) {
	info = read.table(files[1],header=TRUE)
} else {
	info = data.frame(file=files,legend=letters[seq(along=files)],col=seq(along=files))
}

for (i in seq(along=info$file)) {
	cat(". checking file",info$file[i],"\n")
	if (! file.exists(info$file[i])) stop("file not found")

	nd = readLines(info$file[i])

	if (length(grep("SPECTRAL NORMS",nd)) == 0) stop("no spectral norms")

	if (! is.null(lev)) {
		nflevg = getvar("NFLEVG",nd)
		has.levels = getvar("NSPPR",nd) > 0
		if (! has.levels) stop("no levels in norms (NSPPR=0)")
		if (lev > nflevg) stop("level asked for > NFLEVG")
	}
}

sp = lapply(info$file,spnorm,time=time,nmax=nmax,lev=lev)

times = lapply(sp,function(spdf) attr(spdf,"step")*attr(spdf,"tstep"))

noms = unique(unlist(lapply(sp,function(spdf) names(spdf))))
for (i in seq(along=sp)) {
	indc = match(noms,names(sp[[i]]))
	ic = na.omit(indc)
	for (j in attr(ic,"na.action")) sp[[i]][,j] = NA_real_
	sp[[i]][which(! is.na(indc))] = sp[[i]][ic]
	names(sp[[i]]) = noms
}

ip = grep("LOG\\(P/P_hyd\\)|d4",noms,invert=TRUE)
noms[ip] = abbreviate(noms[ip])

titre = paste("Spectral norm of",noms)

hasx11 = capabilities("X11")

if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

sp1 = sp[[1]]
nf = dim(sp1)[2]

if (length(time) == 1) {
	cat("Spectral norms at step",sprintf("%gs (= %.5gh)",time,time/3600),"\n")
	lev = attr(sp1,"level")
	ylim = rev(range(lev))
	nr = 1+(nf-1)%/%4
	nc = min(nf,4)

	if (! hasx11) png("spnormv.png")

	par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	for (i in seq(min(nf,nr*nc))) {
		plot(sp1[,i],lev,type="l",ylim=ylim,lty=1,xlab=noms[i],
			ylab="Level",main=titre[i],col=info$col[1])
		abline(v=0,col="darkgrey",lty=2)
	}

	if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
	invisible(dev.off())

	stop("vertical profile")
}

ttime = unique(sort(unlist(times)))

if (diff(range(ttime)) < 86400) {
	x = pretty(ttime/6)*6
	xlab = "fc time (s)"
} else if (diff(range(ttime)) < 86400*6) {
	ttime = ttime/3600
	times = lapply(times,"/",3600)
	x = pretty(ttime/6)*6
	xlab = "fc time (h)"
} else {
	ttime = ttime/86400
	times = lapply(times,"/",86400)
	x = pretty(ttime)
	xlab = "fc time (days)"
}

time = times[[1]]
xlim = range(ttime)

if (! is.null(lev)) {
	ind3d = which(noms != "SP")
	titre[ind3d] = paste("Spectral norm of",noms[ind3d],"- lev",lev)
}

regs = lapply(sp1,function(y) line(time,y))

if (! hasx11) png("spnorm.png")

par(mfrow=c(3,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
for (i in 1:min(nf,6)) {
	ylim = range(unlist(lapply(sp,"[[",i)),na.rm=TRUE)
	plott(time,sp1[,i],type="l",lty=1,xlim=xlim,ylim=ylim,xlab=xlab,ylab=noms[i],
		main=titre[i],xaxis=x,col=info$col[1],reg=regs[[i]])
	if (length(sp) > 1) {
		if (i == 1) legend("topleft",info$legend,lty=1,col=info$col,cex=1.4*par("cex"))
		for (j in seq(along=sp)[-1]) lines(times[[j]],sp[[j]][,i],col=info$col[j])
	}
}

if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
invisible(dev.off())

if (nf > 6) {
	if (! hasx11) png("spnorm2.png")

	par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	for (i in 5:min(nf,8)) {
		ylim = range(unlist(lapply(sp,"[[",i)),na.rm=TRUE)
		plott(time,sp1[,i],type="l",lty=1,xlim=xlim,ylim=ylim,xlab=xlab,ylab=noms[i],
			main=titre[i],xaxis=x,col=info$col[1],reg=regs[[i]])
		if (length(sp) > 1) {
			legend("topleft",info$legend,lty=1,col=info$col,cex=1.3*par("cex"))
			for (j in seq(along=sp)[-1]) lines(times[[j]],sp[[j]][[i]],col=info$col[j])
		}
	}

	if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
	invisible(dev.off())
}

if (length(files) == 2) {
	spref = spnorm(files[2],attr(sp1,"time"),nmax=nmax,lev=lev)

	spdiff = sp1-spref
	sp0 = sp1$SP[1]
	if (spdiff$SP[1] != 0)
		cat("--> different initial value for SP:",sp0,spref$SP[1],"- diff:",spdiff$SP[1],"\n")

	if (! hasx11) png("spdiff.png")

	titre = paste("Difference in spectral norm of",noms)
	par(mfrow=c(3,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	for (i in 1:min(nf,6)) {
		plott(time,spdiff[,i],type="l",lty=1,xlab=xlab,ylab=noms[i],main=titre[i],xaxis=x,
			col=1,reg=NULL)
		abline(h=0,col="darkgrey")
	}

	if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
	invisible(dev.off())

	if (nf > 6) {
		if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
		if (! hasx11) png("spdiff2.png")

		par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		for (i in 7:min(nf,10)) {
			plott(time,spdiff[,i],type="l",lty=1,xlab=xlab,ylab=noms[i],main=titre[i],xaxis=x,
				col=1,reg=NULL)
			abline(h=0,col="darkgrey")
		}

		if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
		invisible(dev.off())
	}
}

