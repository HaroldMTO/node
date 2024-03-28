library(mfnode)

gpfre = sprintf("%s|DIV",gpfre)

getarg = function(x,args)
{
	ind = grep(sprintf("\\<%s=",x),args)
	if (length(ind) == 0) return(NULL)

	strsplit(sub(sprintf("\\<%s=",x),"",args[ind]),split=":")[[1]]
}

args = commandArgs(trailingOnly=TRUE)

hasx11 = is.null(getarg("png",args)) && capabilities("X11")
ask = hasx11 && interactive()
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

xaxis = data.frame(unit=c(1,60,3600,86400),label=sprintf("fc time (%s)",
	c("s","mn","h","days")),mindiff=c(0,240,14400,6*86400),freq=c(6,1,6,1))

pngd = getarg("png",args)
if (is.null(pngd)) pngd = "."
lev = as.integer(getarg("lev",args))
if (is.null(lev)) lev = 0L
hmin = as.numeric(getarg("hmin",args))
hmax = as.numeric(getarg("hmax",args))
ptype = getarg("type",args)
if (is.null(ptype)) ptype = "spec"
detail = regexpr("gp.+detail",ptype) > 0
spre = getarg("spre",args)
spref = getarg("spref",args)
gpre = getarg("gpre",args)
gpref = getarg("gpref",args)

fnode = grep("=",args,invert=TRUE,value=TRUE)
stopifnot(length(fnode) == 1)

cat("Read file",fnode,"\n")
nd = readLines(fnode)
nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
nflevg = getvar("NFLEVG",nd)
has.levels = getvar(".*NSPPR",nd) > 0

grouplev = NULL
if (length(lev) > 1) {
	grouplev = lev
	lev = seq(nflevg)
} else if (lev == -1) {
	if (! has.levels) {
		cat("--> no level norms\n")
		q("no")
	}

	lev = seq(nflevg)
}

if (! identical(lev,0L)) stopifnot(has.levels)

tstep = getvar("TSTEP",nd)

if (interactive()) browser()

if (ptype == "spec") {
	cat("Parse spectral norms\n")
	if (is.null(spref)) {
		sp1 = spnorm(nd,lev)
	} else {
		sp1 = spnorm(nd,lev,spref)
	}

	if (is.null(sp1)) {
		cat("--> no SP norms\n")
		quit("no")
	}

	step = dimnames(sp1)[[1]]
	ix = grep("^X",step)
	if (length(ix) == length(step)) {
		cat("--> SP norms for STEPX only, quit\n")
		quit("no")
	} else if (length(ix) > 0) {
		cat("--> SP norms present for STEPX, remove",length(ix),"steps\n")
		sp1 = sp1[-ix,,,drop=FALSE]
		step = dimnames(gp1)[[1]]
	}

	istep = as.numeric(gsub("C(\\d+)","\\1.5",step))
	cat(". steps:",head(step[-length(step)]),"...",step[length(step)],"\n")

	spnoms = dimnames(sp1)[[3]]
	spl = list(sp1)
	leg = c("t0",sub("spnorm +","",spre))
	for (i in seq(along=spre)) {
		cat("Parse spectral norms, pattern",spre[i],"\n")
		spi = spnorm(nd,lev,spre[i])
		if (is.null(spi)) {
			cat("--> no norms for pattern:",spre[i],"\n")
			next
		}

		indv = match(spnoms,dimnames(spi)[[3]])
		stopifnot(any(! is.na(indv)))

		stepi = dimnames(spi)[[1]]
		inds = match(step,stepi)
		stopifnot(any(! is.na(inds)))
		spl[[i+1]] = spi[inds,,indv,drop=FALSE]
	}

	leg = leg[which(! sapply(spl,is.null))]
	spl = spl[! sapply(spl,is.null)]

	stepc = step[regexpr("C\\d+",step) > 0]
	if (length(stepc) > 0) {
		cat("--> PC scheme, split P/C steps\n")
		indc = match(stepc,step)
		ip = which(step %in% step[-indc])
		indc = match(paste("C",step[ip],sep=""),step)
		splc = lapply(spl,function(spi) spi[indc,,,drop=FALSE])
		splp = lapply(spl,function(spi) spi[ip,,,drop=FALSE])
		spl = c(splp,splc)
		leg = c(leg,paste(leg,"Cor"))
		istep = istep[ip]
		sp1 = spl[[1]]
	}

	nf = length(spnoms)
	nt = dim(sp1)[1]
	times = tstep*istep

	if (length(lev) > 1) {
		cat("Produce vertical profiles at start/mid/end of time\n")
		indt = c(1,(nt+1)%/%2,nt)
		tt = sprintf("SP norm of %s",spnoms)
		ts = sprintf("t%s",paste(istep[indt],collapse="/"))
		nc = min(nf,2)
		nr = min(length(spl),2)
		nj = nc

		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/spnormv%d.png",pngd,i))

			par(mfcol=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			for (j in seq(min(nf-nj*i,nj))+nj*i) {
				y = t(sp1[indt,,j])
				plotvmean(y,lev,xlab=spnoms[j],main=c(tt[j],ts))

				if (length(spl) == 1) next

				y = sapply(spl,function(x) x[indt[2],,j],simplify="array")
				ts2 = sprintf("t%s",istep[indt[2]])
				plotvmean(y,lev,xlab=spnoms[j],main=c(tt[j],ts2),lty=c(1,2,2),
					col=seq(along=spl),legend=leg)
			}

			if (! hasx11) dev.off()
		}

		if (! is.null(grouplev)) {
			cat("--> groups of levels:",grouplev,"\n")
			indg = 0
			leg = character()
			for (i in seq(along=grouplev)) {
				indi = which(lev <= grouplev[i])
				indg = indi[indi > max(indg)]
				cat(".. level indices:",range(indg),"\n")
				leg[i] = sprintf("Lev %d-%d",min(indg),max(indg))
				spi = sp1[,1,,drop=FALSE]
				spi[] = apply(sp1[,indg,,drop=FALSE],c(1,3),mean)
				spl[[i]] = spi
			}

			sp1 = spl[[1]]
			lev = 0
		}
	}

	if (length(lev) == 1) {
		cat("Produce time-series per variable\n")
		if (dim(sp1)[1] == 1) {
			cat("--> 1 time-step only, quit\n")
			quit("no")
		}

		if (length(spl) > 1) {
			con = file(sprintf("%s/%s.txt",pngd,ptype),"wt")
			for (i in seq(along=spl)[-1]) {
				cat("Norm diff. tendency for patterns:\n",sprintf("'%s'",leg[(i-1):i]),"\n",
					file=con)
				iv = which(! is.na(dimnames(spl[[i]])[[3]]))
				spdiff = spl[[i]][,1,iv,drop=FALSE]-spl[[i-1]][,1,iv,drop=FALSE]
				a = apply(spdiff,3,function(x) coef(line(x))[2])
				s = apply(spdiff,3,function(x) sd(residuals(line(x))))
				cat("\n",file=con)
				evo = rbind(tend=a,sd=s)
				evo = format(evo,digits=3,width=9)
				cat("field",format(dimnames(evo)[[2]],width=9),"\n",sep="\t",file=con)
				write.table(evo,con,quote=FALSE,col.names=FALSE,sep="\t")
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
		spl = lapply(spl,function(x) x[it,,,drop=FALSE])
		sp1 = spl[[1]]

		cat("Spectral norms",spref,"\n")
		titre = paste("SP norm of",spnoms)
		if (lev > 0) titre = paste(titre,"- lev",lev)

		# 2 rows, 2 columns
		nc = 2
		nr = min(1+(nf-1)%/%2,2)
		nj = nr*nc
		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/spnorm%d.png",pngd,i))

			par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				cat(". SP field",spnoms[j],"\n")
				y = sapply(spl,function(x) x[,1,j])
				ymax = max(abs(y),na.rm=TRUE)
				scal = 10^-round(log10(ymax/1.5))
				if (.001 <= scal && scal < 1) scal = 1
				plotmean(ttime,y,titre[j],leg,tunit,xlim=xlim,xlab=xlab,ylab=spnoms[j],
					xaxp=xaxp,scale=scal)
			}
		}
	}
}

if (regexpr("^gp",ptype) > 0) {
	cat("Parse GP norms of type",ptype,"\n")
	if (regexpr("gpgfl",ptype) > 0) {
		gp1 = gpnorm(nd,lev,gpout=gpfre)
	} else {
		gp1 = gpnorm(nd,lev,gpref,gpfre)
	}

	if (is.null(gp1)) {
		cat("--> no GP norms, quit\n")
		quit("no")
	}

	i0 = apply(gp1,1,function(x) all(x==0))
	if (any(i0) && ! all(i0)) {
		cat("--> GP norms all 0 for some steps, removed\n")
		gp1 = gp1[-which(i0),,,,drop=FALSE]
	}

	step = dimnames(gp1)[[1]]
	cat(". steps:",head(step[-length(step)]),"...",step[length(step)],"\n")

	ix = grep("^X",step)
	if (length(ix) == length(step)) {
		cat("--> GP norms for STEPX only, quit\n")
		quit("no")
	} else if (length(ix) > 0) {
		cat("--> GP norms present for STEPX, remove",length(ix),"steps\n")
		gp1 = gp1[-ix,,,,drop=FALSE]
		step = dimnames(gp1)[[1]]
	}

	istep = as.numeric(gsub("C(\\d+)","\\1.5",step))
	gpl = list(gp1)

	gpnoms = dimnames(gp1)[[4]]
	leg = c("t0",sub("gpnorm +","",gpre))
	for (i in seq(along=gpre)) {
		if (regexpr("gpgfl",ptype) > 0) {
			gpi = gpnorm(nd,lev,gpre[i],gpout=gpfre)
		} else {
			gpi = gpnorm(nd,lev,gpre[i],gpfre)
		}

		if (is.null(gpi)) {
			cat("--> no norms for pattern",gpre[i],"\n")
			next
		}

		i0 = apply(gpi,1,function(x) all(x==0))
		if (all(i0)) {
			cat("--> GFL all 0 for all steps, continue\n")
			next
		} else if (any(i0)) {
			cat("--> GFL all 0 for some steps, removed\n")
			gpi = gpi[-which(i0),,,,drop=FALSE]
		}

		indv = match(gpnoms,dimnames(gpi)[[4]])
		stopifnot(any(! is.na(indv)))

		stepi = dimnames(gpi)[[1]]
		ix = grep("^X",stepi)
		if (length(ix) > 0) {
			gpi = gpi[-ix,,,,drop=FALSE]
			stepi = dimnames(gpi)[[1]]
		}

		inds = match(step,stepi)
		stopifnot(any(! is.na(inds)))
		gpl[[i+1]] = gpi[inds,,,indv,drop=FALSE]
	}

	leg = leg[which(! sapply(gpl,is.null))]
	gpl = gpl[! sapply(gpl,is.null)]

	nf = length(gpnoms)
	nt = dim(gp1)[1]

	times = tstep*istep

	if (length(lev) > 1) {
		cat("Produce vertical profiles\n")
		tt = sprintf("GP norm of %s",gpnoms)
		indt = c(1,(nt+1)%/%2,nt)
		ts = sprintf("t%s",paste(istep[indt],collapse="/"))
		nc = 3
		nr = min(length(gpl),2)
		nj = 3-nr

		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/%snormv%d.png",pngd,ptype,i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				cat(". GP field",gpnoms[j],"\n")
				y = aperm(gp1[indt,,,j],c(2,1,3))
				plotvmnx(y,lev,xlab=gpnoms[j],main=c(tt[j],ts))

				if (length(gpl) == 1) next

				y = sapply(gpl,function(x) x[indt[2],,,j],simplify="array")
				y = aperm(y,c(1,3,2))
				il = which(apply(y,2,function(x) any(! is.na(x))))
				ts2 = sprintf("t%d",indt[2])
				plotvmnx(y,lev,xlab=gpnoms[j],main=c(tt[j],ts2),legend=leg[il],lty=2,col=il)
			}
		}

		if (! is.null(grouplev)) {
			cat("--> groups of levels:",grouplev,"\n")
			indg = 0
			leg = character()
			for (i in seq(along=grouplev)) {
				indi = which(lev <= grouplev[i])
				indg = indi[indi > max(indg)]
				cat(".. level indices:",range(indg),"\n")
				leg[i] = sprintf("L %d-%d",min(indg),max(indg))
				gpi = gp1[,1,,,drop=FALSE]
				gpi[,,1,] = apply(gp1[,indg,1,,drop=FALSE],c(1,3,4),mean)
				gpi[,,2,] = apply(gp1[,indg,2,,drop=FALSE],c(1,3,4),min)
				gpi[,,3,] = apply(gp1[,indg,3,,drop=FALSE],c(1,3,4),max)
				gpl[[i]] = gpi
			}

			gp1 = gpl[[1]]
			lev = 0
		}
	}

	if (length(lev) == 1) {
		cat("Produce time-series, level",lev,"\n")
		if (nt == 1) stop("1 time-step only (stop)\n")

		if (length(gpl) > 1) {
			con = file(sprintf("%s/%s.txt",pngd,ptype),"wt")
			for (i in seq(along=gpl)[-1]) {
				cat("Norm diff. tendency for patterns:\n",sprintf("'%s'",leg[(i-1):i]),"\n",
					file=con)
				iv = which(! is.na(dimnames(gpl[[i]])[[4]]))
				gpdiff = gpl[[i]][,1,1,iv,drop=FALSE]-gpl[[i-1]][,1,1,iv,drop=FALSE]
				a = apply(gpdiff,4,function(x) coef(line(x))[2])
				s = apply(gpdiff,4,function(x) sd(residuals(line(x))))
				cat("\n",file=con)
				evo = rbind(tend=a,sd=s)
				evo = format(evo,digits=3,width=9)
				cat("field",format(dimnames(evo)[[2]],width=9),"\n",sep="\t",file=con)
				write.table(evo,con,quote=FALSE,col.names=FALSE,sep="\t")
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
		gpl = lapply(gpl,function(x) x[it,,,,drop=FALSE])
		gp1 = gpl[[1]]

		cat("GP norms",gpref,"\n")
		gpnoms = dimnames(gp1)[[4]]
		titre = paste("GP norm of",gpnoms)
		if (lev > 0) titre = paste(titre,"- lev",lev)

		nf = length(gpnoms)
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
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/%snorm%d.png",pngd,ptype,i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				y = sapply(gpl,function(x) x[,1,,j],simplify="array")
				il = apply(y,3,function(x) any(! is.na(x)))
				ymax = max(abs(y[,1,]),na.rm=TRUE)
				scal = 10^-round(log10(ymax/1.5))
				if (.001 <= scal && scal < 1 || is.infinite(scal)) scal = 1
				cat(". GP field",gpnoms[j],"- scaling:",scal,ymax,"\n")
				plotmean(ttime,y[,1,],main=titre[j],leg[il],tunit,xlim=xlim,xlab=xlab,
					ylab=gpnoms[j],xaxp=xaxp,scale=scal)

				plotmnx(ttime,y,titre[j],xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)

				if (! detail) next

				plotmnx(ttime,y,titre[j],imnx=2,xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)
				plotmnx(ttime,y,titre[j],imnx=3,xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)
			}
		}
	}
}

if (ptype == "surf") {
	surf = gpnorm2D(nd)
	if (length(surf) == 0) {
		cat("--> no GP norms for surface fields\n")
	} else {
		msurf = surf[[1]]
		for (i in seq(along=surf)[-1]) msurf = rbind(msurf,surf[[i]])
		cat("Surface fields:\n")
		print(msurf)
	}
}

