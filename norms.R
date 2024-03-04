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

line2num = function(nd)
{
	lre = regmatches(nd,gregexpr(sprintf("(%s|\\<NaN\\>)",Gnum),nd))
	lre = lapply(lre,function(x) gsub("(\\d+)([-+]\\d+)","\\1E\\2",x))
	sapply(lre,as.numeric)
}

gpnorm2D = function(nd)
{
	ind = grep("^ *NUMFLDS=",nd)
	indo = grep("^ *GPNORM OUTPUT",nd)

	nfg = as.integer(sub(" *NUMFLDS= *(\\d+) .+","\\1",nd[ind]))
	ind = ind[nfg > 0]
	nfg = nfg[nfg > 0]
	surf = list()
	group = character()

	for (i in seq(along=ind)) {
		group[i] = sub("^.+ (\\w+) +- +.+","\\1",nd[ind[i]-1])
		gnames = sub("^ *\\w+( +\\d+)+ +(\\w+(\\.\\w+)?).+","\\2",nd[ind[i]+seq(nfg[i])])

		ii = grep(sprintf("\\<%s\\>",group[i]),nd[indo-1])
		if (length(ii) == 0) {
			group[i] = sub("^ *(\\w+) +.+","\\1",nd[ind[i]+1])
			ii = grep(sprintf("\\<%s\\>",group[i]),nd[indo-1])
			if (length(ii) == 0) {
				cat("--> no GP norms for group",group[i],i,"\n")
				next
			}
		}

		ii = ii[1]
		if (regexpr(", +FIELD +\\d+",nd[indo[ii]-1]) > 0) {
			# nfg lines AVE, every 4 lines (group, GPNORM, AVE, 1)
			indi = indo[ii]+(seq(nfg[i])-1)*4+1
		} else if (regexpr(" \\d+ +FIELDS\\>",nd[indo[ii]-1]) > 0) {
			# nfg lines AVE, every 3 lines (GPNORM, AVE, 1)
			indi = indo[ii]+(seq(nfg[i])-1)*3+1
		} else {
			# nfg lines after GPNORM and AVE
			indi = indo[ii]+seq(nfg[i])+1
		}

		gpn = t(line2num(nd[indi]))
		dimnames(gpn) = list(gnames,c("ave","min","max"))
		surf[[i]] = gpn
	}

	surf[sapply(surf,length) > 0]
}

gpnorm = function(nd,lev,ind,noms)
{
	if (missing(ind)) {
		ind = grep("GPNORM +\\w+.* +AVERAGE",nd)
		indo = grep("GPNORM OUTPUT",nd[ind],invert=TRUE)
		ind = ind[indo]
	}

	if (length(ind) == 0) stop("no GP norms")

	if (missing(noms)) {
		noms = unique(sub(" *GPNORM +(\\w+.+?) +AVERAGE.+","\\1",nd[ind]))
	} else {
		i = grep(sprintf(" *GPNORM +(%s)\\>",paste(noms,collapse="|")),nd[ind])
		ind = ind[i]
   }

	indi = rep(ind,each=length(lev))+lev+1

	gpn = line2num(nd[indi])

	nt = length(gpn)/(3*length(lev)*length(noms))
	if (nt > as.integer(nt)) {
		nn = sapply(noms,function(x) length(grep(sprintf("\\<%s\\>",x),nd[ind])))
		noms = noms[nn == max(nn)]
		i = grep(sprintf(" *GPNORM +(%s)\\>",paste(noms,collapse="|")),nd[ind])
		ind = ind[i]
		indi = rep(ind,each=length(lev))+lev+1
		gpn = line2num(nd[indi])
		nt = length(gpn)/(3*length(lev)*length(noms))
	}

	stopifnot(nt == as.integer(nt))

	noms = unique(sub(" *GPNORM +(\\w+.+?) +AVERAGE.+","\\1",nd[ind]))
	noms[noms == "SURFACE PRESSURE"] = "SURF P"
	noms[noms == "TEMPRATURE"] = "TEMP"
	noms[noms == "U VELOCITY"] = "U VELOC."
	noms[noms == "V VELOCITY"] = "V VELOC."

	dim(gpn) = c(3,length(lev),length(noms),nt)
	gpl = aperm(gpn,c(4,2,1,3))

	dimnames(gpl) = list(NULL,lev,c("ave","min","max"),noms)

	gpl
}

spnorm = function(nd,lev,ind)
{
	if (missing(ind)) {
		ind = grep("SPECTRAL NORMS",nd)
		inds = grep("NORMS AT NSTEP CNT4",nd[ind-1])
		ind = ind[inds]
	}

	spsp = as.numeric(gsub("SPECTRAL NORMS.+ ([-0-9.E+]+|NaN)$","\\1",nd[ind]))

	noms = strsplit(nd[ind[1]+1]," {2,}")[[1]][-1]
	noms[noms == "KINETIC ENERGY"] = "TKE"

	indi = rep(ind,each=length(lev))+lev+2

	spn = line2num(nd[indi])

	nt = length(spn)/(length(lev)*length(noms))
	stopifnot(nt == as.integer(nt))

	dim(spn) = c(length(noms),length(lev),nt)
	spn = aperm(spn,c(3,2,1))

	if (length(lev) == 1) {
		spn = c(spn,spsp)
		noms = c(noms,"SP")
		dim(spn) = c(nt,length(lev),length(noms))
	}

	if (any(regexpr("\\<LNHDYN *= *T",nd) > 0)) {
		if (has.levels) {
			indi = indi+nflevg+2
			indn = ind[1]+nflevg+3
		} else {
			indi = indi+2
			indn = ind[1]+3
		}

		spnh = line2num(nd[indi])

		nomsnh = strsplit(nd[indn]," {2,}")[[1]][-1]
		nomsnh[nomsnh == "LOG(PRE/PREHYD)"] = "LOG(P/P_hyd)"
		nomsnh[nomsnh == "d4 = VERT DIV + X"] = "d4 (= vdiv+X)"

		dim(spnh) = c(length(nomsnh),length(lev),nt)
		spnh = aperm(spnh,c(3,2,1))
		spn = c(spn,spnh)
		noms = c(noms,nomsnh)
		dim(spn) = c(nt,length(lev),length(noms))
	}

	ip = grep("LOG\\(P/P_hyd\\)|d4",noms,invert=TRUE)
	noms[ip] = abbreviate(noms[ip])

	istep = sub(" *NORMS AT NSTEP CNT4( \\(PREDICTOR\\))? +(\\d+)","\\2",nd[ind-1])
	dimnames(spn) = list(istep,lev,noms)

	spn
}

countfield = function(ind,ind2,nl2)
{
	which(diff(ind[seq(ind2[1],ind2[2])]) > nl2)[1]
}

indexpand = function(ind,nf,nl2)
{
	rep(ind,each=nf)+(seq(nf)-1)*nl2
}

args = commandArgs(trailingOnly=TRUE)

hasx11 = is.null(getarg("png",args)) && capabilities("X11")
ask = hasx11 && interactive()
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

xaxis = data.frame(unit=c(1,3600,86400),label=sprintf("fc time (%s)",c("s","h","days")),
	mindiff=c(0,86400,6*86400),freq=c(6,6,1))

pngd = getarg("png",args)
lev = as.integer(getarg("lev",args))
if (is.null(lev)) lev = 0
hmin = as.numeric(getarg("hmin",args))
hmax = as.numeric(getarg("hmax",args))
ptype = getarg("type",args)
if (is.null(ptype)) ptype = "spec"
detail = regexpr("gp.+detail",ptype) > 0
spre = getarg("spre",args)
spref = getarg("spref",args)
if (is.null(spref)) spref = "spnorm t0"
gpre = getarg("gpre",args)
gpref = getarg("gpref",args)
if (is.null(gpref)) {
	gpref = "gpnorm gmvt0"
	if (any(ptype == "gpgfl")) gpref = "gpnorm gflt0"
}
fnode = grep("=",args,invert=TRUE,value=TRUE)
stopifnot(length(fnode) == 1)

cat("Read file",fnode,"\n")
nd = readLines(fnode)
nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
nflevg = getvar("NFLEVG",nd)
has.levels = getvar("NSPPR",nd) > 0

if (length(lev) > 1 || lev != 0) stopifnot(has.levels)
if (length(lev) > 1) {
	grouplev = lev
	lev = seq(nflevg)
} else if (lev == -1) {
	lev = seq(nflevg)
}

nstop = getvar("NSTOP",nd)
tstep = getvar("TSTEP",nd)
icnt4 = grep("^ *START CNT4",nd)[1]

if (interactive()) browser()

if (ptype == "spec") {
	cat("Parse spectral norms\n")
	ind = grep("SPECTRAL NORMS",nd)
	ind = ind[ind > icnt4]
	# ind2: for corrector (one line may be interleaved)
	ind1 = grep(sub("spnorm t0","NORMS AT NSTEP CNT4",spref),nd[ind-1])
	ind2 = grep(sub("spnorm t0","NORMS AT NSTEP CNT4",spref),nd[ind-2])
	if (length(ind2) > 0) ind1 = sort(c(ind1,ind2))
	sp1 = spnorm(nd,lev,ind[ind1[-1]])
	spnoms = dimnames(sp1)[[3]]
	spl = list(sp1)
	leg = c("t0",sub("spnorm +","",spre))
	for (i in seq(along=spre)) {
		cat("Parse spectral norms, pattern",spre[i],"\n")
		ind2 = grep(spre[i],nd[ind-1],ignore.case=TRUE)
		if (length(ind2) == 0) {
			cat("--> no norms for pattern:",spre[i],"\n")
			next
		}

		indv = match(spnoms,dimnames(sp1)[[3]])
		stopifnot(any(! is.na(indv)))
		spi = spnorm(nd,lev,ind[ind2])
		spl[[i+1]] = spi[,,indv,drop=FALSE]
	}

	leg = leg[seq(along=spl)]

	nfrsdi = getvar(".+ NFRSDI",nd)
	nsdits = getvar("NSDITS",nd)
	if (nsdits[1] != 0) message("Warning: possible crash for SP plot (NSDITS != 0)")
	istep = seq(1,nstop,by=nfrsdi)

	nf = length(spnoms)
	nt = dim(sp1)[1]
	if (length(istep) > nt) length(istep) = nt
	times = tstep*istep

	if (length(lev) > 1) {
		cat("Produce vertical profiles at start/mid/end of time\n")
		levels = seq(nflevg)
		indt = c(1,(nt+1)%/%2,nt)
		tt = sprintf("Spectral norm of %s",spnoms)
		ylim = rev(range(levels))
		nc = min(nf,3)
		nr = min(1+(nf-1)%/%nc,3)
		nj = nc*nr

		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/spnormv%d.png",pngd,i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				cat(". SP field",spnoms[j],"\n")
				matplot(t(sp1[indt,,j]),levels,type="l",lty=1:3,col=1,ylim=ylim,
					xlab=spnoms[j],ylab="Level",main=tt[j])
				abline(v=0,col="darkgrey",lty=2)
			}

			if (! hasx11) dev.off()
		}

		if (length(grouplev) > 1) {
			cat("--> groups of levels:",grouplev,"\n")
			indg = 0
			leg = character()
			for (i in seq(along=grouplev)) {
				indi = which(levels <= grouplev[i])
				indg = indi[indi > max(indg)]
				cat(".. level indices:",range(indg),"\n")
				leg[i] = sprintf("L %d-%d",min(indg),max(indg))
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
		if (dim(sp1)[1] == 1) stop("1 time-step only (stop)\n")

		if (any(sapply(spl,function(x) dim(x)[1]) > nt)) {
			cat("--> limiting norms to",nt,"occurrences\n")
			spl = lapply(spl,function(x) x[1:nt,,,drop=FALSE])
		}

		if (length(spl) > 1) {
			con = file(sprintf("%s/%s.txt",pngd,ptype),"wt")
			cat("Norm diff. tendency for patterns:\n",
				paste(sprintf("'%s'",spre),collapse=" "),"\n",file=con)
			for (i in seq(along=spl)[-1]) {
				spdiff = matrix(spl[[i]][,1,]-spl[[i-1]][,1,],ncol=nf)
				iv = which(! is.na(dimnames(spl[[i]])[[3]]))
				a = apply(spdiff[,iv,drop=FALSE],2,function(x) coef(line(x))[2])
				s = apply(spdiff[,iv,drop=FALSE],2,function(x) sd(residuals(line(x))))
				cat("\n",file=con)
				write.table(rbind(tend=a,sd=s),con,quote=FALSE,col.names=FALSE)
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
		spl = lapply(spl,function(x) x[it,,,drop=FALSE])
		sp1 = spl[[1]]

		cat("Spectral norms",spref,"\n")
		titre = paste("Spectral norm of",spnoms)
		if (lev > 0) titre = paste(titre,"- lev",lev)

		# 4 rows max, 2 columns
		nc = 2
		nr = min(1+(nf-1)%/%2,4)
		nj = nr*nc
		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/spnorm%d.png",pngd,i))

			par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				cat(". SP field",spnoms[j],"\n")
				y = sapply(spl,function(x) x[,1,j])
				ymax = max(abs(y))
				scal = 10^-round(log10(ymax/1.5))
				if (.001 <= scal && scal < 1) scal = 1
				col = lty = seq(along=spl)
				lty[-1] = 2
				matplot(ttime,y*scal,type="l",lty=lty,col=col,xlim=xlim,xlab=xlab,
					ylab=spnoms[j],main=titre[j],xaxt="n")
				legend("topleft",leg,col=col,lty=lty,bg="transparent")

				axis(1,x)
				reg = line(ttime,sp1[,1,j]*scal)
				abline(reg,lty=2)
				tend = coef(reg)[2]*86400/tunit*scal
				tt = sprintf("Trend: %+.2e [unit]/day",tend)
				if (scal != 1) tt = sprintf("%s, scaling: *%.0e",tt,scal)
				if (tunit == 86400) {
					tt = sprintf("%s, variation: %.2g",tt,tend*diff(range(ttime)))
				}

				mtext(tt,lty=2,cex=par("cex"))
			}

			if (! hasx11) dev.off()
		}
	}
}

if (regexpr("^gp",ptype) > 0) {
	cat("Parse GP norms of type",ptype,"\n")
	# don't take step 0
	nfrgdi = getvar(".+ NFRGDI",nd)
	ngdits = getvar("NGDITS",nd)
	if (ngdits[1] != 0) message("Warning: possible crash for GP plot (NGDITS != 0)")
	istep = seq(1,nstop,by=nfrgdi)

	gpfre1 = "[UVW] VELOCITY|(SURFACE )?PRESSURE|TEMPERATURE|GRAD[LM]_\\w+|GEOPOTENTIAL"
	gpfre2 = "MOIST AIR SPECIF|ISOBARE CAPACITY|SURFACE DIV|d\\(DIV\\)\\*dP"
	gpfre3 = "(ATND|ADIAB|CTY|(SI)?SL)_\\w+"
	gpfre = paste(gpfre1,gpfre2,gpfre3,sep="|")

	nl2 = 2+has.levels*nflevg

	if (regexpr("gpgfl",ptype) > 0) {
		ind = grep("GPNORM +\\w+.* +AVERAGE",nd)
		ind = ind[ind > icnt4]
		indo = grep(sprintf("GPNORM +(%s|OUTPUT) +AVERAGE",gpfre),nd[ind],invert=TRUE)
		ind = ind[indo]
		ind1 = grep(gpref,nd[ind-1])
		if (length(ind1) == 0) {
			ind1 = grep("NORMS AT NSTEP CNT4",nd[ind-2-nl2])
			if (length(ind1) == 0) ind1 = grep("NORMS AT NSTEP CNT4",nd[ind-3-nl2])
			if (length(ind1) == 0) ind1 = grep("NORMS AT NSTEP CNT4",nd[ind-4-nl2])
		}

		nf = countfield(ind,ind1,nl2)
		indi = indexpand(ind[ind1],nf,nl2)
		gp1 = gpnorm(nd,lev,indi)
		if (dim(gp1)[1] > length(istep)) gp1 = gp1[-1,,,,drop=FALSE]

		gpl = list(gp1)
	} else {
		ind = grep(sprintf("GPNORM +(%s) +AVERAGE",gpfre),nd)
		ind = ind[ind > icnt4]
		ind1 = grep(gpref,nd[ind-1],ignore.case=TRUE)
		nf = countfield(ind,ind1,nl2)
		indi = indexpand(ind[ind1],nf,nl2)
		gp1 = gpnorm(nd,lev,indi)
		if (dim(gp1)[1] > length(istep)) gp1 = gp1[-1,,,,drop=FALSE]
		gpl = list(gp1)
	}

	gpnoms = dimnames(gp1)[[4]]
	leg = c("t0",sub("gpnorm +","",gpre))
	for (i in seq(along=gpre)) {
		ind2 = grep(gpre[i],nd[ind-1],ignore.case=TRUE)
		if (length(ind2) == 0) {
			cat("--> no norms for pattern",gpre[i],"\n")
			next
		}

		nf = countfield(ind,ind2,nl2)
		indi = indexpand(ind[ind2],nf,nl2)
		gpi = gpnorm(nd,lev,indi)
		if (dim(gpi)[1] > length(istep)) gpi = gpi[-1,,,,drop=FALSE]
		indv = match(gpnoms,dimnames(gpi)[[4]])
		stopifnot(any(! is.na(indv)))
		gpl[[i+1]] = gpi[,,,indv,drop=FALSE]
	}

	leg = leg[seq(along=gpl)]

	nf = length(gpnoms)
	nt = dim(gp1)[1]
	if (length(istep) > nt) length(istep) = nt

	times = tstep*istep

	if (length(lev) > 1) {
		cat("Produce vertical profiles\n")
		levels = seq(nflevg)
		indt = c(1,(nt+1)%/%2,nt)
		nc = 3
		nr = 2
		nj = nr
		ylim = rev(range(levels))

		for (i in seq((nf-1)%/%nj+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%s/%snormv%d.png",pngd,ptype,i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nj*i,nj)+nj*i) {
				cat(". GP field",gpnoms[j],"\n")
				tt = sprintf("GP norm of %s",gpnoms[j])
				tt2 = c("ave","min","max")
				for (k in c(2,1,3)) {
					matplot(t(gp1[indt,,k,j]),levels,type="l",lty=1:3,col=1,ylim=ylim,
						xlab=gpnoms[j],ylab="Level",main=c(tt,tt2[k]))
					abline(v=0,col="darkgrey",lty=2)
				}
			}

			if (! hasx11) dev.off()
		}

		if (length(grouplev) > 1) {
			cat("--> groups of levels:",grouplev,"\n")
			indg = 0
			leg = character()
			for (i in seq(along=grouplev)) {
				indi = which(levels <= grouplev[i])
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
			cat("Norm diff. tendency for patterns:\n",
				paste(sprintf("'%s'",gpre),collapse=" "),"\n",file=con)
			for (i in seq(along=gpl)[-1]) {
				gpdiff = matrix(gpl[[i]][,1,1,]-gpl[[i-1]][,1,1,],ncol=dim(gpl[[i]])[4])
				iv = which(! is.na(dimnames(gpl[[i]])[[4]]))
				a = apply(gpdiff[,iv,drop=FALSE],2,function(x) coef(line(x))[2])
				s = apply(gpdiff[,iv,drop=FALSE],2,function(x) sd(residuals(line(x))))
				cat("\n",file=con)
				write.table(rbind(tend=a,sd=s),con,quote=FALSE,col.names=FALSE)
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
				y = sapply(gpl,function(x) x[,1,1,j])
				ymax = max(abs(y))
				scal = 10^-round(log10(ymax/1.5))
				if (.001 <= scal && scal < 1 || is.infinite(scal)) scal = 1
				cat(". GP field",gpnoms[j],"- scaling:",scal,ymax,"\n")
				col = lty = seq(along=gpl)
				lty[-1] = 2
				matplot(ttime,y*scal,type="l",lty=lty,col=col,xlim=xlim,xlab=xlab,
					ylab=gpnoms[j],main=paste(titre[j],"(ave)"),xaxt="n")

				reg = line(ttime,gp1[,1,1,j]*scal)
				abline(reg,lty=2)
				tend = coef(reg)[2]*86400/tunit
				tt = sprintf("Trend: %+.2e [unit]/day",tend)
				if (scal != 1) tt = sprintf("%s, scaling: *%.0e",tt,scal)
				if (tunit == 86400) {
					tt = sprintf("%s, variation: %.2g",tt,tend*diff(range(ttime)))
				}

				legend("topleft",leg,col=col,lty=lty,bg="transparent")
				mtext(tt,cex=1.1*par("cex"))
				axis(1,x)
				if (all(gp1[,1,,j] == 0)) text(sum(range(ttime))/2,.5,"all values = 0")

				y = sapply(gpl,function(x) x[,1,,j])
				matplot(ttime,gp1[,1,,j],type="l",lty=ltymnx,xlab=xlab,ylab=gpnoms[j],
					xlim=xlim,ylim=range(y),main=paste(titre[j],"(ave/min/max)"),col=1,
					xaxt="n")
				axis(1,x)
				for (k in seq(along=gpl)[-1]) {
					if (all(is.na(gpl[[k]][,1,,j]))) next
					matlines(ttime,gpl[[k]][,1,,j],lty=ltymnx,col=k)
				}

				if (! detail) next

				y = sapply(gpl,function(x) x[,1,2,j])
				matplot(ttime,gp1[,1,2,j],type="l",lty=ltymnx[2],xlab=xlab,
					ylab=gpnoms[j],xlim=xlim,ylim=range(y),main=paste(titre[j],"(min)"),col=1,
					xaxt="n")
				axis(1,x)
				for (k in seq(along=gpl)[-1]) {
					if (all(is.na(gpl[[k]][,1,2,j]))) next
					matlines(ttime,gpl[[k]][,1,2,j],lty=ltymnx[2],col=k)
				}

				y = sapply(gpl,function(x) x[,1,3,j])
				matplot(ttime,gp1[,1,3,j],type="l",lty=ltymnx[3],xlab=xlab,
					ylab=gpnoms[j],xlim=xlim,ylim=range(y),main=paste(titre[j],"(max)"),col=1,
					xaxt="n")
				axis(1,x)
				for (k in seq(along=gpl)[-1]) {
					if (all(is.na(gpl[[k]][,1,3,j]))) next
					matlines(ttime,gpl[[k]][,1,3,j],lty=ltymnx[3],col=k)
				}
			}

			if (! hasx11) dev.off()
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
