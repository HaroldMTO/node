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

gpnorm = function(nd,lev,ind)
{
	if (missing(ind)) {
		ind = grep("GPNORM +\\w+.* +AVERAGE",nd)
		indo = grep("GPNORM OUTPUT",nd[ind],invert=TRUE)
		ind = ind[indo]
	}

	if (length(ind) == 0) stop("no GP norms")

	indi = rep(ind,each=length(lev))+lev+1

	gpn = line2num(nd[indi])

	noms = unique(sub(" *GPNORM +(\\w+.+?) +AVERAGE.+","\\1",nd[ind]))
	noms[noms == "SURFACE PRESSURE"] = "SURF P"
	noms[noms == "TEMPRATURE"] = "TEMP"
	noms[noms == "U VELOCITY"] = "U VELOC."
	noms[noms == "V VELOCITY"] = "V VELOC."

	nt = length(gpn)/(3*length(lev)*length(noms))
	stopifnot(nt == as.integer(nt))

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

	istep = sub("NORMS AT NSTEP CNT4( \\(PREDICTOR\\))? +(\\d+)","\\2",nd[ind-1])
	dimnames(spn) = list(istep,lev,noms)

	spn
}

countfield = function(ind,ind2,nl2)
{
	which(diff(ind[seq(ind2[1],ind2[2])]) > nl2)[1]
}

indexpand = function(ind,nf,nl)
{
	rep(ind,each=nf)+(seq(nf)-1)*nl2
}

hasx11 = capabilities("X11")
ask = hasx11 && interactive()
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

xaxis = data.frame(unit=c(1,3600,86400),label=sprintf("fc time (%s)",c("s","h","days")),
	mindiff=c(0,86400,6*86400),freq=c(6,6,1))

args = commandArgs(trailingOnly=TRUE)
lev = as.integer(getarg("lev",args))
hmin = as.numeric(getarg("hmin",args))
hmax = as.numeric(getarg("hmax",args))
splot = getarg("plot",args)
if (is.null(splot)) splot = c("spec","gp")
type = getarg("type",args)
if (is.null(type)) type = "gpgfl"
spre = getarg("spre",args)
spref = getarg("spref",args)
if (is.null(spref)) spref = "spnorm t0"
gpre = getarg("gpre",args)
gpref = getarg("gpref",args)
if (is.null(gpref)) gpref = "gpnorm gmvt0"

fnode = grep("=",args,invert=TRUE,value=TRUE)
stopifnot(length(fnode) == 1)

cat("Read file",fnode,"\n")
nd = readLines(fnode)
nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
nflevg = getvar("NFLEVG",nd)
has.levels = getvar("NSPPR",nd) > 0

if (length(lev) == 0) {
	stopifnot(has.levels)
	lev = seq(nflevg)
} else if (length(lev) > 1) {
	stopifnot(has.levels)
	stopifnot(all(! duplicated(lev)) && all(lev %in% 1:nflevg))
} else if (lev > 0) {
	stopifnot(has.levels)
}

nstop = getvar("NSTOP",nd)
tstep = getvar("TSTEP",nd)
i1 = grep("START CNT4",nd)

if (any(splot == "sp")) {
	ind = grep("SPECTRAL NORMS",nd)
	ind = ind[ind > i1]
	ind1 = grep(sub("spnorm t0","NORMS AT NSTEP CNT4",spref),nd[ind-1])
	sp1 = spnorm(nd,lev,ind[ind1[-1]])
	spnoms = dimnames(sp1)[[3]]
	spl = list(sp1)
	leg = "t0"
	for (i in seq(along=spre)) {
		ind2 = grep(spre[i],nd[ind-1])
		if (length(ind2) == 0) {
			cat("--> no norms for pattern:",spre[i],"\n")
			next
		}

		indv = match(spnoms,dimnames(sp1)[[3]])
		stopifnot(any(! is.na(indv)))
		spi = spnorm(nd,lev,ind[ind2])
		spl[[i+1]] = spi[,,indv,drop=FALSE]
		leg[i+1] = sub("spnorm +","",spre[i])
	}

	nfrsdi = getvar(".+ NFRSDI",nd)
	istep = seq(1,nstop,by=nfrsdi)

	nt = dim(sp1)[1]
	if (length(istep) > nt) length(istep) = nt
	times = tstep*istep

	if (length(lev) > 1) {
		cat("Norms at step",sprintf("%gs (= %.5gh)",times[nt],times[nt]/3600),"\n")
		nf = length(spnoms)
		tt = sprintf("Spectral norm of %s - step %d",spnoms,istep[nt])
		ylim = rev(range(lev))
		nc = min(nf,3)
		nr = min(1+(nf-1)%/%nc,3)
		np = nc*nr

		for (i in seq((nf-1)%/%np+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("spnormv%d.png",i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			for (j in seq(min(nf-np*i,np)+np*i)) {
				plot(sp1[nt,,j],lev,type="l",ylim=ylim,xlab=spnoms[j],ylab="Level",main=tt[j])
				abline(v=0,col="darkgrey",lty=2)
			}
		}
	} else {
		if (dim(sp1)[1] == 1) stop("1 time-step only (stop)\n")

		cat("Norm diff. tendency for patterns:\n",
			paste(sprintf("'%s'",spre),collapse=" "),"\n")
		for (i in seq(along=spl)[-1]) {
			spdiff = spl[[i]][,1,]-spl[[i-1]][,1,]
			iv = which(! is.na(dimnames(spdiff)[[2]]))
			a = apply(spdiff[,iv],2,function(x) coef(line(x))[2])
			s = apply(spdiff[,iv],2,function(x) sd(residuals(line(x))))
			cat("\n")
			print(rbind(tend=a,sd=s))
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
		nf = length(spnoms)
		nr = min(1+(nf-1)%/%2,4)
		np = 2*nr
		for (i in seq((nf-1)%/%np+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("spnorm%d.png",i))

			par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-np*i,np)+np*i) {
				y = sapply(spl,function(x) x[,1,j])
				col = lty = seq(along=spl)
				lty[-1] = 2
				matplot(ttime,y,type="l",lty=lty,col=col,xlim=xlim,xlab=xlab,
					ylab=spnoms[j],main=titre[j],xaxt="n")
				legend("topleft",leg,col=col,lty=lty)

				axis(1,x)
				reg = line(ttime,sp1[,1,j])
				abline(reg,lty=2)
				tend = coef(reg)[2]*86400/tunit
				tt = sprintf("Tend: %+.3e [unit]/day",tend)
				if (tunit == 86400) {
					tt = sprintf("%s, variation: %.2g",tt,tend*diff(range(ttime)))
				}

				mtext(tt,lty=2,cex=par("cex"))
			}
		}
	}
}

if (any(splot == "gp")) {
	gpfre1 = "[UVW] VELOCITY|(SURFACE )?PRESSURE|TEMPERATURE|GRAD[LM]_\\w+|GEOPOTENTIAL"
	gpfre2 = "MOIST AIR SPECIF|ISOBARE CAPACITY|SURFACE DIV|d\\(DIV\\)\\*dP"
	gpfre3 = "(ATND|ADIAB|CTY|SISL)_\\w+"
	gpfre = paste(gpfre1,gpfre2,gpfre3,sep="|")

	if (type == "gpgfl") {
		ind = grep("GPNORM +\\w+.* +AVERAGE",nd)
		ind = ind[ind > i1]
		indo = grep(sprintf("GPNORM +(%s|OUTPUT) +AVERAGE",gpfre),nd[ind],invert=TRUE)
		gp1 = gpnorm(nd,lev,ind[indo])
		gp1 = gp1[-1,,,,drop=FALSE]
		gpl = list(gp1)
		leg = "t0"
	} else {
		nl2 = 2+has.levels*nflevg
		ind = grep(sprintf("GPNORM +(%s) +AVERAGE",gpfre),nd)
		ind = ind[ind > i1]
		ind1 = grep(gpref,nd[ind-1],ignore.case=TRUE)
		nf = countfield(ind,ind1,nl2)
		indi = indexpand(ind[ind1],nf,nl2)
		gp1 = gpnorm(nd,lev,indi)
		gpl = list(gp1)
		gpnoms = dimnames(gp1)[[4]]
		leg = sub("gpnorm +","",gpref)
		for (i in seq(along=gpre)) {
			ind2 = grep(gpre[i],nd[ind-1],ignore.case=TRUE)
			if (length(ind2) == 0) {
				cat("--> no norms for pattern",gpre[i],"\n")
				next
			}

			nf = countfield(ind,ind2,nl2)
			indi = indexpand(ind[ind2],nf,nl2)
			gpi = gpnorm(nd,lev,indi)
			indv = match(gpnoms,dimnames(gpi)[[4]])
			stopifnot(any(! is.na(indv)))
			gpl[[i+1]] = gpi[,,,indv,drop=FALSE]
			leg[i+1] = sub("gpnorm +","",gpre[i])
		}

		gpl = gpl[! sapply(gpl,is.null)]
	}

	nfrgdi = getvar(".+ NFRGDI",nd)
	istep = seq(1,nstop,by=nfrgdi)
	nt = dim(gp1)[1]
	if (length(istep) > nt) length(istep) = nt

	times = tstep*istep

	if (length(lev) > 1) {
		nf = length(gpnoms)
		tt = sprintf("%s - step %d",c("Average","Minimum","Maximum"),istep[nt])
		nc = 3
		nr = max(2,min(1+(nf-1)%/%nc,3))
		np = nc*nr

		for (i in seq((nf-1)%/%nr+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%snormv%d.png",type,i))

			par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in seq(min(nf-nr*i,nr)+nr*i)) {
				for (k in c(2,1,3)) {
					plot(gp1[nt,,k,j],lev,type="l",ylim=ylim,xlab=gpnoms[j],ylab="Level",
						main=tt[k])
					abline(v=0,col="darkgrey",lty=2)
				}
			}
		}
	} else {
		if (nt == 1) stop("1 time-step only (stop)\n")

		cat("Norm diff. tendency for patterns:\n",
			paste(sprintf("'%s'",gpre),collapse=" "),"\n")
		for (i in seq(along=gpl)[-1]) {
			gpdiff = gpl[[i]][,1,1,]-gpl[[i-1]][,1,1,]
			iv = which(! is.na(dimnames(gpdiff)[[2]]))
			a = apply(gpdiff[,iv],2,function(x) coef(line(x))[2])
			s = apply(gpdiff[,iv],2,function(x) sd(residuals(line(x))))
			cat("\n")
			print(rbind(tend=a,sd=s))
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

		# 3 rows max, 2 columns
		nf = length(gpnoms)
		nr = max(2,min(3,nf))
		for (i in seq((nf-1)%/%nr+1)-1) {
			if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
			if (! hasx11) png(sprintf("%snorm%d.png",type,i))

			par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (j in 1:min(nf-nr*i,nr)+nr*i) {
				y = sapply(gpl,function(x) x[,1,1,j])
				col = lty = seq(along=gpl)
				lty[-1] = 2
				matplot(ttime,y,type="l",lty=lty,col=col,xlim=xlim,xlab=xlab,ylab=gpnoms[j],
					main=paste(titre[j],"average"),xaxt="n")

				reg = line(ttime,gp1[,1,1,j])
				abline(reg,lty=2)
				tend = coef(reg)[2]*86400/tunit
				tt = sprintf("Tend: %+.3e [unit]/day",tend)
				if (tunit == 86400) {
					tt = sprintf("%s, variation: %.2g",tt,tend*diff(range(ttime)))
				}

				legend("topleft",leg,col=col,lty=lty)
				mtext(tt,lty=2,cex=par("cex"))
				axis(1,x)
				if (all(gp1[,1,,j] == 0)) text(sum(range(ttime))/2,.5,"all values = 0")

				matplot(ttime,gp1[,1,,j],type="l",lty=c(1,3,3),xlab=xlab,ylab=gpnoms[j],
					xlim=xlim,main=paste(titre[j],"(ave/min/max)"),col=1,xaxt="n")
				axis(1,x)
				for (k in seq(along=gpl)[-1]) {
					if (all(is.na(gpl[[k]][,1,,j]))) next
					matlines(ttime,gpl[[k]][,1,,j],lty=1:3,col=k)
				}
			}
		}
	}
}

if (any(splot == "surf")) {
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
