Gnum = "-?[0-9]*\\.[0-9]+[-+0-9E]+[0-9]"

getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(\\d+).*",var,sep)
	unique(as.integer(gsub(re,"\\1",grep(re,nd,value=TRUE))))
}

timestep = function(nd)
{
	s = grep("^ *TSTEP += *[-+]?[0-9][-+0-9E.]+",nd,value=TRUE)
	stopifnot(length(s) == 1)

	as.numeric(sub("TSTEP += *([-+]?[0-9][-+0-9E.]+) .+","\\1",s))
}

gpnorm2D = function(nd)
{
	nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
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

		gpre = regmatches(nd[indi],gregexpr(sprintf("(%s|NaN)",Gnum),nd[indi]))
		gpre = lapply(gpre,function(x) gsub("(\\d+)(\\-\\d+)","\\1E\\2",x))
		gpn = t(sapply(gpre,as.numeric))
		dimnames(gpn) = list(gnames,c("ave","min","max"))
		surf[[i]] = gpn
	}

	surf[sapply(surf,length) > 0]
}

gpnorm = function(nd,lev)
{
	nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
	ind = grep("GPNORM +\\w+.* +AVERAGE",nd)
	indo = grep("GPNORM OUTPUT",nd)
	ind = ind[! ind %in% indo]
	if (length(ind) == 0) return(NULL)

	has.levels = regexpr(sprintf("^ *1 +%s",Gnum),nd[ind[1]+2]) > 0 &&
		regexpr(sprintf("^ *\\d+ +%s",Gnum),nd[ind[1]+nflevg+1]) > 0

	if (length(lev) > 1) {
		stopifnot(has.levels)
		stopifnot(all(! duplicated(lev)) && all(lev %in% 1:nflevg))
	} else if (lev > 0) {
		stopifnot(has.levels)
	}

	indi = rep(ind,each=length(lev))+lev+1

	gpre = regmatches(nd[indi],gregexpr(sprintf("(%s|NaN)",Gnum),nd[indi]))
	gpre = lapply(gpre,function(x) gsub("(\\d+)(\\-\\d+)","\\1E\\2",x))
	gpn = sapply(gpre,as.numeric)
	gpn[is.na(gpn)] = 0

	noms = unique(sub(" +GPNORM +((\\w|\\.)+) +AVERAGE.+","\\1",nd[ind]))
	nt = length(gpn)/(3*length(lev)*length(noms))
	stopifnot(nt == as.integer(nt))

	dim(gpn) = c(3,length(lev),length(noms),nt)
	gpl = aperm(gpn,c(4,2,1,3))

	dimnames(gpl) = list(NULL,lev,c("ave","min","max"),noms)

	gpl
}

args = commandArgs(trailingOnly=TRUE)
lev = NULL
ind = grep("lev=\\d+",args)
if (length(ind) > 0) {
	lev = strsplit(sub("lev=","",args[ind[1]]),split=":")[[1]]
	lev = as.integer(lev)
}

fnode = grep("=",args,invert=TRUE,value=TRUE)
stopifnot(length(fnode) == 1)

hasx11 = capabilities("X11")
ask = hasx11 && interactive()
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

cat("Read GP norm in file",fnode,"\n")
nd = readLines(fnode)
nflevg = getvar("NFLEVG",nd)

if (is.null(lev)) lev = seq(nflevg)

gp1 = gpnorm(nd,lev)
if (is.null(gp1)) stop("no GP norms")

nstop = getvar("NSTOP",nd)
nfrgdi = getvar(".+ NFRGDI",nd)
tstep = timestep(nd)

istep = seq(0,nstop,by=nfrgdi)
if (length(istep) > dim(gp1)[1]) length(istep) = dim(gp1)[1]
times = tstep*istep

noms = dimnames(gp1)[[4]]

nf = length(noms)
nt = length(times)

if (length(lev) > 1) {
	cat("GP norms at step",sprintf("%gs (= %.5gh)",times[nt],times[nt]/3600),"\n")
	tt = sprintf("%s - step %d",c("Average","Minimum","Maximum"),istep[nt])
	ylim = rev(range(lev))
	nr = min(1+(nf-1)%/%3,3)

	for (i in seq((nf-1)%/%nr+1)-1) {
		if (! hasx11) png(sprintf("gpnormv%d.png",i))

		par(mfrow=c(nr,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in seq(min(nf-nr*i,nr)+nr*i)) {
			for (k in c(2,1,3)) {
				plot(gp1[nt,,k,j],lev,type="l",ylim=ylim,xlab=noms[j],ylab="Level",main=tt[k])
				abline(v=0,col="darkgrey",lty=2)
			}
		}

		if (ask) invisible(readline("Press enter to continue"))
	}
} else {
	if (dim(gp1)[1] == 1) stop("1 time-step only (stop)\n")

	if (diff(range(times)) < 86400) {
		tunit = 1
		xlab = "fc time (s)"
	} else if (diff(range(times)) < 86400*6) {
		tunit = 3600
		xlab = "fc time (h)"
	} else {
		tunit = 86400
		xlab = "fc time (days)"
	}

	ttime = times/tunit
	if (tunit == 86400) {
		x = pretty(ttime,7)
	} else {
		x = pretty(ttime/6,7)*6
	}

	if (lev == 0) {
		titre = paste("GP norm of",noms)
	} else {
		titre = paste("GP norm of",noms,"- lev",lev)
	}

	# 3 rows, 6 plots per page
	nr = 3
	for (i in seq((nf-1)%/%nr+1)-1) {
		if (! hasx11) png(sprintf("gpnorm%d.png",i))

		par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			plot(ttime,gp1[,1,1,j],type="l",xlab=xlab,ylab=noms[j],
				main=paste(titre[j],"average"),xaxt="n")
			reg = line(ttime,gp1[,1,1,j])
			abline(reg,lty=2)
			mtext(sprintf("Tend: %+.3g [unit]/day",coef(reg)[2]*86400/tunit),lty=2,
				cex=par("cex"))
			axis(1,x)
			if (all(gp1[,1,,j] == 0)) text((1+dim(gp1)[1])/2,.7,"all values = 0")
			matplot(ttime,gp1[,1,,j],type="l",lty=c(1,2,2),xlab=xlab,ylab=noms[j],
				main=paste(titre[j],"(ave/min/max)"),col=1,xaxt="n")
			axis(1,x)
		}

		if (ask) invisible(readline("Press enter to continue"))
	}
}

surf = gpnorm2D(nd)
if (length(surf) == 0) {
	cat("--> no GP norms for surface fields\n")
} else {
	msurf = surf[[1]]
	for (i in seq(along=surf)[-1]) msurf = rbind(msurf,surf[[i]])
	cat("Surface fields:\n")
	print(msurf)
}
