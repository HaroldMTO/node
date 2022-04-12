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

	tstep = as.numeric(sub("TSTEP += *([-+]?[0-9][-+0-9E.]+) .+","\\1",s))
}

spnorm = function(nd,lev)
{
	nd = grep("^ *$",nd,value=TRUE,invert=TRUE)
	ind = grep("SPECTRAL NORMS",nd)

	inds = grep("NORMS AT NSTEP CNT4",nd[ind-1])
	ind = ind[inds]

	spsp = as.numeric(gsub("SPECTRAL NORMS.+ ([-0-9.E+]+|NaN)$","\\1",nd[ind]))

	has.levels = regexpr(sprintf("^ *1 +%s",Gnum),nd[ind[1]+3]) > 0 &&
		regexpr(sprintf("^ *\\d+ +%s",Gnum),nd[ind[1]+nflevg+1]) > 0

	noms = strsplit(nd[ind[1]+1]," {2,}")[[1]][-1]
	noms[noms == "KINETIC ENERGY"] = "TKE"

	if (length(lev) > 1) {
		stopifnot(has.levels)
		stopifnot(all(! duplicated(lev)) && all(lev %in% 1:nflevg))
	} else if (lev > 0) {
		stopifnot(has.levels)
	}

	indi = rep(ind,each=length(lev))+lev+2

	spre = regmatches(nd[indi],gregexpr(sprintf("(%s|NaN)",Gnum),nd[indi]))
	spre = lapply(spre,function(x) gsub("(\\d+)(\\-\\d+)","\\1E\\2",x))
	spn = sapply(spre,as.numeric)

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

		spre = regmatches(nd[indi],gregexpr(sprintf("(%s|NaN)",Gnum),nd[indi]))
		spre = lapply(spre,function(x) gsub("(\\d+)(\\-\\d+)","\\1E\\2",x))
		spnh = sapply(spre,as.numeric)

		nomsnh = strsplit(nd[indn]," {2,}")[[1]][-1]
		nomsnh[nomsnh == "LOG(PRE/PREHYD)"] = "LOG(P/P_hyd)"
		nomsnh[nomsnh == "d4 = VERT DIV + X"] = "d4 (= vdiv+X)"

		dim(spnh) = c(length(nomsnh),length(lev),nt)
		spnh = aperm(spnh,c(3,2,1))
		spn = c(spn,spnh)
		noms = c(noms,nomsnh)
		dim(spn) = c(nt,length(lev),length(noms))
	}

	spn[is.na(spn)] = 0

	ip = grep("LOG\\(P/P_hyd\\)|d4",noms,invert=TRUE)
	noms[ip] = abbreviate(noms[ip])

	istep = sub("NORMS AT NSTEP CNT4( \\(PREDICTOR\\))? +(\\d+)","\\2",nd[ind-1])
	dimnames(spn) = list(istep,lev,noms)

	spn
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
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")
ask = hasx11 && interactive()

cat("Read SP norm in file",fnode,"\n")
nd = readLines(fnode)

nflevg = getvar("NFLEVG",nd)
if (is.null(lev)) lev = seq(nflevg)

nstop = getvar("NSTOP",nd)
nfrsdi = getvar(".+ NFRSDI",nd)
tstep = timestep(nd)

sp1 = spnorm(nd,lev)

times = as.integer(dimnames(sp1)[[1]])*tstep

noms = dimnames(sp1)[[3]]
titre = paste("Spectral norm of",noms)

nf = length(noms)
nt = length(times)

if (length(lev) > 1) {
	cat("Spectral norms at step",sprintf("%gs (= %.5gh)",times[nt],times[nt]/3600),"\n")
	ylim = rev(range(lev))
	nc = min(nf,3)
	nr = min(1+(nf-1)%/%nc,3)

	np = nc*nr
	for (i in seq((nf-1)%/%np+1)-1) {
		if (! hasx11) png(sprintf("spnormv%d.png",i))

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		for (j in seq(min(nf-np*i,np)+np*i)) {
			plot(sp1[nt,,j],lev,type="l",ylim=ylim,lty=1,xlab=noms[j],
				ylab="Level",main=titre[j])
			abline(v=0,col="darkgrey",lty=2)
		}

		if (ask) invisible(readline("Press enter to continue"))
	}
} else {
	if (dim(sp1)[1] == 1) stop("1 time-step only (stop)\n")

	if (diff(range(times)) < 600) {
		tunit = 1
		xlab = "fc time (s)"
	} else if (diff(range(times)) < 14400) {
		tunit = 60
		xlab = "fc time (mn)"
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
		titre = paste("Spectral norm of",noms)
	} else {
		titre = paste("Spectral norm of",noms,"- lev",lev)
	}

	# 4 rows max, np plots per page
	nr = min(1+(nf-1)%/%2,4)
	np = 2*nr
	for (i in seq((nf-1)%/%np+1)-1) {
		if (! hasx11) png(sprintf("spnorm%d.png",i))

		par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-np*i,np)+np*i) {
			plot(ttime,sp1[,1,j],type="l",xlab=xlab,ylab=noms[j],main=titre[j],xaxt="n")
			axis(1,x)
			reg = line(ttime,sp1[,1,j])
			abline(reg,lty=2)
			tend = coef(reg)[2]*86400/tunit
			if (FALSE && tunit == 86400) {
				tt = sprintf("Tend: %+.3g [unit]/day, variation: %.2g",tend,
					tend*diff(range(ttime)))
			} else {
				tt = sprintf("Tend: %+.3g [unit]/day",tend)
			}

			mtext(tt,lty=2,cex=par("cex"))
		}

		if (ask) invisible(readline("Press enter to continue"))
	}
}
