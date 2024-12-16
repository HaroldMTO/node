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
detail = regexpr("gp.+detail",ptype) > 0
gpre = getarg("gpre",args)
gpref = getarg("gpref",args)
leg = getarg("leg",args)

fnode = grep("=",args,invert=TRUE,value=TRUE)

cat("Read file",fnode[1],"\n")
nd = readLines(fnode[1])
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
nstop = getvar("NSTOP",nd)

if (interactive()) browser()

cat("Parse GP norms of type",ptype,"\n")
if (is.null(gpref)) {
	gpref = "NORMS AT (START|NSTEP|END) CNT4"
	ind = grep(gpref,nd)
	if (length(ind) == 0) gpref = ""
}

if (regexpr("gfl",ptype) > 0) {
	gp1 = gpnorm(nd,lev,gpref,gpout=gpfre)
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
if (nstop == 0 && dim(gp1)[1] > 1) cat("--> steps are events of the job\n")

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
times = tstep*istep

gpl = list(gp1)

gpnoms = dimnames(gp1)[[4]]

if (is.null(gpre) && length(fnode) == 1) {
	leg = "GP"
	nt = length(times)
	nv = length(gpnoms)
	cat("Select times and variables among:
times:",head(times[-nt]),"...",times[nt],"
vars:",head(gpnoms[-nv]),"...",gpnoms[nv],"\n")
	sp1 = spnorm(nd,lev,abbrev=FALSE)
	if (! is.null(sp1)) {
		indt = match(times,dimnames(sp1)[[1]])
		indv = match(gpnoms,dimnames(sp1)[[3]])
		if (any(! is.na(indt)) && any(! is.na(indv))) {
			leg = c(leg,"SP")
			sp1 = sp1[indt,,indv,drop=FALSE]
			gpl = c(gpl,list(sp1))
		} else {
			cat("--> SP times:",dimnames(sp1)[[1]],"\n")
			cat("--> SP vars:",dimnames(sp1)[[3]],"\n")
		}
	}

	if (length(leg) == 1) {
		cat("--> no mixed norms, quit\n")
		quit("no")
	}
} else {
	if (length(fnode) > 1) {
		if (is.null(leg)) leg = sub(".*node\\.?","",fnode,ignore.case=TRUE)
		gpre = rep(gpref,length(fnode)-1)
	} else if (is.null(leg)) {
		leg = c("t0",sub("gpnorm +(gmv|gfl)?","",gpre))
	}

	for (i in seq(along=gpre)) {
		if (length(fnode) > 1) {
			cat("Read file",fnode[i+1],"\n")
			nd = readLines(fnode[i+1])
		}

		cat("Parse GP norms, pattern",gpre[i],"\n")
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
		if (all(is.na(inds))) next

		gpl[[i+1]] = gpi[inds,,,indv,drop=FALSE]
	}
}

leg = leg[which(! sapply(gpl,is.null))]
gpl = gpl[! sapply(gpl,is.null)]

nf = length(gpnoms)
nt = dim(gp1)[1]

if (length(lev) > 1) {
	cat("Produce vertical profiles\n")
	tt = sprintf("GP norm of %s",gpnoms)
	indt = c(1,(nt+1)%/%2,nt)
	ts = sprintf("t%s",paste(istep[indt],collapse="/"))
	nc = 3
	nr = min(length(gpl),2)
	nj = 1

	for (i in seq((nf-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/%snormv%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			cat(". GP field",gpnoms[j],"in",ficpng,"\n")
			y = aperm(gp1[indt,,,j],c(2,1,3))
			plotvmnx(y,lev,xlab=gpnoms[j],main=c(tt[j],ts))

			if (nr == 1) next

			y = sapply(gpl,function(x) x[indt[2],,,j],simplify="array")
			y = aperm(y,c(1,3,2))
			il = which(apply(y,2,function(x) any(! is.na(x))))
			ts2 = sprintf("t%d",indt[2])
			plotvmnx(y[,il,,drop=FALSE],lev,xlab=gpnoms[j],main=c(tt[j],ts2),legend=leg[il],
				lty=2,col=il)
		}

		dev.off()
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
	if (nt == 1) {
		cat("1 time-step only, quit\n")
		quit("no")
	}

	if (length(gpl) > 1) {
		con = file(sprintf("%s/%s.txt",pngd,ptype),"wt")
		for (i in seq(along=gpl)[-1]) {
			cat("\nNorm diff. tendency for patterns:",sprintf("'%s'",leg[(i-1):i]),"\n",
				file=con)
			iv = which(! is.na(dimnames(gpl[[i]])[[4]]))
			gpdiff = gpl[[i]][,1,1,iv,drop=FALSE]-gpl[[i-1]][,1,1,iv,drop=FALSE]
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
		ficpng = sprintf("%s/%snorm%d.png",pngd,ptype,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nj*i,nj)+nj*i) {
			y = sapply(gpl,function(x) x[,1,,j],simplify="array")
			il = which(apply(y,3,function(x) any(! is.na(x))))
			scal = 1/scale10(y[,1,])
			if (! is.finite(scal) || .001 <= scal && scal < 1) scal = 1
			plotmean(ttime,y[,1,il],main=titre[j],leg[il],tunit,xlim=xlim,
				xlab=xlab,ylab=gpnoms[j],xaxp=xaxp,scale=scal,col=il)

			plotmnx(ttime,y,titre[j],xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)

			if (! detail) next

			plotmnx(ttime,y,titre[j],imnx=2,xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)
			plotmnx(ttime,y,titre[j],imnx=3,xlim=xlim,xlab=xlab,ylab=gpnoms[j],xaxp=xaxp)
		}

		dev.off()
	}
}
