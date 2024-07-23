library(mfnode)

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
spre = getarg("spre",args)
spref = getarg("spref",args)
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

cat("Parse spectral norms\n")
if (is.null(spref)) {
	spref = "NORMS AT (START|NSTEP|END) CNT4"
	ind = grep(spref,nd)
	if (length(ind) == 0) spref = ""
}

sp1 = spnorm(nd,lev,spref)
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
if (nstop == 0 && dim(sp1)[1] > 1) cat("--> steps are events of the job\n")

spnoms = dimnames(sp1)[[3]]

if (length(fnode) > 1) {
	if (is.null(leg)) leg = sub(".*node\\.?","",fnode,ignore.case=TRUE)
	spre = rep(spref,length(fnode)-1)
} else if (is.null(leg)) {
	leg = c("t0",sub("spnorm +","",spre))
}

spl = list(sp1)
for (i in seq(along=spre)) {
	if (length(fnode) > 1) {
		cat("Read file",fnode[i+1],"\n")
		nd = readLines(fnode[i+1])
	}

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
		ficpng = sprintf("%s/spnormv%d.png",pngd,i)
		png(ficpng)

		par(mfcol=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in seq(min(nf-nj*i,nj))+nj*i) {
			y = t(sp1[indt,,j])
			plotvmean(y,lev,type="o",pch=".",cex=1.1,xlab=spnoms[j],main=c(tt[j],ts))

			if (nr == 1) next

			y = sapply(spl,function(x) x[indt[2],,j],simplify="array")
			ts2 = sprintf("t%s",istep[indt[2]])
			plotvmean(y,lev,type="o",pch=".",cex=1.1,xlab=spnoms[j],main=c(tt[j],ts2),
				lty=c(1,2,2),col=seq(along=spl),legend=leg)
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
		con = file(sprintf("%s/spec.txt",pngd),"wt")
		for (i in seq(along=spl)[-1]) {
			cat("\nNorm diff. tendency for patterns:",sprintf("'%s'",leg[(i-1):i]),"\n",
				file=con)
			iv = which(! is.na(dimnames(spl[[i]])[[3]]))
			spdiff = spl[[i]][,1,iv,drop=FALSE]-spl[[i-1]][,1,iv,drop=FALSE]
			a = apply(spdiff,3,function(x) coef(line(x))[2])
			s = apply(spdiff,3,function(x) sd(residuals(line(x))))
			evo = rbind(tend=a,sd=s)
			evo = format(evo,digits=3,width=12)
			cat("field",sprintf("% 12.12s",dimnames(evo)[[2]]),"\n",file=con)
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
		ficpng = sprintf("%s/spnorm%d.png",pngd,i)
		png(ficpng)

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

		dev.off()
	}
}
