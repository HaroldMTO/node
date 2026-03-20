library(mfnode)

args = commandArgs(trailingOnly=TRUE)
largs = strsplit(args,split="=")
cargs = lapply(largs,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(largs,function(x) x[1])

if ("png" %in% names(cargs)) {
	cat("--> sending plots to PNG files in",cargs$png,"\n")
	stopifnot(file.info(cargs$png)$isdir)
} else if (! capabilities("X11")) {
	stop("no X11 device, please pass option png")
} else {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
}

cat("Search files - path/pattern:",cargs$path,cargs$patt,"\n")
fnode = dir(cargs$path,cargs$patt,full.names=TRUE,recursive=TRUE)
if (length(fnode) == 0) {
	cat("--> no file found\n")
	q("no")
}

dd = as.Date(gsub(".*(\\<\\d{2}/20\\d{6})/.+","\\1",fnode),"%H/%Y%m%d")
res = as.numeric(gsub(".*(\\<\\d{2})/20\\d{6}/.+","\\1",fnode))
if ("res" %in% names(cargs)) {
	cat("--> limit base times to",cargs$res,"\n")
	ind = res %in% as.numeric(cargs$res)
	fnode = fnode[ind]
	dd = dd[ind]
	res = res[ind]
}

dd = dd+res/24
if ("start" %in% names(cargs)) {
	cat("--> limit time to after",cargs$start,"\n")
	ind = as.POSIXct(dd) >= as.POSIXct(cargs$start)
	fnode = fnode[ind]
	dd = dd[ind]
}

if ("end" %in% names(cargs)) {
	cat("--> limit time to before",cargs$end,"\n")
	ind = as.POSIXct(dd) <= as.POSIXct(cargs$end)
	fnode = fnode[ind]
	dd = dd[ind]
}

ind = order(dd)
fnode = fnode[ind]
dd = dd[ind]
nf = length(fnode)
cat("--> found",nf,"log files - date range:",as.character(range(dd)),"\n")

nomj = c("JO","JB","JC","JP")
jacs = array(dim=c(nf,3,4))
nomg = c("Ritz1","Ritz2","gradient","J")
grads = array(dim=c(nf,3,4))
sps = gps = fps = fpsps = NULL
rts = rtx = numeric(nf)
ij2 = 6
ig2 = 10

if (interactive()) browser()

spref = "NORMS AT (START|NSTEP|END) CNT4"

for (i in seq(nf)) {
	nd = readLines(fnode[i],skipNul=TRUE)
	ind = grep("GREPCOST - ITER",nd)
	if (length(ind) > 1) {
		jacob = t(matrix(numlines(nd[ind]),ncol=length(ind)))
		jac = sub("GREPCOST.+SIM,(.+JCVARBC).+","\\1",nd[ind[1]])
		jac = trimws(strsplit(jac,split=",")[[1]])
		stopifnot(all(nomj %in% jac))
		dimnames(jacob)[[2]] = jac
		stopifnot(dim(jacob)[1] > ij2)
		jacs[i,,] = jacob[c(1,ij2,dim(jacob)[1]),nomj]
	}

	ind = grep("ritz values",nd,ignore.case=TRUE)
	if (length(ind) > 2) {
		ritz = t(matrix(numlines(nd[ind[-1]]),nrow=2))
		grad = numlines(grep("estimated reduction in norm",nd,value=TRUE,ignore.case=TRUE))
		quad = numlines(grep("estimated quadratic cost",nd,value=TRUE,ignore.case=TRUE))
		stopifnot(dim(ritz)[1] > ig2)
		stopifnot(length(grad) > ig2)
		stopifnot(length(quad) > ig2)
		grads[i,,1:2] = ritz[c(1,ig2,dim(ritz)[1]),]
		grads[i,,3] = grad[c(1,ig2,length(grad))]
		grads[i,,4] = quad[c(1,ig2,length(quad))]
	}

	spn = spnorm(nd,0,spref,abbrev=FALSE)
	if (! is.null(spn)) {
		if (is.null(sps)) {
			sps = array(NA_real_,dim=c(nf,dim(spn)[-1]),
				dimnames=c(list(seq(nf)),dimnames(spn)[-1]))
		}

		sps[i,,] = spn[1,,]
	}

	gpn = gpnorm(nd,0,spref)
	if (! is.null(gpn)) {
		if (is.null(gps)) {
			gps = array(NA_real_,dim=c(nf,dim(gpn)[-1]),
				dimnames=c(list(seq(nf)),dimnames(gpn)[-1]))
		}

		gps[i,,,] = gpn[1,,,]
	}

	nd = sub("\\<Savg(\\w+)","S000\\1",nd)

	fpn = fpgpnorm(nd,0,quiet=TRUE)
	if (! is.null(fpn)) {
		if (is.null(fps)) {
			fps = array(NA_real_,dim=c(nf,dim(fpn)[-1]),
				dimnames=c(list(seq(nf)),dimnames(fpn)[-1]))
		}

		fps[i,,,] = fpn[1,,,]
	}

	fpspn = fpspnorm(nd,0)
	if (! is.null(fpspn)) {
		if (is.null(fpsps)) {
			fpsps = array(NA_real_,dim=c(nf,dim(fpspn)[-1]),
				dimnames=c(list(seq(nf)),dimnames(fpspn)[-1]))
		}

		fpsps[i,] = fpspn[1,]
	}

	rt = runtime(nd)
	if (! is.null(rt)) {
		rts[i] = as.numeric(rt$wall[dim(rt)[1]]-attr(rt,"start"),units="secs") %% 86400
		rtx[i] = as.numeric(max(rt$dwall),units="secs")
	}
}

if ("save" %in% names(cargs)) {
	cat("Save GP/SP norms in",cargs$save,"\n")
	save(dd,sps,gps,fps,fpsps,file=cargs$save)
}

if (any(! is.na(rts))) {
	png(sprintf("%s/rts.png",cargs$png))
	par(mfrow=c(2,1),mgp=c(2,1,0))
	plot(dd,rts,type="o",main="Elapsed time",pch=20,xlab="Date",ylab="Time (s)",xaxt="n")
	axis.Date(1,dd)

	plot(dd,rtx,type="o",main="Max wall-time of steps",pch=20,xlab="Date",ylab="Time (s)",
		xaxt="n")
	axis.Date(1,dd)
	dev.off()
}

if (any(! is.na(jacs))) {
	cat("Jacobian values\n")
	png(sprintf("%s/jacobts.png",cargs$png))
	tt = "Iter 1, 6 and last one"
	par(mfrow=c(2,2),mgp=c(2,1,0))
	for (i in seq(along=nomj)) {
		matplot(dd,jacs[,,i],type="o",main=c(nomj[i],tt),xlab="Date",ylab=nomj[i],lty=1:3,
			col=1,pch=20,cex=.7,xaxt="n")
		axis.Date(1,dd)
		#legend("topleft",legend=c("1st","last"),col=1:2,lty=1,pch=20,bg="transparent")
	}

	dev.off()
}

if (any(! is.na(grads))) {
	cat("Ritz values\n")
	png(sprintf("%s/gradts.png",cargs$png))
	tt = "Iter 1, 10 and last one"
	par(mfrow=c(2,2),mgp=c(2,1,0))
	for (i in seq(along=nomg)) {
		matplot(dd,grads[,,i],type="o",main=c(nomg[i],tt),xlab="Date",ylab=nomg[i],lty=1:3,
			col=1,pch=20,cex=.7,xaxt="n")
		axis.Date(1,dd)
	}

	dev.off()
}

if (! is.null(sps)) {
	cat("Spectral norms\n")
	spnoms = dimnames(sps)[[3]]
	tt = paste("SP norm of",spnoms)
	nt = dim(sps)[1]
	nv = length(spnoms)
	# 2 rows, 2 columns
	nc = 2
	nr = min(1+(nv-1)%/%2,2)
	nj = nr*nc

	for (i in seq((nv-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/spnormts%d.png",cargs$png,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nv-nj*i,nj)+nj*i) {
			cat(". SP field",spnoms[j],"\n")
			plot(dd,sps[,,j],type="o",main=tt[j],xlab="Date",ylab=spnoms[j],lty=1,
				col=1,pch=20,cex=.7,xaxt="n")
			axis.Date(1,dd)
		}

		dev.off()
	}
}

if (! is.null(gps)) {
	cat("GP norms\n")
	gpnoms = dimnames(gps)[[4]]
	tt = paste("GP norm of",gpnoms)
	nt = dim(gps)[1]
	nv = length(gpnoms)
	# 2 rows, 2 columns
	nc = 2
	#nr = min(1+(nv-1)%/%2,2)
	nr = 2
	nj = 2

	for (i in seq((nv-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/gpnormts%d.png",cargs$png,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nv-nj*i,nj)+nj*i) {
			cat(". GP field",gpnoms[j],"\n")
			matplot(dd,gps[,,1,j],type="o",main=tt[j],xlab="Date",ylab=gpnoms[j],lty=1:3,
				col=1,pch=20,cex=.7,xaxt="n")
			axis.Date(1,dd)

			matplot(dd,gps[,,,j],type="l",main=tt[j],xlab="Date",ylab=gpnoms[j],lty=c(1,3,3),
				col=1,cex=.7,xaxt="n")
			axis.Date(1,dd)
		}

		dev.off()
	}
}

if (! is.null(fpsps)) {
	cat("FP SP norms\n")
	fpnoms = dimnames(fpsps)[[2]]
	tt = paste("FP SP norm of",fpnoms)
	nt = dim(fpsps)[1]
	nv = length(fpnoms)
	# 2 rows, 2 columns
	nc = 2
	nr = min(1+(nv-1)%/%2,2)
	nj = nr*nc

	for (i in seq((nv-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/fpspnormts%d.png",cargs$png,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nv-nj*i,nj)+nj*i) {
			cat(". FP field",fpnoms[j],":",range(fpsps[,j],na.rm=TRUE),"\n")
			ina = which(is.na(fpsps[,j]))
			if (length(ina) == dim(fpsps)[1]) next

			plot(dd,fpsps[,j],type="o",main=tt[j],xlab="Date",ylab=fpnoms[j],lty=1,
				col=1,pch=20,cex=.7,xaxt="n")
			axis.Date(1,dd)
		}

		dev.off()
	}
}

if (! is.null(fps)) {
	cat("FP norms\n")
	fpnoms = dimnames(fps)[[4]]
	tt = paste("FP norm of",fpnoms)
	nt = dim(fps)[1]
	nv = length(fpnoms)
	# 2 rows, 2 columns
	nc = 2
	#nr = min(1+(nv-1)%/%2,2)
	nr = 2
	nj = 2

	for (i in seq((nv-1)%/%nj+1)-1) {
		ficpng = sprintf("%s/fpnormts%d.png",cargs$png,i)
		png(ficpng)

		par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nv-nj*i,nj)+nj*i) {
			cat(". FP field",fpnoms[j],"\n")
			matplot(dd,fps[,,1,j],type="o",main=tt[j],xlab="Date",ylab=fpnoms[j],lty=1:3,
				col=1,pch=20,cex=.7,xaxt="n")
			axis.Date(1,dd)

			matplot(dd,fps[,,,j],type="l",main=tt[j],xlab="Date",ylab=fpnoms[j],lty=c(1,3,3),
				col=1,cex=.7,xaxt="n")
			axis.Date(1,dd)
		}

		dev.off()
	}
}
