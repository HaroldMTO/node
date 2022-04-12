getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(\\d+).*",var,sep)
	unique(as.integer(gsub(re,"\\1",grep(re,nd,value=TRUE))))
}

longend = function(nd,ndglg)
{
	ij = grep("\\( *JGL,NLOENG *\\)",nd)
	is = grep("Set up transforms",nd)

	ind = seq(ij+1,is-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *\\-?\\d+ +\\d+\\)",nd[ind])))
	nloeng = as.integer(gsub("\\( *\\-?\\d+ +(\\d+)\\)","\\1",s))
   stopifnot(length(nloeng)%/%2 == (length(nloeng)+1)%/%2)

   off = (length(nloeng)-ndglg)/2
	nloeng[off+seq(ndglg%/%2)]
}

stdatm = function(nd,nlev)
{
	ia = grep(" +pressure +temperature +height +density",nd,ignore.case=TRUE)
	ind = ia+seq(nlev)
	ire = regexec(" *\\d+ +(\\d+\\.\\d+) +(\\d+\\.\\d+) +(\\d+\\.\\d+) +(\\-?\\d+\\.\\d+)",nd[ind])
	P = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	T = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	Z = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	rho = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))

	data.frame(P=P,T=T,Z=Z,rho=rho)
}

abh = function(nd,nlev)
{
	ih = grep("A and B at half levels",nd)
	ind = ih+1+seq(nlev+1)
	ire = regexec(" *\\d+ +(\\d+\\.\\d+) +(\\-?\\d+\\.\\d+) +(\\-?\\d+\\.\\d+)",nd[ind])
	alh = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))

	data.frame(Ah=ah,Bh=bh,alpha=alh)
}

getgem = function(nd)
{
	ig = grep("Printings in SUGEM_NAML",nd,ignore.case=TRUE)
	mucen = as.numeric(sub(".* RMUCEN *= *(\\-?\\d+\\.[-+0-9eE]+).*","\\1",nd[ig+1]))
	locen = as.numeric(sub(".* RLOCEN *= *(\\-?\\d+\\.[-+0-9eE]+).*","\\1",nd[ig+1]))
	stret = as.numeric(sub(".* RSTRET *= *(\\-?\\d+\\.[-+0-9eE]+).*","\\1",nd[ig+1]))
	nlginc = as.numeric(sub(".* RNLGINC *= *(\\-?\\d+\\.[-+0-9eE]+).*","\\1",nd[ig+3]))

	list(mucen=mucen,locen=locen,stret=stret,nlginc=nlginc)
}

procmap = function(nd)
{
	la = lapply(regmatches(nd,gregexpr("\\d+",nd)),as.integer)
	ia = sapply(la,"[",1)

	# in case of FPOS job with change of resolution, there are prints for each geom
	ind = which(diff(ia) < 0)
	if (length(ind) > 0) {
		nd = nd[seq(ind[1])]
		la = lapply(regmatches(nd,gregexpr("\\d+",nd)),as.integer)
		ia = sapply(la,"[",1)
	}

	seta = tapply(la,ia,simplify2array)
	seta = lapply(seta,function(x) t(x[-(1:2),,drop=FALSE]))
}

latmap = function(sta,end,offlat,nlat)
{
	nx = max(end)
	ny = dim(sta)[1]
	z = array(NA_real_,dim=c(nx,ny))

	for (ib in seq(dim(sta)[2])) {
		for (j in seq(ny)) {
			i = seq(sta[j,ib],end[j,ib])
			z[i,ny+1-j] = ib
		}
	}

	list(x=seq(nx),y=offlat+seq(ny),z=z)
}

xoff = function(x)
{
	x[x > 0] = x[x > 0]-1
	x
}

plotmap = function(sta,onl)
{
	mp1 = sapply(sta,dim)
	nlat = sum(mp1[1,])
	end = mapply("+",sta,onl,SIMPLIFY=FALSE)
	end = lapply(end,xoff)
	nlon = max(unlist(end))

	nlon1 = max(end[[1]][1,])
	ngp = sum(sapply(onl,sum))
	ngp1 = sum(onl[[1]][,1])

	isplit = sapply(sta,function(x) any(x[1,] == 0) || x[1,1] > 1)
	nlat = nlat-sum(isplit)
	offlat = nlat-(cumsum(mp1[1,])-cumsum(isplit))
	titre = c("MPI tasks and their grid-points",
		sprintf("#gp: %d %d - nlat: %d - ndlon min/max: %d %d",ngp1,ngp,nlat,nlon1,nlon))
	lmap = latmap(sta[[1]],end[[1]],offlat[1],nlat)
	image(lmap,xlim=c(1,nlon),ylim=c(1,nlat),col=2+(seq(mp1[2,1])-1)%%7,
		main=titre,xlab="Longitude index",ylab="Latitude index",yaxt="n")
	y = pretty(seq(nlat))
	axis(2,at=nlat+1-y,labels=y)

	if (length(sta) > 1) {
		for (ia in seq(along=sta)[-1]) {
			lmap = latmap(sta[[ia]],end[[ia]],offlat[ia],nlat)
			image(lmap,add=TRUE,col=2+(2*(ia-1)+seq(mp1[2,ia])-1)%%7)
		}
	}

	j = mp1[1,]%/%2
	x = unlist(lapply(seq(along=j),function(i) sta[[i]][j[i],]+onl[[i]][j[i],]%/%2))
	y = offlat+j
	y = unlist(lapply(seq(along=y),function(i) rep(y[i],each=mp1[2,i])))
	procs = seq(sum(mp1[2,]))
	if (max(mp1[2,]) > 15) {
		text(x,y,procs,cex=.5)
	} else {
		text(x,y,procs,cex=.7)
	}

	c(length(procs),length(sta),nlat,nlon,ngp1,ngp)
}

runtime = function(nd)
{
   indw = grep("STEP +\\d+ +H=.+\\+CPU=",nd)
	if (length(indw) == 0) return(NULL)

   walls = as.difftime(gsub("^ *([[:digit:]:]+) .+","\\1",nd[indw]),units="secs")
   cpus = as.difftime(as.numeric(gsub(".+\\+CPU= *","",nd[indw])),units="secs")

	dwalls = diff(walls)

	# in case of change of date (time goes to 00:00)
	ind = which(dwalls < 0)
	for (i in ind) dwalls[-(1:i)] = dwalls[-(1:i)]+86400
	rt = data.frame(wall=walls,dwall=c(0,dwalls),cpu=cpus)

	i1 = grep("TIME OF START *=",nd)
	attr(rt,"start") = as.difftime(gsub("^ *TIME OF START *= *","",nd[i1]),units="secs")

	rt
}

equilon = function(ndgnh,ndlon)
{
	lats = 90*(1-(seq(ndgnh)-.5)/ndgnh)
	round(ndlon*cos(pi/180*lats))
}

plotnlon = function(nlong,nlon90,nlon45)
{
   lty = c(1,2,2)
   nlat = length(nlong)
	dlon = abs(nlong-nlon90)
	ix = which.max(dlon)
   titre = c("Nb of grid-points per Gaussian latitude",
      sprintf("nlat: %d - max diff from equi0: %d",nlat,dlon[ix]))
	matplot(data.frame(nlong,nlon90,nlon45),type="l",lty=lty,col=1,
		xlab="Latitude index",ylab="Nb of points",main=titre)
	legend("topleft",c("actual grid def","equi0","equi45"),lty=lty)
   text(c(1,nlat/2,nlat),nlong[c(1,nlat/2,nlat)],nlong[c(1,nlat/2,nlat)],pos=3,cex=.75)
   text(nlat,nlon45[nlat],nlon45[nlat],pos=3,cex=.75,offset=0.2)
	text(ix,nlong[ix],sprintf("diff: %d",dlon[ix]),pos=3,cex=.75)
	#ind = which(nlon-nlon90 != 0)
	#x = seq(along=nlon)
	#arrows(x[ind],nlon[ind],x[ind],nlon90[ind],.1,90)
}

dumpGem = function(con,gem,nsmax,nsttyp,nhtyp)
{
	cat("NFPMAX =",nsmax,"\nNFPTTYP =",nsttyp,"\nNFPHTYP =",nhtyp,"\n",file=con)
	cat("FPMUCEN =",gem$mucen,"\nFPLOCEN =",gem$locen,"\nFPSTRET =",gem$stret,
		"\nFPNLGINC =",gem$nlginc,"\n",file=con)
}

dumpLon = function(con,nloeng,nlat)
{
	# nloeng has supplementary lats at both poles
	i1 = (length(nloeng)-nlat)/2
	cat("NFPRGRI(:) =",file=con)
	write(nloeng[i1+seq(nlat%/%2)],con,ncolumns=10,sep=",")
}

dumpAB = function(con,ab)
{
	cat("NFPLEV =",dim(ab)[1]-1,"\n",file=con)
	cat("FPVALH(0:) =",file=con)
	write(ab$Ah,con,ncolumns=5,sep=",")
	cat("FPVBH(0:) =",file=con)
	write(ab$Bh,con,ncolumns=5,sep=",")
}

args = commandArgs(TRUE)
if (length(args) == 0) args = "NODE.001_01"

nd = readLines(args[1])

nproc = getvar("NPROC",nd)
nprgpns = getvar("NPRGPNS",nd)
ndglg = getvar("NDGLG",nd)
ndgnh = ndglg%/%2
ndlon = getvar("NDLON",nd)
ig = grep("NGPTOTG",nd)
ngptot = as.integer(strsplit(gsub("^ *","",nd[ig+1])," +")[[1]])
nlong = longend(nd,ndglg)

cat("Write geometry namelist for FPOS jobs\n")
gem = getgem(nd)
nsmax= getvar("NSMAX",nd)
nmsmax= getvar("NMSMAX",nd)
nsttyp = getvar("NSTTYP",nd)
nhtyp = getvar("NHTYP",nd)

nlev = getvar("NFLEVG",nd)
itropo = getvar("SUSTA: CLOSEST FULL LEVEL",nd,":")
ab = abh(nd,nlev)
std = stdatm(nd,nlev)
itropt = which(std$T == std$T[itropo-1])[1]

eta = ab$alpha+ab$Bh

con = file("fp.txt","w")
cat("&NAMFPD\nNLAT =",ndglg,"\nNLON =",ndlon,"\n/\n",file=con)
cat("&NAMFPG\n",file=con)
dumpGem(con,gem,nsmax,nsttyp,nhtyp)
dumpLon(con,nlong,ndgnh)
dumpAB(con,ab)
cat("/\n",file=con)
close(con)

nfp = readLines("fp.txt")
ind = grep("[&*/\\] *nam\\w+",nfp,invert=TRUE,ignore.case=TRUE)
nfp[ind] = sub("([.0-9]) *$","\\1,",nfp[ind])
writeLines(nfp,"fp.txt")

hasx11 = capabilities("X11")
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

if (! hasx11) png("levels.png")
tt = c("Vertical levels in hybrid coordinate",
	sprintf("%d levels - tropopause: ~%d-%d",nlev,itropt,itropo))
matplot(cbind(ab[c("Bh","alpha")],eta=eta),type="o",lty=1,pch="|",main=tt,
	xlab="Level",ylab=expression(eta))
legend("topleft",c("Bh","Ah/Pref",expression(eta)),lty=1,pch="|",col=1:3,inset=.01)
if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
invisible(dev.off())

if (! hasx11) png("stdatm.png")
par(mfrow=c(1,3))
ttstd = "Standard atmosphere"
plot(std$P/100,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Pressure"),
	xlab="Pressure (hPa)",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
plot(std$T-273.15,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Temperature (°C)"),
	xlab="Temperature",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
plot(std$rho,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Density of air"),
	xlab="Density (-)",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
invisible(dev.off())

if (length(unique(nlong)) == 1) {
	cat("--> regular Gaussian grid\n")
} else {
	stopifnot(ndglg%/%2 == (ndglg+1)%/%2)
	nlon90 = equilon(ndgnh,ndlon)
	n45 = nlong[length(nlong)%/%2]
	nlon45 = equilon(ndgnh,n45*sqrt(2))

	if (! hasx11) png("ndlon.png")
	plotnlon(nlong,nlon90,nlon45)
	if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
	invisible(dev.off())
}

cat("Grid-point mapping wrt MPI tasks\n")
s = grep("SETA=.+ LAT=.+ NSTA=",nd,value=TRUE)
sta = procmap(s)
s = grep("SETA=.+ LAT=.+ D%NONL=",nd,value=TRUE)
onl = procmap(s)

if (! hasx11) png("procmap.png")
ncomp = plotmap(sta,onl)
if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
invisible(dev.off())

ngp = min(ngptot)
ndim = c(nproc,nprgpns,ndglg,ndlon,min(ngptot),max(ngptot))
cat("Values from log file (nproc/nprgpns/ndglg/ndlon/ngptot/ngptotg):\n",ndim,"\n")
if (any(ncomp != ndim)) cat("--> computed values differ from log file:\n",ncomp,"\n")

if (any(regexpr("LSLAG *= *T(RUE)?",nd) > 0)) {
	naslb1 = getvar("\\w.+ YDSL%NASLB1",nd)
	nslrpt = getvar("SLRSET: +NSLRPT",nd)
	nslspt = getvar("SLRSET: +NSLSPT",nd)
	islwide = getvar("ISLWIDE",nd)
	cat("SL halo/total sizes:",nslspt,"/",nslrpt,"(send/recv) +",ngp,"+ pad =",naslb1,"\n")
	cat("SL halo width:",islwide,"lats and longs\n")
	cat("Core ratio in SL:",round(ngp/naslb1,3)*100,"%\n")
	isl = grep("NSLCOMM",nd,ignore.case=TRUE)
	icomm = as.integer(strsplit(gsub("^ +","",nd[isl+1])," +")[[1]])
	ip = grep("MYPROC",nd,ignore.case=TRUE)
	cat("SLCOMM for MPI task",sub(" *MYPROC += +(\\d+).*","\\1",nd[ip[1]]),":",icomm,"\n")
}

cat("Run-time information\n")
tt = runtime(nd)
if (is.null(tt)) stop("no time-steps")

nts = dim(tt)[1]

t0 = as.numeric(tt$wall[1]-attr(tt,"start"),units="secs") %% 86400
tint = as.numeric(tt$wall[nts]-tt$wall[1],units="secs") %% 86400
total = as.numeric(tt$wall[nts]-attr(tt,"start"),units="secs") %% 86400
sstt = sprintf("setup+step0, forecast, total: %gs, %gs, %gs",t0,tint,total)
cat(sstt,"\n")

if (! hasx11) png("runtime.png")
par(mfrow=c(2,1))
its = seq(nts)-1
plot(its,tt$dwall,type="h",main=c("Wall-time",sstt),xlab="Time-step",ylab="Time (s)")
plot(its,tt$cpu,type="h",main="CPU-time",xlab="Time-step",ylab="Time (s)")
if (hasx11 && interactive()) invisible(readline("Press enter to continue"))
invisible(dev.off())
