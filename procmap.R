Gnum = "-?(\\d+\\.\\d*|\\d*\\.\\d+([eE]?[-+]?\\d+)?\\>)"
Gint = "-?[0-9]+\\>"

getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(%s|%s).*",var,sep,Gint,Gnum)
	unique(as.numeric(gsub(re,"\\1",grep(re,nd,value=TRUE))))
}

longend = function(nd,ndglg)
{
	ij = grep("\\( *JGL,NLOENG *\\)",nd)
	is = grep("Set up transforms",nd)

	ind = seq(ij+1,is-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
	nloeng = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))
   stopifnot(length(nloeng)%%2 == 0)

   off = (length(nloeng)-ndglg)/2
	nloeng[off+seq(ndglg%/%2)]
}

intlines = function(nd)
{
	as.integer(unlist(regmatches(nd,gregexpr(Gint,nd))))
}

numlines = function(nd)
{
	as.numeric(unlist(regmatches(nd,gregexpr(Gnum,nd))))
}

spec = function(nd,ndglg)
{
	i1 = grep("^ *NUMPP\\>",nd)
	i2 = grep("NUMBER OF THREADS",nd)
	ind = seq(i1+1,i2-1)
	numpp = intlines(nd[ind])

	i1 = grep("^ *NPROCM",nd)
	i2 = grep("^ *NFRSTLAT",nd)
	ind = seq(i1+1,i2-1)
	nprocm = intlines(nd[ind])

	i1 = grep("^ *NALLMS",nd)
	i2 = grep("^ *NPTRMS",nd)
	ind = seq(i1+1,i2-1)
	nallms = intlines(nd[ind])

	i1 = grep("^ *MYLEVS",nd)
	i2 = grep("^ *NUMLL *$",nd)
	ind = seq(i1+1,i2-1)
	mylevs = intlines(nd[ind])

	i1 = grep("^ *NBSETLEV",nd)
	i2 = grep("^ *MYLATS",nd)
	ind = seq(i1+1,i2-1)
	nbsetlev = intlines(nd[ind])

	i1 = grep("^ *YDLAP%MYMS",nd)
	i2 = grep("^ *NASM0 *$",nd)
	ind = seq(i1+1,i2-1)
	myms = intlines(nd[ind])

	ij = grep("\\( *JGL,NMENG *\\)",nd)
	is = grep("\\( *JM,NDGLU *\\)",nd)
	ind = seq(ij+1,is-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
	nmeng = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))
	off = (length(nmeng)-ndglg)/2
	nmeng = nmeng[off+seq(ndglg%/%2)]

	ij = grep("\\( *JM,NDGLU *\\)",nd)
	is = grep("Set up distributed",nd)
	ind = seq(ij+1,is-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
	ndglu = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))

	i1 = grep("^ *\\(JGL,NLOEN,NMEN\\)",nd)
	i2 = grep("^ *ARRAY +NSTAGP +ALLOCATED",nd)
	ind = seq(i1+1,i2-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+ +\\d+\\)",nd[ind])))
	nmen = as.integer(gsub("\\( *\\d+ +\\d+ +(\\d+)\\)","\\1",s))

	i1 = grep("^ *EIGEN-VALUES OF THE LAPLACIAN",nd)
	i2 = grep("^ *EIGEN-VALUES OF ITS INVERSE",nd)
	i3 = grep("^ *YDLAP%NASM0G",nd)
	ind = seq(i1+1,i2-1)
	rlapdi = numlines(nd[ind])

	ind = seq(i2+1,i3-1)
	rlapin = numlines(nd[ind])

	list(numpp=numpp,nprocm=nprocm,nallms=nallms,mylevs=mylevs,nbsetlev=nbsetlev,myms=myms,
		nmeng=nmeng,ndglu=ndglu,nmen=nmen,rlapdi=rlapdi,rlapin=rlapin)
}

stdatm = function(nd,nflevg)
{
	ia = grep(" +pressure +temperature +height +density",nd,ignore.case=TRUE)
	ind = ia+seq(nflevg)
	snum = "-?\\d+\\.\\d+"
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s) +(%s)",snum,snum,snum,snum),nd[ind])
	P = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	T = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	Z = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	rho = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))

	data.frame(P=P,T=T,Z=Z,rho=rho)
}

abh = function(nd,nflevg)
{
	ih = grep("A and B at half levels",nd)
	ind = ih+1+seq(nflevg+1)
	snum = "-?\\d+\\.\\d+"
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",snum,snum,snum),nd[ind])
	alh = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	data.frame(Ah=ah,Bh=bh,alpha=alh)
}

silev = function(nd,nflevg)
{
	il = grep("JLEV *= \\d+ +SITLAF *=",nd)
	if (length(il) > 0) {
		ire = regexec(sprintf(" *\\d+ +SITLAF += +(%s) +SIDPHI += +(%s)",Gnum,Gnum),nd[il])
	} else {
		il = grep("Level +SITLAF +SIDPHI",nd)
		stopifnot(length(il) == 1)
		il = seq(il+1,il+nflevg)
		ire = regexec(sprintf(" *\\d+ +(%s) +(%s)",Gnum,Gnum),nd[il])
	}

	sitlaf = as.numeric(sapply(regmatches(nd[il],ire),"[",2))
	sidphi = as.numeric(sapply(regmatches(nd[il],ire),"[",5))

	il = grep("(VERTICAL|Level) +.+ +EIGENVALUES",nd)
	il = seq(il+1,il+nflevg)
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",Gnum,Gnum,Gnum),nd[il])
	sivp = as.numeric(sapply(regmatches(nd[il],ire),"[",8))

	i1 = grep("KNSHD *:",nd)
	i2 = grep("^ *SUHDF",nd)
	il = seq(i1+1,i2-1)
	knshd = intlines(nd[il])

	i1 = grep("PDILEV",nd)
	i2 = grep("SUHDVPN",nd)
	il = seq(i1+1,i2-1)
	pdi = numlines(nd[il])

	i1 = grep("SURCORDI",nd)
	i2 = grep("Set up relaxation",nd,ignore.case=TRUE)
	ind = seq(i1+1,i2-1)
	ic = grep("RCORDI",nd[ind])
	noms = sub("^ *(RCORDI\\w+).*","\\1",nd[ind[ic]])
	stopifnot(all(noms == sprintf("RCORDI%s",c("T","H","F"))))

	il = ind[seq(ic[1]+1,ic[2]-1)]
	rcordit = numlines(nd[il])
	il = ind[seq(ic[2]+1,ic[3]-1)]
	rcordih = numlines(nd[il])[-1]
	il = ind[-(1:ic[3])]
	rcordif = numlines(nd[il])

	data.frame(sivp=sivp,sitlaf=sitlaf,sidphi=sidphi,pdi=pdi,knshd=knshd,rcordit=rcordit,
		rcordif=rcordif,rcordih=rcordih)
}

cuico = function(nd,nflevg)
{
	snum = "-?\\d+\\.\\d+"

	ind = grep("JLEV +VCUICO\\>(\\(1)?",nd)
	if (length(ind) == 1) {
		indi = ind+seq(nflevg-3)
		ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s) +(%s)",snum,snum,snum,snum),nd[indi])
		vintw = matrix(as.numeric(sapply(regmatches(nd[indi],ire),"[",2:5)),nrow=4)
		vintw = t(vintw)
	} else {
		vintw = matrix(nrow=nflevg-3,ncol=4)
		for (i in seq(along=ind)) {
			indi = ind[i]+seq(nflevg-3)
			ire = regexec(sprintf(" *\\d+ +(%s)",snum),nd[indi])
			vintw[,i] = as.numeric(sapply(regmatches(nd[indi],ire),"[",2))
		}
	}

	vintw
}

getgem = function(nd)
{
	ig = grep("Printings in SUGEM_NAML",nd,ignore.case=TRUE)
	mucen = as.numeric(sub(sprintf(".* RMUCEN *= *(%s).*",Gnum),"\\1",nd[ig+1]))
	locen = as.numeric(sub(sprintf(".* RLOCEN *= *(%s).*",Gnum),"\\1",nd[ig+1]))
	stret = as.numeric(sub(sprintf(".* RSTRET *= *(%s).*",Gnum),"\\1",nd[ig+1]))
	nlginc = as.numeric(sub(sprintf(".* RNLGINC *= *(%s).*",Gnum),"\\1",nd[ig+3]))

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
		cols = c("black","brown3","deepskyblue","turquoise","hotpink2","olivedrab1",
			"darkgrey","wheat")
		cols = list(1:8,cols)

		for (ia in seq(along=sta)[-1]) {
			lmap = latmap(sta[[ia]],end[[ia]],offlat[ia],nlat)
			image(lmap,add=TRUE,col=cols[[1+(ia-1)%%2]][2+(2*(ia-1)+seq(mp1[2,ia])-1)%%7])
		}
	}

	j = mp1[1,]%/%2
	#rug(offlat,.01,side=4)
	mtext(mp1[2,],2,at=offlat+j,adj=1.2,las=2,cex=.7)
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

equilon = function(ndgnh,ndlon,nadd=0)
{
	lats = 90*(1-(seq(ndgnh)-.5)/ndgnh)
	if (nadd > 0) {
		nlat = ndgnh+nadd
		lats2 = 90*(1-(seq(nlat)-.5)/nlat)
		lats = lats2[-(1:nadd)]
	}

	round(ndlon*cos(pi/180*lats))
}

plotnlon = function(nlong,nlon90,nlon45)
{
   nlat = length(nlong)
   titre = c("Nb of grid-points per Gaussian latitude",
      sprintf("actual grid + equilong curves - #lat: %d",nlat))
	matplot(data.frame(nlong,nlon90,nlon45),type="l",lty=1:3,col=1:3,
		xlab="Latitude index",ylab="Nb of points",main=titre)
	tt = paste(c("actual grid","equi0","equi45"),c(sum(nlong),sum(nlon90),sum(nlon45)))
	legend("topleft",tt,lty=1:3,col=1:3)
	il = c(1,nlat)
   text(il,nlong[il],nlong[il],pos=3,cex=.8)
   text(il,nlon45[il],nlon45[il],pos=3,cex=.8,offset=0.2,col=3)
   text(1,nlon90[1],nlon90[1],pos=3,cex=.8,offset=0.2,col=2)
	dlon = abs(nlong-nlon90)
	ix = which.max(dlon)
	text(ix,nlong[ix],sprintf("diff max: %+d",dlon[ix]),pos=3,cex=.9)
	arrows(ix,nlong[ix],ix,nlon90[ix],.05,90,3)
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

hasx11 = length(grep("png",args)) == 0 && capabilities("X11")
ask = hasx11 && interactive()
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")

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

nflevg = getvar("NFLEVG",nd)
ab = abh(nd,nflevg)

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

std = stdatm(nd,nflevg)
itropo = getvar("SUSTA: CLOSEST FULL LEVEL",nd,":")
itropt = which(std$T == std$T[itropo-1])[1]

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("eta.png")
tt = c("Vertical hybrid coordinate",
	sprintf("%d levels - tropopause: ~%d-%d",nflevg,itropt,itropo))
matplot(cbind(ab[c("Bh","alpha")],eta=eta),type="o",lty=1,pch="|",main=tt,
	xlab="Level",ylab=expression(eta))
legend("topleft",c("Bh","Ah/Pref",expression(eta)),lty=1,pch="|",col=1:3,inset=.01)
abline(h=0,col="darkgrey")
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("levels.png")
op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ss = "Layer interfaces"
ylim = c(nflevg,0)
plot(eta,0:nflevg,type="o",lty=1,pch="-",main=c("Hybrid coordinate",ss),
	xlab=expression(eta),ylab="Level",ylim=ylim,cex=1.5,yaxs="i")
abline(h=c(itropt,itropo),lty=2)
plot(ab$Bh,0:nflevg,type="o",lty=1,pch="-",main=c("Coefficient B",ss),
	xlab="Coef B",ylab="Level",ylim=ylim,cex=1.5,yaxs="i")
abline(h=c(itropt,itropo),lty=2)
plot(ab$Ah,0:nflevg,type="o",lty=1,pch="-",main=c("Coefficient A/Pref",ss),
	xlab="alpha (=A/Pref)",ylab="Level",ylim=ylim,cex=1.5,yaxs="i")
abline(h=c(itropt,itropo),lty=2)
par(op)
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("stdatm.png")
op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ttstd = "Standard atmosphere"
plot(std$P/100,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Pressure"),
	xlab="Pressure (hPa)",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
plot(std$T-273.15,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Temperature (°C)"),
	xlab="Temperature",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
abline(v=0,lty=1)
plot(std$rho,std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Density of air"),
	xlab="Density (-)",ylab="Z (mgp)",cex=1.5)
abline(h=c(0,std$Z[c(itropt,itropo)]),lty=2)
par(op)
if (! hasx11) invisible(dev.off())

cat("Vertical SI system and spectral horizontal diffusion\n")
si = silev(nd,nflevg)

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("sipre.png")
op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ylim = c(nflevg,1)
plot(si$sivp,1:nflevg,type="o",lty=1,pch="-",main="SIVP: vert. modes (= freq.)",
	xlab="Mode",ylab="Level",ylim=ylim,cex=1.5,log="x")
plot(si$sitlaf/100,1:nflevg,type="o",lty=1,pch="-",main="SITLAF: d(ln(P))/ln(P)",
	xlab="Pressure (hPa)",ylab="Level",ylim=ylim,cex=1.5)
plot(si$sidphi,1:nflevg,type="o",lty=1,pch="-",main="SIDPHI: diff. of geopotential",
	xlab="Geopotential",ylab="Level",ylim=ylim,cex=1.5)
par(op)
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("sihd.png")
op = par(mfrow=c(2,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
plot(si$pdi/100,1:nflevg,type="l",lty=1,main="PDILEV: 1+7.5*(3-log10(P))",
	xlab="PDILEV (hPa)",ylab="Level",ylim=ylim,cex=1.5)
plot(si$knshd,1:nflevg,type="l",lty=1,main="KNSHD",xlab="KNSHD",ylab="Level",ylim=ylim,
	cex=1.5)
plot(si$rcordit,1:nflevg,type="l",lty=1,main="RCORDIT (tropo)",xlab="RCORDIT",
	ylab="Level",ylim=ylim,cex=1.5)
plot(si$rcordih,1:nflevg,type="l",lty=1,main="RCORDIH",xlab="RCORDIH",ylab="Level",
	ylim=ylim,cex=1.5)
plot(si$rcordif,1:nflevg,type="l",lty=1,main="RCORDIF",xlab="RCORDIF",ylab="Level",
	ylim=ylim,cex=1.5)
par(op)
if (! hasx11) invisible(dev.off())

cat("Spectral and vertical partitionning\n")
nprtrw = getvar("NPRTRW",nd)
nprtrn = getvar("NPRTRN",nd)
nprtrns = getvar("NPRTRNS",nd)
nprtrv = getvar("NPRTRV",nd)

sp = spec(nd,ndglg)
nm = length(sp$ndglu)-1

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("specgp.png")
op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ss = sprintf("nmeng: %d... %d",min(sp$nmeng),max(sp$nmeng))
plot(sp$nmeng,type="l",main=c("Wave cut-off per latitude",ss),xlab="Latitude index",
	ylab="Nb of waves",xaxt="n")
axis(1,pretty(seq(along=sp$nmeng)/8,8)*8)
ss = sprintf("ndglu: %d... %d",min(sp$ndglu),max(sp$ndglu))
plot(0:nm,sp$ndglu,type="l",main=c("Nb of longitudes per wave",ss),
	xlab="Wave index 'jm'",ylab="Nb of longitudes",xaxt="n")
axis(1,pretty((seq(along=sp$ndglu)-1)/8,8)*8)
par(op)
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("specproc.png")
op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
pr = unlist(lapply(seq(along=sp$numpp),function(i) rep(i,each=sp$numpp[i])))
plot(0:nm,sp$nprocm,type="h",main=c("W-set and waves","'nprocm'"),xlab="Wave index 'jm'",
	ylab="W-set index")
plot(sp$nallms,type="h",main=c("Waves and W-set","'nallms'"),
	xlab="W-set index (nb of waves 'numpp')",ylab="Wave index 'jm'",lwd=1.5,
	col=(pr-1)%%8+1,xaxt="n")
x = cumsum(sp$numpp)
axis(1,c(x[1]/2,(x[-1]+x[-length(x)])/2),sprintf("%d (%d)",seq(sp$numpp),sp$numpp))
par(op)
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("vset.png")
plot(sp$nbsetlev,1:nflevg,type="p",ylim=ylim,main="V-set and levels",xlab="V-set",
	ylab="Level",pch="-")
if (! hasx11) invisible(dev.off())

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("speclap.png")
op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
plot(sp$rlapdi,type="l",main=c("Eigen-values of the Laplacian","'rlapdi'"),
	xlab="Wave index 'jm'",ylab="Eigen-value")
plot(sp$rlapin,type="l",main=c("Eigen-values of inverse of Laplacian","'rlapin'"),
	xlab="Wave index 'jm'",ylab="Eigen-value")
plot(sp$rlapin,type="l",xlim=c(1,min(ndglg,20)),main="First Eigen-values 'rlapin'",
	xlab="Wave index 'jm'",ylab="Eigen-value")
par(op)
if (! hasx11) invisible(dev.off())

cat("Vertical cubic weights (SL)\n")
vintw = cuico(nd,nflevg)

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("vintw.png")
op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
nl3 = nflevg-3
ntop = min(nl3,5)
nmid = max(ntop,1+nl3/2)
matplot(abs(vintw[,2:4]),1:nl3,type="o",lty=1,pch="-",
	main="Cubic weight for vert. interp.",xlab="abs(Weight)",ylab="Level",ylim=c(nl3,1),
	log="x")
legend("bottomright",sprintf("w(,%d)",2:4),lty=1,col=1:3)
matplot(abs(vintw[1:ntop,2:4]),1:ntop,type="o",lty=1,pch="-",
	main="Weight at top levels",xlab="abs(Weight)",ylab="Level",ylim=c(ntop,1))
matplot(abs(vintw[ntop:nmid,2:4]),ntop:nmid,type="o",lty=1,pch="-",
	main="Weight at mid levels",xlab="abs(Weight)",ylab="Level",ylim=c(nmid,ntop))
matplot(abs(vintw[nmid:nl3,2:4]),nmid:nl3,type="o",lty=1,pch="-",
	main="Weight at bottom",xlab="abs(Weight)",ylab="Level",ylim=c(nl3,nmid))
par(op)
if (! hasx11) invisible(dev.off())

if (length(unique(nlong)) == 1) {
	cat("--> regular Gaussian grid\n")
} else {
	stopifnot(ndglg%%2 == 0)
	nlon90 = equilon(ndgnh,ndlon,1)
	n45 = nlong[length(nlong)%/%2]
	nlon45 = equilon(ndgnh,n45*sqrt(2),1)

	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
	if (! hasx11) png("ndlon.png")
	cat(par("mfrow"),"\n")
	par(mfrow=c(1,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	plotnlon(nlong,nlon90,nlon45)
	if (! hasx11) invisible(dev.off())
}

cat("Grid-point mapping wrt MPI tasks\n")
s = grep("SETA=.+ LAT=.+ NSTA=",nd,value=TRUE)
sta = procmap(s)
s = grep("SETA=.+ LAT=.+ D%NONL=",nd,value=TRUE)
onl = procmap(s)

if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
if (! hasx11) png("procmap.png")
ncomp = plotmap(sta,onl)

ngp = min(ngptot)
ndim = c(nproc,nprgpns,ndglg,ndlon,min(ngptot),max(ngptot))
cat("Values from log file (nproc/nprgpns/ndglg/ndlon/ngptot/ngptotg):\n",ndim,"\n")
if (any(ncomp != ndim)) cat("--> computed values differ from log file:\n",ncomp,"\n")

if (any(regexpr("LSLAG *= *T(RUE)?",nd) > 0)) {
	naslb1 = getvar("\\w.+ YDSL%NASLB1",nd)
	nslrpt = getvar("SLRSET: +NSLRPT",nd)
	nslspt = getvar("SLRSET: +NSLSPT",nd)
	islwide = getvar("ISLWIDE",nd)
	cat("SL halo comms (send/recv):",nslspt,"/",nslrpt,"points\n")
	cat("SL total size (naslb1):",naslb1,"points\n")
	cat("SL halo width:",islwide,"lats and longs\n")
	cat("Core ratio in SL:",round(ngp/naslb1,3)*100,"%\n")
	isl = grep("NSLCOMM",nd,ignore.case=TRUE)
	icomm = as.integer(strsplit(gsub("^ +","",nd[isl+1])," +")[[1]])
	ip = grep("MYPROC",nd,ignore.case=TRUE)
	cat("SL comms for MPI task",sub(" *MYPROC += +(\\d+).*","\\1",nd[ip[1]]),":",icomm,"\n")
}

cat("Run-time information\n")
nstop = getvar("NSTOP",nd)
tstep = getvar("TSTEP",nd)
tt = runtime(nd)
if (! is.null(tt)) {
	nts = dim(tt)[1]

	t0 = as.numeric(tt$wall[1]-attr(tt,"start"),units="secs") %% 86400
	tint = as.numeric(tt$wall[nts]-tt$wall[1],units="secs") %% 86400
	total = as.numeric(tt$wall[nts]-attr(tt,"start"),units="secs") %% 86400
	sstt = sprintf("setup+step0, forecast, total: %gs, %gs, %gs",t0,tint,total)
	cat(sstt,"\n")

	if (! hasx11) png("runtime.png")
	pas = c(1,2,5,10,20,40,50)
	nt = c(0,100,200,500,1000,2000)
	it = findInterval(nts,nt)

	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
	op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	its = seq(1,nts,by=pas[it])
	plot(its-1,tt$dwall[its],type="h",main=c("Wall-time",sstt),xlab="Time-step",
		ylab="Time (s)")
	plot(its-1,tt$cpu[its],type="h",main="CPU-time",xlab="Time-step",ylab="Time (s)")
	par(op)
	if (! hasx11) invisible(dev.off())
}

if (nsttyp == 2) {
	cat("--> tilted grid\n")
	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
	if (! hasx11) png("pole.png")
	library(maps)
	xp = 180/pi*gem$locen
	yp = 180/pi*asin(gem$mucen)
	xlim = xp+c(-40,40)
	ylim = yp+c(-20,20)
	map("world",xlim=xlim,ylim=ylim)
	points(xp,yp,pch="+",col="red")
	text(xp,yp,sprintf("pole (lat/long): %.3g %.3g",yp,xp),pos=3,col="red")
	if (! hasx11) invisible(dev.off())
}
