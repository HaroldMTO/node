Gnum = "-?(\\d+\\.\\d*|\\d*\\.\\d+([eE]?[-+]?\\d+)?\\>)"
Gint = "-?[0-9]+\\>"

getvar = function(var,nd,sep="=")
{
	re = sprintf("^ *\\<%s *%s *(%s|%s).*",var,sep,Gint,Gnum)
	unique(as.numeric(gsub(re,"\\1",grep(re,nd,value=TRUE))))
}

getgem = function(nd)
{
	#ig = grep("Printings in SUGEM_NAML",nd,ignore.case=TRUE)
	ig = grep("RMUCEN.+ RLOCEN.+RSTRET",nd,ignore.case=TRUE)
	mucen = as.numeric(sub(sprintf(".* RMUCEN *= *(%s).*",Gnum),"\\1",nd[ig]))
	locen = as.numeric(sub(sprintf(".* RLOCEN *= *(%s).*",Gnum),"\\1",nd[ig]))
	stret = as.numeric(sub(sprintf(".* RSTRET *= *(%s).*",Gnum),"\\1",nd[ig]))
	nlginc = as.numeric(sub(sprintf(".* RNLGINC *= *(%s).*",Gnum),"\\1",nd[ig+2]))

	list(mucen=mucen,locen=locen,stret=stret,nlginc=nlginc)
}

longend = function(nd,ndglg)
{
	ij = grep("\\( *JGL,NLOENG *\\)",nd)
	if (length(ij) == 0) {
		ij = grep("\\( *JGL,NLOENG,NMENG *\\)",nd)
		if (length(ij) == 0) return(NULL)

		is = grep("\\(JM,NDGLU\\)",nd)
		ind = seq(ij+1,is-1)
		s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+ +\\d+\\)",nd[ind])))
		nloeng = as.integer(gsub("\\( *-?\\d+ +(\\d+) +\\d+\\)","\\1",s))
	} else {
		is = grep("Set up transforms",nd)
		ind = seq(ij+1,is-1)
		s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
		nloeng = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))
	}

   stopifnot(length(nloeng)%%2 == 0)

   off = (length(nloeng)-ndglg)/2
	nloeng[off+seq(ndglg%/%2)]
}

wavend = function(nd,ndglg)
{
	ij = grep("\\( *JGL,NMENG *\\)",nd)
	if (length(ij) == 0) {
		ij = grep("\\( *JGL,NLOENG,NMENG *\\)",nd)
		is = grep("\\(JM,NDGLU\\)",nd)
		ind = seq(ij+1,is-1)
		s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+ +\\d+\\)",nd[ind])))
		nmeng = as.integer(gsub("\\( *-?\\d+ +\\d+ +(\\d+)\\)","\\1",s))
	} else {
		is = grep("\\( *JM,NDGLU *\\)",nd)
		ind = seq(ij+1,is-1)
		s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
		nmeng = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))
	}

	off = (length(nmeng)-ndglg)/2
	nmeng[off+seq(ndglg%/%2)]
}

specnb = function(nd)
{
	ij = grep("\\( *JM,NDGLU *\\)",nd)
	is = grep("Set up distributed",nd)
	ind = seq(ij+1,is-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+\\)",nd[ind])))
	ndglu = as.integer(gsub("\\( *-?\\d+ +(\\d+)\\)","\\1",s))

	ndglu
}

wavenb = function(nd)
{
	i1 = grep("^ *\\(JGL,NLOEN,NMEN\\)",nd)
	i2 = grep("^ *ARRAY +NSTAGP +ALLOCATED",nd)
	ind = seq(i1+1,i2-1)
	s = unlist(regmatches(nd[ind],gregexpr("\\( *-?\\d+ +\\d+ +\\d+\\)",nd[ind])))
	nmen = as.integer(gsub("\\( *\\d+ +\\d+ +(\\d+)\\)","\\1",s))

	nmen
}

intlines = function(nd)
{
	as.integer(unlist(regmatches(nd,gregexpr(Gint,nd))))
}

numlines = function(nd)
{
	as.numeric(unlist(regmatches(nd,gregexpr(Gnum,nd))))
}

specdis = function(nd)
{
	i1 = grep("^ *NUMPP\\>",nd)
	i2 = grep("MAXIMUM NUMBER OF THREADS",nd)
	ind = seq(i1+1,i2-1)
	numpp = intlines(nd[ind])

	i1 = grep("^ *NPROCM",nd)
	i2 = grep("^ *NFRSTLAT",nd)
	ind = seq(i1+1,i2-1)
	nprocm = intlines(nd[ind])

	i1 = grep("^ *NALLMS",nd)
	i2 = grep("^ *NPTRMS",nd)
	ind = seq(i1+1,i2-1)
	nds = gsub(" \\*{3,}"," 0",nd[ind])
	nallms = intlines(nds)

	i1 = grep("^ *MYLEVS",nd)
	i2 = grep("^ *NUMLL *$",nd)
	ind = seq(i1+1,i2-1)
	mylevs = intlines(nd[ind])

	i1 = grep("^ *NBSETLEV",nd)
	i2 = grep("^ *MYLATS",nd)
	ind = seq(i1+1,i2-1)
	nbsetlev = intlines(nd[ind])

	i1 = grep("^ *(YDLAP%)?MYMS",nd)
	i2 = grep("^ *(NASM0|YDLAP%NASN0|NALLMS) *$",nd)
	i2 = i2[i2 > i1]
	ind = seq(i1+1,i2[1]-1)
	myms = intlines(nd[ind])

	i1 = grep("^ *EIGEN-VALUES OF THE LAPLACIAN",nd)
	i2 = grep("^ *EIGEN-VALUES OF ITS INVERSE",nd)
	i3 = grep("^ *((YDLAP%)?NASM0G|YDLEP%NESM0G)",nd)
	ind = seq(i1+1,i2-1)
	rlapdi = numlines(nd[ind])

	ind = seq(i2+1,i3-1)
	rlapin = numlines(nd[ind])

	list(numpp=numpp,nprocm=nprocm,nallms=nallms,mylevs=mylevs,nbsetlev=nbsetlev,myms=myms,
		rlapdi=rlapdi,rlapin=rlapin)
}

stdatm = function(nd,nflevg)
{
	snum = "-?\\d+\\.\\d+"

	ia = grep(" +pressure +temperature +height +density",nd,ignore.case=TRUE)
	if (length(ia) == 0) {
		ia = grep("standard +atmosphere +height",nd,ignore.case=TRUE)
		ind = ia+seq(nflevg)
		ire = regexec(sprintf(" *\\d+ +(%s)",snum),nd[ind])
		Z = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
		return(data.frame(Z=Z))
	}

	ind = ia+seq(nflevg)
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s) +(%s)",snum,snum,snum,snum),nd[ind])
	P = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	T = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	Z = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	rho = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))

	data.frame(P=P,T=T,Z=Z,rho=rho)
}

abh = function(nd,nflevg)
{
	ih = grep("A and B (at half levels|on half layers)",nd)
	ind = ih+1+seq(nflevg+1)
	snum = "-?\\d+\\.\\d+"
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",snum,snum,snum),nd[ind])
	alh = as.numeric(sapply(regmatches(nd[ind],ire),"[",2))
	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",4))
	data.frame(Ah=ah,Bh=bh,alpha=alh)
}

abhfp = function(nd,nfplev)
{
	ih = grep("Set up F-post processing, vertical geometry",nd)
	ind = ih+1+seq(nfplev+1)
	snum = "-?\\d+\\.\\d+"
	ire = regexec(sprintf(" *FPVALH\\(( *\\d+|\\*)+\\) *= *(%s) +FPVBH\\(( *\\d+|\\*+)\\) *= *(%s)",
		snum,snum),nd[ind])
	ah = as.numeric(sapply(regmatches(nd[ind],ire),"[",3))
	bh = as.numeric(sapply(regmatches(nd[ind],ire),"[",5))
	data.frame(Ah=ah,Bh=bh)
}

vfe = function(nd,nflevg,oper)
{
	i1 = grep(sprintf("VFE operator %s",oper),nd)
	if (length(i1) == 0) return(NULL)

	i2 = grep("VFE operator|A and B (at half levels|on half layers)",nd)
	il = seq(i1+1,min(i2[i2 > i1])-1)
	roper = numlines(nd[il])
	stopifnot(length(roper)/nflevg == length(roper)%/%nflevg)
	t(matrix(roper,ncol=nflevg))
}

silev = function(nd,nflevg)
{
	il = grep("^( *JLEV *=)? *\\d+ +SITLAF *=",nd)
	if (length(il) > 0) {
		ire = regexec(sprintf(" *\\d+ +SITLAF *= +(%s) +SIDPHI *= +(%s)",Gnum,Gnum),nd[il])
	} else {
		il = grep("Level +SITLAF +SIDPHI",nd)
		stopifnot(length(il) == 1)
		il = seq(il+1,il+nflevg)
		ire = regexec(sprintf(" *\\d+ +(%s) +(%s)",Gnum,Gnum),nd[il])
	}

	sitlaf = as.numeric(sapply(regmatches(nd[il],ire),"[",2))
	sidphi = as.numeric(sapply(regmatches(nd[il],ire),"[",5))

	i1 = grep("(VERTICAL|Level) +.+ +EIGENVALUES",nd)
	if (length(i1) == 1) {
		il = seq(i1+1,i1+nflevg)
		ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",Gnum,Gnum,Gnum),nd[il])
		sivp = as.numeric(sapply(regmatches(nd[il],ire),"[",8))
	} else {
		sivp = rep(0,nflevg)
	}

	i1 = grep("\\<PDILEV",nd)
	if (length(i1) == 1) {
		i2 = grep("SUE?HDVPN",nd)[1]
		il = seq(i1+1,i2-1)
		pdi = numlines(nd[il])
		pdis = rep(0,nflevg)
	} else if (length(i1) > 1) {
		il = seq(i1[1]+1,i1[2]-1)
		pdi = numlines(nd[il])
		i1 = grep("PDILEVS\\>",nd)
		i2 = grep("SUE?HDVPN",nd)[1]
		il = seq(i1+1,i2-1)
		pdis = numlines(nd[il])
	} else {
		pdi = rep(0,nflevg)
		pdis = rep(0,nflevg)
	}

	i1 = grep("KNSHD *:",nd)
	i2 = grep("^ *SUHDF",nd)
	if (length(i1) == 1 && length(i2) == 1) {
		il = seq(i1+1,i2-1)
		knshd = intlines(nd[il])
	} else {
		knshd = rep(0,nflevg)
	}

	data.frame(sivp=sivp,sitlaf=sitlaf,sidphi=sidphi,pdi=pdi,pdis=pdis,knshd=knshd)
}

sicor = function(nd,nflevg)
{
	i1 = grep("SURCORDI",nd)
	if (length(i1) == 0) return(NULL)

	i2 = grep("Set up relaxation",nd,ignore.case=TRUE)
	if (length(i2) == 0) i2 = grep("NSLDIMK *=",nd,ignore.case=TRUE)
	ind = seq(i1+1,i2-1)
	ic = grep("\\<RCORDI",nd[ind])
	noms = sub("^ *(RCORDI\\w+).*","\\1",nd[ind[ic]])
	stopifnot(all(noms == sprintf("RCORDI%s",c("T","H","F"))))

	il = ind[seq(ic[1]+1,ic[2]-1)]
	rcordit = numlines(nd[il])
	il = ind[seq(ic[2]+1,ic[3]-1)]
	rcordih = numlines(nd[il])[-1]
	il = ind[-(1:ic[3])]
	rcordif = numlines(nd[il])

	data.frame(rcordit=rcordit,rcordif=rcordif,rcordih=rcordih)
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

weno = function(nd,nflevg)
{
	snum = "-?\\d+\\.\\d+"

	ind = grep("GAMMA_WENO ",nd)
	if (length(ind) == 0) return(NULL)

	indi = ind+seq(nflevg-4)
	ire = regexec(sprintf(" *\\d+ +(%s) +(%s) +(%s)",snum,snum,snum),nd[indi])
	gamma = matrix(as.numeric(sapply(regmatches(nd[indi],ire),"[",2:4)),nrow=3)
	t(gamma)
}

mesodrag = function(nd,nflevg)
{
	snum = "-?\\d+\\.\\d+"

	id = grep("PROFIL VERTICAL DE DRAG MESO",nd)
	ind = seq(id+2,id+1+nflevg)
	m = matrix(numlines(nd[ind]),nr=6)
	gwd = as.data.frame(t(m))
	names(gwd) = c("u","uday","t","tday","q","qday")

	gwd
}

varqc = function(nd)
{
	indo = grep("^ *VAR *= *\\d+",nd)
	iv = as.integer(sub("VAR *= *(\\d+).+","\\1",nd[indo]))
	qc = vector(length(unique(iv)),mode="list")
	for (i in iv) {
		ind = which(iv == i)
		notvar = intlines(sub(".+ NOTVAR *=","",nd[indo[ind[1]]]))
		v = numlines(gsub("\\*+","-9999.",nd[indo[ind[-1]]]))
		v[v==-9999] = NA
		qc[[i]] = matrix(c(v,notvar),nc=length(ind))
	}

	qc
}

jotable = function(nd)
{
	iobst = grep("Obstype +\\d+ +=+",nd,ignore.case=TRUE)
	ijog = grep("Jo Global",nd)
	ndo = nd[iobst[1]:ijog[1]]
	#ijoh = grep("Jo_Costfunction",ndo)
	ijoh = grep("Codetype +\\d+ +=+",ndo)
	njo = length(ijoh)
	ijot = grep("^ +\\w+ +\\d+( +\\d+\\.\\d+){3}",ndo)
	lj = strsplit(ndo[ijot],split=" +")
	jot = t(sapply(lj,function(x) as.numeric(x[3:6])))
	jot = as.data.frame(jot)
	jot = cbind(sapply(lj,"[",2),jot)

	names(jot) = c("Variable","DataCount","Jo_Costfunction","Jo/n","ObsErr")
	code = integer(dim(jot)[[1]])

	ijoh = c(ijoh,ijog[1])
	for (i in seq(njo)) {
		indi = ijot > ijoh[i] & ijot < ijoh[i+1]
		code[indi] = sub(" +Codetype +(\\d+) +.+","\\1",ndo[ijoh[i]])
	}

	code = as.integer(code)
	jot = cbind(Codetype=code,jot)

	jot
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
	mtext(mp1[2,],2,line=1,at=offlat+j,adj=0,las=2,cex=.7)
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
   walls = as.difftime(gsub("^ *([[:digit:]:]+) .+","\\1",nd[indw]),units="secs")
   cpus = as.difftime(as.numeric(gsub(".+\\+CPU= *","",nd[indw])),units="secs")

	dwalls = diff(walls)

	# in case of change of date (time goes to 00:00)
	ind = which(dwalls < 0)
	for (i in ind) dwalls[-(1:i)] = dwalls[-(1:i)]+86400

	# small escalating over steps within 1s
	i1 = 1
	ind = which(c(dwalls,dwalls[length(dwalls)]+1) > 0)
	for (i in ind) {
		if (i > i1) {
			n = i-i1+1
			dt = seq(0,1,length.out=n+1)[-(n+1)]
			walls[i1:i] = walls[i1]+dt
		}

		i1 = i+1
	}

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

pngalt = function(...)
{
	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))
	if (! hasx11) png(...)
}

pngoff = function(op)
{
	if (! hasx11) {
		invisible(dev.off())
	} else if (! missing(op)) {
		par(op)
	}
}

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

hasx11 = ! "png" %in% names(cargs) && capabilities("X11")
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")
ask = hasx11 && interactive()
if (! "png" %in% names(cargs)) cargs$png = "."

nd = readLines(cargs$ficin)

cat("Grid information\n")
nproc = getvar("NPROC",nd)
nprgpns = getvar(".*\\<NPRGPNS",nd)
ndglg = getvar("NDGLG",nd)
ndgnh = ndglg%/%2
ndlon = getvar("NDLON",nd)
ig = grep("NGPTOTG",nd)
ngptot = as.integer(strsplit(gsub("^ *","",nd[ig+1])," +")[[1]])

cat("Grid type and truncature\n")
gem = getgem(nd)
nsmax = getvar("NSMAX",nd)
nmsmax = getvar("NSMAX.+NMSMAX",nd)
nsttyp = getvar("NSTTYP",nd)
if (length(nsttyp) == 0) nsttyp = getvar(".+ NSTTYP",nd)
nhtyp = getvar("NHTYP",nd)
if (length(nhtyp) == 0) nhtyp = getvar(".+ NHTYP",nd)
nlong = longend(nd,ndglg)
if (is.null(nlong)) {
	cat("--> nlong guessed with NDLON/NDGLG\n")
	nlong = rep(ndlon,ndglg)
}

nflevg = getvar("NFLEVG",nd)

cat("Standard atmosphere\n")
std = stdatm(nd,nflevg)
pngalt(sprintf("%s/stdatm.png",cargs$png))
ttstd = "Standard atmosphere"
if (dim(std)[2] == 1) {
	itropo = itropt = NA_integer_
	op = par(mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	plot(seq(nflevg),std$Z,type="o",lty=1,pch="-",main=c(ttstd,"Height of level"),
		xlab="Level",ylab="Height above MSL (m)")
} else {
	itropo = getvar("SUSTA: CLOSEST FULL LEVEL",nd,":")
	itropt = which(std$T == std$T[itropo-1])[1]
	op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
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
}

pngoff(op)

ab = abh(nd,nflevg)
eta = ab$alpha+ab$Bh

cat("Vertical coordinate (1)\n")
pngalt(sprintf("%s/levels.png",cargs$png))
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
pngoff(op)

cat("Vertical coordinate (2)\n")
pngalt(sprintf("%s/eta.png",cargs$png))
tt = c("Vertical hybrid coordinate",
	sprintf("%d levels - tropopause: ~%d-%d",nflevg,itropt,itropo))
matplot(cbind(ab[c("Bh","alpha")],eta=eta),type="o",lty=1,pch="|",main=tt,
	xlab="Level",ylab=expression(eta))
legend("topleft",c("Bh","Ah/Pref",expression(eta)),lty=1,pch="|",col=1:3,inset=.01)
abline(h=0,col="darkgrey")
pngoff()

cat("Write geometry namelist for FPOS jobs\n")
con = file(sprintf("%s/fp.txt",cargs$png),"w")
cat("&NAMFPD\nNLAT =",ndglg,"\nNLON =",ndlon,"\n/\n",file=con)
cat("&NAMFPG\n",file=con)
dumpGem(con,gem,nsmax,nsttyp,nhtyp)
dumpLon(con,nlong,ndgnh)
dumpAB(con,ab)
cat("/\n",file=con)
close(con)

nfp = readLines(sprintf("%s/fp.txt",cargs$png))
ind = grep("[&*/\\] *nam\\w+",nfp,invert=TRUE,ignore.case=TRUE)
nfp[ind] = sub("([.0-9]) *$","\\1,",nfp[ind])
writeLines(nfp,sprintf("%s/fp.txt",cargs$png))
nfplev = getvar("NFPLEV",nd)
if (length(nfplev) == 1) {
	cat("FP levels:",nfplev,"\n")
	abfp = abhfp(nd,nfplev)
	pngalt(sprintf("%s/fplevels.png",cargs$png))
	op = par(mfrow=c(1,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	ss = "Layer interfaces"
	ylim = c(nfplev,0)
	plot(abfp$Ah,0:nfplev,type="o",lty=1,pch="-",main=c("Coefficient A",ss),
		xlab="Coef A",ylab="Level",ylim=ylim,cex=1.5,yaxs="i")
	plot(abfp$Bh,0:nfplev,type="o",lty=1,pch="-",main=c("Coefficient B",ss),
		xlab="Coef B",ylab="Level",ylim=ylim,cex=1.5,yaxs="i")
	pngoff(op)
}

cat("Vertical scheme\n")
rinte = vfe(nd,nflevg,"RINTE")
pngalt(sprintf("%s/vfeint.png",cargs$png))
if (is.null(rinte)) {
	plot(1,xlab="",ylab="",pch="",xaxt="n",yaxt="n")
	text(1,1,"no VFE array RINTE")
} else {
	ilev = c(1,nflevg%/%2,nflevg,nflevg+1)
	op = par(mfrow=c(1,4),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),pch="-")
	for (i in c(1,nflevg%/%2,nflevg)) {
		s = sum(rinte[,i])
		plot(rinte[,i],seq(nflevg),type="o",xlab="Rinte_l",ylab="Level",
			main=c(sprintf("VFE Sum_l=1,%d",i),sprintf("sum(Rinte): %.3g",s)),
			ylim=c(nflevg,1))
		abline(h=i,lty=2,col="grey")
		abline(v=0,col="grey")
	}

	s = sum(rinte[,nflevg+1])
	plot(rinte[,nflevg+1],seq(nflevg),type="o",xlab="Rinte_l",ylab="Level",
		main=c("Integral",sprintf("sum(Rinte*1): %.3g",s)),ylim=c(nflevg,1))
	abline(v=0,col="grey")
}
pngoff(op)

rderi = vfe(nd,nflevg,"RDERI")
pngalt(sprintf("%s/vfeder.png",cargs$png))
if (is.null(rderi)) {
	plot(1,xlab="",ylab="",pch="",xaxt="n",yaxt="n")
	text(1,1,"no VFE array RDERI")
} else {
	ilev = c(1,nflevg%/%2,nflevg,nflevg+1)
	op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),pch="-")
	for (i in c(1,nflevg%/%2,nflevg)) {
		s = sum(rderi[,i])
		plot(rderi[,i],seq(nflevg),type="o",xlab="Rderi_l",ylab="Level",
			main=c(sprintf("VFE Sum_l=1,%d",i),sprintf("sum(Rderi): %.3g",s)),
			ylim=c(nflevg,1))
		abline(v=0,col="grey")
	}
}
pngoff(op)

cat("Vertical SI system and spectral horizontal diffusion\n")
si = silev(nd,nflevg)

pngalt(sprintf("%s/sipre.png",cargs$png))
op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ylim = c(nflevg,1)
plot(si$sivp,1:nflevg,type="o",lty=1,pch="-",main="SIVP: vert. modes (= freq.)",
	xlab="Mode",ylab="Level",ylim=ylim,cex=1.5)
if (all(si$sivp == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
plot(si$sitlaf/100,1:nflevg,type="o",lty=1,pch="-",main="SITLAF: d(ln(P))/ln(P)",
	xlab="Pressure (hPa)",ylab="Level",ylim=ylim,cex=1.5)
plot(si$sidphi,1:nflevg,type="o",lty=1,pch="-",main="SIDPHI: diff. of geopotential",
	xlab="Geopotential",ylab="Level",ylim=ylim,cex=1.5)
pngoff(op)

pngalt(sprintf("%s/sihd.png",cargs$png))
op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
plot(si$pdi/100,1:nflevg,type="l",lty=1,main="PDILEV: 1+7.5*(3-log10(P))",
	xlab="PDILEV (hPa)",ylab="Level",ylim=ylim,cex=1.5)
if (all(si$pdi == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
plot(si$pdis/100,1:nflevg,type="l",lty=1,main="PDILEVS: cf PDILEV)",
	xlab="PDILEVS (hPa)",ylab="Level",ylim=ylim,cex=1.5)
plot(si$knshd,1:nflevg,type="l",lty=1,main="KNSHD",xlab="KNSHD",ylab="Level",ylim=ylim,
	cex=1.5)
if (all(si$knshd == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
pngoff(op)

cor = sicor(nd,nflevg)
pngalt(sprintf("%s/sicor.png",cargs$png))
if (is.null(cor)) {
	op = par(xaxt="n",yaxt="n")
	plot(1,xlab="",ylab="",col=0)
	text(1,1,"no CORDI")
} else {
	op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	plot(cor$rcordit,1:nflevg,type="l",lty=1,main="RCORDIT (tropo)",xlab="RCORDIT",
		ylab="Level",ylim=ylim,cex=1.5)
	if (all(cor$rcordit == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
	plot(cor$rcordih,1:nflevg,type="l",lty=1,main="RCORDIH",xlab="RCORDIH",ylab="Level",
		ylim=ylim,cex=1.5)
	if (all(cor$rcordih == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
	plot(cor$rcordif,1:nflevg,type="l",lty=1,main="RCORDIF",xlab="RCORDIF",ylab="Level",
		ylim=ylim,cex=1.5)
	if (all(cor$rcordif == 0)) text(0,nflevg/2,"not decoded",.5,col=2)
}
pngoff(op)

cat("Spectral and vertical dimensions\n")
nprtrw = getvar("NPRTRW",nd)
if (length(nprtrw) == 0) {
	nprtrw = getvar(".*\\<NPRTRW",nd)
	if (length(nprtrw) == 0) {
		i1 = grep("^.+ NPRTRW *= *$",nd)
		nprtrw = getvar(".*\\<NPRTRW",paste(nd[i1+0:1],collapse=""))
	}
	nprtrn = getvar(".*\\<NPRTRN",nd)
	nprtrns = getvar(".*\\<NPRTRNS",nd)
	nprtrv = getvar(".*\\<NPRTRV",nd)
} else {
	nprtrn = getvar("NPRTRN",nd)
	nprtrns = getvar("NPRTRNS",nd)
	nprtrv = getvar("NPRTRV",nd)
}

cat("Spectral partitionning\n")
if (getvar("NPRINTLEV",nd) > 0) {
	nmeng = try(wavend(nd,ndglg))
	ndglu = try(specnb(nd))
	#nmen = wavenb(nd)

	nm = length(ndglu)-1

	pngalt(sprintf("%s/specgp.png",cargs$png))
	op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	ss = sprintf("nmeng: %d... %d",min(nmeng),max(nmeng))
	plot(nmeng,type="l",main=c("Wave cut-off per latitude",ss),xlab="Latitude index",
		ylab="Nb of waves",xaxt="n")
	axis(1,pretty(seq(along=nmeng)/8,8)*8)
	ss = sprintf("ndglu: %d... %d",min(ndglu),max(ndglu))
	plot(0:nm,ndglu,type="l",main=c("Nb of longitudes per wave",ss),
		xlab="Wave index 'jm'",ylab="Nb of longitudes",xaxt="n")
	axis(1,pretty((seq(along=ndglu)-1)/8,8)*8)
	pngoff(op)
}

cat("Spectral partitionning, part 2\n")
if (getvar("NPRINTLEV",nd) > 0) {
	sp = specdis(nd)
	pngalt(sprintf("%s/specproc.png",cargs$png))
	op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	pr = unlist(lapply(seq(along=sp$numpp),function(i) rep(i,each=sp$numpp[i])))
	np = length(sp$nprocm)-1
	plot(0:np,sp$nprocm,type="h",main=c("W-set and waves","'nprocm'"),xlab="Wave index 'jm'",
		ylab="W-set index")
	if (any(sp$nprocm == 0)) text(np/2,max(sp$nprocm)/2,"NPROCM not fully decoded",.5,col=2)
	plot(sp$nallms,type="h",main=c("Waves and W-set","'nallms'"),
		xlab="W-set index (nb of waves 'numpp')",ylab="Wave index 'jm'",lwd=1.5,
		col=(pr-1)%%8+1,xaxt="n")
	x = cumsum(sp$numpp)
	axis(1,c(x[1]/2,(x[-1]+x[-length(x)])/2),sprintf("%d (%d)",seq(sp$numpp),sp$numpp))
	if (any(duplicated(sp$nallms))) {
		stopifnot(all(! duplicated(sp$nallms[sp$nallms!=0])))
		text(nm/2,nm/2,"NALLMS not fully decoded",.5,col=2)
	}
	pngoff(op)

	pngalt(sprintf("%s/vset.png",cargs$png))
	plot(sp$nbsetlev,1:nflevg,type="p",ylim=ylim,main="V-set and levels",xlab="V-set",
		ylab="Level",pch="-")
	pngoff()

	pngalt(sprintf("%s/speclap.png",cargs$png))
	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	plot(sp$rlapdi,type="l",main=c("Eigen-values of the Laplacian","'rlapdi'"),
		xlab="Wave index 'jm'",ylab="Eigen-value")
	plot(sp$rlapin,type="l",main=c("Eigen-values of inverse of Laplacian","'rlapin'"),
		xlab="Wave index 'jm'",ylab="Eigen-value")
	plot(sp$rlapin,type="l",xlim=c(1,min(ndglg,20)),main="First Eigen-values 'rlapin'",
		xlab="Wave index 'jm'",ylab="Eigen-value")
	pngoff(op)
}

cat("Vertical cubic weights (SL)\n")
vintw = cuico(nd,nflevg)

pngalt(sprintf("%s/vintw.png",cargs$png))
if (all(is.na(vintw))) {
	plot(1,xaxt="n",yaxt="n",xlab="",ylab="",pch="")
	text(1,1,"no vertical cubic weights (VINTW)")
} else {
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
}
pngoff(op)

cat("Vertical WENO coefficients\n")
gamma = weno(nd,nflevg)
pngalt(sprintf("%s/weno.png",cargs$png))
if (is.null(gamma)) {
	plot(1,xaxt="n",yaxt="n",xlab="",ylab="",pch="")
	text(1,1,"no gamma weights (WENO)")
} else {
	op = par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

	for (i in 1:3) {
		tt = sprintf("WENO weights %d",i)
		plot(gamma[,i],2:nl3,type="o",lty=1,pch="-",main=tt,xlab="Gamma",ylab="Level",ylim=c(nl3,2))
	}
}
pngoff(op)

if (length(unique(nlong)) == 1) {
	cat("--> regular Gaussian grid\n")
} else {
	stopifnot(ndglg%%2 == 0)
	nlon90 = equilon(ndgnh,ndlon,1)
	n45 = nlong[length(nlong)%/%2]
	nlon45 = equilon(ndgnh,n45*sqrt(2),1)

	pngalt(sprintf("%s/ndlon.png",cargs$png))
	cat(par("mfrow"),"\n")
	op = par(mfrow=c(1,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	plotnlon(nlong,nlon90,nlon45)
	pngoff(op)
}

cat("Grid-point mapping wrt MPI tasks\n")
ndim = c(nproc,nprgpns,ndglg,ndlon,min(ngptot),max(ngptot))
cat("Values from log file (nproc/nprgpns/ndglg/ndlon/ngptot/ngptotg):\n",ndim,"\n")

s = grep("SETA=.+ LAT=.+ NSTA=",nd,value=TRUE)
pngalt(sprintf("%s/procmap.png",cargs$png))
if (length(s) == 0) {
	plot(1,xaxt="n",yaxt="n",xlab="",ylab="",pch="")
	text(1,1,"no map of MPI tasks")
} else {
	sta = procmap(s)
	s = grep("SETA=.+ LAT=.+ (D%)?NONL=",nd,value=TRUE)
	onl = procmap(s)

	ncomp = plotmap(sta,onl)

	if (any(ncomp != ndim)) cat("--> computed values differ from log file:\n",ncomp,"\n")
}
pngoff()

if (any(regexpr("LSLAG *= *T(RUE)?",nd) > 0)) {
	ngp = min(ngptot)
	cat("SL scheme information\n")
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

cat("Vertical mesoscale drag\n")
gwd = mesodrag(nd,nflevg)
pngalt(sprintf("%s/mesodrag.png",cargs$png))
op = par(mfrow=c(1,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
ylim = c(nflevg,1)
plot(gwd$u,1:nflevg,type="o",lty=1,pch="-",main="Mesoscale drag",
	xlab="Wind speed",ylab="Level",ylim=ylim,cex=1.5)
abline(v=0,col="grey")
plot(gwd$t,1:nflevg,type="o",lty=1,pch="-",main="Mesoscale drag",
	xlab="Temperature",ylab="Level",ylim=ylim,cex=1.5)
abline(v=0,col="grey")
pngoff()

cat("Values of Var QC\n")
if (! is.null(cargs$varqc) && length(grep("OBSERVATION TYPE:",nd)) > 0) {
	qc = varqc(nd)
	qc1 = qc[[1]]
	types = seq(dim(qc1)[1])
	pngalt(sprintf("%s/varqc.png",cargs$png))
	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	barplot(qc1[,1],names.arg=types,main="",xlab="Obs type",ylab="RAQC")
	barplot(qc1[,2],names.arg=types,main="",xlab="Obs type",ylab="RLQC")
	#barplot(t(qc1[,3:5]),names.arg=types,beside=TRUE,main="",xlab="Obs type",ylab="RBGQC(1:3)")
	matplot(qc1[,3:5],main="",xlab="Obs type",ylab="RBGQC(1:3)",type="o",lty=1,col=1:3,pch=20)
	barplot(qc1[,6],names.arg=types,main="",xlab="Obs type",ylab="NOTVAR")
	pngoff()
}

cat("Jo tables\n")
if (length(grep("JOT-sname",nd)) > 0) {
	jot = jotable(nd)

	jog = by(jot[-(1:2)],jot$Codetype,colMeans)
	jod = jog[[1]]
	for (jo in jog[-1]) jod = rbind(jod,jo,deparse.level=0)
	jod = cbind(data.frame(Codetype=dimnames(jog)[[1]]),jod)

	cat("Means of Jo by codetype:\n")
	print(jod)

	jog = by(jot[-(1:2)],jot$Variable,colMeans)
	jod = jog[[1]]
	for (jo in jog[-1]) jod = rbind(jod,jo,deparse.level=0)
	jod = cbind(data.frame(Variable=dimnames(jog)[[1]]),jod)

	cat("Means of Jo by variable:\n")
	print(jod)
}

cat("Run-time information\n")
nstop = getvar("NSTOP",nd)
tstep = getvar("TSTEP",nd)
if (length(grep("STEP +\\d+ +H=.+\\+CPU=",nd)) > 0) {
	tt = runtime(nd)
	nts = dim(tt)[1]

	t0 = round(as.numeric(tt$wall[1]-attr(tt,"start"),units="secs") %% 86400,3)
	tint = round(as.numeric(tt$wall[nts]-tt$wall[1],units="secs") %% 86400,3)
	total = round(as.numeric(tt$wall[nts]-attr(tt,"start"),units="secs") %% 86400,3)
	sstt = sprintf("setup+step0, forecast, total: %gs, %gs, %gs",t0,tint,total)
	cat(sstt,"\n")

	pas = c(1,2,3,5,6,10,15,20,40,50)
	nt = c(0,200,400,600,1200,2000,3000,4000,8000,10000)
	it = findInterval(nts,nt)

	pngalt(sprintf("%s/runtime.png",cargs$png))
	op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	its = seq(1,nts,by=pas[it])
	dwall = c(0,diff(tt$wall))
	plot(its-1,dwall[its],type="h",main=c("Wall-time",sstt),xlab="Time-step",
		ylab="Time (s)")
	plot(its-1,tt$cpu[its],type="h",main="CPU-time",xlab="Time-step",ylab="Time (s)")
	pngoff(op)
}

if (nsttyp == 2) {
	cat("--> tilted grid\n")
	library(maps)
	xp = 180/pi*gem$locen
	yp = 180/pi*asin(gem$mucen)
	xlim = xp+c(-40,40)
	ylim = yp+c(-20,20)
	pngalt(sprintf("%s/pole.png",cargs$png))
	map("world",xlim=xlim,ylim=ylim)
	points(xp,yp,pch="+",col="red")
	text(xp,yp,sprintf("pole (lat/long): %.3g %.3g",yp,xp),pos=3,col="red")
	pngoff()
}
