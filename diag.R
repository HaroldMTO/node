#R --args fic=ICMSHARPE+0000 date=20210329/00 path=/scratch/work/petithommeh/oper/arpege/4dvarfr/OPER
#R --args fic=PFARPE0000+2400 date=20210329/00 path=/scratch/work/petithommeh/oper/arpege/4dvarfr/OPER
# fic=pfarpe.txt level=levels.txt params=uv.txt
diags = "~/util/diags"

library(maps)

Gvp0 = 101325

getFrame = function(file)
{
	system(sprintf("epy_dump.py %s -f frame -o %s",file,Gficbin))
	con = file(Gficbin,"rb")
	dims = readBin(con,"integer",size=8,n=5)
	if (dims[5] > 0) stop("LAM grid not supported")

	nlong = readBin(con,"integer",size=8,n=dims[1])

	if (dims[2] > 0) nwave = readBin(con,"integer",size=8,n=dims[2])

	base = readBin(con,"numeric",n=1)
	base = as.POSIXct(base,origin="1970-01-01")
	step = readBin(con,"integer",size=8,n=1)
	if (regexpr(".+\\+0*([0-9]+)",file) > 0) {
		ech = as.integer(gsub(".+\\+0*([0-9]+)","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	lats = readBin(con,"numeric",n=dims[4])
	longs = readBin(con,"numeric",n=dims[4])

	Ai = readBin(con,"numeric",n=dims[3])
	Bi = readBin(con,"numeric",n=dims[3])

	close(con)

	nlat = length(nlong)
	theta = 90-180*(seq(nlat)-.5)/nlat

	frame = list(nlat=nlat,nwave=dims[2],nlevel=dims[3],npdg=dims[4],nlong=nlong,
		theta=theta,lat=lats,long=longs,A=Ai,B=Bi,base=base,step=step)

	frame
}

getGPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin))
	con = file(Gficbin,"rb")
	nf = readBin(con,"integer",size=8,n=1)
	stopifnot(length(nf) == 1)

	data = matrix(nrow=frame$npdg,ncol=nf)
	for (j in seq(nf)) data[,j] = readBin(con,"numeric",n=frame$npdg)
	nf2 = readBin(con,"integer",size=8,n=1)
	close(con)

	stopifnot(nf == nf2)
	data
}

getSPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin))
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",size=8)
	data = readBin(con,"numeric",n=nsp2)
	close(con)

	stopifnot(nsp2 == frame$nwave*(frame$nwave-1))

	matrix(data,nrow=frame$nwave)
}

degrade = function(frame,ilat,nlat=2,nlon=1)
{
	if (is.null(ilat)) {
		stopifnot(frame$nlat%/%2*2 == frame$nlat)
		ilatn = seq(1,frame$nlat/2,nlat)
		ilat = c(2*ilatn-1,rev(frame$nlat-2*(ilatn-1)))
		nlon = 2
	}

	stopifnot(all(frame$nlong%/%nlon*nlon == frame$nlong))

	clats = c(0,cumsum(frame$nlong))
	lats = datatoGauss(frame$lat,frame)
	longs = datatoGauss(frame$long,frame)

	stopifnot(all(ilat %in% seq(frame$nlat)) && identical(ilat,sort(ilat)))
	# save
	nlong = frame$nlong

	# update
	frame$nlat = length(ilat)
	frame$nlong = frame$nlong[ilat]%/%nlon
	frame$npdg = sum(frame$nlong)
	frame$theta = frame$theta[ilat]

	frame$lat = frame$long = numeric(frame$npdg)
	frame$ind = integer(frame$npdg)

	ip = 0
	for (i in ilat) {
		ind = seq(1,nlong[i],nlon)
		np = length(ind)
		frame$lat[ip+1:np] = lats[ind,i]
		frame$long[ip+1:np] = longs[ind,i]
		frame$ind[ip+1:np] = clats[i]+ind
		ip = ip+np
	}

	frame
}

toGauss = function(frame)
{
	longs = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (i in seq(frame$nlat)) {
		nl = frame$nlong[i]
		longs[seq(nl),i] = 360*(seq(nl)-1)/nl
	}

	list(lats=frame$theta,longs=longs)
}

datatoGauss = function(data,frame)
{
	clats = c(0,cumsum(frame$nlong))
	datag = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (ilat in seq(frame$nlat)) {
		ip = clats[ilat]
		np = frame$nlong[ilat]
		datag[1:np,ilat] = data[(ip+1):(ip+np)]
	}

	datag
}

equalize = function(mapf,nc=5,offset=2)
{
	# nc equally numbered groups of latitudes (mapf: normalized grid-point density)
	mapc = quantile(mapf,prob=seq(nc-1)/nc)

	# mean of mapf for each group
	im = findInterval(mapf,mapc)+1
	mapm = tapply(mapf,im,mean)

	# index of last element of each group
	indl = c(match(seq(2,nc),im)-1,length(mapf))

	# vector of lat indices
	ilat = integer()
	i1 = 1
	for (i in seq(nc)) {
		# equalizing density (offset to limit n)
		n = max(as.integer((max(mapm)+offset)/(mapm[i]+offset)),1)
		cat("--> lat group/from/to/by:",i,i1,indl[i],n,"\n")
		ilat = c(ilat,seq(i1,indl[i],by=n))
		i1 = indl[i]+1
	}

	ilat
}

domain = function(frame,ind,dlat,dlon)
{
	stopifnot(length(ind) == 1)
	stopifnot(length(dlat) == 2 && length(dlon) == 2)

	xlim = frame$long[ind]+dlon
	ylim = frame$lat[ind]+dlat

	list(xlim=xlim,ylim=ylim)
}

interp = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	ip1 = clato[ilato]+ilono1
	ip2 = ip1%%clato[ilato+1]+1
	d1 = datao[ip1,]+e1*(datao[ip2,]-datao[ip1,])
	ip1 = clato[ilato+1]+ilono2
	ip2 = ip1%%clato[ilato+2]+1
	d2 = datao[ip1,]+e2*(datao[ip2,]-datao[ip1,])
	d1+e0*(d2-d1)
}

interpGauss = function(datao,framo,frame)
{
	# only for Gauss grids nested inside
	stopifnot(length(framo$nlong) >= length(frame$nlong))
	stopifnot(max(framo$nlong) >= max(frame$nlong))

	grid = toGauss(frame)
	grido = toGauss(framo)

	data = matrix(nrow=frame$npdg,ncol=dim(datao)[2])

	clats = c(0,cumsum(frame$nlong))
	clato = c(0,cumsum(framo$nlong))
	dlat = 180/framo$nlat
	dlon = 360/framo$nlong

	for (ip in seq(frame$npdg)) {
		ilat = which.max(ip <= clats[-1])
		ilon = ip-clats[ilat]
		e0 = (90-dlat/2-grid$lats[ilat])/dlat
		ilato = floor(e0)+1
		e0 = e0-(ilato-1)
		stopifnot(0 < ilato && ilato < length(framo$nlong))
		stopifnot(0 <= e0 && e0 < 1)

		e1 = grid$longs[ilon,ilat]/dlon[ilato]
		e2 = grid$longs[ilon,ilat]/dlon[ilato+1]
		ilono1 = floor(e1)+1
		ilono2 = floor(e2)+1
		stopifnot(0 < ilono1 && ilono1 <= framo$nlong[ilato])
		stopifnot(0 < ilono2 && ilono2 <= framo$nlong[ilato+1])
		e1 = e1-(ilono1-1)
		e2 = e2-(ilono2-1)
		stopifnot(0 <= e1 && e1 < 1)
		stopifnot(0 <= e2 && e2 < 1)
		data[ip,] = interp(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
	}

	data
}

xdiff = function(data,frame)
{
	grid = toGauss(frame)
	dlon = cos(grid$lats*pi/180)*360/frame$nlong
	gm = dlon[frame$nlat%/%2]/dlon

	clats = c(0,cumsum(frame$nlong))

	for (i in seq(frame$nlat)) {
		ip = clats[i]+seq(frame$nlong[i])

		d1 = gm[i]*(data[ip[1],]-data[ip[length(ip)],])
		for (j in seq(dim(data)[2])) data[ip[-1],j] = gm[i]*diff(data[ip,j])
		data[ip[1],] = d1
	}

	data
}

dilat = function(frame)
{
	clats = c(0,cumsum(frame$nlong))

	ddm = numeric(frame$nlat)

	for (i in seq(frame$nlat)) {
		ip = clats[i]+seq(frame$nlong[i])

		ddm[i] = mean(sqrt(diff(frame$lat[ip])^2+diff(frame$long[ip])^2))
	}

	ddm/max(ddm[1:frame$nlat/2])
}

lat45 = function(data,frame)
{
	grid = toGauss(frame)
	dlon = cos(grid$lats*pi/180)*360/frame$nlong
	gm = dlon[frame$nlat%/%2]/dlon

	clats = c(0,cumsum(frame$nlong))

	i = frame$nlat%/%4
	ip = clats[i]+seq(frame$nlong[i])
	data[ip,,drop=FALSE]
}

interpAB = function(datao,framo,frame)
{
	eta = frame$A/Gvp0+frame$B
	etao = framo$A/Gvp0+framo$B

	ind = findInterval(eta,etao)
	e = (eta-etao[ind])/(etao[ind+1]-etao[ind])
	stopifnot(all(0 < ind & ind < length(etao)))
	stopifnot(all(0 <= e & e < 1))

	data = matrix(nrow=dim(datao)[1],ncol=frame$nlevel)

	for (i in seq(frame$nlevel)) {
		d1 = datao[,ind[i]]
		d2 = datao[,ind[i]+1]
		data[,i] = d1+e[i]*(d2-d1)
	}

	data
}

vgrad = function(data,frame)
{
	eta = frame$A/Gvp0+frame$B
	deta = diff(eta[frame$ilev])
	data = apply(data,1,function(x) diff(x)/deta)

	t(data)
}

vmean = function(data,m=median,decreasing=TRUE)
{
	np = dim(data)[1]
	nlev = dim(data)[2]

	indf = apply(data,1,order,decreasing=decreasing)
	mm = matrix(c(1:np,indf[,1]),ncol=2)
	ddm01 = data[mm]
	if (nlev < 4) return(ddm01)

	mm = matrix(c(1:np,indf[,1]),ncol=2)
	ddm02 = data[mm]

	dd = data
	ip = logical(np)

	for (im in seq(nlev%/%4)) {
		for (i in seq(nlev-2*im)+im) dd[,i] = apply(dd,1,function(x) m(x[c(i-im,i+im)]))
		#indf = apply(dd,1,order,decreasing=decreasing)
		mm = matrix(c(1:np,indf[,1]),ncol=2)
		ddm1 = dd[mm]
		mm = matrix(c(1:np,indf[,2]),ncol=2)
		ddm2 = dd[mm]

		ip = ip | abs(ddm1-ddm01) < abs(ddm2-ddm02)
		ddm1[ip] = ddm2[ip]
		if (all(ip)) break

		iip = ! ip
		ddm01[iip] = ddm1[iip]
		ddm02[iip] = ddm2[iip]
	}

	ddm1
}

inDomain = function(frame,domain)
{
	ilat = domain$ylim[1] <= frame$lat & frame$lat <= domain$ylim[2]

	xlim = (domain$xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= frame$long & frame$long <= xlim[2]
	} else {
		ilon = frame$long <= xlim[2] | frame$long >= xlim[1]
	}

	ilat & ilon
}

subFrame = function(frame,ind)
{
	list(lat=frame$lat[ind],long=frame$lon[ind])
}

zoom = function(data,frame,domain)
{
	mask = inDomain(frame,domain)
	ind = which(mask)

	list(lat=frame$lat[ind],long=frame$lon[ind],data=data[ind,,drop=FALSE])
}

section = function(data,frame,long,lat=c(-90,90))
{
	stopifnot(0 <= long && long < 360)
	stopifnot(all(-90 <= lat & lat <= 90))

	grid = toGauss(frame)
	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	if (length(lat) == 1) {
		dlat = 180/frame$nlat
		ilat = which.max(grid$lats <= lat)
		e0 = (90-dlat/2-grid$lats[ilat])/dlat-ilat+1
		ind = c(ilat,ilat+1)
	} else {
		ind = which(min(lat) <= grid$lats & grid$lats <= max(lat))
	}

	data1 = matrix(nrow=length(ind),ncol=dim(data)[2])
	lats = longs = numeric(length(ind))

	for (i in seq(along=ind)) {
		ilat = ind[i]
		e = long/dlon[ilat]
		ilon = floor(e)+1
		e = e-(ilon-1)
		ip1 = clats[ilat]+ilon
		ip2 = ip1%%clats[ilat+1]+1
		data1[i,] = (1-e)*data[ip1,]+e*data[ip2,]
		lats[i] = (1-e)*frame$lat[ip1]+e*frame$lat[ip2]
		longs[i] = frame$long[ip1]
	}

	if (length(lat) == 1) data1 = data1[1,]+e0*(data1[2,]-data1[1,])

	list(lat=grid$lats[ind],lats=lats,longs=longs,data=data1)
}

sectionf = function(data,frame,lat,long=c(0,360))
{
	long = unique(long)
	stopifnot(-90 <= lat && lat < 90)
	stopifnot(length(long) == 2)
	stopifnot(all(0 <= long & long <= 360))

	grid = toGauss(frame)
	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	ilat = max(which(grid$lats >= lat))

	if (long[1] < long[2]) {
		ind = which(long[1] <= grid$longs[,ilat] & grid$longs[,ilat] <= long[2])
	} else {
		ind = c(which(long[1] <= grid$longs[,ilat] & grid$longs[,ilat] <= 360),
			which(0 <= grid$longs[,ilat] & grid$longs[,ilat] <= long[2]))
	}

	data1 = matrix(nrow=length(ind),ncol=dim(data)[2])

	for (j in seq(dim(data)[2])) {
		datag = datatoGauss(data[,1],frame)
		data1[,j] = datag[ind,ilat]
	}

	lats = longs = numeric(length(ind))

	e = (lat-grid$lats[ilat])/diff(grid$lats[ilat+0:1])

	ip = clats[ilat]+ind
	lats = frame$lat[ip]
	longs = frame$long[ip]

	list(long=grid$longs[ind,ilat],lats=lats,longs=longs,data=data1)
}

mapxy = function(xy,main=NULL,axes=TRUE,frame.plot=TRUE,cex.axis=.9,mar=rep(3,4)+.1,...)
{
	l = map(xlim=xy$xlim,ylim=xy$ylim,mar=mar,...)

	title(main)

	if (axes) {
		axis(1,cex.axis=cex.axis)
		axis(2,cex.axis=cex.axis)
	}

	if (frame.plot) box()

	l
}

plotxy = function(x,y,data,breaks="Sturges",palette=terrain.colors,pch=15,legend.title="",
	ppi=54,cex.axis=.9,...)
{
	u = par("usr")
	f = par("fin")

	if (ppi > 144) stop("ppi > 144\n")

	npmax = prod(f*ppi)
	if (npmax < .Machine$integer.max && length(data) > as.integer(npmax)) {
		npmax = as.integer(npmax)
		cat("--> reducing plot from",length(data),"to",npmax,"points\n")
		ind = seq(1,length(data),length.out=npmax)

		x = x[ind]
		y = y[ind]
		data = data[ind]
	}

	h = hist(data,breaks,plot=FALSE)
	br = h$breaks

	ind = findInterval(data,br)
	cols = palette(length(br))

	points(x,y,col=cols[ind],pch=pch,...)

	if (! is.na(legend.title)) {
		DOPfill.legend(levels=br,col=cols,title=legend.title,cex.axis=cex.axis)
	}
}

plotxy2 = function(x,y,zx,zy,breaks="Sturges",colvec=TRUE,length=.05,angle=10,ppi=6,
	legend.title="",cex.axis=.9,...)
{
	f = par("fin")

	if (ppi > 144) stop("ppi > 144\n")

	npmax = prod(f*ppi)
	if (npmax < .Machine$integer.max && length(data) > as.integer(npmax)) {
		npmax = as.integer(npmax)
		cat("--> reducing plot from",length(data),"to",npmax,"points\n")
		ind = seq(1,length(x),length.out=npmax)

		x = x[ind]
		y = y[ind]
		zx = zx[ind]
		zy = zy[ind]
	}

	ffz = sqrt(zx^2+zy^2)

	u = par("usr")
	ux = diff(u[1:2])
	uy = diff(u[3:4])

	# scale of wind vectors: a fraction of a "small" diagonal
	ffu = sqrt(ux^2+uy^2/4)/(2*sqrt(length(zx)))

	if (colvec) {
		fx = zx/ffz*ffu
		fy = zy/ffz*ffu

		h = hist(ffz,plot=FALSE)
		br = h$breaks

		ic = findInterval(ffz,br)
		palette = rainbow(length(br),start=2/3,end=0)
		cols = palette[ic]
	} else {
		ff1 = ffu*scale(ffz,center=FALSE)
		fx = zx/ffz*ff1
		fy = zy/ffz*ff1
		cols = 1
	}

	x2 = x+fx
	y2 = y+fy

	ffi = sqrt((fx/ux*f[1])^2+(fy/uy*f[2])^2)

	ind = which(ffi > 1.e-3)
	arrows(x[ind],y[ind],x2[ind],y2[ind],length,angle,col=cols,...)

	if (colvec && ! is.na(legend.title)) {
		DOPfill.legend(levels=br,col=palette,title=legend.title,cex.axis=cex.axis)
	}
}

DOPfill.legend = function(x,y,width,height,levels,col,title,border=NA,cex.axis,...)
{
	u = par("usr")
	p = par("plt")
	if (missing(width)) width = (1-p[2])/5*diff(u[1:2])/diff(p[1:2])
	if (missing(x)) x = u[2]+width/5

	nl = length(levels)
	if (missing(height)) height = diff(u[3:4])
	dy = height/(nl-1)
	if (missing(y)) y = u[3]
	ybas = y + dy*(seq(nl-1)-1)
	yhaut = ybas + dy

	op = par(xpd=TRUE)
	rect(x,ybas,x+width,yhaut,col=col,border=border)
	if (! missing(title)) text(x+width/2,y+height,title,adj=c(.5,0))
	par(op)

	op = par(las=2,yaxt="s")
	axis(4,at=c(ybas[1],yhaut),labels=levels,tick=FALSE,pos=x,cex.axis=cex.axis,...)
	par(op)
}

prettydiff = function(x,...)
{
	lab = pretty(x,...)
	n = length(lab)
	i = n%/%8
	j = i
	if (n%/%2*2 == n && j > 0) j = j-1
	lab[-(n%/%2+seq(-j,i))]
}

revpos = function(x,at)
{
	# for any xi in [a,b], xi = a+(xi-a), xi' = b-(xi-a) = a+b-xi
	sum(range(x))-at
}

plotv = function(x,y,z,ylim=NULL,xaxt="s",yaxt="s",xrev=FALSE,...)
{
	if (! is.null(ylim)) {
		indy = which(ylim[1] <= y & y <= ylim[2])
		y = y[indy]
		indy = rev(indy)
	} else {
		indy = rev(seq(along=y))
	}

	if (xrev) {
		contour(rev(x),y,z[,indy],xaxt="n",yaxt="n",...)
	} else {
		contour(x,y,z[,indy],xaxt="n",yaxt="n",...)
	}

	if (xaxt != "n") {
		lab = pretty(x/30)*30

		if (xrev) {
			axis(1,at=revpos(x,lab),lab)
		} else {
			axis(1,at=lab)
		}
	}

	if (yaxt != "n") {
		lab = pretty(y)
		axis(2,revpos(y,lab),lab)
	}
}

plotf = function(y,z,ylim=NULL,plot.axes=TRUE,...)
{
   if (is.null(ylim)) {
      indy = rev(seq(along=y))
   } else {
      indy = which(ylim[1] <= y & y <= ylim[2])
      y = y[indy]
      indy = rev(indy)
   }

	filled.contour(y=y,z=z[,indy],plot.axes=plot.axes,...)

	if (plot.axes) {
		lab = pretty(y)
		axis(2,revpos(y,lab),lab)
   }
}

plotz = function(x,y,type="o",pch="+",...)
{
	plot(x,y,type=type,pch=pch,...)
}

getField = function(fic,param,symbol,frame,frlow=frame,interp=FALSE)
{
	stopifnot(! is.null(frlow$ilev))

	fsave = sprintf("%s/e%d/%s.RData",dirname(fic),frame$step/3600,symbol)
	if (file.exists(fsave)) {
		ilev = 0
		load(fsave)
		if (dim(data)[1] == frlow$npdg && identical(ilev,frlow$ilev)) return(data)

		cat("--> differing geometry, read data again\n")
	}

	data = getGPFields(fic,param,frame)

	if (frame$npdg != frlow$npdg) {
		if (interp) {
			data = interpGauss(data,frame,frlow)
		} else {
			stopifnot("ind" %in% names(frlow))
			data = data[frlow$ind,,drop=FALSE]
		}
	}

	if (dim(data)[2] > 1 && frame$nlevel != frlow$nlevel)
		data = interpAB(data,frame,frlow)

	ilev = frlow$ilev
	data = data[,ilev,drop=FALSE]

	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	save("data","ilev",file=fsave)

	data
}

mapdom = function(dom,frame,data,...)
{
	l = mapxy(dom,...)
	plotxy(frame$long,frame$lat,data)
	lines(l)
}

mapex = function(data,frame,desc)
{
	xy = zoom(data,frame,monde)
	if (! hasx11) png(sprintf("world_%s.png",desc$symbol))
	mapdom(monde,xy,xy$data[,1],main=desc$longname)
	if (hasx11 && ask) invisible(readline("Press enter to continue"))

	if (FALSE) {
	xy = zoom(data,frame,france)
	if (! hasx11) png(sprintf("france_%s.png",desc$symbol))
	mapdom(france,xy,xy$data[,1],main=desc$longname)
	if (hasx11 && ask) invisible(readline("Press enter to continue"))
	}

	xy = zoom(data,frame,europe)
	if (! hasx11) png(sprintf("euro_%s.png",desc$symbol))
	mapdom(europe,xy,xy$data[,1],main=sprintf("Field %s",desc$longname))
	if (hasx11 && ask) invisible(readline("Press enter to continue"))

	#xy = zoom(data,frame,hima)
	#if (! hasx11) png(sprintf("hima_%s.png",desc$symbol))
	#mapdom(hima,xy,xy$data[,1],main=desc$longname)
	#if (hasx11 && ask) invisible(readline("Press enter to continue"))
	#dev.off()

	#xz = sectionf(data,frame,30,c(65,110))
	#xzdf = as.data.frame(xz)
	#write.table(xzdf,"hima30N.txt",quote=FALSE,row.names=FALSE)

	#yz = section(data,frame,80,hima$ylim)
	#yzdf = as.data.frame(yz)
	#write.table(yzdf,"hima80E.txt",quote=FALSE,row.names=FALSE)
}

Gficbin = tempfile(fileext=".bin")
Gxlab="<- Equator  ...  //  ...  France  ...  //  ...  North pole ->"

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

if (file.exists(cargs$fic) && file.info(cargs$fic)$isdir) {
	fics = dir(cargs$fic,sprintf("\\.%s",cargs$ext),full.names=TRUE)
} else if (regexpr("\\.txt$",cargs$fic) > 0) {
	fics = scan(cargs$fic,"character",quiet=TRUE)
} else {
	fics = cargs$fic
}

stopifnot(length(fics) > 0 && all(file.exists(fics)))

if ("date" %in% names(cargs)) {
	cbase = strsplit(cargs$date,split="/")[[1]]
	base = as.POSIXct(cbase[1],format="%Y%m%d")+as.integer(cbase[2])
	cat("--> base date:",base,"\n")
}

ilev = 0
if ("level" %in% names(cargs)) {
	if (length(cargs$level) == 1 && file.exists(cargs$level)) {
		ilev = scan(cargs$level,what=integer())
	} else {
		ilev = as.integer(cargs$level)
	}

	cat("--> level indices/number (0: all levels):",ilev,"\n")
}

#graph = "graph" %in% names(cargs) || length(fics) == 1

monde = list(xlim=c(-170,170),ylim=c(-85,85))
france = list(xlim=c(-6,11),ylim=c(40,53))
europe = list(xlim=c(-15,32),ylim=c(35,59))
asie = list(xlim=c(60,120),ylim=c(20,55))
hima = list(xlim=c(65,110),ylim=c(25,50))

desc = read.table(sprintf("%s/params.txt",diags),header=TRUE)

if ("params" %in% names(cargs)) {
	if (length(cargs$params) == 1 && file.exists(cargs$params)) {
		params = scan(cargs$params,what=character())
		ind = match(params,desc$faname)
	} else {
		ind = match(cargs$params,desc$faname)
	}

	stopifnot(all(! is.na(ind)))
	desc = desc[ind,]
}

npar = dim(desc)[1]

if (! "path" %in% names(cargs)) cargs$path="."

ficana = "ARPE.0000.000"
#ficana[inds] = "analysis.surf-arpege.tl1798-c22.fa"
#ficana[! inds] = "analysis.atm-arpege.tl1798-c22.fa"

dstats = array(NA_real_,c(length(fics),5,npar))

hasx11 = capabilities("X11")
ask = hasx11 && interactive()

if (interactive()) browser()

for (i in seq(along=fics)) {
	cat("\nfile",basename(fics[i]),":\n")

	frame = getFrame(fics[i])
	cat("--> base/step:",as.character(frame$base,"%F %R %Z"),"+",frame$step/3600,"h\n")
	if (frame$nlat > 800) {
		frlow = degrade(frame,2,2)
	} else {
		frlow = frame
	}

	mapf = dilat(frlow)
	#ilateq = equalize(mapf,offset=3)
	#if (length(ilateq) < frlow$nlat) frlow = degrade(frlow,ilateq)

	if (! hasx11) png("stretch_%s.png")
	plot(mapf,type="l",main=c("Stretching coefficient by latitude","normalized at Equator"),
		xlab="Latitude index",ylab="Stretching coef.",log="y")
	#rug(ilateq,.02)
	#points(ilateq,mapf[ilateq],pch=20,cex=.5)
	if (ask) invisible(readline("Press enter to continue"))

	eta = frlow$A/Gvp0+frlow$B
	stopifnot(all(abs(ilev) <= length(eta)))
	nlev = length(eta)
	if (length(ilev) == 1 && ilev < 0 && nlev > -ilev) {
		# selection
		e = seq(min(eta),max(eta),length.out=-ilev)
		frlow$ilev = sapply(e,function(x) which.min(abs(x-eta)))
		indi = frlow$ilev
	} else if (length(ilev) < nlev && all(ilev > 0)) {
		frlow$ilev = ilev
		indi = ilev
	} else {
		frlow$ilev = seq(nlev)
		indi = frlow$ilev
	}

	etai = eta[indi]
	yeta = rev(range(eta))

	date = frame$base+frame$step
	#date = base+frame$step
	ficref = sprintf("%s/%s/%s",cargs$path,strftime(date,"%H/%Y%m%d"),ficana)
	if (! file.exists(ficref))
		ficref = sprintf("%s/%s/%s",cargs$path,strftime(date,"%H/%Y%m%d"),fics[i])

	datax = datay = dataz = NULL
	dataox = dataoy = dataoz = NULL

	for (j in seq(npar)) {
		if (desc$ltype[j] == "-") {
			patt = desc$faname[j]
		} else if (length(ilev) > 1) {
			patt = sprintf("%s*%s",desc$ltype[j],desc$faname[j])
		} else if (desc$ltype[j] == "S") {
			patt = sprintf("S%03d%s",ilev,desc$faname[j])
		} else {
			cat("--> level type not 'S' for reading specified level\n")
			next
		}

		cat(". get fields matching pattern",patt,"\n")
		data = getField(fics[i],patt,desc$symbol[j],frame,frlow)

		if (! hasx11) png(sprintf("hist_%s.png",desc$symbol[j]))
		hist(data,col="whitesmoke",main=sprintf("Distribution of %s",desc$longname[j]))
		if (ask) invisible(readline("Press enter to continue"))

		nl = dim(data)[2]
		dstats[i,,j] = quantile(data[,nl],c(0,.1,.5,.9,1),names=FALSE)

		ddx = xdiff(data,frlow)

		cat("Statistics of zonal difference:\n")
		print(summary(as.vector(ddx)))

		if (! hasx11) png(sprintf("spec45_%s.png",desc$symbol[j]))
		d45 = lat45(data,frlow)
		spectrum(d45[,nl],c(3,3),main="Smoothed spectral density of lat ~45")
		if (ask) invisible(readline("Press enter to continue"))

		if (nl == 1) {
			mapex(data,frlow,desc[j,])
		} else {
			ix = arrayInd(which.max(abs(data)),dim(data))
			il = ix[1,2]
			ix = ix[1,1]
			cat("Level of max value:",il,"/",nl,"\n")
			mapex(data[,il,drop=F],frlow,desc[j,])

			cat("Vertical profile at max value\n")
			dom = domain(frlow,ix,c(-20,20),c(-30,30))
			xy = zoom(data,frlow,dom)
			print(summary(as.vector(xy$data)))

			if (! hasx11) png(sprintf("maxprof_%s.png",desc$symbol[j]))
			tt = "Vertical profile at max value"
			plotz(data[ix,],etai,ylim=yeta,main=c(tt,desc$longname[j]),
				xlab=desc$longname[j],ylab=expression(eta))
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

			cat("Mean vertical profile\n")
			datam = apply(data,2,mean)

			if (! hasx11) png(sprintf("meanprof_%s.png",desc$symbol[j]))
			tt = sprintf("Field %s",desc$longname[j])
			plotz(datam,etai,ylim=yeta,main=c("Mean vertical profile",tt),
				xlab=desc$longname[j],ylab=expression(eta))
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

			if (FALSE) {
			ddz = vgrad(data,frlow)
			ddz = t(apply(data,1,diff))

			cat("Statistics of vertical gradient:\n")
			print(summary(as.vector(ddz)))
			ix = arrayInd(which.max(abs(data)),dim(data))
			il = ix[1,2]
			ix = ix[1,1]
			plotz(ddz[ix,],etai[-1],ylim=yeta,main=c("Mean vertical profile",tt),
				xlab=desc$longname[j],ylab=expression(eta))
			}

			if (desc$ltype[j] == "S") {
				yzs = section(data,frlow,180,c(0,90))
				yzn = section(data,frlow,0,c(0,90))
				ind = order(yzs$lat)

				if (! hasx11) png(sprintf("section_%s.png",desc$symbol[j]),720)
				par(mfrow=c(1,2),xaxs="i",yaxs="i")
				tt = sprintf("Cross-section of %s",desc$longname[j])
				plotv(yzs$lat[ind],etai,yzs$data[ind,],main=c(tt,"tilted meridian 180"),
					xlab="tilted lat",ylab="eta")
				plotv(yzn$lat,etai,yzn$data,main=c(tt,"tilted meridian 0"),
					xlab="tilted lat",ylab="eta",xrev=TRUE)
				if (hasx11 && ask) invisible(readline("Press enter to continue"))
			}
		}

		if (desc$save[j] == "X") {
			datax = data
		} else if (desc$save[j] == "Y") {
			datay = data
		} else if (desc$save[j] == "Z") {
			dataz = data
		}

		if (! file.exists(ficref)) next

		cat(". get ref field",desc$faname[j],"\n")
		framo = getFrame(ficref)
		datao = getField(ficref,patt,desc$symbol[j],framo,frlow,interp=TRUE)

		if (desc$save[j] == "X") {
			dataox = datao
		} else if (desc$save[j] == "Y") {
			dataoy = datao
		} else if (desc$save[j] == "Z") {
			dataoz = datao
		}

		ddiff = data-datao

		cat("Statistics around diff min/max value:\n")
		i1 = arrayInd(which.min(ddiff),dim(ddiff))[1,1]
		i2 = arrayInd(which.max(ddiff),dim(ddiff))[1,1]
		dom = domain(frlow,i1,c(-20,20),c(-30,30))
		xy = zoom(ddiff,frlow,dom)
		print(summary(as.vector(xy$data)))
		drms = sqrt(mean(xy$data^2))
		cat("RMSE at local min:",drms,"\n")

		if (! hasx11) png(sprintf("diffmin_%s.png",desc$symbol[j]))
		l = mapxy(dom,main=desc$longname[j])
		plotxy(frlow$long,frlow$lat,ddiff[,1])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))

		dom = domain(frlow,i2,c(-20,20),c(-30,30))
		xy = zoom(ddiff,frlow,dom)
		print(summary(as.vector(xy$data)))
		drms = sqrt(mean(xy$data^2))
		cat("RMSE at local max:",drms,"\n")

		if (! hasx11) png(sprintf("diffmin_%s.png",desc$symbol[j]))
		l = mapxy(dom,main=desc$longname[j])
		plotxy(frlow$long,frlow$lat,ddiff[,1])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))

		cat("Statistics of forecast error:\n")
		print(summary(as.vector(ddiff)))
		drms = sqrt(mean(ddiff^2))
		cat("Global bias/RMSE:",mean(ddiff),drms,"\n")

		if (! hasx11) png(sprintf("histdiff_%s.png",desc$symbol[j]))
		tt = sprintf("Distribution of difference from analysis of %s",desc$longname[j])
		hist(ddiff,main=tt)
		if (ask) invisible(readline("Press enter to continue"))

		if (dim(data)[2] == 1) {
			mapex(ddiff,frlow,desc[j,])
		} else if (desc$ltype[j] == "S") {
			if (! hasx11) png(sprintf("profile_%s.png",desc$symbol[j]))
			par(mfrow=c(1,2))
			tt = sprintf("Field %s",desc$longname[j])
			plotz(data[i1,],etai,ylim=yeta,main=c(tt,"Profile at point 'min of error'"),
				xlab=desc$longname[j],ylab=expression(eta))
			plotz(data[i2,],etai,ylim=yeta,main=c(tt,"Profile at point 'max of error'"),
				xlab=desc$longname[j],ylab=expression(eta))
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

         yzs = section(ddiff,frlow,180,c(0,90))
         yzn = section(ddiff,frlow,0,c(0,90))
         ind = order(yzs$lat)

         if (! hasx11) png(sprintf("sectiondiff_%s.png",desc$symbol[j]),720)
         par(mfrow=c(1,2),xaxs="i",yaxs="i")
			tt = sprintf("Cross-section diff of %s",desc$longname[j])
         plotv(yzs$lat[ind],etai,yzs$data[ind,],main=c(tt,"tilted meridian 180"),
            xlab="tilted lat",ylab="eta",nlevels=7)
         plotv(yzn$lat,etai,yzn$data,main=c(tt,"tilted meridian 0"),
            xlab="tilted lat",ylab="eta",nlevels=7,xrev=TRUE)
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

			ns = length(yzs$lat)
			dataf = yzn$data
			dataf[1:ns,] = yzs$data[ind,]
			#levels = prettydiff(dataf,7)

			if (! hasx11) png(sprintf("sectionvdiff_%s.png",desc$symbol[j]),720)
			x = seq(0,1,length.out=dim(dataf)[1])
			plotv(x,etai,dataf,nlevels=7,main=tt,xlab=Gxlab,ylab="eta",xaxs="i",yaxs="i")
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

			if (! hasx11) png(sprintf("sectionfill_%s.png",desc$symbol[j]),720)
			plotf(etai,dataf,nlevels=7,color.palette=cm.colors,main=tt,xlab=Gxlab,ylab="eta")
			if (hasx11 && ask) invisible(readline("Press enter to continue"))

			dbias = apply(ddiff,2,mean)
			drms = apply(ddiff,2,function(x) sqrt(mean(x^2)))
			cat("Global bias/RMSE by level:\n")
			print(summary(dbias))
			print(summary(drms))
		}
	}

	if (! is.null(datax) && ! is.null(datay)) {
		ff = sqrt(datax^2+datay^2)
		hist(ff,col="whitesmoke",main="Distribution of wind speed")
		if (ask) invisible(readline("Press enter to continue"))

		print(summary(as.vector(ff)))
		x = zoom(datax,frlow,france)
		y = zoom(datay,frlow,france)
		if (! hasx11) png(sprintf("france_uv.png"))
		l = mapxy(france,main="u/v wind components")
		plotxy2(x$long,x$lat,x$data[,nl],y$data[,nl])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))

		x = zoom(datax,frlow,europe)
		y = zoom(datay,frlow,europe)
		if (! hasx11) png(sprintf("euro_uv.png"))
		l = mapxy(europe,main="u/v wind components")
		plotxy2(x$long,x$lat,x$data[,nl],y$data[,nl])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))
	}

	if (! is.null(dataox) && ! is.null(dataoy)) {
		ffo = sqrt(dataox^2+dataoy^2)
		hist(ffo,col="whitesmoke",main="Distribution of wind speed")
		if (ask) invisible(readline("Press enter to continue"))

		print(summary(as.vector(ffo)))
		x = zoom(dataox,frlow,france)
		y = zoom(dataoy,frlow,france)
		if (! hasx11) png(sprintf("france_uvo.png"))
		l = mapxy(france,main="u/v wind components")
		plotxy2(x$long,x$lat,x$data[,nl],y$data[,nl])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))

		x = zoom(dataox,frlow,europe)
		y = zoom(dataoy,frlow,europe)
		if (! hasx11) png(sprintf("euro_uvo.png"))
		l = mapxy(europe,main="u/v wind components")
		plotxy2(x$long,x$lat,x$data[,nl],y$data[,nl])
		lines(l)
		if (ask) invisible(readline("Press enter to continue"))
	}

	if (! is.null(datax) && ! is.null(datay) && ! is.null(dataox) && ! is.null(dataoy)) {
		ffd = ff-ffo
		print(summary(as.vector(ffd)))
		drms = sqrt(mean(ffd^2))
		cat("Global bias/RMSE:",mean(ffd),drms,"\n")

		if (! hasx11) png("histdiff_ff.png")
		hist(ffd,col="whitesmoke",main="Distribution of difference from ref of wind speed")
		if (ask) invisible(readline("Press enter to continue"))

		xy = zoom(ffd,frlow,europe)
		if (! hasx11) png(sprintf("euro_ffd.png"))
		mapdom(europe,xy,xy$data[,1],main="Difference of wind speed")
		if (ask) invisible(readline("Press enter to continue"))
	}
}

if (length(fics) > 1) {
	if (! hasx11) png("dstats.png")
	tt = sprintf("Forecast of %s",desc$longname)
	for (j in seq(along=desc$faname)) {
		matplot(dstats[,,j],type="o",lty=1,col=1,pch=c("n","1","5","9","x"),
			xlab="Forecast time",ylab=desc$symbol[j],main=tt[j])
		if (hasx11 && ask) invisible(readline("Press enter to continue"))
	}
}
