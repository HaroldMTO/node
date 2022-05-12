#R --args fic=ICMSHARPE+0000 date=20210329/00 path=/scratch/work/petithommeh/oper/arpege/4dvarfr/OPER
#R --args fic=PFARPE0000+2400 date=20210329/00 path=/scratch/work/petithommeh/oper/arpege/4dvarfr/OPER
# fic=pfarpe.txt level=levels.txt params=uv.txt
diags = "~/util/diags"

library(maps)

Gvp0 = 101325

getFrame = function(file)
{
	system(sprintf("epy_dump.py %s -f frame -o %s",file,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	dims = readBin(con,"integer",size=8,n=5)

	# global or LAM
	if (dims[5] == 0) {
		nlong = readBin(con,"integer",size=8,n=dims[1])
		nwave = readBin(con,"integer",size=8,n=dims[2])
		ngp = dims[4]
		nl = dims[3]
	} else {
		ngp = prod(dims[1:2])
		nl = dims[5]
	}

	base = readBin(con,"numeric",n=1)
	base = as.POSIXct(base,origin="1970-01-01")
	step = readBin(con,"integer",size=8,n=1)
	if (regexpr(".+\\+0*([0-9]+)",file) > 0) {
		ech = as.integer(gsub(".+\\+0*([0-9]+)","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	lats = readBin(con,"numeric",n=ngp)
	longs = readBin(con,"numeric",n=ngp)

	Ai = readBin(con,"numeric",n=nl)
	Bi = readBin(con,"numeric",n=nl)

	close(con)

	if (dims[5] == 0) {
		nlat = length(nlong)
		theta = 90-180*(seq(nlat)-.5)/nlat

		frame = list(nlat=nlat,nwave=dims[2],nlevel=nl,npdg=ngp,nlong=nlong,
			theta=theta,lat=lats,long=longs,A=Ai,B=Bi,base=base,step=step,lam=FALSE)
	} else {
		frame = list(nlat=dims[1],nlong=dims[2],nwavex=dims[3],nwavey=dims[4],nlevel=nl,
			npdg=ngp,lat=lats,long=longs,A=Ai,B=Bi,base=base,step=step,lam=TRUE)
	}

	frame
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

	if (dim(data)[2] > 1) {
		if (frame$nlevel != frlow$nlevel) data = interpAB(data,frame,frlow)

		ilev = frlow$ilev
		data = data[,ilev,drop=FALSE]
	}

	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	save("data","ilev",file=fsave)

	data
}

getGPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
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
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",size=8)
	data = readBin(con,"numeric",n=nsp2)
	close(con)

	stopifnot(nsp2 == frame$nwave*(frame$nwave-1))

	matrix(data,nrow=frame$nwave)
}

degrade = function(frame,ilat,nlat,nlon=2)
{
	if (frame$lam) {
		stopifnot(length(frame$nlong) == 1)

		ilat = seq(1,frame$nlat,nlat)
		ilon = seq(1,frame$nlong,nlon)
		frame$ind = rep((ilat-1)*frame$nlong,each=length(ilon))+ilon
		frame$lat = frame$lat[frame$ind]
		frame$long = frame$long[frame$ind]
		frame$nlat = length(ilat)
		frame$nlong = length(ilon)
		frame$npdg = length(frame$ind)
	} else {
		stopifnot(all(frame$nlong%/%nlon*nlon == frame$nlong))

		if (missing(ilat)) {
			stopifnot(frame$nlat%/%2*2 == frame$nlat)

			ilatn = seq(1,frame$nlat/2,nlat)
			ilat = c(ilatn,rev(frame$nlat-ilatn+1))
		}

		stopifnot(all(ilat %in% seq(frame$nlat)) && identical(ilat,sort(ilat)))

		clats = c(0,cumsum(frame$nlong))
		lats = datatoGauss(frame$lat,frame)
		longs = datatoGauss(frame$long,frame)

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

latPole = function(frame)
{
	n1 = frame$nlong[1]

	if (all(duplicated(frame$lat[1:n1])[-1])) {
		latcen = 90
	} else if (diff(range(frame$long[1:n1])) < 180) {
		# lat of Pole is so that Pole is out of 1st lat circle
		latcen = mean(frame$lat[c(1,n1/2+1)])
	} else {
		# lat of Pole is close to geo North Pole and is inside 1st lat circle
		# latcen = (latx+2*(90-latx)+latn)/2 = 90+(latn-latx)/2 = 90-(latx-latn)/2
		latcen = 90-abs(diff(frame$lat[c(1,n1/2+1)]))/2
	}

	latcen
}

lonPole = function(frame)
{
	n1 = frame$nlong[1]
	frame$long[n1/2+1]
}

compass = function(frame)
{
	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	l = m = numeric(frame$npdg)

	off = 0
	for (ilat in seq(frame$nlat)) {
		nlon = frame$nlong[ilat]
		cslon = 2*pi*seq(0,nlon-1)/nlon
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		gemu = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
		rcoslat = 1/sqrt(1-gemu^2)
		l[off+1:nlon] = -sqm2*sin(cslon)*rcoslat
		m[off+1:nlon] = ca1*rcoslat
		off = off+nlon
	}

	data.frame(l=l,m=m)
}

toU = function(nord,x,y)
{
	x*nord$m-y*nord$l
}

toV = function(nord,x,y)
{
	y*nord$m+x*nord$l
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

toDomain = function(frame,ind,dlat,dlon)
{
	if (! missing(ind)) {
		stopifnot(length(ind) == 1)
		stopifnot(length(dlat) == 2 && length(dlon) == 2)
		stopifnot(diff(dlat) > 0 && diff(dlon) > 0)

		xlim = frame$long[ind]+dlon
		ylim = frame$lat[ind]+dlat
	} else {
		xlim = longRange(frame$long)
		ylim = range(frame$lat)
	}

	list(xlim=xlim,ylim=ylim)
}

longRange = function(long)
{
	dlong1 = range(long%%360)
	dlong2 = range((long+180)%%360-180)

	if (diff(dlong1) > 180 && diff(dlong2) > 180)
		stop("longitudes span more than 180dg (case not supported)")

	if (diff(dlong1) <= diff(dlong2)) {
		return(dlong1)
	} else {
		return(dlong2)
	}
}

inDomain = function(frame,domain)
{
	ilat = domain$ylim[1] <= frame$lat & frame$lat <= domain$ylim[2]

	# case xlim=[0,0] (whole globe) is in alternative
	xlim = (domain$xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= frame$long & frame$long <= xlim[2]
	} else {
		ilon = frame$long <= xlim[2] | frame$long >= xlim[1]
	}

	ilat & ilon
}

area = function(dom1,dom2)
{
	diff(dom1$xlim)*diff(dom1$ylim)/(diff(dom2$xlim)*diff(dom2$ylim))
}

select = function(frame,ind)
{
	list(lat=frame$lat[ind],long=frame$long[ind])
}

zoom = function(data,frame,domain)
{
	mask = inDomain(frame,domain)
	ind = which(mask)

	lam = diff(domain$xlim) < 360
	list(lat=frame$lat[ind],long=frame$long[ind],lam=lam,data=data[ind,,drop=FALSE])
}

section = function(data,frame,long,lat=c(-90,90))
{
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

	# frame longs belong to [-180,180[
	long = (long+180)%%360-180

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

sectiongeo = function(data,frame,long=c(0,360),lat=c(-90,90))
{
	stopifnot(all(-90 <= lat & lat <= 90))

	# frame longs belong to [-180,180[
	if (length(long) == 1) {
		stopifnot(diff(lat) > 0)
		xf = frame$lat
		yf = frame$long
		x = lat
		y = (long+180)%%360-180
	} else if (length(lat) == 1) {
		stopifnot(diff(long) > 0)
		xf = frame$long
		yf = frame$lat
		x = (long+180)%%360-180
		y = lat
	} else {
		stop("lat or long must be of length 1")
	}

	clats = c(0,cumsum(frame$nlong))

	data1 = array(dim=c(2,frame$nlat,dim(data)[2]))
	x1 = array(dim=c(2,frame$nlat))

	for (ilat in seq(frame$nlat)) {
		off = clats[ilat]
		nlon = frame$nlong[ilat]
		yfn = yf[off+1:nlon]
		yi = y
		if (length(long) == 1) {
			# 1 more point in longitude (1st one) since frame is a Gaussian global grid
			yfn = c(yfn,yfn[1])
			il = which(abs(diff(yfn)) > 180)
			stopifnot(length(il) <= 2)
			for (i in rev(il)) yfn[-(1:i)] = yfn[-(1:i)]-360*sign(diff(yfn)[i])
			if (yi < min(yfn)) yi = yi+360

			stopifnot(all(diff(yfn) > -180))
		}

		if (FALSE && yfn[nlon]%%360 == yi) {
			data1[1,ilat,] = data[off+nlon,]
			x1[1,ilat] = xf[off+nlon]
			next
		}

		ind = yfn[-nlon] <= yi & yi < yfn[-1] | yfn[-1] <= yi & yi < yfn[-nlon]
		if (all(! ind)) next

		xfn = xf[off+1:nlon]
		i1 = which(ind)
		stopifnot(length(i1) <= 2)
		i2 = i1+1
		i2[i1 == nlon] = 1
		e = (y-yfn[i1])/(yfn[i2]-yfn[i1])
		stopifnot(all(0 <= e & e < 1))

		for (i in seq(along=i1)) {
			data1[i,ilat,] = (1-e[i])*data[off+i1[i],]+e[i]*data[off+i2[i],]
			x1[i,ilat] = (1-e[i])*xfn[i1[i]]+e[i]*xfn[i2[i]]
		}
	}

	x1 = as.vector(x1)
	if (all(is.na(x1))) stop("no lat/long crossing given parameter")

	ll = sort(x1,index.return=TRUE,na.last=TRUE)
	ind = ll$ix[1:length(na.omit(x1))]
	if (diff(x) > 0) {
		ii = which(x[1] <= x1[ind] & x1[ind] <= x[2])
	} else {
		ii = which(x[1] <= x1[ind] | x1[ind] <= x[2])
	}

	ind = ind[ii]

	dim(data1) = c(2*frame$nlat,dim(data)[2])

	if (length(long) == 1) {
		list(lats=x1[ind],long=long,data=data1[ind,])
	} else {
		list(longs=x1[ind],lat=lat,data=data1[ind,])
	}
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
		datag = datatoGauss(data[,j],frame)
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

plotxy = function(x,y,data,breaks="Sturges",palette=terrain.colors,pch=15,ppi=54,
	legend.title="",cex.axis=.9,...)
{
	u = par("usr")
	f = par("fin")

	if (ppi > 144) stop("ppi > 144\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(data) > npmax) {
		cat("--> reducing plot from",length(data),"to",npmax,"points\n")
		ind = seq(1,length(data),length.out=npmax)
		dn = diff(ind[1:2])%/%3
		if (dn > 2) {
			ind[seq(2,length(ind)-1,by=3)] = ind[seq(2,length(ind)-1,by=3)]+dn%/%2
			ind[seq(3,length(ind),by=3)] = ind[seq(3,length(ind),by=3)]-dn%/%2
		}

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

plotxy2 = function(x,y,zx,zy,breaks="Sturges",colvec=TRUE,length=.05,angle=15,ppi=6,
	legend.title="",cex.axis=.9,...)
{
	f = par("fin")

	if (ppi > 96) stop("ppi > 96\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(zx) > npmax) {
		cat("--> reducing plot from",length(zx),"to",npmax,"points\n")
		ind = seq(1,length(x),length.out=npmax)
		dn = diff(ind[1:2])%/%3
		if (dn > 2) {
			ind[seq(2,length(ind)-1,by=3)] = ind[seq(2,length(ind)-1,by=3)]+dn%/%2
			ind[seq(3,length(ind),by=3)] = ind[seq(3,length(ind),by=3)]-dn%/%2
		}

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
		cols = rep(1,length(ff1))
	}

	x2 = x+fx
	y2 = y+fy

	ffi = sqrt((fx/ux*f[1])^2+(fy/uy*f[2])^2)

	ind = ffi > 1.e-3
	arrows(x[ind],y[ind],x2[ind],y2[ind],length,angle,col=cols[ind],...)
	if (any(! ind)) points(x[! ind],y[! ind],pch=1,col=1)

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

mapdom = function(dom,frame,data,breaks="Sturges",palette=terrain.colors,pch=15,
	ppi=54,...)
{
	l = mapxy(dom,...)
	plotxy(frame$long,frame$lat,data,breaks,palette,pch,ppi)
	lines(l)
}

mapdom2 = function(dom,frame,datax,datay,breaks="Sturges",colvec=TRUE,length=.05,
	angle=15,ppi=6,...)
{
	l = mapxy(dom,...)
	plotxy2(frame$long,frame$lat,datax,datay,breaks,colvec,length,angle,ppi)
	lines(l)
}

mapex = function(data,frame,desc,dom=character())
{
	dlat = diff(range(frame$lat))

	if ("monde" %in% dom || dlat > 100 || diff(range(frame$lon)) > 200) {
		cat("plot domain world\n")
		xy = zoom(data,frame,monde)
		pngalt(sprintf("%s/mapworld_%s.png",pngd,desc$symbol))
		mapdom(monde,xy,xy$data[,1],main=desc$longname)
		pngoff()
	}

	if ("europe" %in% dom || length(which(inDomain(frame,europe))) > .7*frame$npdg) {
		cat("plot domain Europe\n")
		xy = zoom(data,frame,europe)
		pngalt(sprintf("%s/mapeuro_%s.png",pngd,desc$symbol))
		mapdom(europe,xy,xy$data[,1],main=sprintf("Field %s",desc$longname))
		pngoff()
	}

	if ("hima" %in% dom || length(which(inDomain(frame,hima))) > 100) {
		xy = zoom(data,frame,hima)
		cat("plot domain Himalaya\n")
		pngalt(sprintf("%s/maphima_%s.png",pngd,desc$symbol))
		mapdom(hima,xy,xy$data[,1],main=desc$longname)
		pngoff()
	}

	if ("france" %in% dom || length(which(inDomain(frame,france))) > 100) {
		xy = zoom(data,frame,france)
		cat("plot domain France\n")
		pngalt(sprintf("%s/mapfrance_%s.png",pngd,desc$symbol))
		mapdom(france,xy,xy$data[,1],main=desc$longname)
		pngoff()
	}

	#xz = sectionf(data,frame,30,c(65,110))
	#xzdf = as.data.frame(xz)
	#write.table(xzdf,"hima30N.txt",quote=FALSE,row.names=FALSE)

	#yz = section(data,frame,80,hima$ylim)
	#yzdf = as.data.frame(yz)
	#write.table(yzdf,"hima80E.txt",quote=FALSE,row.names=FALSE)
}

mapexi = function(data,frame,desc,doms,np=100,...)
{
	nl = dim(data)[2]
	dmean = dsd = drms = array(dim=c(nl,length(doms)))

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (names(doms)[i] == "monde") {
			# no domain, no zoom
			domin = c(doms[[i]],list(lam=FALSE))
			xy = list(lat=frame$lat,long=frame$long,data=data)
		} else {
			pts = select(frame,ind)
			domin = toDomain(pts)
			if (area(domin,doms[[i]]) < .1) next

			xy = zoom(data,frame,doms[[i]])
		}

		cat(". domain:",names(doms)[i],"\n")

		pngalt(sprintf("%s/map%s_%s.png",pngd,names(doms)[i],desc$symbol))
		if (nl == 1) {
			op = par(mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			mapdom(doms[[i]],xy,xy$data[,1],main=desc$longname,...)
		} else {
			op = par(mfcol=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			tt = c(desc$longname,sprintf("level %d",frame$ilev[1]))
			mapdom(doms[[i]],xy,xy$data[,1],main=tt,...)
			tt = c(desc$longname,sprintf("level %d",frame$ilev[nl]))
			mapdom(doms[[i]],xy,xy$data[,nl],main=tt,...)
			yz = sectiongeo(data,frame,long=mean(domin$xlim),lat=domin$ylim)
			tt = c(desc$longname,sprintf("S->N section, longitude %g",
				round(mean(domin$xlim),1)))
			plotv(yz$lats,etai,yz$data,main=tt,xlab="Latitude",ylab="eta")
			xz = sectiongeo(data,frame,lat=mean(domin$ylim),long=domin$xlim)
			tt = c(desc$longname,sprintf("W->E section, latitude %g",
				round(mean(domin$ylim),1)))
			plotv(xz$longs,etai,xz$data,main=tt,xlab="Longitude",ylab="eta")
		}

		pngoff(op)

		pngalt(sprintf("%s/hist%s_%s.png",pngd,names(doms)[i],desc$symbol))
		op = par(mfrow=c(min(nl,2),1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		hist(xy$data,col="whitesmoke",main=sprintf("Distribution of %s",desc$longname))

		if (nl > 1) {
			tt = c(sprintf("Distribution of %s",desc$longname),
				sprintf("level %d",frame$ilev[nl]))
			hist(xy$data[,nl],col="whitesmoke",main=tt)
		}

		pngoff(op)

		dmean[,i] = apply(xy$data,2,mean)
		dsd[,i] = apply(xy$data,2,sd)
		drms[,i] = apply(xy$data,2,function(x) sqrt(mean(x^2)))
	}

	list(mean=dmean,sd=dsd,rms=drms)
}

mapex2 = function(u,v,frame,doms,np=100,...)
{
	nl = dim(u)[2]
	dmean = dsd = drms = array(dim=c(nl,length(doms)))

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (names(doms)[i] == "monde") {
			# no domain, no zoom
			domin = c(doms[[i]],list(lam=FALSE))
			xy = list(lat=frame$lat,long=frame$long,data=data)
		} else {
			pts = select(frame,ind)
			domin = toDomain(pts)
			if (area(domin,doms[[i]]) < .1) next

			xy = zoom(data,frame,doms[[i]])
		}

		cat(". domain:",names(doms)[i],"\n")
		ff = sqrt(u^2+v^2)
		x = zoom(u,frame,doms[[i]])
		y = zoom(v,frame,doms[[i]])
		ffd = sqrt(x$data^2+y$data^2)

		pngalt(sprintf("%s/map%s_ff.png",pngd,names(doms)[i]))
		if (nl == 1) {
			op = par(mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			tt = "u/v wind"
			mapdom2(doms[[i]],ff,main=tt,...)
			mapdom2(doms[[i]],x,x$data,y$data,main=tt,...)
		} else {
			op = par(mfrow=c(min(nl,2),1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			tt = c("wind speed",sprintf("level %d",frame$ilev[1]))
			mapdom(doms[[i]],x,ffd[,1],main=tt,...)
			tt = c("wind speed",sprintf("level %d",frame$ilev[nl]))
			mapdom(doms[[i]],x,ffd[,nl],main=tt,...)
			yz = sectiongeo(ff,frame,long=mean(domin$xlim),lat=domin$ylim)
			tt = c("wind speed",sprintf("S->N section, longitude %g",
				round(mean(domin$xlim),1)))
			plotv(yz$lats,etai,yz$data,main=tt,xlab="Latitude",ylab="eta")
			xz = sectiongeo(ff,frame,lat=mean(domin$ylim),long=domin$xlim)
			tt = c("wind speed",sprintf("W->E section, latitude %g",
				round(mean(domin$ylim),1)))
			plotv(xz$longs,etai,xz$data,main=tt,xlab="Longitude",ylab="eta")
		}

		pngoff(op)

		pngalt(sprintf("%s/hist%s_ff.png",pngd,names(doms)[i]))
		op = par(mfrow=c(min(nl,2),1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		hist(ffd,col="whitesmoke",main="Distribution of wind speed")

		if (nl > 1) {
			pngalt(sprintf("%s/histn%s_ff.png",pngd,names(doms)[i]))
			tt = c("Distribution of wind speed",sprintf("level %d",frame$ilev[nl]))
			hist(ff[,nl],col="whitesmoke",main=tt)
		}

		pngoff(op)

		pngalt(sprintf("%s/map%s_uv.png",pngd,names(doms)[i]))
		op = par(mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		mapdom2(doms[[i]],x,x$data[,nl],y$data[,nl],main=tt,...)
		pngoff(op)

		dmean[,i] = apply(ff,2,mean)
		dsd[,i] = apply(ff,2,sd)
		drms[,i] = apply(ff,2,function(x) sqrt(mean(x^2)))
	}

	list(mean=dmean,sd=dsd,rms=drms)
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


Gficbin = tempfile(fileext=".bin")
Gxlab="<- Equator  ...  //  ...  France  ...  //  ...  North pole ->"

monde = list(xlim=c(-170,170),ylim=c(-85,85))
france = list(xlim=c(-6,11),ylim=c(40,53))
europe = list(xlim=c(-15,30),ylim=c(35,60))
asie = list(xlim=c(60,120),ylim=c(20,55))
hima = list(xlim=c(65,110),ylim=c(25,55))
amnord = list(xlim=c(-140,-50),ylim=c(25,70))
amsud = list(xlim=c(-85,-30),ylim=c(-60,10))
eurat = list(xlim=c(-30,35),ylim=c(20,70))
tropics = list(xlim=c(-50,80),ylim=c(-30,30))
austral = list(xlim=c(110,180),ylim=c(-50,-10))
poly = list(xlim=c(-155,-135),ylim=c(-25,-5))
doms = list(monde=monde,france=france,europe=europe,tropics=tropics,hima=hima,
	austral=austral,poly=poly)
ams = list(amnord=amnord,amsud=amsud,eurat=eurat,asie=asie)
stopifnot(all(sapply(doms,function(dom) diff(dom$xlim) > 0 && diff(dom$ylim) > 0)))
stopifnot(all(sapply(ams,function(dom) diff(dom$xlim) > 0 && diff(dom$ylim) > 0)))

desc = read.table(sprintf("%s/config/params.txt",diags),header=TRUE)

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

if (! "fic" %in% names(cargs)) {
	fics = scan("config/file.txt","character",quiet=TRUE)
} else if (regexpr("\\.txt$",cargs$fic) > 0) {
	fics = scan(cargs$fic,"character",quiet=TRUE)
} else {
	fics = cargs$fic
}

stopifnot(length(fics) > 0 && all(file.exists(fics)))

if (! "params" %in% names(cargs)) {
	params = scan("config/param.txt",what=character(),quiet=TRUE)
} else if (length(cargs$params) == 1 && file.exists(cargs$params)) {
	params = scan(cargs$params,what=character(),quiet=TRUE)
} else {
	params = cargs$params
}

ind = match(params,desc$faname)
if (any(is.na(ind))) stop("unknown parameters, see config/params.txt")

desc = desc[ind,]
npar = dim(desc)[1]

ilev = 0
if ("level" %in% names(cargs)) {
	if (length(cargs$level) == 1 && file.exists(cargs$level)) {
		ilev = scan(cargs$level,what=integer(),quiet=TRUE)
	} else {
		ilev = as.integer(cargs$level)
	}
} else if (file.exists("config/level.txt")) {
	ilev = scan("config/level.txt",what=integer(),quiet=TRUE)
}

cat("--> level indices/number (0: all levels):",ilev,"\n")

pngd = "."
if ("png" %in% names(cargs)) pngd = cargs$png

path = "."
if ("path" %in% names(cargs)) path = cargs$path

if ("date" %in% names(cargs)) {
	cbase = strsplit(cargs$date,split="/")[[1]]
	base = as.POSIXct(cbase[1],format="%Y%m%d")+as.integer(cbase[2])*3600
	cat("--> base date:",format(base),"\n")
}

ficana = "ARPE.0000.000"
#ficana[inds] = "analysis.surf-arpege.tl1798-c22.fa"
#ficana[! inds] = "analysis.atm-arpege.tl1798-c22.fa"

dstats = dstatd = array(NA_real_,c(length(fics),5,npar))
lstat = lstatd = lapply(seq(dim(desc)[1]),function(x) list())
lstatff = lapply(seq(dim(desc)[1]),function(x) list())

hasx11 = ! "png" %in% names(cargs) && capabilities("X11")
ask = ! "noask" %in% names(cargs) && hasx11 && interactive()
graph2 = FALSE

#if (interactive()) browser()

for (i in seq(along=fics)) {
	cat("File",basename(fics[i]),":\n")

	frame = getFrame(fics[i])
	cat("--> base/step:",as.character(frame$base,"%F %R %Z"),"+",frame$step/3600,"h\n")
	if (frame$nlat > 400) {
		cat("--> halving initial resolution:",frame$nlat,max(frame$nlong),"\n")
		frlow = degrade(frame,nlat=2,nlon=2)
	} else {
		frlow = frame

		if (! frame$lam) {
			mapf = dilat(frlow)
			ilateq = equalize(mapf,offset=3)
			if (length(ilateq) < frlow$nlat) frlow = degrade(frlow,ilateq)

			pngalt(sprintf("%s/stretch_%d.png",pngd,i))
			plot(mapf,type="l",main=c("Stretching coefficient by latitude",
				"normalized at Equator"),xlab="Latitude index",ylab="Stretching coef.",
				log="y")
			#points(ilateq,mapf[ilateq],pch=20,cex=.5)
			pngoff()
		}
	}

	eta = frlow$A/Gvp0+frlow$B
	stopifnot(all(abs(ilev) <= length(eta)))

	nlev = length(eta)
	if (length(ilev) == 1 && ilev < 0) {
		stopifnot(-ilev < nlev)

		cat("--> uniform eta selection of levels:",-ilev,"among",nlev,"\n")
		e = seq(min(eta),max(eta),length.out=-ilev)
		frlow$ilev = sapply(e,function(x) which.min(abs(x-eta)))
	} else if (length(ilev) < nlev) {
		stopifnot(all(ilev %in% seq(nlev)))

		cat("--> selection of levels:",ilev,"\n")
		frlow$ilev = ilev
	} else {
		frlow$ilev = seq(nlev)
	}

	etai = eta[frlow$ilev]
	yeta = rev(range(eta))

	date = frame$base+frame$step
	#date = base+frame$step
	ficref = sprintf("%s/%s/%s",path,strftime(date,"%H/%Y%m%d"),ficana)
	if (! file.exists(ficref))
		ficref = sprintf("%s/%s/%s",path,strftime(date,"%H/%Y%m%d"),fics[i])

	datax = datay = dataz = NULL
	dataox = dataoy = dataoz = NULL

	for (j in seq(npar)) {
		if (desc$ltype[j] == "-") {
			patt = desc$faname[j]
		} else if (length(frlow$ilev) > 1) {
			patt = sprintf("%s*%s",desc$ltype[j],desc$faname[j])
		} else if (desc$ltype[j] == "S") {
			patt = sprintf("S%03d%s",ilev,desc$faname[j])
		} else {
			cat("--> level type not 'S' for reading specified level\n")
			next
		}

		cat(". get fields matching pattern",patt,"\n")
		data = getField(fics[i],patt,desc$symbol[j],frame,frlow)

		nl = dim(data)[2]
		dstats[i,,j] = quantile(data[,nl],c(0,.1,.5,.9,1),names=FALSE)

		lstat[[j]][[i]] = mapexi(data,frlow,desc[j,],doms)

		if (FALSE) {
		if (! frame$lam) {
			ddx = xdiff(data,frlow)

			cat("Statistics of zonal difference:\n")
			print(summary(as.vector(ddx)))

			pngalt(sprintf("%s/spec45_%s.png",pngd,desc$symbol[j]))
			d45 = lat45(data,frlow)
			spectrum(d45[,nl],c(3,3),main="Smoothed spectral density of lat ~45")
			pngoff()
		}

		if (graph2 && nl > 2) {
			ix = arrayInd(which.max(abs(data)),dim(data))
			il = ix[1,2]
			ix = ix[1,1]
			cat("Level of max value:",il,"/",nl,"\n")
			mapex(data[,il,drop=F],frlow,desc[j,],names(doms))

			cat("Vertical profile at max value\n")
			dom = toDomain(frlow,ix,c(-20,20),c(-30,30))
			xy = zoom(data,frlow,dom)
			print(summary(as.vector(xy$data)))

			pngalt(sprintf("%s/profmax_%s.png",pngd,desc$symbol[j]))
			tt = "Vertical profile at max value"
			plotz(data[ix,],etai,ylim=yeta,main=c(tt,desc$longname[j]),
				xlab=desc$longname[j],ylab=expression(eta))
			pngoff()

			cat("Mean vertical profile\n")
			datam = apply(data,2,mean)

			pngalt(sprintf("%s/profmean_%s.png",pngd,desc$symbol[j]))
			tt = sprintf("Field %s",desc$longname[j])
			plotz(datam,etai,ylim=yeta,main=c("Mean vertical profile",tt),
				xlab=desc$longname[j],ylab=expression(eta))
			pngoff()

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

		dstatd[i,,j] = quantile(ddiff[,nl],c(0,.1,.5,.9,1),names=FALSE)

		cat("Statistics of forecast error:\n")
		print(summary(as.vector(ddiff)))
		drms = sqrt(mean(ddiff^2))
		cat("Global bias/RMSE:",mean(ddiff),drms,"\n")

		pngalt(sprintf("%s/histdiff_%s.png",pngd,desc$symbol[j]))
		tt = sprintf("Distribution of difference from analysis of %s",desc$longname[j])
		hist(ddiff,main=tt)
		pngoff()

		if (graph2) {
			cat("Statistics around diff min/max value:\n")
			i1 = arrayInd(which.min(ddiff),dim(ddiff))[1,1]
			i2 = arrayInd(which.max(ddiff),dim(ddiff))[1,1]
			dom = toDomain(frlow,i1,c(-20,20),c(-30,30))
			xy = zoom(ddiff,frlow,dom)
			print(summary(as.vector(xy$data)))
			drms = sqrt(mean(xy$data^2))
			cat("RMSE at local min:",drms,"\n")

			pngalt(sprintf("%s/diffmin_%s.png",pngd,desc$symbol[j]))
			mapdom(dom,frlow,ddiff[,nl],main=desc$longname[j])
			pngoff()

			dom = toDomain(frlow,i2,c(-20,20),c(-30,30))
			xy = zoom(ddiff,frlow,dom)
			print(summary(as.vector(xy$data)))
			drms = sqrt(mean(xy$data^2))
			cat("RMSE at local max:",drms,"\n")

			pngalt(sprintf("%s/diffmin_%s.png",pngd,desc$symbol[j]))
			mapdom(dom,frlow,ddiff[,nl],main=desc$longname[j])
			pngoff()
		}

		lstatd[[j]][[i]] = mapexi(ddiff,frlow,desc[j,],doms)

		if (FALSE) {
		if (graph2 && nl > 2) {
			pngalt(sprintf("%s/profile_%s.png",pngd,desc$symbol[j]))
			par(mfrow=c(1,2))
			tt = sprintf("Field %s",desc$longname[j])
			plotz(data[i1,],etai,ylim=yeta,main=c(tt,"Profile at point 'min of error'"),
				xlab=desc$longname[j],ylab=expression(eta))
			plotz(data[i2,],etai,ylim=yeta,main=c(tt,"Profile at point 'max of error'"),
				xlab=desc$longname[j],ylab=expression(eta))
			pngoff()

			ns = length(yzs$lat)
			dataf = yzn$data
			dataf[1:ns,] = yzs$data[ind,]
			#levels = prettydiff(dataf,7)

			pngalt(sprintf("%s/vsecdiff_%s.png",pngd,desc$symbol[j]),720)
			x = seq(0,1,length.out=dim(dataf)[1])
			plotv(x,etai,dataf,nlevels=7,main=tt,xlab=Gxlab,ylab="eta",xaxs="i",yaxs="i")
			pngoff()

			pngalt(sprintf("%s/vsecfill_%s.png",pngd,desc$symbol[j]),720)
			plotf(etai,dataf,nlevels=7,color.palette=cm.colors,main=tt,xlab=Gxlab,ylab="eta")
			pngoff()
		}

		dbias = apply(ddiff,2,mean)
		drms = apply(ddiff,2,function(x) sqrt(mean(x^2)))
		cat("Global bias/RMSE by level:\n")
		print(summary(dbias))
		print(summary(drms))
		}
	}

	if (! is.null(datax) && ! is.null(datay)) {
		cat(". compound field wind speed\n")
		ff = sqrt(datax^2+datay^2)
		pngalt(sprintf("hist_ff.png"))
		hist(ff,col="whitesmoke",main="Distribution of wind speed")
		pngoff()

		print(summary(as.vector(ff)))

		if (frame$lam || latPole(frame) == 90) {
			u = datax
			v = datay
		} else {
			cat("rotating wind to geo North\n")
			nord = compass(frame)
			if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
			u = toU(nord,datax,datay)
			v = toV(nord,datax,datay)

			pole = list(xlim=lonPole(frame)+c(-5,5),ylim=latPole(frame)+c(-5,5))
			x = zoom(u,frlow,pole)
			y = zoom(v,frlow,pole)
			pngalt(sprintf("mappole_uv.png"))
			mapdom2(pole,x,x$data[,nl],y$data[,nl],main="u/v wind components")
			pngoff()
		}

		lstatff[[i]] = mapex2(u,v,frlow,doms)
	}

	if (! is.null(dataox) && ! is.null(dataoy)) {
		ffo = sqrt(uo^2+vo^2)
		pngalt(sprintf("hist_ffo.png"))
		hist(ffo,col="whitesmoke",main="Distribution of wind speed")
		pngoff()

		print(summary(as.vector(ffo)))

		if (! framo$lam && latPole(framo) < 90) {
			cat("rotating wind to geo North\n")
			nord = compass(framo)
			if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
			uo = toU(nord,dataox,dataoy)
			vo = toV(nord,dataox,dataoy)

			for (id in seq(along=doms)) {
				if (length(which(inDomain(framo,doms[[id]]))) <= 100) next

				cat("plot domain",names(doms)[id],"\n")
				x = zoom(uo,frlow,doms[[id]])
				y = zoom(vo,frlow,doms[[id]])
				pngalt(sprintf("%s/map%s_uvo.png",pngd,names(doms)[id]))
				mapdom2(doms[[id]],x,x$data[,nl],y$data[,nl],main="u/v wind components")
				pngoff()
			}
		}
	}

	if (! is.null(datax) && ! is.null(datay) && ! is.null(dataox) && ! is.null(dataoy)) {
		ffd = ff-ffo
		print(summary(as.vector(ffd)))
		drms = sqrt(mean(ffd^2))
		cat("Global bias/RMSE:",mean(ffd),drms,"\n")

		pngalt(sprinf("histdiff_ff.png",pngd))
		hist(ffd,col="whitesmoke",main="Distribution of difference from ref of wind speed")
		pngoff()

		xy = zoom(ffd,frlow,europe)
		pngalt(sprintf("%s/mapeuro_ffd.png",pngd))
		mapdom(europe,xy,xy$data[,1],main="Difference of wind speed")
		pngoff()
	}
}

if (length(fics) > 1) {
	tt = sprintf("Forecast of %s",desc$longname)

	nf = dim(desc)
	nr = 3

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/stats%d.png",pngd,i))
		op = par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			lines(ftime,dstats[,3,j],lty=1,col=1,ylim=range(dstats[,,j]))
			boxplot(dstatd[,,j],col="bisque",xlab="Lead time",ylab=desc$symbol[j],
				main=tt[j],at=ftime)
		}

		pngoff(op)
	}
}

if (any(! is.na(dstatd))) {
	tt = sprintf("Error of %s",desc$longname)

	nf = dim(desc)
	nr = 3

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/err%d.png",pngd,i))
		op = par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			lines(dstatd[,3,j],lty=1,col=1,ylim=range(dstatd[,,j]))
			boxplot(dstatd[,,j],col="bisque",xlab="Lead time",ylab=desc$symbol[j],
				main=tt[j],at=ftime)
		}

		pngoff(op)
	}

	ndom = dim(lstatd[[1]][[1]]$mean)[2]
	nr = 3

	for (j in seq(along=lstatd)) {
		m = sapply(lstatd[[j]],function(x) simplify2array(x$mean))
		s = sapply(lstatd[[j]],function(x) simplify2array(x$sd))
		r = sapply(lstatd[[j]],function(x) simplify2array(x$rms))
		dim(m) = dim(s) = dim(r) = c(length(frlow$ilev),ndom,length(fics))

		tt = sprintf("Score of %s",desc$longname[j])

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/score%s_%s.png",pngd,names(doms)[i],desc$shortname[j]))
			op = par(mfrow=c(nr,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				tt[2] = sprintf("domain %s",names(doms)[i])
				ds = data.frame(r[1,id,],m[1,id,],s[1,id,])
				matplot(ds,type="o",lty=1:3,col=1,pch=20,xlab="Lead time",ylab=desc$symbol[j],
					main=tt)
				legend("topleft",c("RMS","std dev","bias"),col=1,lty=1:3,pch=20,
					bg="transparent")

				tt[2] = sprintf("domain %s, level %d",names(doms)[i],frlow$ilev[nl])
				ds = data.frame(r[nl,id,],m[nl,id,],s[nl,id,])
				matplot(ds,type="o",lty=1:3,col=1,pch=20,xlab="Lead time",ylab=desc$symbol[j],
					main=tt)
				legend("topleft",c("RMS","std dev","bias"),col=1,lty=1:3,pch=20,
					bg="transparent")
			}

			pngoff(op)
		}
	}
}
