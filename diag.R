diags = "~/util/diags"

library(maps)

Gvp0 = 101325

readDom = function(filename,header=TRUE,...)
{
	df = read.table(filename,header,...)

	stopifnot(all(df$east > df$west & df$north > df$south))
	stopifnot(all(-90 < df$north & df$north <= 90))
	stopifnot(all(-90 <= df$south & df$south < 90))

	doms = lapply(seq(dim(df)[1]),
		function(i) list(xlim=c(df$west[i],df$east[i]),ylim=c(df$south[i],df$north[i])))

	names(doms) = df$name
	doms
}

getGrid = function(con)
{
	nx = readBin(con,"integer",1,endian="swap")
	type = readBin(con,"integer",nx/4,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	# int is 4 bytes, num is nfp bytes (= 4 or 8)
	nfp = type[1]

	# 0: global (Gauss) grid, 1: LAM grid
	if (type[2] == 0) {
		nx = readBin(con,"integer",1,endian="swap")
		gem = readBin(con,"numeric",nx/nfp,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		ndgnh = nx/4
		nloeng = readBin(con,"integer",ndgnh,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		ndglg = 2*ndgnh
		nloeng[(ndgnh+1):ndglg] = rev(nloeng)

		grid = list(nfp=nfp,mucen=gem[1],locen=gem[2],rstret=gem[3],nlat=ndglg,nlong=nloeng,
			npdg=sum(nloeng),gauss=TRUE,lam=FALSE)
	} else if (type[2] == 1) {
		ndglg = type[3]
		ndlon = type[4]
		grid = list(nfp=nfp,nlat=ndglg,nlong=rep(ndlon,ndglg),npdg=ndglg*ndlon,
			gauss=FALSE,lam=TRUE)
	} else {
		stop("unknown grid type")
	}

	nx = readBin(con,"integer",1,endian="swap")
	nl1 = nx/nfp
	Ah = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))
	nx = readBin(con,"integer",1,endian="swap")
	Bh = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	grid$eta = ((Ah[-1]+Ah[-nl1])/Gvp0+(Bh[-1]+Bh[-nl1]))/2
	grid$nlevel = nl1-1

	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/4 == 2)
	base = readBin(con,"integer",2,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/nfp == 1)
	grid$step = readBin(con,"numeric",1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	grid$base = as.POSIXct(as.Date(as.character(base[1]),"%Y%m%d"))+base[2]

	grid
}

is.equal = function(grid,frame)
{
	if (grid$npdg != frame$npdg || grid$nlevel != frame$nlevel ||
		grid$nlat != frame$nlat) return(FALSE)
	if (! identical(grid$nlong,frame$nlong)) return(FALSE)

	if (grid$gauss) {
		if (! isTRUE(all.equal(grid$mucen,sin(pi/180*latPole(frame)))) ||
			! isTRUE(all.equal(grid$locen,pi/180*lonPole(frame)))) return(FALSE)
	}

	grid$gauss == frame$gauss
}

toFrame = function(grid,frame)
{
	grid = c(grid,frame[c("lat","long","base")])
	if (grid$gauss) grid$theta = frame$theta

	grid
}

getVars = function(con,grid,vars)
{
	lvar = vector("list",length(vars))

	while (TRUE) {
		nx = readBin(con,"integer",1,endian="swap")
		if (length(nx) == 0) break

		s = readChar(con,nx)
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		stopifnot(nx == 4)
		nl = readBin(con,"integer",1,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		stopifnot(nx/grid$nfp == nl*grid$npdg)

		i = match(sub("^\\.?(\\w+) +$","\\1",s),sub("^\\.?","",vars))
		if (is.na(i)) {
			seek(con,nl*grid$npdg*grid$nfp+4,"current")
			next
		}

		lvar[[i]] = readBin(con,"numeric",nl*grid$npdg,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))
		dim(lvar[[i]]) = c(grid$npdg,nl)
		if (all(! sapply(lvar,is.null))) break
	}

	lvar
}

getFrame = function(file)
{
	system(sprintf("epy_dump.py %s -f frame -o %s",file,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	dims = readBin(con,"integer",5,size=8)

	# global (Gauss grid) or LAM (Cartesian grid)
	if (dims[5] == 0) {
		nlong = readBin(con,"integer",dims[1],size=8)
		nwave = readBin(con,"integer",dims[2],size=8)
		ngp = dims[4]
		nl = dims[3]
	} else {
		ngp = prod(dims[1:2])
		nl = dims[5]
	}

	base = readBin(con,"numeric",1)
	base = as.POSIXct(base,origin="1970-01-01")
	step = readBin(con,"integer",1,size=8)
	if (regexpr(".+\\+0*([0-9]+)",file) > 0) {
		ech = as.integer(gsub(".+\\+0*([0-9]+)","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	lats = readBin(con,"numeric",ngp)
	longs = readBin(con,"numeric",ngp)

	Ah = readBin(con,"numeric",nl+1)
	Bh = readBin(con,"numeric",nl+1)
	eta = ((Ah[-1]+Ah[-(nl+1)])/Gvp0+(Bh[-1]+Bh[-(nl+1)]))/2

	close(con)

	if (dims[5] == 0) {
		nlat = length(nlong)
		theta = 90-180*(seq(nlat)-.5)/nlat

		frame = list(nlat=nlat,nwave=dims[2],nlevel=nl,npdg=ngp,nlong=nlong,theta=theta,
			lat=lats,long=longs,eta=eta,base=base,step=step,gauss=TRUE,lam=FALSE)
	} else {
		nlong = rep(dims[2],dims[1])
		frame = list(nlat=dims[1],nlong=nlong,nwavex=dims[3],nwavey=dims[4],nlevel=nl,
			npdg=ngp,lat=lats,long=longs,eta=eta,base=base,step=step,gauss=FALSE,lam=TRUE)
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

		cat("--> different grid or levels, read data again\n")
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

	ilev = frlow$ilev

	if (dim(data)[2] > 1) {
		if (frame$nlevel != frlow$nlevel) data = interpAB(data,frame,frlow)

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
	nl = readBin(con,"integer",1,size=8)
	stopifnot(length(nl) == 1)

	data = matrix(nrow=frame$npdg,ncol=nl)
	for (j in seq(nl)) data[,j] = readBin(con,"numeric",frame$npdg)
	stopifnot(nl == readBin(con,"integer",1,size=8))
	close(con)

	data
}

getSPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",1,size=8)
	data = readBin(con,"numeric",nsp2)
	close(con)

	stopifnot(nsp2 == frame$nwave*(frame$nwave-1))

	matrix(data,nrow=frame$nwave)
}

degrade = function(frame,ilat,nlat,nlon=2)
{
	stopifnot(all(frame$nlong%%nlon == 0))

	if (frame$gauss) {
		if (missing(ilat)) {
			stopifnot(frame$nlat%%2 == 0)

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
	} else {
		stopifnot(all(duplicated(frame$nlong)[-1]))
		ndlon = frame$nlong[1]

		ilat = seq(1,frame$nlat,nlat)
		ilon = seq(1,ndlon,nlon)
		frame$ind = rep((ilat-1)*ndlon,each=length(ilon))+ilon
		frame$lat = frame$lat[frame$ind]
		frame$long = frame$long[frame$ind]
		frame$nlat = length(ilat)
		frame$nlong = rep(length(ilon),frame$nlat)
		frame$npdg = length(frame$ind)
	}

	frame
}

csLat = function(frame)
{
	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	lat = numeric(frame$nlat)

	ind = frame$lat[cumsum(frame$nlong)-frame$nlong/2]
	cslon = pi

	for (ilat in seq(frame$nlat)) {
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		lat[ilat] = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
	}

	lat
}

csLong = function(frame)
{
	longs = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (i in seq(frame$nlat)) {
		nl = frame$nlong[i]
		longs[seq(nl),i] = 360*(seq(nl)-1)/nl
	}

	longs
}

datatoGauss = function(data,frame)
{
	clats = c(0,cumsum(frame$nlong))
	datag = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (ilat in seq(frame$nlat)) {
		ip = clats[ilat]
		np = frame$nlong[ilat]
		datag[1:np,ilat] = data[ip+1:np]
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
	x*nord$m+y*nord$l
}

toV = function(nord,x,y)
{
	y*nord$m-x*nord$l
}

ucs = function(nord,x,y)
{
	x*nord$m-y*nord$l
}

vcs = function(nord,x,y)
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

interpGauss = function(datao,framo,frame)
{
	stopifnot(frame$gauss)

	# only for (Gauss) grids nested inside
	stopifnot(length(framo$nlong) >= length(frame$nlong))
	stopifnot(max(framo$nlong) >= max(frame$nlong))

	longs = csLong(frame)

	data = matrix(nrow=frame$npdg,ncol=dim(datao)[2])

	clats = c(0,cumsum(frame$nlong))
	clato = c(0,cumsum(framo$nlong))
	dlat = 180/framo$nlat
	dlon = 360/framo$nlong

	for (ip in seq(frame$npdg)) {
		ilat = which.max(ip <= clats[-1])
		ilon = ip-clats[ilat]
		e0 = (90-dlat/2-frame$theta[ilat])/dlat
		ilato = floor(e0)+1
		e0 = e0-(ilato-1)
		stopifnot(0 < ilato && ilato < length(framo$nlong))
		stopifnot(0 <= e0 && e0 < 1)

		e1 = longs[ilon,ilat]/dlon[ilato]
		e2 = longs[ilon,ilat]/dlon[ilato+1]
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

xdiff = function(data,frame)
{
	stopifnot(frame$gauss)

	dlon = cos(frame$theta*pi/180)*360/frame$nlong
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
	stopifnot(frame$gauss)

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
	stopifnot(frame$gauss)

	dlon = cos(frame$theta*pi/180)*360/frame$nlong
	gm = dlon[frame$nlat%/%2]/dlon

	clats = c(0,cumsum(frame$nlong))

	i = frame$nlat%/%4
	ip = clats[i]+seq(frame$nlong[i])
	data[ip,,drop=FALSE]
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

makeDomain = function(frame,ind,dlat,dlon)
{
	stopifnot(length(ind) == 1)
	stopifnot(length(dlat) == 2 && length(dlon) == 2)
	stopifnot(diff(dlat) > 0 && diff(dlon) > 0)

	xlim = frame$long[ind]+dlon
	ylim = frame$lat[ind]+dlat
}

getDomain = function(frame)
{
	xlim = longRange(frame$long)
	ylim = range(frame$lat)

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

inDomain = function(points,domain)
{
	ilat = domain$ylim[1] <= points$lat & points$lat <= domain$ylim[2]

	# case xlim=[0,0] (whole globe) is in alternative
	xlim = (domain$xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= points$long & points$long <= xlim[2]
	} else {
		ilon = points$long <= xlim[2] | points$long >= xlim[1]
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

	lam = frame$lam || diff(domain$xlim) < 360
	list(lat=frame$lat[ind],long=frame$long[ind],lam=lam,data=data[ind,,drop=FALSE])
}

section = function(data,frame,long,lat=c(-90,90))
{
	stopifnot(frame$gauss)
	stopifnot(all(-90 <= lat & lat <= 90))

	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	if (length(lat) == 1) {
		dlat = 180/frame$nlat
		ilat = which.max(frame$theta <= lat)
		e0 = (90-dlat/2-frame$theta[ilat])/dlat-ilat+1
		ind = c(ilat,ilat+1)
	} else {
		ind = which(min(lat) <= frame$theta & frame$theta <= max(lat))
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

	list(lat=frame$theta[ind],lats=lats,longs=longs,data=data1)
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
		if (frame$gauss) {
			# 1 more point in longitude (1st one) since frame is a Gaussian global grid
			yfn = c(yfn,yfn[1])
			il = which(abs(diff(yfn)) > 180)
			stopifnot(length(il) <= 2)
			for (i in rev(il)) yfn[-(1:i)] = yfn[-(1:i)]-360*sign(diff(yfn)[i])
			if (yi < min(yfn)) yi = yi+360

			stopifnot(all(diff(yfn) > -180))
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
	stopifnot(frame$gauss)

	long = unique(long)
	stopifnot(-90 <= lat && lat < 90)
	stopifnot(length(long) == 2)
	stopifnot(all(0 <= long & long <= 360))

	longs = csLong(frame)
	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	ilat = max(which(frame$theta >= lat))

	if (long[1] < long[2]) {
		ind = which(long[1] <= longs[,ilat] & longs[,ilat] <= long[2])
	} else {
		ind = c(which(long[1] <= longs[,ilat] & longs[,ilat] <= 360),
			which(0 <= longs[,ilat] & longs[,ilat] <= long[2]))
	}

	data1 = matrix(nrow=length(ind),ncol=dim(data)[2])

	for (j in seq(dim(data)[2])) {
		datag = datatoGauss(data[,j],frame)
		data1[,j] = datag[ind,ilat]
	}

	lats = longs = numeric(length(ind))

	e = (lat-frame$theta[ilat])/diff(frame$theta[ilat+0:1])

	ip = clats[ilat]+ind

	list(long=longs[ind,ilat],lats=frame$lat[ip],longs=frame$long[ip],data=data1)
}

mapxy = function(xy,axes=TRUE,frame.plot=TRUE,...)
{
	l = map(xlim=xy$xlim,ylim=xy$ylim,...)

	if (axes) {
		axis(1)
		axis(2)
	}

	if (frame.plot) box()

	l
}

plotxy = function(x,y,data,breaks="Sturges",palette=terrain.colors,pch=15,ppi=50,
	cex.leg=.8,...)
{
	u = par("usr")
	f = par("fin")

	if (ppi > 144) stop("ppi > 144\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(data) > npmax) {
		cat("--> reducing xy plot from",length(data),"to",npmax,"points\n")
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
	if (is.character(breaks) && length(br) > 8) {
		h = hist(data,breaks=8,plot=FALSE)
		br = h$breaks
	}

	ind = findInterval(data,br)
	cols = palette(length(br))

	points(x,y,col=cols[ind],pch=pch,...)

	levels = sprintf("% .3g",br)
	DOPfill.legend(levels,col=cols,cex.axis=cex.leg)
}

plotxy2 = function(x,y,zx,zy,breaks="Sturges",colvec=TRUE,length=.05,angle=15,ppi=8,
	cex.leg=.8,...)
{
	f = par("fin")

	if (ppi > 96) stop("ppi > 96\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(zx) > npmax) {
		cat("--> reducing xy2 plot from",length(zx),"to",npmax,"points\n")
		ind = as.integer(seq(1,length(x),length.out=npmax))
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

	# scale of wind vectors: a fraction of the smallest side
	asp = min(ux,uy)/max(ux,uy)
	ffu = asp*sqrt((ux^2+uy^2)/length(zx))

	if (colvec) {
		fx = zx/ffz*ffu
		fy = zy/ffz*ffu

		h = hist(ffz,breaks,plot=FALSE)
		br = h$breaks
		if (is.character(breaks) && length(br) > 8) {
			h = hist(ffz,breaks=8,plot=FALSE)
			br = h$breaks
		}

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

	if (colvec) {
		levels = sprintf("% .3g",br)
		DOPfill.legend(levels,col=palette,cex.axis=cex.leg)
	}
}

DOPfill.legend = function(levels,col,...)
{
	u = par("usr")
	p = par("plt")

	width = (1-p[2])/6*diff(u[1:2])/diff(p[1:2])
	x = u[2]+width/3

	nl = length(levels)
	height = diff(u[3:4])
	dy = height/(nl-1)
	y = u[3]
	ybas = y + dy*(seq(nl-1)-1)
	yhaut = ybas + dy

	rect(x,ybas,x+width,yhaut,col=col,border=NA,xpd=TRUE)

	op = par(las=2,yaxt="s")
	axis(4,c(ybas[1],yhaut),levels,tick=FALSE,pos=x,mgp=c(1,.5,0),...)
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

plotv0 = function(x,y,z,ylim=NULL,xaxt="s",yaxt="s",...)
{
	if (! is.null(ylim)) {
		indy = which(ylim[1] <= y & y <= ylim[2])
		y = y[indy]
		indy = rev(indy)
	} else {
		indy = rev(seq(along=y))
	}

	contour(x,y,z[,indy],xaxt="n",yaxt="n",...)

	if (xaxt != "n") {
		lab = pretty(x/30)*30
		axis(1,at=lab)
	}

	if (yaxt != "n") {
		lab = pretty(y)
		axis(2,revpos(y,lab),lab)
	}
}

plotv = function(x,y,z,nlevels=8,ylim=rev(range(y,finite=TRUE)),xaxs="i",yaxs="i",...)
{
	contour(x,y,z,nlevels=nlevels,ylim=ylim,xaxs=xaxs,yaxs=yaxs,...)
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

mapdom = function(dom,points,data,axes=TRUE,main=NULL,mar=par("mar"),...)
{
	l = mapxy(dom,axes=axes,mar=mar)
	plotxy(points$long,points$lat,data,...)
	lines(l)
	title(main)
}

mapdom2 = function(dom,points,datax,datay,axes=TRUE,main=NULL,mar=par("mar"),...)
{
	l = mapxy(dom,axes=axes,mar=mar)
	plotxy2(points$long,points$lat,datax,datay,...)
	lines(l)
	title(main)
}

mapexi = function(data,frame,desc,doms,np=4,prefix="",graph=TRUE,...)
{
	nl = dim(data)[2]
	qv = array(dim=c(101,nl,length(doms)))

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (frame$gauss && names(doms)[i] == "monde") {
			# no domain, no zoom
			domin = c(doms[[i]],list(lam=FALSE))
			xy = list(lat=frame$lat,long=frame$long,data=data)
		} else {
			pts = select(frame,ind)
			domin = getDomain(pts)
			if (area(domin,doms[[i]]) < .1) next

			xy = zoom(data,frame,doms[[i]])
		}

		qv[,,i] = apply(xy$data,2,quantile,prob=seq(0,100)/100)

		if (! graph) next

		nom = names(doms)[i]
		tt = desc$longname
		mar = c(1,1,2,5)+.1
		cat(".. domain:",nom,"\n")

		pngalt(sprintf("%s/map%s%s_%s.png",pngd,prefix,nom,desc$symbol))
		if (nl == 1) {
			op = par(mar=mar,mgp=c(2,.75,0))
			mapdom(doms[[i]],xy,xy$data[,1],main=tt,...)
		} else {
			long = round(mean(domin$xlim),1)
			lat = round(mean(domin$ylim),1)

			op = par(mfcol=c(2,2),mar=mar,mgp=c(2,.75,0))

			for (il in c(1,nl)) {
				tt[2] = sprintf("level %d",frame$ilev[il])
				mapdom(doms[[i]],xy,xy$data[,il],main=tt,mar=mar,...)
				abline(h=lat,col="darkgrey",lty=2)
				abline(v=long,col="darkgrey",lty=2)
			}

			yz = sectiongeo(data,frame,long=mean(domin$xlim),lat=domin$ylim)
			tt[2] = sprintf("S->N section, longitude %g",long)
			plotv(yz$lats,etai,yz$data,main=tt,xlab="Latitude",ylab="eta")
			xz = sectiongeo(data,frame,lat=mean(domin$ylim),long=domin$xlim)
			tt[2] = sprintf("W->E section, latitude %g",lat)
			plotv(xz$longs,etai,xz$data,main=tt,xlab="Longitude",ylab="eta")
		}

		pngoff(op)

		pngalt(sprintf("%s/hist%s%s_%s.png",pngd,prefix,nom,desc$symbol))
		nr = min(nl,3)
		if (nl == 1) {
			op = par(mfrow=c(nr,1),mar=c(3,3,3,1)+.1)
			hist(xy$data[,1],col="whitesmoke",main=tt,xlab=desc$longname)
		} else {
			op = par(mfrow=c(nr,2),mar=c(3,3,3,1)+.1)
			indl = c(1,1+nl%/%2,nl)
			for (il in indl[1:nr]) {
				mapdom(doms[[i]],xy,xy$data[,il],mar=c(0,0,0,5)+.1,axes=FALSE,...)
				abline(h=lat,col="darkgrey",lty=2)
				abline(v=long,col="darkgrey",lty=2)
				tt[2] = sprintf("level %d",frame$ilev[il])
				hist(xy$data[,il],col="whitesmoke",main=tt,xlab=desc$longname)
			}
		}

		pngoff(op)
	}

	qv
}

mapex2 = function(u,v,frame,doms,np=4,prefix="",graph=TRUE,...)
{
	nl = dim(u)[2]
	qv = array(dim=c(101,nl,length(doms)))

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (names(doms)[i] == "monde") {
			# no domain, no zoom
			domin = c(doms[[i]],list(lam=FALSE))
			x = list(lat=frame$lat,long=frame$long,data=u)
			y = list(lat=frame$lat,long=frame$long,data=v)
		} else {
			pts = select(frame,ind)
			domin = getDomain(pts)
			x = zoom(u,frame,doms[[i]])
			y = zoom(v,frame,doms[[i]])
		}

		ffdom = sqrt(x$data^2+y$data^2)
		qv[,,i] = apply(ff,2,quantile,prob=seq(0,100)/100)

		if (! graph || area(domin,doms[[i]]) < .1) next

		nom = names(doms)[i]
		tt = "u/v wind"
		cat(".. domain:",nom,"\n")
		ff = sqrt(u^2+v^2)
		mar = c(1,1,2,5)+.1

		pngalt(sprintf("%s/map%s%s_ff.png",pngd,prefix,nom))
		if (nl == 1) {
			op = par(mar=mar,mgp=c(2,.75,0))
			mapdom2(doms[[i]],x,x$data,y$data,main=tt,...)
		} else {
			long = round(mean(domin$xlim),1)
			lat = round(mean(domin$ylim),1)

			op = par(mfcol=c(2,2),mar=mar,mgp=c(2,.75,0))

			for (il in c(1,nl)) {
				tt[2] = sprintf("level %d",frame$ilev[il])
				mapdom2(doms[[i]],x,x$data[,il],y$data[,il],main=tt,mar=mar,...)
				abline(h=lat,col="darkgrey",lty=2)
				abline(v=long,col="darkgrey",lty=2)
			}

			yz = sectiongeo(ff,frame,long=mean(domin$xlim),lat=domin$ylim)
			tt[2] = sprintf("S->N section, longitude %g",long)
			plotv(yz$lats,etai,yz$data,main=tt,xlab="Latitude",ylab="eta")
			xz = sectiongeo(ff,frame,lat=mean(domin$ylim),long=domin$xlim)
			tt[2] = sprintf("W->E section, latitude %g",lat)
			plotv(xz$longs,etai,xz$data,main=tt,xlab="Longitude",ylab="eta")
		}

		pngoff(op)

		pngalt(sprintf("%s/hist%s%s_ff.png",pngd,prefix,nom))
		tt = "wind speed"
		nr = min(nl,3)
		if (nl == 1) {
			op = par(mfrow=c(nr,1),mar=c(3,3,3,1)+.1)
			hist(ffdom[,il],col="whitesmoke",main=tt,xlab="wind speed")
		} else {
			op = par(mfrow=c(nr,2),mar=c(3,3,3,1)+.1)
			indl = c(1,1+nl%/%2,nl)
			for (il in indl[1:nr]) {
				mapdom2(doms[[i]],x,x$data[,il],y$data[,il],axes=FALSE,mar=mar,...)
				abline(h=lat,col="darkgrey",lty=2)
				abline(v=long,col="darkgrey",lty=2)
				tt[2] = sprintf("level %d",frame$ilev[il])
				hist(ffdom[,il],col="whitesmoke",main=tt,xlab="wind speed")
			}
		}

		pngoff(op)
	}

	qv
}

rms = function(x,...)
{
	c(sqrt(mean(x^2,...)),mean(x,...),sd(x,...))
}

plotb = function(x,y,...)
{
	plot(x,y[3,],type="l",ylim=range(y,na.rm=TRUE),...)
	w = diff(range(x))/dim(y)[1]
	boxplot(y,range=0,col="bisque",lty=1,pars=list(boxwex=.15*w),add=TRUE,at=x,xaxt="n",
		yaxt="n")
}

plotbv = function(x,y,...)
{
	plot(x[3,],y,type="l",xlim=range(x,na.rm=TRUE),...)
	w = diff(range(y))/dim(x)[1]
	boxplot(x,range=0,col="bisque",lty=1,pars=list(boxwex=.15*w),horizontal=TRUE,add=TRUE,
		at=y,xaxt="n",yaxt="n")
}

matplott = function(x,y,...)
{
	matplot(x,y,type="o",lty=1:3,col=1,pch="+",...)
	abline(h=0,col="grey")
	legend("topleft",c("RMS","bias","std dev"),col=1,lty=1:3,pch="+",bg="transparent")
}

matplotv = function(x,y,...)
{
	matplot(x,y,type="o",lty=1:3,col=1,pch="-",...)
	abline(v=0,col="grey")
	legend("topleft",c("RMS","bias","std dev"),col=1,lty=1:3,pch="+",bg="transparent")
}

pngalt = function(...)
{
	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))

	if (! hasx11) {
		stopifnot(is.null(dev.list()))
		png(...)
	}
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

if (file.exists("config/domain.txt")) {
	doms = readDom("config/domain.txt")
} else {
	doms = readDom(sprintf("%s/config/domain.txt",diags))
}

descall = read.table(sprintf("%s/config/params.txt",diags),header=TRUE)

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

if ("fic" %in% names(cargs)) {
	dff = read.table(cargs$fic,header=TRUE)
} else {
	dff = read.table("config/file.txt",header=TRUE)
}

dates = NA

if (file.exists("config/date.txt")) {
	dates = scan("config/date.txt",quiet=TRUE)
	res = as.integer(dates[1])
	dates = dates[-1]
} else if ("date" %in% names(cargs)) {
	cdate = strsplit(cargs$date,split="/")[[1]]
	dates = as.integer(cdate[1])
	res = as.integer(cdate[2])
}

if (! "params" %in% names(cargs)) {
	params = scan("config/param.txt",what=character(),quiet=TRUE,comment.char="#")
} else if (length(cargs$params) == 1 && file.exists(cargs$params)) {
	params = scan(cargs$params,what=character(),quiet=TRUE,comment.char="#")
} else {
	params = cargs$params
}

ind = match(params,descall$faname)
if (any(is.na(ind))) stop("unknown parameters, see config/params.txt")

desc = descall[ind,]
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

ref = "."
if ("ref" %in% names(cargs)) ref = cargs$ref

hasx11 = ! "png" %in% names(cargs) && capabilities("X11")
ask = ! "noask" %in% names(cargs) && hasx11 && interactive()
graph2 = FALSE

dstats = dstatd = array(NA_real_,c(length(dff$file),5,length(dates),npar))
lstatff = lstatffd = lapply(dff$file,function(x) lapply(dates,list))
lstat = lstatd = lapply(seq(dim(desc)[1]),function(x) lstatff)
htime = integer(length(dff$file))
tstep = rep(NA,length(dff$file))

if (interactive()) browser()

for (id in seq(along=dates)) {
	if (is.na(dates[id])) {
		fics = dff$file
		ficsave = sprintf("%s/diag.Rdata",dirname(fics[1]))
	} else {
		fics = sprintf("%02d/%d/%s",res,dates[id],dff$file)
		ficsave = sprintf("%s/diag_%d%02d.Rdata",dirname(fics[1]),dates[id],res)
		cat("\nReading files for base:",dates[id],"\n")
	}

	if (any(! file.exists(fics))) {
		cat("-->",length(which(! file.exists(fics))),"files missing:\n",
			head(fics[! file.exists(fics)]),"...\n")
		stop("Error")
	}

	isbin = regexpr("\\.bin$",fics) > 0

	if (file.exists(ficsave) && all(! isbin)) {
		cat("--> load",ficsave,"and skip maps\n")
		load(ficsave)
		next
	}

	for (i in seq(along=fics)) {
		cat("\nFile",basename(fics[i]),":\n")

		if (isbin[i]) {
			con = file(fics[i],"rb")
			grid = getGrid(con)
			stopifnot(exists("frame") && is.equal(grid,frame))
			frame = toFrame(grid,frame)
			ldata = getVars(con,grid,desc$faname)
			close(con)
		} else {
			frame = getFrame(fics[i])
		}

		if (! is.na(dates[id]) && as.character(frame$base,"%Y%m%d") != dates[id]) {
			message("dates mismatch: ",frame$base," != ",dates[id])
			stop("date mismatch")
		}

		cat("--> Gauss grid:",frame$gauss,"\n")
		if (frame$step %% 60 != 0) {
			s = sprintf("step: +%gs - graph: %s\n",frame$step,dff$graph[i])
		} else if (frame$step %% 3600 != 0) {
			s = sprintf("step: +%gmn - graph: %s\n",frame$step/60,dff$graph[i])
		} else {
			s = sprintf("step: +%gh - graph: %s\n",frame$step/3600,dff$graph[i])
		}

		if (is.na(tstep[i])) {
			cat("-->",s)
			tstep[i] = s
		} else {
			stopifnot(s == tstep[i])
		}

		gg = id == 1 && dff$graph[i]
		if (gg && frame$gauss) {
			mapf = dilat(frame)
			ilateq = equalize(mapf,offset=3)
			if (FALSE && length(ilateq) < frame$nlat) frlow = degrade(frame,ilateq)

			pngalt(sprintf("%s/stretch%d.png",pngd,i))
			tt = c("Stretching coefficient","normalized at Equator")
			plot(mapf,type="l",main=tt,xlab="Latitude index",ylab="Stretching coef.",log="y")
			#points(ilateq,mapf[ilateq],pch=20,cex=.5)
			pngoff()
		}

		if (frame$nlat > 800) {
			cat("--> halving initial resolution:",frame$nlat,max(frame$nlong),"\n")
			frlow = degrade(frame,nlat=2,nlon=2)
		} else {
			frlow = frame
		}

		stopifnot(all(abs(ilev) <= length(frlow$eta)))

		nlev = length(frlow$eta)
		if (length(ilev) == 1 && ilev < 0) {
			stopifnot(-ilev < nlev)

			cat("--> uniform eta selection of levels:",-ilev,"among",nlev,"\n")
			e = seq(min(frlow$eta),max(frlow$eta),length.out=-ilev)
			frlow$ilev = sapply(e,function(x) which.min(abs(x-frlow$eta)))
		} else if (length(ilev) < nlev && ilev > 0) {
			stopifnot(all(ilev %in% seq(nlev)))

			cat("--> selection of levels:",ilev,"\n")
			frlow$ilev = ilev
		} else {
			frlow$ilev = seq(nlev)
		}

		etai = frlow$eta[frlow$ilev]
		yeta = rev(range(frlow$eta))

		date = frame$base+frame$step
		htime[i] = frame$step
		ficref = sprintf("%s/%s/%s",ref,strftime(date,"%H/%Y%m%d"),dff$ref[i])
		if (! file.exists(ficref)) {
			cat("--> no ref file",ficref,"\n")
		} else {
			framo = getFrame(ficref)
		}

		datax = datay = dataz = NULL
		dataox = dataoy = dataoz = NULL

		for (j in seq(npar)) {
			if (desc$ltype[j] == "-") {
				patt = desc$faname[j]
			} else if (length(frlow$ilev) > 1) {
				patt = sprintf("%s\\\\d+%s",desc$ltype[j],desc$faname[j])
			} else if (desc$ltype[j] == "S") {
				patt = sprintf("S%03d%s",ilev,desc$faname[j])
			} else {
				cat("--> level type not 'S' for reading specified level\n")
				next
			}

			ss = desc$symbol[j]
			nom = desc$longname[j]

			cat(". get fields matching pattern",patt,"\n")
			if (regexpr("^\\.",desc$faname[j]) > 0) {
				if (! isbin[i]) next
				data = ldata[[j]]
				if (is.null(data)) next
				if (dim(data)[2] > 1 && dim(data)[2] != length(frlow$ilev)) {
					data = data[,frlow$ilev,drop=FALSE]
				}
			} else {
				if (isbin[i]) next
				data = getField(fics[i],patt,ss,frame,frlow)
			}

			nl = dim(data)[2]

			# primary checks
			stopifnot(all(! is.infinite(data)))
			stopifnot(all(apply(data,2,function(x) length(unique(na.omit(x))) > 1)))
			stopifnot(all(apply(data,2,function(x) any(! is.na(x)))))

			dstats[i,,id,j] = quantile(data[,nl],c(0,.1,.5,.9,1),names=FALSE)

			lstat[[j]][[i]][[id]] = mapexi(data,frlow,desc[j,],doms,prefix=i,graph=gg)

			if (gg && frame$gauss) {
				cat("Statistics of zonal difference:\n")
				ddx = xdiff(data,frlow)
				print(summary(as.vector(ddx)))

				if (FALSE && any(ddx != 0)) {
					pngalt(sprintf("%s/spec45_%s.png",pngd,ss))
					d45 = lat45(data,frlow)
					spectrum(d45[,nl],c(3,3),main="Smoothed spectral density of lat ~45")
					pngoff()
				}
			}

			if (graph2 && nl > 2) {
				cat("Mean vertical profile\n")
				datam = apply(data,2,mean)

				pngalt(sprintf("%s/profmean_%s.png",pngd,ss))
				tt = c("Mean vertical profile",sprintf("Field %s",nom))
				plotz(datam,etai,ylim=yeta,main=tt,xlab=nom,ylab="eta")
				pngoff()

				ddz = vgrad(data,frlow)
				ddz = t(apply(data,1,diff))

				cat("Statistics of vertical gradient:\n")
				print(summary(as.vector(ddz)))
				ix = arrayInd(which.max(abs(data)),dim(data))
				il = ix[1,2]
				ix = ix[1,1]
				plotz(ddz[ix,],etai[-1],ylim=yeta,main=tt,xlab=nom,ylab="eta")
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
			if (regexpr("^\\.",desc$faname[j]) > 0) {
				if (! isbin[i]) next
				datao = ldata[[j]]
				if (is.null(datao)) next
				if (dim(datao)[2] > 1 && dim(datao)[2] != length(frlow$ilev)) {
					datao = datao[,frlow$ilev,drop=FALSE]
				}
			} else {
				if (isbin[i]) next
				datao = getField(ficref,patt,ss,framo,frlow)
			}

			if (desc$save[j] == "X") {
				dataox = datao
			} else if (desc$save[j] == "Y") {
				dataoy = datao
			} else if (desc$save[j] == "Z") {
				dataoz = datao
			}

			ddiff = data-datao

			dstatd[i,,id,j] = quantile(ddiff[,nl],c(0,.1,.5,.9,1),names=FALSE)

			cat("Statistics of forecast error:\n")
			print(summary(as.vector(ddiff)))
			drms = sqrt(mean(ddiff^2))
			cat("Global bias/RMSE:",mean(ddiff),drms,"\n")

			cat("Global bias/RMSE by level:\n")
			dbias = apply(ddiff,2,mean)
			drms = apply(ddiff,2,function(x) sqrt(mean(x^2)))
			print(summary(dbias))
			print(summary(drms))

			lstatd[[j]][[i]][[id]] = mapexi(ddiff,frlow,desc[j,],doms,
				prefix=sprintf("diff%d",i),graph=gg)

			if (graph2) {
				cat("Statistics around diff min/max value:\n")
				i1 = arrayInd(which.min(ddiff),dim(ddiff))[1,1]
				i2 = arrayInd(which.max(ddiff),dim(ddiff))[1,1]
				dom = makeDomain(frlow,i1,c(-20,20),c(-30,30))
				xy = zoom(ddiff,frlow,dom)
				print(summary(as.vector(xy$data)))
				drms = sqrt(mean(xy$data^2))
				cat("RMSE at local min:",drms,"\n")

				pngalt(sprintf("%s/diffmin_%s.png",pngd,ss))
				mapdom(dom,frlow,ddiff[,nl],main=nom)
				pngoff()

				dom = makeDomain(frlow,i2,c(-20,20),c(-30,30))
				xy = zoom(ddiff,frlow,dom)
				print(summary(as.vector(xy$data)))
				drms = sqrt(mean(xy$data^2))
				cat("RMSE at local max:",drms,"\n")

				pngalt(sprintf("%s/diffmin_%s.png",pngd,ss))
				mapdom(dom,frlow,ddiff[,nl],main=nom)
				pngoff()

				if (nl > 2) {
					pngalt(sprintf("%s/profile_%s.png",pngd,ss))
					par(mfrow=c(1,2))
					tt = sprintf("Field %s",nom)
					tt[2] = "Profile at point 'min of error'"
					plotz(data[i1,],etai,ylim=yeta,main=tt,xlab=nom,ylab="eta")
					tt[2] = "Profile at point 'max of error'"
					plotz(data[i2,],etai,ylim=yeta,main=tt,xlab=nom,ylab="eta")
					pngoff()

					ns = length(yzs$lat)
					dataf = yzn$data
					dataf[1:ns,] = yzs$data[ind,]
					#levels = prettydiff(dataf,7)

					pngalt(sprintf("%s/vsecdiff_%s.png",pngd,ss),720)
					x = seq(0,1,length.out=dim(dataf)[1])
					xlab = "<- Equator  ...  //  ...  France  ...  //  ...  North pole ->"
					plotv(x,etai,dataf,nlevels=7,main=tt[1],xlab=xlab,ylab="eta")
					pngoff()

					pngalt(sprintf("%s/vsecfill_%s.png",pngd,ss),720)
					plotf(etai,dataf,nlevels=7,color.palette=cm.colors,main=tt[1],xlab=xlab,
						ylab="eta")
					pngoff()
				}
			}
		}

		if (! is.null(datax) && ! is.null(datay)) {
			cat(". compound field wind speed\n")
			ff = sqrt(datax^2+datay^2)
			print(summary(as.vector(ff)))

			if (gg && frame$gauss && latPole(frame) < 90) {
				cat("rotating wind to geo North\n")
				nord = compass(frame)
				if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
				u = ucs(nord,datax,datay)
				v = vcs(nord,datax,datay)

				pole = list(xlim=lonPole(frame)+c(-6,6),ylim=latPole(frame)+c(-5,5))
				x = zoom(u,frlow,pole)
				y = zoom(v,frlow,pole)
				pngalt(sprintf("%s/map%dpole_ff.png",pngd,i))
				mapdom2(pole,x,x$data[,1],y$data[,1],main="u/v wind components")
				pngoff()
			} else {
				u = datax
				v = datay
			}

			lstatff[[i]][[id]] = mapex2(u,v,frlow,doms,prefix=i)
		}

		if (! is.null(dataox) && ! is.null(dataoy)) {
			cat(". compound field wind speed\n")
			ffo = sqrt(dataox^2+dataoy^2)
			print(summary(as.vector(ffo)))

			if (framo$gauss && latPole(framo) < 90) {
				cat("rotating wind to geo North\n")
				nord = compass(framo)
				if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
				uo = ucs(nord,dataox,dataoy)
				vo = vcs(nord,dataox,dataoy)

				pole = list(xlim=lonPole(frame)+c(-6,6),ylim=latPole(frame)+c(-5,5))
				x = zoom(dataox,frlow,pole)
				y = zoom(dataoy,frlow,pole)
				pngalt(sprintf("%s/mapdiff%dpole_ff.png",pngd,i))
				mapdom2(pole,x,x$data[,nl],y$data[,nl],main="u/v wind components")
				pngoff()
			}

			ffd = ff-ffo
			print(summary(as.vector(ffd)))
			drms = sqrt(mean(ffd^2))
			cat("Global bias/RMSE:",mean(ffd),drms,"\n")

			cat("Global bias/RMSE by level:\n")
			dbias = apply(ffd,2,mean)
			drms = apply(ffd,2,function(x) sqrt(mean(x^2)))
			print(summary(dbias))
			print(summary(drms))

			descff = data.frame(longname="wind speed",symbol="ff")
			lstatffd[[i]][[id]] = mapexi(ffd,frlow,descff,doms,prefix=sprintf("diff%d",i))
		}
	}

	save("lstat","lstatd","lstatff","lstatffd","dstats","dstatd","tstep","htime",
		file=ficsave)
}

if (! is.null(cargs$png)) cat(tstep,file=sprintf("%s/steps.txt",cargs$png))

labs = sprintf("Lead time (%s)",c("s","mn","h","days"))
xaxis = data.frame(unit=c(1,60,3600,86400),label=labs,mindiff=5*c(0,60,3600,86400),
	cunit=c("s","mn","h","d"))
iu = findInterval(diff(range(htime)),xaxis$mindiff)
tunit = xaxis$unit[iu]
cunit = xaxis$cunit[iu]
xlabt = xaxis$label[iu]
htime = htime/tunit
tfreq = 1
if (tunit == 3600 && diff(range(htime)) > 12) tfreq = 3
ht6 = pretty(htime/tfreq)*tfreq

if (length(fics) > 1) {
	cat("Summary statistics of forecasts (whole domain, whole atm.)\n")
	tt = sprintf("Forecast of %s",desc$longname)

	nf = dim(desc)[1]
	nr = min(3,nf)

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/stats%d.png",pngd,i))
		op = par(mfrow=c(nr,1),mar=c(3,3,3,1)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			ds = apply(dstats[,,,j,drop=FALSE],1:2,mean,na=TRUE)
			plotb(htime,t(ds),main=tt[j],xlab=xlabt,ylab=desc$symbol[j],xaxt="n")
			axis(1,ht6)
		}

		pngoff(op)
	}

	ndom = length(doms)
	nt = length(lstat[[1]])
	if (nt < length(htime)) {
		length(htime) = nt
		length(ht6) = nt
	}

	cat("Statistics of forecasts over",ndom,"domains\n")
	longname = desc$longname
	symbol = desc$symbol

	if (any(! is.na(lstatff[[1]][[1]][[1]]))) {
		lstat = c(lstat,list(lstatff))
		longname = c(desc$longname,"wind speed")
		symbol = c(desc$symbol,"ff")
	}

	for (j in seq(along=lstat)) {
		ss = symbol[j]
		nom = longname[j]

		# complexity because of bin files
		nl = dim(lstat[[j]][[1]][[1]])[2]
		if (is.null(nl)) nl = dim(lstat[[j]][[2]][[1]])[2]
		qv100 = array(NA_real_,dim=c(101,nl,ndom,length(dates),length(fics)))
		for (i in seq(along=fics)) {
			for (id in seq(along=dates)) {
				if (is.list(lstat[[j]][[i]][[id]])) next
				qv100[,,,id,i] = simplify2array(lstat[[j]][[i]][[id]])
			}
		}

		qv = apply(qv100,c(2,3,5),quantile,prob=c(0,2,5,8,10)/10,na.rm=TRUE)

		tt = sprintf("Statistics of %s",nom)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = c(1,1+nl%/%2,nl)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/stat%d_%s.png",pngd,i,ss))
			op = par(mfcol=c(nr,nc),mar=c(3,3,3,1)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				for (il in indl[1:nr]) {
					tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il])
					plotb(htime,qv[,il,id,],main=tt,xlab=xlabt,ylab=ss,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff(op)
		}

		if (nl == 1) next

		etai = frlow$eta[frlow$ilev]
		yeta = rev(range(frlow$eta))
		indt = which(apply(qv,4,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[seq(1,length(indt),length.out=min(length(indt),3))]
		nr = min(2,ndom)
		nc = min(3,length(indt))

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/statv%d_%s.png",pngd,i,ss))
			op = par(mfrow=c(nr,nc),mar=c(3,3,3,0)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				for (it in indt[1:nc]) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],cunit)
					plotbv(qv[,,id,it],etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
				}
			}

			pngoff(op)
		}
	}
}

if (any(! is.na(dstatd))) {
	cat("Summary statistics of forecast error (whole domain, whole atm.)\n")
	tt = sprintf("Error of %s",desc$longname)

	nf = dim(desc)[1]
	nr = min(3,nf)

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/dstats%d.png",pngd,i))
		op = par(mfrow=c(nr,2),mar=c(3,3,3,1)+.1,mgp=c(2,.75,0))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			ds = apply(dstatd[,,,j,drop=FALSE],1:2,mean,na.rm=TRUE)
			plotb(htime,t(ds),main=tt[j],xlab=xlabt,ylab=ss,xaxt="n")
			axis(1,ht6)
		}

		pngoff(op)
	}

	ndom = length(doms)
	nt = length(lstatd[[1]])
	if (nt < length(htime)) {
		length(htime) = nt
		length(ht6) = nt
	}

	cat("Scores of forecasts over",ndom,"domains\n")
	longname = desc$longname
	symbol = desc$symbol

	if (any(! is.na(lstatffd[[1]][[1]][[1]]))) {
		lstatd = c(lstatd,list(lstatffd))
		longname = c(desc$longname,"wind speed")
		symbol = c(desc$symbol,"ff")
	}

	for (j in seq(along=lstatd)) {
		ss = symbol[j]
		nom = longname[j]

		# complexity because of bin files
		nl = dim(lstatd[[j]][[1]][[1]])[2]
		if (is.null(nl)) nl = dim(lstatd[[j]][[2]][[1]])[2]
		qv100 = array(NA_real_,dim=c(101,nl,ndom,length(dates),length(fics)))
		for (i in seq(along=fics)) {
			for (id in seq(along=dates)) {
				if (is.list(lstatd[[j]][[i]][[id]])) next
				qv100[,,,id,i] = simplify2array(lstatd[[j]][[i]][[id]])
			}
		}

		qv = apply(qv100,c(2,3,5),quantile,prob=c(0,2,5,8,10)/10,na.rm=TRUE)
		rmsv = apply(qv100,c(2,3,5),rms,na.rm=TRUE)

		tt = sprintf("Error of %s",nom)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = c(1,1+nl%/%2,nl)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/err%d_%s.png",pngd,i,ss))
			op = par(mfcol=c(nr,nc),mar=c(3,3,3,1)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				for (il in indl[1:nr]) {
					tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il])
					plotb(htime,qv[,il,id,],main=tt,xlab=xlabt,ylab=ss,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff(op)
		}

		tt = sprintf("Score of %s",nom)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/score%d_%s.png",pngd,i,ss))
			op = par(mfcol=c(nr,nc),mar=c(3,3,3,1)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				for (il in indl) {
					tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il])
					matplott(htime,t(rmsv[,il,id,]),xlab=xlabt,ylab=ss,main=tt,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff(op)
		}

		if (nl == 1) next

		etai = frlow$eta[frlow$ilev]
		yeta = rev(range(frlow$eta))
		indt = which(apply(qv,4,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[seq(1,length(indt),length.out=min(length(indt),3))]
		nr = min(2,ndom)
		nc = min(3,length(indt))

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/scorev%d_%s.png",pngd,i,ss))
			op = par(mfrow=c(nr,nc),mar=c(3,3,3,0)+.1,mgp=c(2,.75,0))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				for (it in indt[1:nc]) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],cunit)
					matplotv(t(rmsv[,,id,it]),etai,xlab=ss,ylab="eta",main=tt,ylim=yeta)
				}
			}

			pngoff(op)
		}
	}
}

