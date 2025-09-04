library(mfnode)

varqc = function(nd)
{
	if (length(grep("RAQC=",nd)) == 0 || length(grep("RLQC=",nd)) == 0) return(NULL)

	indo = grep("^ *VAR *= *\\d+",nd)
	iv = as.integer(sub("VAR *= *(\\d+).+","\\1",nd[indo]))
	qc = vector(length(unique(iv)),mode="list")
	for (i in unique(iv)) {
		ind = which(iv == i)
		notvar = intlines(sub(".+ NOTVAR *=","",nd[indo[ind[1]]]))
		v = numlines(gsub("\\*+","-9999.",nd[indo[ind[-1]]]))
		v[v==-9999] = NA
		qc[[i]] = matrix(c(v,notvar),nc=length(ind))
	}

	iq = sapply(qc,function(x) length(na.omit(x[,3:6])) > 0)
	qc = qc[iq]

	qc
}

qcheck = function(qc)
{
	ra = sapply(qc,"[",,1)
	rl = sapply(qc,"[",,2)
	bg = sapply(qc,"[",,3:5,simplify="array")

	list(ra=ra,rl=rl,bg=bg)
}

jcdfi = function(nd,nflevg)
{
	idfi = grep("Jc-DFI Diagnostic Table",nd)
	if (length(idfi) != 2) return(NULL)

	jc = t(matrix(numlines(nd[seq(idfi[1]+1,idfi[2]-2)]),ncol=nflevg))
	jc
}

zonalcoef1 = function(nd,ind)
{
	noms = strsplit(gsub("^ *","",nd[ind+1]),split=" +")[[1]]
	lc = strsplit(gsub("^ *","",nd[ind+1+seq(nsmax+1)]),split=" +")
	mc = t(sapply(lc,as.numeric))
	n = mc[,1]
	mc = mc[,-1]
	dimnames(mc) = list(n=n,var=noms[-1])
	mc
}

zonalcoef = function(nd)
{
	ind = grep("SUQNORM *- printing the coefficients of the norm",nd)
	if (length(ind) == 0) return(NULL)

	if (length(ind) == 1) {
		nsmax = getvar("NSMAX",nd)
		mc = zonalcoef1(nd,ind)
		return(list(mc))
	} else {
		ins = grep("COMPUTE NORM FOR RESOLUTION",nd)
		stopifnot(length(ind) == length(ins))
		lmc = list()
		for (i in seq(along=ind)) {
			nsmax = intlines(nd[ins[i]])
			lmc[[i]] = zonalcoef1(nd,ind)
		}

		return(lmc)
	}
}

gom = function(nd)
{
	ind = grep("GOM variable +RMS",nd)
	if (length(ind) == 0) return(NULL)

	i1 = grep("Norm of interpolation operator .+ - start",nd)
	i2 = grep("Norm of interpolation operator .+ - end",nd)
	stopifnot(length(i1) == length(i2))

	lgom = list()

	for (i in seq(along=i1)) {
		ind = i1[i]:i2[i]
		ig = grep(sprintf(" *\\w+ +%s *$",Gnum),nd[ind])
		lg = strsplit(gsub("^ *","",nd[ind[ig]]),split=" +")
		gomvar = sapply(lg,"[",1)
		gomrms = as.numeric(sapply(lg,"[",2))

		it = grep("Observation type: +\\d+ +Number of obs:",nd[ind])

		type = gsub("Observation type: +(\\d+) +Number of obs: *(\\d+)","\\1",nd[ind[it]])
		nb = gsub("Observation type: +(\\d+) +Number of obs: *(\\d+)","\\2",nd[ind[it]])

		nvar = length(gomrms)/length(it)
		stopifnot(nvar*length(it) == length(gomrms))

		#cat(". group",i,":",nvar,"variables,",length(it),"types\n")
		gom = matrix(gomrms,ncol=length(it),dimnames=list(var=gomvar[1:nvar],type=type))
		attr(gom,"nobs") = as.integer(nb)
		lgom[[i]] = gom
	}

	lgom
}

jotable = function(nd)
{
	if (length(grep("JOT-sname",nd)) == 0) return(NULL)

	iobst = grep("Obstype +\\d+ +=+",nd,ignore.case=TRUE)
	ijog = grep("Jo Global",nd)
	ndo = nd[iobst[1]:ijog[1]]
	ijoh = grep("Codetype +\\d+ +=+",ndo)
	ijot = grep("^ +\\w+ +\\d+( +\\d+\\.\\d+){3}",ndo)
	lj = strsplit(ndo[ijot],split=" +")
	jot = t(sapply(lj,function(x) as.numeric(x[3:6])))
	jot = as.data.frame(jot)
	names(jot) = c("DataCount","Jo_Costfunction","Jo/n","ObsErr")
	vars = sapply(lj,"[",2)

	indi = findInterval(ijot,ijoh)
	code = sub(" +Codetype +(\\d+) +.+","\\1",ndo[ijoh[indi]])

	jot = cbind(Codetype=as.integer(code),Variable=vars,jot)

	jot
}

joMean = function(jot,by=c("Codetype","Variable"))
{
	if (by == "Codetype") {
		jog = by(jot[-(1:2)],jot$Codetype,colMeans)
		jod = jog[[1]]
		for (jo in jog[-1]) jod = rbind(jod,jo,deparse.level=0)
		jod = cbind(data.frame(Codetype=dimnames(jog)[[1]]),jod)
	} else if (by == "Variable") {
		jog = by(jot[-(1:2)],jot$Variable,colMeans)
		jod = jog[[1]]
		for (jo in jog[-1]) jod = rbind(jod,jo,deparse.level=0)
		jod = cbind(data.frame(Variable=dimnames(jog)[[1]]),jod)
	}

	jod
}

scatqc = function(nd)
{
	i1 = grep("CELL NO. +N.OBS",nd)
	i2 = grep("TOTAL OF +\\d+ +ASCAT +SCATT",nd)
	if (length(i1) == 0 || length(i2) == 0) return(NULL)

	stopifnot(length(i1) == length(i2))

	lqc = vector("list",length(i1))

	for (i in seq(along=i1)) {
		ind = seq(i1[i]+1,i2[i]-1)
		l = strsplit(sub("^ +","",nd[ind]),split=" +")
		m = sapply(l,as.character)
		df = as.data.frame(t(m))
		names(df) = strsplit(gsub("^ +CELL NO.","cell",nd[i1[i]]),split=" +")[[1]]
		lqc[[i]] = df
	}

	lqc
}

screenstat = function(nd)
{
	i1 = grep("OB.TYP +REPORTS +ACTIVE +PASSIVE",nd)
	if (length(i1) == 0) return(NULL)

	i2 = grep("^ +-{30,}",nd)
	ii = grep("EVENT SUMMARY OF",nd)
	i2 = i2[i1[1] < i2 & i2 < ii[1]]

	stopifnot(length(i1) == length(i2))

	noms = strsplit(gsub("^ *","",nd[i1[1]]),split=" +")[[1]]

	lst = vector("list",length(i1))

	for (i in seq(along=i1)) {
		ind = seq(i1[i]+1,i2[i]-1)
		l = strsplit(sub("^ +","",nd[ind]),split=" +")
		m = sapply(l,as.integer)
		nums = m[1,]
		m = m[-1,]
		dimnames(m) = list(noms[-1],nums)
		lst[[i]] = m
	}

	ii = grep("STATUS SUMMARY OF",nd)
	names(lst) = gsub("^ *STATUS SUMMARY OF +(\\w+).*","\\1",nd[ii])
	lst
}

screenobs = function(nd)
{
	i1 = grep("NUMBER OF +(.+?) +IN DIFFERENT OBSERVATION TYPES",nd)
	if (length(i1) == 0) return(NULL)

	i2 = grep("^ +-{30,}",nd)
	ii = grep("FINISH SCREENING OF OBSERVATIONS",nd)
	i2 = i2[i1[1] < i2 & i2 < ii]

	stopifnot(length(i1) == length(i2))

	noms = gsub("^ *NUMBER OF +(.+?) +IN DIFFERENT OBSERVATION TYPES.*","\\1",nd[i1])

	lno = vector("list",length(i1))

	for (i in seq(along=i1)) {
		ind = seq(i1[i]+2,i2[i]-1)
		ii = grep("VARIAB",nd[ind])
		stopifnot(length(ii) == 1)
		l = strsplit(sub("^ +","",nd[ind[-(1:ii)]]),split=" +")
		m = sapply(l,as.integer)
		ii = grep("OB\\.TYP",nd[ind])
		stopifnot(length(ii) == 1)
		nums = m[1,]
		m = m[-1,]
		dimnames(m) = list(intlines(nd[ind[ii]]),nums)
		lno[[i]] = m
	}

	names(lno) = noms
	lno
}

canawagons = function(nd)
{
	i1 = grep("nombre total de wagons",nd)
	i2 = grep("Elimination des SYNOP sur relief",nd)
	if (length(i1) == 0 || length(i2) == 0) return(NULL)

	stopifnot(length(i1) == length(i2))

	it = grep("Type d'observations numero",nd)

	# list of histograms per obs. type (there should be 2 types)
	lht = vector("list",2)

	for (i in seq(along=i1)) {
		cat(". step",i,"\n")
		iit = it[i1[i] < it & it < i2[i]]
		ind = seq(i1[i],i2[i])
		indi = grep("wagons rejetes",nd[ind])

		lh = list()
		n = 0
		for (j in indi) {
			l = strsplit(sub("^ +","",nd[ind[(j-5):(j-2)]]),split=" +")
			m = sapply(l,as.numeric)

			k = which(iit < ind[j])
			ik = iit[k[length(k)]]
			attr(m,"type") = as.integer(sub("Type d'observations numero","",nd[ik]))
			param = gsub("^ +| +$","",nd[ind[j]-7])
			attr(m,"param") = param
			n = n+1
			lh[[n]] = m
			names(lh)[n] = sub(" +metres?","m",param)
		}

		types = as.integer(sub("Type d'observations numero","",nd[iit]))
		if (i == 1) {
			type1 = types
		} else {
			stopifnot(identical(types,type1))
		}

		# in lht, reverse types (j) and occurrences (i)
		for (j in seq(along=types)) {
			mh = simplify2array(lh[sapply(lh,attr,"type") == types[j]])
			if (i == 1) lht[[j]] = list()
			lht[[j]][[i]] = mh
		}

		if (i == 1) names(lht) = types
	}

	lh = lapply(lht,simplify2array)
}

jacobian = function(nd)
{
	ind = grep("GREPCOST - ITER",nd)
	if (length(ind) == 0) return(NULL)

	jacob = t(matrix(numlines(nd[ind]),ncol=length(ind)))
	jac = sub("GREPCOST.+SIM,(.+JCVARBC).+","\\1",nd[ind[1]])
	jac = trimws(strsplit(jac,split=",")[[1]])
	dimnames(jacob)[[2]] = jac

	iter = as.integer(sub(".+JCVARBC +(\\d+) .+","\\1",nd[ind]))
	if (iter[length(iter)] > 40) iter[length(iter)] = 40
	dimnames(jacob)[[1]] = iter

	jacob
}

ritzvalues = function(nd)
{
	ind = grep("ritz values",nd,ignore.case=TRUE)
	if (length(ind) == 0) return(NULL)

	ritz = t(matrix(numlines(nd[ind[-1]]),nrow=2))
	ritz
}

costgrad = function(nd)
{
	grad = numlines(grep("estimated reduction in norm",nd,value=TRUE,ignore.case=TRUE))
	quad = numlines(grep("estimated quadratic cost",nd,value=TRUE,ignore.case=TRUE))

	list(grad=grad,quad=quad)
}

bmat = function(nd,nflevg)
{
	i1 = grep("correlations *\\(\\*100\\)",nd)
	if (length(i1) == 0) return(NULL)

	ii = grep("Diagnostics of the horizontal correlations",nd)

	nf = length(i1)
	if (nflevg > 100) {
		ind = seq(i1[nf]+5,ii-1)-i1[nf]
	} else {
		ind = seq(i1[nf]+3,ii-1)-i1[nf]
	}

	lcor = list()

	nl = length(ind)
	for (i in seq(along=i1)) {
		indi = i1[i]+ind
		cor = sapply(indi,function(j) gsub("(...)"," \\1",substring(nd[j],2)))
		mb = intlines(cor)
		m = matrix(mb,ncol=nflevg)
		lcor[[i]] = m/100
	}

	noms = gsub("^ *(\\w.+\\>) +correlations.+","\\1",nd[i1])
	ips = grep("T, *Ps",noms)
	if (length(ips) == 1) {
		lcor[[nf+1]] = lcor[[ips]][nflevg+1,]
		lcor[[ips]] = lcor[[ips]][-(nflevg+1),]
		noms[ips] = "T"
		noms = c(noms,"Ps")
	}

	names(lcor) = noms

	lcor
}

bcor = function(nd)
{
	i1 = grep("VARIABLE +\\d+ +CORREL LENGTH SCALES",nd)
	if (length(i1) == 0) return(NULL)

	lb = vector("list",length(i1))
	noms = gsub("^ *(\\w.+\\>) +correlations.+","\\1",nd[i1])

	ii = grep("Diagnostics of the horizontal correlations",nd)

	nf = length(i1)
	m = intlines(nd[seq(i1[nf]+3,ii-1)])
	for (i in seq(nf-1)) {
		ind = seq(i1[i]+1,i1[i+1]-1)
		lb[[i]] = numlines(nd[ind])
	}

	lb[[nf]] = numlines(nd[i1[nf]+1])
	lb
}

jberr = function(nd,nflevg)
{
	i1 = grep("Vor +unbal Div",nd)
	if (length(i1) == 0) return(NULL)

	ii = grep("Q +unbal O3",nd)
	stopifnot(length(ii) == 1)

	noms = unlist(strsplit(sub("^ +","",nd[i1:ii]),split="  +"))
	jb = numlines(sub("^ +\\d+","",nd[ii+seq(nflevg)]))
	jb = matrix(jb,nrow=nflevg,byrow=TRUE,dimnames=list(1:nflevg,noms))

	i0 = apply(jb,2,function(x) any(x != -999))
	jb[,i0]
}

args = commandArgs(trailingOnly=TRUE)
largs = strsplit(args,split="=")
cargs = lapply(largs,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(largs,function(x) x[1])

if (! "png" %in% names(cargs) && capabilities("X11")) {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
} else {
	cat("--> no X11 device, sending plots to PNG files\n")
}

if (! "png" %in% names(cargs)) cargs$png = "."

if (interactive()) browser()

fnode = grep("=",args,invert=TRUE,value=TRUE)
leg = NULL
nd = readLines(fnode[1],skipNul=TRUE)
if (length(fnode) > 1) {
	library(mfnode)
	nds = lapply(fnode[-1],readLines,skipNul=TRUE)
	leg = gsub(".*node","",fnode,ignore.case=TRUE)
}

if ("leg" %in% names(cargs)) leg = cargs$leg

nflevg = getvar("NFLEVG",nd)

cat("Values of Var QC\n")
qc = varqc(nd)
if (! is.null(qc)) {
	qc1 = qc[[1]]
	types = seq(dim(qc1)[1])
	q = qcheck(qc)

	png(sprintf("%s/varqc.png",cargs$png))
	par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	barplot(t(q$ra),names.arg=types,main="RAQC, variable 1",xlab="Obs type",ylab="RAQC")
	barplot(t(q$rl),names.arg=types,main="RLQC, variable 1",xlab="Obs type",ylab="RLQC")
	nc = dim(qc1)[2]
	#barplot(q$raqc,names.arg=types,beside=TRUE,main="",xlab="Obs type",ylab="RBGQC(1:3)")
	matplot(q$bg[,,1],main="Variable 1",xlab="Obs type",ylab="RBGQC(1:3)",type="o",lty=1,
		pch=20)
	#barplot(qc1[,nc],names.arg=types,main="Variable 1",xlab="Obs type",ylab="NOTVAR")
	ix = which.max(apply(q$bg[,1,,drop=FALSE],1,max,na.rm=TRUE))
	sstt = sprintf("max: %g (var: %d)",max(q$bg[ix,1,],na.rm=TRUE),ix)
	boxplot(q$bg[,1,],outline=FALSE,main=c("RBGQC(1)",sstt),xlab="Variable",
		ylab="RBG among obs. types")
	dev.off()
}

cat("Jo tables\n")
jot = jotable(nd)
if (! is.null(jot)) {
	jod = joMean(jot,"Codetype")
	con = file(sprintf("%s/jo.txt",cargs$png),open="w+")
	cat("<pre>Means of Jo by codetype:\n",file=con)
	jod[,-1] = signif(jod[,-1],5)
	write.table(format(jod,width=12,justify="right"),file=con,quote=FALSE,row.names=FALSE,
		col.names=sprintf("% 12.12s",names(jod)))

	jod = joMean(jot,"Variable")
	cat("\nMeans of Jo by variable:\n",file=con)
	jod[,-1] = signif(jod[,-1],5)
	write.table(format(jod,width=12,justify="right"),file=con,quote=FALSE,row.names=FALSE,
		col.names=sprintf("% 12.12s",names(jod)))
	cat("</pre>\n",file=con)
	close(con)
}

cat("GOM values\n")
lgom = gom(nd)
if (! is.null(lgom)) {
	gomvar = unique(unlist(sapply(lgom,function(x) dimnames(x)[[1]])))
	nvar = sapply(lgom,function(x) dim(x)[1])
	nt = sapply(lgom,function(x) dim(x)[2])
	stopifnot(all(duplicated(nt)[-1]))
	nvarx = max(nvar)
	lgomx = lgom[nvar == nvarx]
	gomvar = dimnames(lgomx[[1]])[[1]]
	mgom = simplify2array(lgomx)
	cat("Global means of GOM variables:\n")
	print(apply(mgom,1,mean))

	ngom = t(sapply(lgomx,function(x) attr(x,"nobs")))

	png(sprintf("%s/gomt.png",cargs$png))
	par(mfrow=c(3,3))
	for (i in seq(9)) {
		m = t(mgom[i,,])
		tt = sprintf("GOM variable %s",gomvar[i])
		matplot(m,type="b",lty=1,pch=20,main=c(tt,"per obs type"),xlab="Time-step",
			ylab="RMS")
	}

	dev.off()

	nvarx = min(nvar)
	lgomx = lgom[nvar == nvarx]
	gomvar = dimnames(lgomx[[1]])[[1]]
	mgom = simplify2array(lgomx)
	cat("Global means of GOM variables:\n")
	print(apply(mgom,1,mean))

	png(sprintf("%s/gom.png",cargs$png))
	par(mfrow=c(3,2))
	for (i in seq(min(nvarx,6))) {
		m = t(mgom[i,,])
		tt = sprintf("GOM variable %s",gomvar[i])
		matplot(m,type="b",lty=1,pch=20,main=c(tt,"per obs type"),xlab="Time-step",
			ylab="RMS")
	}

	dev.off()
}

cat("Jc DFI\n")
jc = jcdfi(nd,nflevg)
if (! is.null(jc)) {
	if (length(fnode) > 1) {
		jcs = lapply(nds,jcdfi,nflevg)
		jcs = simplify2array(c(list(jc),jcs))
	}

	png(sprintf("%s/jcdfi.png",cargs$png))
	ic = which(apply(jc,2,function(x) any(x > 0)))
	nr = nc = min(2,length(ic))
	op = par(mfrow=c(nr,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	if (length(fnode) > 1) {
		cols = seq(dim(jcs)[2])
		for (i in ic) {
			plotvmean(jcs[,i,],1:nflevg,main="Jc after DFI",col=cols,lty=1,xlab="Jc",
				ylab="Level",legend=leg)
		}
	} else {
		for (i in ic) {
			plot(jc[,i],1:nflevg,type="l",xlab="Jc",ylab="Level",ylim=c(nflevg,1))
		}
	}

	dev.off()
}

cat("Jacobian values\n")
jacob = jacobian(nd)
if (! is.null(jacob)) {
	jacobs = list(jacob)
	if (length(fnode) > 1) {
		jacobs = lapply(nds,jacobian)
		jacobs = simplify2array(c(list(jacob),jacobs))
	}

	iter = as.integer(dimnames(jacob)[[1]])
	jac = dimnames(jacob)[[2]]

	png(sprintf("%s/jacob.png",cargs$png))
	par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),pch=20)

	if (length(fnode) > 1) {
		for (i in 1:min(dim(jacob)[2],4)) {
			m = matrix(jacobs[,i,],nrow=length(iter))
			matplot(iter,m,type="o",main=jac[i],xlab="Iteration",ylab=jac[i],lty=1)
			legend("topleft",leg,lty=1,pch=20,col=seq(along=leg))
		}
	} else {
		for (i in 1:min(dim(jacob)[2],4)) {
			plot(iter,jacob[,i],type="o",main=jac[i],xlab="Iteration",ylab=jac[i])
		}
	}

	dev.off()
}

cat("Ritz values and norm of gradient\n")
ritz = ritzvalues(nd)
if (! is.null(ritz)) {
	cost = costgrad(nd)

	if (length(fnode) > 1) {
		rits = lapply(nds,ritzvalues)
		rits = simplify2array(c(list(ritz),rits))
		costs = lapply(nds,costgrad)
		grads = sapply(c(list(cost),costs),"[[","grad",simplify="array")
		quads = sapply(c(list(cost),costs),"[[","quad",simplify="array")
	}

	png(sprintf("%s/grad.png",cargs$png))
	xlab = "Iteration"
	iter = seq(dim(ritz)[1])-1
	iterg = seq(along=cost$grad)-1
	iterq = seq(along=cost$quad)-1

	par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),pch=20)
	if (length(fnode) > 1) {
		plotmean(iter,rits[,1,],leg,main="Ritz value 1",xlab=xlab,ylab="Ritz 1")
		plotmean(iter,rits[,2,],leg,main="Ritz value 2",xlab=xlab,ylab="Ritz 2")
		plotmean(iterg,grads,leg,main="Norm of gradient",xlab=xlab,ylab="Gradient")
		plotmean(iterq,quads,leg,main="Quadratic cost J",xlab=xlab,ylab="Cost function")
	} else {
		plot(iter,ritz[,1],type="o",main="Ritz value 1",xlab=xlab,ylab="Ritz 1")
		plot(iter,ritz[,2],type="o",main="Ritz value 2",xlab=xlab,ylab="Ritz 2")
		plot(iterg,cost$grad,type="o",main="Norm of gradient",xlab=xlab,ylab="Gradient")
		plot(iterq,cost$quad,type="o",main="Quadratic cost J",xlab=xlab,
			ylab="Cost function")
	}

	dev.off()
}

cat("B matrix\n")
lcor = bmat(nd,nflevg)
if (! is.null(lcor)) {
	png(sprintf("%s/bmat.png",cargs$png))
	par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	for (i in 1:4) {
		tt = sprintf("Correlation matrix of %s",names(lcor)[i])
		sstt = sprintf("range: %s",paste(range(lcor[[i]]),collapse=" "))
		image(1:nflevg,1:nflevg,lcor[[i]],col=hcl.colors(10,"BlueRed",alpha=.8),
			main=c(tt,sstt),xlab="Level",ylab="Level")
		contour(1:nflevg,1:nflevg,lcor[[i]],5,labcex=.7,add=TRUE)
	}

	dev.off()

	png(sprintf("%s/bmatps.png",cargs$png))
	plot(lcor[[5]],1:nflevg,type="l",ylim=c(nflevg,1),main="Correlation for Ps",
		xlab="Correlation",ylab="Level")
	abline(v=0,col="darkgrey")

	dev.off()
}

cat("Correlation length scales\n")
lcor = bcor(nd)
if (! is.null(lcor)) {
	ml = simplify2array(lcor[1:5])/1000
	png(sprintf("%s/corlen.png",cargs$png))
	tt = c("Correlation lengths of CV (/1000)",sprintf("Cor. length for SP: %g",lcor[[6]]/1000))
	matplot(ml,1:nflevg,type="l",lty=1,ylim=c(nflevg,1),main=tt,
		xlab="Cor. length",ylab="Level")
	legend("topleft",sprintf("Var %d",1:5),lty=1,col=1:5)
	dev.off()
}

cat("Jb standard errors\n")
jb = jberr(nd,nflevg)
if (! is.null(jb)) {
	tt = sprintf("Jb stderr of %s",dimnames(jb)[[2]])
	for (i in 1:3) {
		png(sprintf("%s/jberr%d.png",cargs$png,i))
		nc = min(3,dim(jb)[2]-3*(i-1))
		par(mfrow=c(1,nc),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		for (j in 3*(i-1)+seq(nc)) {
			plot(jb[,j],1:nflevg,type="l",lty=1,ylim=c(nflevg,1),main=tt[j],
				xlab="Standard error",ylab="Level")
			abline(v=0,col="grey")
		}

		dev.off()
	}
}

cat("CANARI statistics\n")
lh = canawagons(nd)
if (! is.null(lh)) {
	types = names(lh)
	if (length(fnode) > 1) {
		lhs = lapply(nds,canawagons)
		stopifnot(all(sapply(lhs,names) == types))
	}

	for (i in seq(along=lh)) {
		h = lh[[i]]
		if (length(fnode) > 1) stopifnot(all(sapply(lhs,function(x) all(dim(x) == dim(h)))))
		vars = dimnames(h)[[3]]

		cat(". type:",types[i],"\n")
		cat(". vars:",vars,"\n")
		for (j in seq(dim(h)[4])) {
			ficpng = sprintf("%s/residu%d_%s.png",cargs$png,j,types[i])
			png(ficpng)

			nr = min(3,length(vars))
			nc = (length(vars)-1)%/%nr+1
			par(mfrow=c(nr,nc),mgp=c(2,1,0))

			for (k in seq(along=vars)) {
				hi = h[,,k,j]
				tt = c(vars[k],sprintf("type: %s",types[i]))
				plot(hi[,1],hi[,2],type="s",main=tt,xlab="Residual",ylab="Frequency")
				if (length(fnode) > 1) {
					for (l in seq(along=lhs)) {
						hs = lhs[[l]][[i]][,,k,j]
						lines(hs[,1],hs[,2],type="s",col=l+1)
					}
				}
			}

			dev.off()
		}
	}

	ind = grep("Renseignements fournis par CAIDGU",nd)
	field = sub("^ *Renseignements.+ champ +(\\w+) +.+","\\1",nd[ind])
	type = sub("^ *Renseignements.+ champ \\w+ (\\w+.*?) +sur.+","\\1",nd[ind])
	i2 = grep("Impression des diagnostics d'utilisation",nd)
	indv = grep("^([^|]+\\|){9}",nd[1:i2])
	indf = findInterval(indv,ind)
	n = as.numeric(sapply(strsplit(nd[indv],"\\|"),"[",2))
	m = as.numeric(sapply(strsplit(nd[indv],"\\|"),"[",3))
	s = as.numeric(sapply(strsplit(nd[indv],"\\|"),"[",4))
	df = data.frame(level=n,mean=m,"std-dev"=s)
	rownames(df) = format(sprintf("%s (%s)",field,type)[indf],width=15)
	con = file(sprintf("%s/canari.txt",cargs$png),open="w+")
	cat("Level, mean and standard-deviation of fields (CANARI diag):\n",file=con)
	cat(format(c("Statistics","Level","Mean","Std-dev."),width=15),"\n",file=con)
	write.table(format(df,width=15),con,quote=FALSE,col.names=FALSE)
	close(con)
}

cat("Run-time information\n")
rt = runtime(nd)
rts = NULL
if (! is.null(rt)) {
	nts = dim(rt)[1]

	t0 = round(as.numeric(rt$wall[1]-attr(rt,"start"),units="secs") %% 86400,3)
	tint = round(as.numeric(rt$wall[nts]-rt$wall[1],units="secs") %% 86400,3)
	total = round(as.numeric(rt$wall[nts]-attr(rt,"start"),units="secs") %% 86400,3)
	tt = sprintf("setup+step0, forecast, total: %gs, %gs, %gs",t0,tint,total)
	cat(tt,"\n")

	if (nts > 1) {
		if (length(fnode) > 1) {
			rts = lapply(nds,runtime)
			#art = function(rt) sapply(rt,as.numeric,simplify="array")
			#rts = sapply(c(list(rt),rts),art,simplify="array")
		}

		# radiation, if any
		s = grep("L[EM]PHYS *=",nd,value=TRUE)
		lphys = any(as.logical(sub(".*\\<L[EM]PHYS *= *([TF]).+","\\1",s)))
		itr = NULL
		if (lphys) itr = attr(rt,"tsrad")

		# hours (for IO)
		tstep = getvar("TSTEP",nd)
		nth = 3600/tstep
		if (as.integer(nth) != nth) nth = 3*3600/tstep
		if (as.integer(nth) != nth) nth = 6*3600/tstep
		if (as.integer(nth) == nth) {
			ith1 = attr(rt,"tshours")
			ith = ith1
			if (length(ith) > 30) ith = seq(1,nts,by=3*nth)
			if (length(ith) > 30) ith = seq(1,nts,by=6*nth)
			if (length(ith) > 30) ith = seq(1,nts,by=12*nth)
			if (length(ith) > 30) ith = seq(1,nts,by=24*nth)
			if (length(ith) > 30) ith = seq(1,nts,by=48*nth)
		}

		if (nts < 100) {
			its = seq(nts)
		} else {
			pas = c(1,2,3,5,7,11,17,23,41)
			nt = c(0,200,500,800,1500,2500,4000,8000,10000)
			it = findInterval(nts,nt)
			its = unique(sort(c(ith1,seq(1,nts,by=pas[it]))))
		}

		if (nts > 1) {
			tsm = signif(tint/(nts-1),2)
			cat("Mean time of forecast steps:",tsm,"(s)\n")
			if (length(itr) > 2) {
				tintr = sum(rt$dwall[itr[-1]])
				trm = signif(tintr/(length(itr)-1),2)
				tsm = signif((tint-tintr)/(nts-length(itr)),2)
				cat("Mean time of non-radiation steps:",tsm,"(s)\n")
				cat("Mean time of radiation steps:",trm,"(s)\n")
			}
		}

		png(sprintf("%s/runtime.png",cargs$png))
		par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(its-1,rt$dwall[its],type="h",main=c("Wall-time",tt),xlab="Time-step",
			ylab="Time (s)",col="grey10",xaxs="i",yaxs="i")
		if (nth > 1) {
			points(ith-1,rt$dwall[ith],type="h",col="blue1")
			mtext((ith-1)*tstep/3600,3,-.8,at=ith-1,cex=.75,col="blue1")
		}

		if (length(itr) > 1) {
			points(itr-1,rt$dwall[itr],type="h",col="red3")
			mtext(c(tsm,trm),2,at=c(tsm,trm),cex=.75,col=c("black","red3"),las=1)
		}

		start = attr(rt,"start")
		xlab = sprintf("Time since start (=%s) (s)",attr(rt,"start-time"))

		y = seq(dim(rt)[1])-1
		reg = line(rt$wall-start,y)
		tt = sprintf("mean pace: %g steps/s (%d steps in %gs)",signif(coef(reg)[2],2),
			nts,round(total))
		plot(rt$wall-start,y,type="l",main=c("Forecast time",tt),xlab=xlab,
			ylab="Time-step",xaxt="n",xaxs="i",yaxs="i")
		xs = pretty((rt$wall-start)/60)*60
		axis(1,xs)
		abline(reg,lty=2,col="grey")
		ith0 = ith[-1]
		if (length(ith0) > 10) ith0 = ith[seq(1,length(ith),by=2)[-1]]
		if (length(ith0) > 0) {
			h = (ith0-1)*tstep/3600
			mtext(h,2,at=ith0-1,cex=.75,col="blue1",adj=0,las=2)
			segments(rt$wall[ith0]-start,0,y1=ith0-1,col="grey",lty=2)
			segments(0,ith0-1,rt$wall[ith0]-start,col="grey",lty=2)
		}

		for (i in seq(along=rts)) {
			starti = attr(rts[[i]],"start")
			if (start < starti && starti < start+total) {
				lines(rts[[i]]$wall-start,seq(dim(rts[[i]])[1])-1,col=i+1)
			} else {
				lines(rts[[i]]$wall-starti,seq(dim(rts[[i]])[1])-1,col=i+1)
			}
		}

		if (! is.null(leg)) legend("topleft",leg,col=seq(along=leg),lty=1,inset=c(.05,.02))

		dev.off()

		its = seq(nts)
		h = 12
		its = its[tstep/3600*(its-1) <= h]
		itr = itr[tstep/3600*(itr-1) <= h]
		if (length(its) > 200) {
			h = 6
			its = its[tstep/3600*(its-1) <= h]
			itr = itr[tstep/3600*(itr-1) <= h]
		}

		if (length(its) > 1 && length(its) < nts%/%2) {
			nth = 3600/tstep
			if (as.integer(nth) != nth) nth = 3*3600/tstep
			if (as.integer(nth) == nth) {
				ith = ith1
				if (length(ith) > 30) ith = seq(1,nts,by=3*nth)
			}

			ith = ith[tstep/3600*(ith-1) <= 12]

			png(sprintf("%s/runtimez.png",cargs$png))
			par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
			tt = sprintf("Wall-time - zoom [0,%dh]",as.integer(h))
			plot(its-1,rt$dwall[its],type="h",main=tt,
				xlab="Time-step",ylab="Time (s)",col="grey10")
			if (length(ith) > 1) mtext((ith-1)*tstep/3600,3,-.8,at=ith-1,cex=.8,col="blue1")

			if (length(itr) > 1) points(itr-1,rt$dwall[itr],type="h",col="red3")

			plot(its-1,rt$cpu[its],type="h",main="CPU-time - zoom",xlab="Time-step",
				ylab="Time (s)")

			dev.off()
		}
	}
}
