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

fnode = grep("=",args,invert=TRUE,value=TRUE)
nd = readLines(fnode[1])
if (length(fnode) > 1) {
	library(mfnode)
	nds = lapply(fnode[-1],readLines)
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

	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0),pch=20)
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
			art = function(rt) sapply(rt,as.numeric,simplify="array")
			rts = sapply(c(list(rt),rts),art,simplify="array")
		}

		pas = c(1,2,3,5,6,10,15,20,40,50)
		nt = c(0,200,400,600,1200,2000,3000,4000,8000,10000)
		it = findInterval(nts,nt)

		png(sprintf("%s/runtime.png",cargs$png))
		par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		its = seq(1,nts,by=pas[it])
		if (length(fnode) > 1) {
			matplot(its-1,rts[its,2,],type="h",lty=1,main=c("Wall-time",tt),xlab="Time-step",
				ylab="Time (s)")
			legend("topleft",leg,lty=1)
			matplot(its-1,rts[its,3,],type="h",lty=1,main="CPU-time",xlab="Time-step",
				ylab="Time (s)")
			legend("topleft",leg,lty=1)
		} else {
			plot(its-1,rt$dwall[its],type="h",main=c("Wall-time",tt),xlab="Time-step",
				ylab="Time (s)")
			plot(its-1,rt$cpu[its],type="h",main="CPU-time",xlab="Time-step",ylab="Time (s)")
		}

		dev.off()
	}
}
