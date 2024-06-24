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

	llh = vector("list",length(i1))

	for (i in seq(along=i1)) {
		cat(". step",i,"\n")
		iit = it[i1[i] < it & it < i2[i]]
		ind = seq(i1[i],i2[i])
		indi = grep("wagons rejetes",nd[ind])

		lh = lapply(indi,function(j) readhist(nd[ind[(j-5):(j-2)]]))

		lh = list()
		n = 0
		for (j in indi) {
			l = strsplit(sub("^ +","",nd[ind[(j-5):(j-2)]]),split=" +")
			m = sapply(l,as.numeric)

			k = which(iit < ind[j])
			ik = iit[k[length(k)]]
			param = gsub("^ +| +$","",nd[ind[j]-7])
			attr(m,"param") = sub(" +metres?","m",param)
			attr(m,"type") = as.integer(sub("Type d'observations numero","",nd[ik]))
			n = n+1
			lh[[n]] = m
			names(lh)[n] = param
		}

		llh[[i]] = lh
	}

	llh
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

	dwalls = c(0,diff(walls))
	rt = data.frame(wall=walls,dwall=dwalls,cpu=cpus)

	i1 = grep("TIME OF START *=",nd)
	attr(rt,"start") = as.difftime(gsub("^ *TIME OF START *= *","",nd[i1]),units="secs")

	rt
}

readhist = function(nd)
{
	l = strsplit(sub("^ +","",nd),split=" +")
	m = sapply(l,as.numeric)
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
	noms = gsub(".*node","",fnode,ignore.case=TRUE)
}

nflevg = getvar("NFLEVG",nd)

cat("Values of Var QC\n")
qc = varqc(nd)
if (! is.null(qc)) {
	qc1 = qc[[1]]
	types = seq(dim(qc1)[1])
	q = qcheck(qc)

	png(sprintf("%s/varqc.png",cargs$png))
	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	barplot(t(q$ra),names.arg=types,main="RAQC, variable 1",xlab="Obs type",ylab="RAQC")
	barplot(t(q$rl),names.arg=types,main="RLQC, variable 1",xlab="Obs type",ylab="RLQC")
	nc = dim(qc1)[2]
	#barplot(q$raqc,names.arg=types,beside=TRUE,main="",xlab="Obs type",ylab="RBGQC(1:3)")
	matplot(q$bg[,,1],main="Variable 1",xlab="Obs type",ylab="RBGQC(1:3)",type="o",lty=1,pch=20)
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
			plotvmean(jcs[,i,],1:nflevg,main="Jc after DFI",col=cols,lty=1,xlab="Jc",ylab="Level",
				legend=noms)
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
	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	if (length(fnode) > 1) {
		for (i in 1:min(dim(jacob)[2],4)) {
			plotmean(iter,jacobs[,i,],jac[i],noms,xlab="Iteration",ylab=jac[i],lty=1,pch=20)
		}
	} else {
		for (i in 1:min(dim(jacob)[2],4)) {
			plot(iter,jacob[,i],type="o",main=jac[i],xlab="Iteration",ylab=jac[i],pch=20)
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

	op = par(mfrow=c(2,2),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	if (length(fnode) > 1) {
		plotmean(iter,rits[,1,],noms,main="Ritz value 1",xlab=xlab,ylab="Ritz 1",pch=20)
		plotmean(iter,rits[,2,],noms,main="Ritz value 2",xlab=xlab,ylab="Ritz 2",pch=20)
		plotmean(iterg,grads,noms,main="Norm of gradient",xlab=xlab,ylab="Gradient",pch=20)
		plotmean(iterq,quads,noms,main="Quadratic cost J",xlab=xlab,ylab="Cost function",
			pch=20)
	} else {
		plot(iter,ritz[,1],type="o",main="Ritz value 1",xlab=xlab,ylab="Ritz 1",pch=20)
		plot(iter,ritz[,2],type="o",main="Ritz value 2",xlab=xlab,ylab="Ritz 2",pch=20)
		plot(iterg,cost$grad,type="o",main="Norm of gradient",xlab=xlab,ylab="Gradient",
			pch=20)
		plot(iterq,cost$quad,type="o",main="Quadratic cost J",xlab=xlab,
			ylab="Cost function",pch=20)
	}

	dev.off()
}

cat("CANARI statistics\n")
llh = canawagons(nd)
if (! is.null(llh)) {
	for (i in seq(along=llh)) {
		lh = llh[[i]]
		types = sapply(lh,attr,"type")

		cat(". types:",unique(types),"\n")
		for (t in unique(types)) {
			ficpng = sprintf("%s/residu%d_%s.png",cargs$png,i,t)
			png(ficpng)

			indt = which(types == t)
			nr = min(3,length(indt))
			nc = (length(indt)-1)%/%nr+1
			par(mfrow=c(nr,nc),mgp=c(2,1,0))

			for (h in lh[indt]) {
				tt = c(attr(h,"param"),sprintf("type: %s",attr(h,"type")))
				plot(h[,1],h[,2],type="s",main=tt,xlab="Residual",ylab="Frequency")
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
nstop = getvar("NSTOP",nd)
tstep = getvar("TSTEP",nd)
if (length(grep("STEP +\\d+ +H=.+\\+CPU=",nd)) > 0) {
	rt = runtime(nd)
	nts = dim(rt)[1]

	t0 = round(as.numeric(rt$wall[1]-attr(rt,"start"),units="secs") %% 86400,3)
	tint = round(as.numeric(rt$wall[nts]-rt$wall[1],units="secs") %% 86400,3)
	total = round(as.numeric(rt$wall[nts]-attr(rt,"start"),units="secs") %% 86400,3)
	tt = sprintf("setup+step0, forecast, total: %gs, %gs, %gs",t0,tint,total)
	cat(tt,"\n")

	if (length(fnode) > 1) {
		rts = lapply(nds,runtime)
		rts = c(list(rt),rts)
		rts = sapply(rts,function(rt) sapply(rt,as.numeric,simplify="array"),simplify="array")
	}

	pas = c(1,2,3,5,6,10,15,20,40,50)
	nt = c(0,200,400,600,1200,2000,3000,4000,8000,10000)
	it = findInterval(nts,nt)

	png(sprintf("%s/runtime.png",cargs$png))
	op = par(mfrow=c(2,1),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
	its = seq(1,nts,by=pas[it])
	if (length(fnode) > 1) {
		matplot(its-1,rts[its,2,],type="h",lty=1,main=c("Wall-time",tt),xlab="Time-step",
			ylab="Time (s)")
		legend("topleft",noms,lty=1)
		matplot(its-1,rts[its,3,],type="h",lty=1,main="CPU-time",xlab="Time-step",
			ylab="Time (s)")
		legend("topleft",noms,lty=1)
	} else {
		plot(its-1,rt$dwall[its],type="h",main=c("Wall-time",tt),xlab="Time-step",
			ylab="Time (s)")
		plot(its-1,rt$cpu[its],type="h",main="CPU-time",xlab="Time-step",ylab="Time (s)")
	}
	dev.off()
}
