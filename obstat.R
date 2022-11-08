readDef = function(nd)
{
	linfo = vector("list",9)

	info = c("statkind","surfaces","types","params","items","instrument","flagfilter")
	names(linfo) = info
	ind = c(2,4:9)
	for (i in seq(along=info)) {
		val = sub(sprintf("^ *%s *= *([^']+) *",info[i]),"\\1",nd[ind[i]])
		linfo[[i]] = as.integer(strsplit(sub(" +\\(.+\\)","",val)," +")[[1]])
	}

	sinfo = c("comment","areaNSEW")
	ind = c(1,3)
	ni = 7
	names(linfo)[-(1:ni)] = sinfo
	for (i in seq(along=sinfo)) {
		val = sub(sprintf("^ *%s *= *'?([^']+)'? *",sinfo[i]),"\\1",nd[ind[i]])
		linfo[[ni+i]] = sub("(.+) +\\(.+\\)","\\1",val)
	}

	stopifnot(nzchar(linfo$comment))
	stopifnot(nzchar(linfo$areaNSEW))

	if (length(nd) == length(linfo)) {
		return(linfo)
	} else if (length(nd) == length(linfo)+1) {
		binfo = "nbbin"
		ind = 10
	} else if (length(nd) == length(linfo)+3) {
		binfo = c("sizebin","refval","nbbin")
		ind = 10:12
	} else if (length(nd) == length(linfo)+4) {
		binfo = c("nbbin","coorditem1","sizebin1?","refval1?")
		ind = 10:13
	} else {
		stop("unknown statdef")
	}

	ni = length(linfo)
	for (i in seq(along=binfo)) {
		val = sub(sprintf("^ *%s *= *'?([^']+)'? *",binfo[i]),"\\1",nd[ind[i]])
		linfo[[ni+i]] = as.numeric(val)
	}

	names(linfo)[-(1:ni)] = binfo
	linfo$nbbin = as.integer(linfo$nbbin)
	linfo
}

readItem = function(nd,stat)
{
	item = sub("^ *#item *: *(.+) *","\\1",nd[1])
	item = sub(" +$","",item)

	ind = 1

	if (stat == 2) {
		stopifnot(regexpr("^ *# *pop +min +max",nd[2]) > 0)
		val = as.numeric(strsplit(nd[3]," +")[[1]])
		ind = 1:3
	}

	con = file()
	writeLines(nd[-ind],con)
	df = read.table(con,header=TRUE,comment.char="")
	close(con)
	names(df)[1] = sub("X\\.","",names(df)[1])
	list(item=item,data=df)
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

args = commandArgs(trailingOnly=TRUE)

hasx11 = length(args) < 2 && capabilities("X11")
if (! hasx11) cat("--> no X11 device, sending plots to PNG files\n")
ask = hasx11 && interactive()

fobs = args[1]
pngd = "obstat"
if (length(args) == 2) pngd = args[2]

cat("Read file",fobs,"\n")
nd = readLines(fobs)

ind = grep("^ *BEGIN STATDEF",nd)+1
inde = grep("^ *END STATDEF",nd)-1
stopifnot(length(ind) == length(inde))

ii = which(ind > inde)
if (length(ii) > 0) cat("-->",length(ii),"nul def\n")

# for next "BEGIN STATDEF"
ind = c(ind,length(nd))

defs = vector("list",length(inde))

if (interactive()) browser()

for (i in seq(along=inde)) {
	def = readDef(nd[ind[i]:inde[i]])
	s = paste(def$params,collapse=":")
	if (length(def$params) != 1)
		cat(". def:",def$comment,def$areaNSEW,"- params/kind:",s,def$statkind,"\n")

	# items
	ii = (inde[i]+1):(ind[i+1]-1)
	indi = grep("^ *BEGIN STATITEM",nd[ii])+ii[1]
	indie = grep("^ *END STATITEM",nd[ii])+ii[1]-2
	inul = which(indi > indie)
	if (length(inul) > 0) {
		#cat("--> remove",length(inul),"nul items\n")
		def$items = def$items[-inul]
	}

	tt = paste(def$comment,def$areaNSEW)
	if (any(def$instrument != 999)) {
		tt[2] = sprintf("instrument(s): %s",paste(def$instrument,collapse=" "))
	}

	if (length(def$items) == 0) {
		#cat("--> no items\n")
		def$items = list()
		defs[[i]] = def
		next
	} else if (length(def$items) > 6) {
		#cat("--> limiting items to 6 out of",length(def$items),"\n")
		length(def$items) = 6
	}

	items = list()
	for (j in seq(along=def$items)) {
		items[[j]] = readItem(nd[indi[j]:indie[j]],def$statkind)
	}

	def$items = items

	defs[[i]] = def

	s = gsub(" ","-",def$comment)
	if (regexpr("-",s) < 0) s = sprintf("%s-X",s)
	ficpng = sprintf("%s/%s_%s.png",pngd,s,gsub("\\.","",def$areaNSEW))
	#cat(". file",ficpng,"\n")
	pngalt(ficpng)

	nc = max(1,length(def$items)%/%3)
	nr = length(def$items)%/%nc
	par(mfrow=c(nr,nc),mgp=c(2,1,0))

	if (def$statkind == 2) {
		for (j in seq(along=def$items)) {
			item = items[[j]]
			barplot(item$data$population,col="grey98",space=0,main=tt,xlab=item$item,
				ylab="population")
			axis(1)
		}
	} else {
		for (j in seq(along=def$items)) {
			item = items[[j]]
			ylim = range(item$data[,1])
			if (names(item$data)[1] == "Pressure") ylim = rev(ylim)
			matplot(item$data[,-(1:2)],item$data[,1],type="o",main=tt,xlab=item$item,
				ylab=names(item$data)[1],ylim=ylim,lty=c(0,0,2,1),pch=c("-","+",NA,NA),
				col=1)
			#abline(h=0)
			abline(v=0,col="darkgrey",lty=2)
		}
	}

	pngoff()
}

q("no")

comm = sapply(defs,"[[","comment")
area = sapply(defs,"[[","areaNSEW")
tt = paste(comm,area)

it = which(duplicated(tt))
if (length(it) > 0) {
	tt1 = unique(tt[it])
	idup = integer()

	for (i in seq(along=tt1)) {
		ind = which(tt == tt1[i])
		items = defs[[ind[1]]]$items
		ii = which(sapply(ind[-1],function(j) identical(items,defs[[j]]$items)))
		if (length(ii) > 0) {
			#cat("--> duplicates for index",ind[1],":",ind[ii+1],"\n")
			idup = c(idup,ind[ii+1])
		}
	}
}

if (length(idup) > 0) {
	#cat("--> duplicates:",idup,"\n")
	defs = defs[-idup]
}
