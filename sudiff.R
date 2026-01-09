library(mfnode)

Glog = "\\.?(F(ALSE)?|T(RUE)?)\\>\\.?"

splitkeyval = function(nd)
{
	rekey = "\\<[A-Z]\\w+\\>"
	reval = sprintf("%s|%s|%s",Gnum,Gint,Glog)

	re = sprintf("^\\s*(%s)[ /]+(%s)[ /]+(%s) *[:=] *(%s) +(%s) +(%s) *$",rekey,rekey,
		rekey,reval,reval,reval)
	ind = grep(re,nd)
	if (length(ind) > 0) {
		l = regmatches(nd[ind],regexec(re,nd[ind]))
		for (i in seq(along=ind)) {
			s = l[[i]]
			stopifnot(all(nzchar(s[c(2:5,11,17)])))
			nd[ind[i]] = paste(s[2],"=",s[5],s[3],"=",s[11],s[4],"=",s[17])
		}
	}

	re = sprintf("^\\s*(%s)[ /]+(%s) *[:=] *(%s) +(%s) *$",rekey,rekey,reval,reval)
	ind = grep(re,nd)
	if (length(ind) > 0) {
		l = regmatches(nd[ind],regexec(re,nd[ind]))
		for (i in seq(along=ind)) {
			s = l[[i]]
			stopifnot(all(nzchar(s[c(2:4,10)])))
			nd[ind[i]] = paste(s[2],"=",s[4],s[3],"=",s[10])
		}
	}

	nd
}

getkeys = function(nd,logical=FALSE)
{
	if (logical) {
		re = sprintf("\\s*(\\<[A-Z]\\w+) *[:=] *%s",Glog)
	} else {
		re = sprintf("\\s*(\\<[A-Z]\\w+) *[:=] *(%s|%s)",Gnum,Gint)
	}

	ind = grep("end of setup",nd,ignore.case=TRUE)
	if (length(ind) > 0) nd = nd[1:ind[1]]
	s = paste(nd,collapse="\n")
	keyval = sub("^\\s*","",regmatches(s,gregexpr(re,s))[[1]])
	keyval = gsub(" *[:=] *"," = ",keyval)
	lval = strsplit(keyval," *= *")
	key = sapply(lval,"[",1)
	val = sapply(lval,"[",2)
	if (! logical) val = as.numeric(val)
	data.frame(key=key,val=val,s=keyval)
}

args = commandArgs(TRUE)
nd = readLines(args[1],skipNul=TRUE)
nd = splitkeyval(nd)
df = getkeys(nd)
dfl = getkeys(nd,TRUE)

nd = readLines(args[2],skipNul=TRUE)
nd = splitkeyval(nd)
df1 = getkeys(nd)
dfl1 = getkeys(nd,TRUE)

dfu = df[! duplicated(df$key),]
df1u = df1[! duplicated(df1$key),]
dflu = dfl[! duplicated(dfl$key),]
dfl1u = dfl1[! duplicated(dfl1$key),]
cat("-->",dim(dfu)[1],"/",dim(df)[1],"unique keys in ref\n")
cat("-->",dim(dflu)[1],"/",dim(dfl)[1],"unique lkeys in ref\n")
cat("-->",dim(df1u)[1],"/",dim(df1)[1],"unique keys in new\n")
cat("-->",dim(dfl1u)[1],"/",dim(dfl1)[1],"unique lkeys in new\n")
ind = match(dfu$key,df1u$key)
indl = match(dflu$key,dfl1u$key)
ina = is.na(ind)
ilna = is.na(indl)

out = "."
if (length(args) > 2) out = args[3]
if (! file.exists(out)) dir.create(out,recursive=TRUE)

fk = sprintf("%s/%s.refkeys",out,basename(args[1]))
fk1 = sprintf("%s/%s.keys",out,basename(args[2]))
con = file(fk,"w+")
con1 = file(fk1,"w+")

cat("Write",length(ind[!ina]),"+",length(indl[!ilna]),"keys in common\n")
if (any(! ina)) {
	cat("# keys in common\n",file=con)
	write.table(dfu[which(!ina),1:2],con,quote=FALSE,row.names=FALSE,col.names=FALSE)
	cat("# keys in common\n",file=con1)
	write.table(df1u[na.omit(ind),1:2],con1,quote=FALSE,row.names=FALSE,col.names=FALSE)
}

if (any(! ilna)) {
	cat("# lkeys in common\n",file=con)
	write.table(dflu[which(!ilna),1:2],con,quote=FALSE,row.names=FALSE,col.names=FALSE)
	cat("# lkeys in common\n",file=con1)
	write.table(dfl1u[na.omit(indl),1:2],con1,quote=FALSE,row.names=FALSE,col.names=FALSE)
}

cat("Write",length(ind[ina]),"+",length(indl[ilna]),"keys in ref only\n")
if (any(ina)) {
	cat("# keys in ref only\n",file=con)
	write.table(dfu[which(ina),1:2],con,quote=FALSE,row.names=FALSE,col.names=FALSE)

	cat("# keys in new only\n",file=con1)
	if (length(ind) == 0 || all(ina)) {
		cat("Write keys in new only (all)\n")
		write.table(df1u[,1:2],con1,quote=FALSE,row.names=FALSE,col.names=FALSE)
	} else {
		cat("Write",dim(df1u)[1]-length(which(!ina)),"keys in new only\n")
		write.table(df1u[-na.omit(ind),1:2],con1,quote=FALSE,row.names=FALSE,
			col.names=FALSE)
	}
}

if (any(ilna)) {
	cat("# lkeys in ref only\n",file=con)
	write.table(dflu[which(ilna),1:2],con,quote=FALSE,row.names=FALSE,
		col.names=FALSE)

	cat("# lkeys in new only\n",file=con1)
	if (length(indl) == 0 || all(ilna)) {
		cat("Write lkeys in new only (all)\n")
		write.table(dfl1u[,1:2],con1,quote=FALSE,row.names=FALSE,col.names=FALSE)
	} else {
		cat("Write",dim(dfl1u)[1]-length(which(!ilna)),"lkeys in new only\n")
		write.table(dfl1u[-na.omit(indl),1:2],con1,quote=FALSE,row.names=FALSE,
			col.names=FALSE)
	}
}

close(con)
close(con1)

