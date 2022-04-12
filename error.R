args = commandArgs(trailingOnly=TRUE)

fic = args[1]
nl = as.integer(args[2])

nd = readLines(fic)

isnew = length(grep("(Vertical|Horizontal) +location of error",nd)) == 0

if (isnew) {
	indb = grep("^.+? out of physical bounds|^ *wind too strong at|velocity out at",nd)
	stopifnot(length(indb) == 1)

	isToaci = regexpr("Toaci",nd[indb]) > 0
	bb = regmatches(nd[indb],gregexpr("([-+]?[0-9]+([-+0-9.E]+[0-9]\\>)?|NaN)",nd[indb]))
	bb = as.numeric(bb[[1]])
	ilev = as.integer(bb[3])
}

ind = grep("^ *Vertical profile of",nd)

hasx11 = capabilities("X11")
ii = 0

for (i in seq(along=ind)) {
	if (regexpr("T qv",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","t","q","z"),nrows=nl,skip=ind[i])

		if (! hasx11) png("errortq.png")

		par(mfrow=c(1,2))
		plot(vv$t-273.15,vv$z,type="o",pch="+",main="Profile of T",xlab="T (°C)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		if (isToaci) {
			points(bb[4]-273.15+c(bb[5],0,bb[6]),rep(vv$z[ilev],3),type="o",pch="+",col=2,
				lty=2)
		}

		plot(vv$q*1000,vv$z,type="o",pch="+",main="Profile of Qv",xlab="Qv (g/kg)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		if (! isToaci) points(bb[3:4],rep(vv$z[ilev],2),type="o",pch="+",col=4,lty=2)
	} else if (regexpr("T\\>(.+)? qv",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("t","q","z","p","lev"),nrows=nl,skip=ind[i])

		if (! hasx11) png("errortq.png")

		par(mfrow=c(1,2))
		plot(vv$t-273.15,vv$z,type="o",pch="+",main="Profile of T",xlab="T (°C)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)

		plot(vv$q*1000,vv$z,type="o",pch="+",main="Profile of Qv",xlab="Qv (g/kg)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
	} else if (regexpr("U V T",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","u","v","t","spd","svd","nhx","z"),
			nrows=nl,skip=ind[i])

		ii = ii+1
		if (! hasx11) png(sprintf("errorgmv1_%d.png",ii))

		par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$u,vv$z,type="o",pch="+",main="Profile of U",xlab="U (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$v,vv$z,type="o",pch="+",main="Profile of V",xlab="V (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$t-273.15,vv$z,type="o",pch="+",main="Profile of T",xlab="T (°C)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)

		if (interactive()) readline("Press enter")

		if (! hasx11) {
			invisible(dev.off())
			png(sprintf("errorgmv2_%d.png",ii))
		}

		par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$spd,vv$z,type="o",pch="+",main="Profile of SPD",xlab="SPD (Pa)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$svd,vv$z,type="o",pch="+",main="Profile of SVD",xlab="SVD (m/s^2)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$nhx,vv$z,type="o",pch="+",main="Profile of NHX",xlab="NHX (m/s^2)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
	} else if (regexpr("NHPRE GW",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","nhpre","gw","rdphi","nhx","z"),
			nrows=nl,skip=ind[i])

		if (! hasx11) png("errornh.png")

		par(mfrow=c(1,4),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$nhpre,vv$z,type="o",pch="+",main="Profile of nhpre",xlab="nhpre (-)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$gw,vv$z,type="o",pch="+",main="Profile of Gw",xlab="Gw (?)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$rdphi,vv$z,type="o",pch="+",main="Profile of 1/dphi",xlab="rdphi (?)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$nhx,vv$z,type="o",pch="+",main="Profile of NHX",xlab="NHX (m/s^2)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
	} else if (regexpr("U V W Z",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","u","v","w","z"),nrows=nl,skip=ind[i])

		if (! hasx11) png("erroruv.png")

		par(mfrow=c(1,3),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$u,vv$z,type="o",pch="+",main="Profile of U",xlab="U (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$v,vv$z,type="o",pch="+",main="Profile of V",xlab="V (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$w,vv$z,type="o",pch="+",main="Profile of W",xlab="W (/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
	} else if (regexpr("Ut Vt Wt",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","u","v","w","etao","etaf","z"),nrows=nl,
			skip=ind[i])

		if (! hasx11) png("erroreta.png")

		par(mfrow=c(1,4),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$u,vv$z,type="o",pch="+",main="Profile of U",xlab="U (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$v,vv$z,type="o",pch="+",main="Profile of V",xlab="V (m/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$w,vv$z,type="o",pch="+",main="Profile of W",xlab="W (/s)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		matplot(vv[c("etaf","etao")],vv$z,type="o",col=c("black","darkgrey"),pch=c("+","o"),
			lty=c(1,0),main="Profile of Vdisp",xlab="eta (-)",ylab="Z (mgp)",log="y")
		segments(vv$etao,vv$z,vv$etaf,vv$z)
	} else if (regexpr("level cslat cslon",nd[ind[i]]) > 0) {
		vv = read.table(fic,col.name=c("lev","levO","latO","lonO","phiO","LatO","z"),
			nrows=nl,skip=ind[i])

		if (! hasx11) png("errorO.png")

		par(mfrow=c(1,5),mar=c(3,3,3,2)+.1,mgp=c(2,.75,0))
		plot(vv$levO,vv$z,type="p",pch="+",main="Level at O",xlab="level",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$latO,vv$z,type="o",pch="+",main="CS latitude",xlab="lat (°)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$lonO,vv$z,type="o",pch="+",main="CS longitude",xlab="long (°)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$phiO,vv$z,type="o",pch="+",main="angle OCF (phi)",xlab="phi (°)",
			ylab="Z (mgp)",log="y")
		abline(v=0,col="darkgrey",lty=2)
		plot(vv$LatO,vv$z,type="o",pch="+",main="GS latitude",xlab="Lat (°)",ylab="Z (mgp)",
			log="y")
		abline(v=0,col="darkgrey",lty=2)
	} else {
		stop("unknown profile type")
	}

	if (interactive()) readline("Press enter")

	if (! hasx11) invisible(dev.off())
}
