#!/bin/bash

diag=~/util/diag

set -e

usage()
{
	printf "
Description:
	Produce HTML output (stdout or file) related to various plots (SP/GP norms) from \
model forecast

Usage:
	plot.sh FILE [-o HTML] [-h]

Options:
	FILE: text file containing lines made of NODE file names and their description
	HTML: output file containing HTML content (text and links to images)
	-h: print this help and exit normally

Details:
	NODE file names must begin with 'node'. Description is what follows this name, \
separated with ':'.
	Plots are produced in directories named after NODE files, prefix 'node' being deleted.
	An HTML output is produced, containing text and images from these NODE files.
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

fin=""
fout=""
while [ $# -ne 0 ]
do
	case $1 in
	-o)
		fout=$2
		shift
		;;
	-h)
		usage
		exit
		;;
	*)
		fin=$1
		;;
	esac

	shift
done

if [ -z "$fin" -o -z "$fout" ]
then
	printf "Error: input option missing
fin: '$fin'
fout: '$fout'
" >&2
	exit 1
fi

file $fin | grep -q text

temp=$(mktemp -d -t plotXXX)

trap 'rm -r $temp' 0

if ! grep -Ei '^ *node\w+ *: *\w' $fin > $temp/info.txt
then
	echo "Error: pattern 'node...: ...' not found in $fin" >&2
	exit 1
fi

loc=$(dirname $fin)
xpdir=$(cd $loc > /dev/null && pwd)

cd $temp > /dev/null

while read -a ff
do
	echo "Line read: '${ff[*]}'"
	fic=$(echo ${ff[0]} | sed -re 's/ *://')
	file $xpdir/$fic | grep -qE "(ASCII|UTF-8 Unicode) text" || continue

	tt=$(echo "${ff[*]}" | sed -re 's/.+: +//')

	echo $fic | grep -qEi '\<node\w' || continue
	grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $xpdir/$fic || continue

	dd=$(echo $fic | sed -re 's:^node::i')
	mkdir -p $xpdir/$dd

	R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=sp \
		spre="spnorm .*t1\>:spnorm .*t1si"
	convert Rplots.pdf spnorm.png

	if grep -qE "gpnorm gmvt0" $xpdir/$fic
	then
		R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpgmv \
			gpre="gpnorm gmvt1 sl:gpnorm gmvt1 lag"
		convert Rplots.pdf gpgmvnorm.png
	fi

	if grep -qE "gpnorm adiab" $xpdir/$fic
	then
		R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpgmv \
			gpref="gpnorm adiab"
		convert Rplots.pdf gpadiabnorm.png
	fi

	if grep -qE "gpnorm zb2" $xpdir/$fic
	then
		R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpgmv \
			gpref="gpnorm zb2 cpg" gpre="gpnorm zb2 sl"
		convert Rplots.pdf gpsinorm.png
	fi

	R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp
	convert Rplots.pdf gpgflnorm.png

	mv *.png $xpdir/$dd

	for pre in sp gpgmv gpgfl gpadiab gpsi
	do
		{
			echo "<tr>"

			n=0
			for ficp in $(ls -1 $xpdir/$dd | grep -E "${pre}norm.*.png")
			do
				if [ $n -eq 2 ]
				then
					echo -e "</tr>\n<tr>\n"
					n=0
				fi

				n=$((n+1))
				printf "\t<td><img src=\"%s\"/></td>\n" $loc/$dd/$ficp
			done

			echo "</tr>"
		} > $pre.html
	done

	[ $(wc -l $pre.html | awk '{print $1}') -le 2 ] && echo > $pre.html

	date=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
	res=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
	base=$(printf "%s %dh" $date $((res/3600)))
	sed -re "s:TAG NODE:$fic:" -e "s:TAG BASE:$base:" -e "s:TAG DESC:$tt:" \
		-e '/TAG SP/r sp.html' -e '/TAG GPGMV/r gpgmv.html' -e '/TAG GPGFL/r gpgfl.html' \
		-e '/TAG GPADIAB/r gpadiab.html' -e '/TAG GPSI/r gpsi.html' \
		$diag/img.html >> img.html
done < info.txt

sed -re '/TAG IMG/r img.html' $diag/norms.html > out.html

cd $OLDPWD
mv $temp/out.html $fout
