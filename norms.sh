#!/bin/sh

diag=~/util/diag

set -e

usage()
{
	printf "
Description:
	Print HTML output related to various plots (SP/GP norms) from model forecast

Usage:
	plot.sh FILE [-h]

Options:
	FILE: text file containing lines NODE file names and their description
	-h: print this help and exit normally

Details:
	NODE file names must begin with 'node'. Description is what follows this name, \
separated with ':'.
	Plots are produced in directories named after NODE files, prefix 'node' being deleted.
	An HTML output is produced, containing text and images from these NODE files.
"
}

if [ $# -eq 0 ] || echo $* | grep -qE '(^| )\-h\>'
then
	usage
	exit
fi

file $1 | grep -q text

temp=$(mktemp -d -t plotXXX)

trap 'rm -r $temp' 0

grep -E '^ *node\w+ *: *\w' $1 > $temp/info.txt
loc=$(dirname $1)
xpdir=$(cd $loc > /dev/null && pwd)

cd $temp > /dev/null

while read -a ff
do
	#echo "Line read: '${ff[*]}'"
	fic=$(echo ${ff[0]} | sed -re 's/ *://')
	file $xpdir/$fic | grep -qE "(ASCII|UTF-8 Unicode) text" || continue

	tt=$(echo "${ff[*]}" | sed -re 's/.+: +//')

	echo $fic | grep -qE '\<node\w' || continue
	grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $xpdir/$fic || continue

	dd=$(echo $fic | sed -re 's:^node::')
	mkdir -p $xpdir/$dd

	R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=sp \
		spre="spnorm spect1\>:spnorm spect1si" > sp.txt
	convert Rplots.pdf spnorm.png

	if grep -qE "gpnorm gmvt0" $xpdir/$fic
	then
		R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpgmv \
			gpre="gpnorm gmvt1 sl:gpnorm gmvt1 lag" > gmv.txt
		convert Rplots.pdf gpgmvnorm.png
	fi

	R --slave -f $diag/norms.R --args $xpdir/$fic lev=0 plot=gp > gfl.txt
	convert Rplots.pdf gpgflnorm.png

	mv *.png $xpdir/$dd

	for pre in sp gpgmv gpgfl
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
		$diag/img.html >> img.html
done < info.txt

sed -re '/TAG IMG/r img.html' $diag/plots.html
