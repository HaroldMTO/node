#!/bin/sh

diags=~/util/diags

set -e

usage()
{
	printf "
Description:
	Produce an HTML file related to plots of SP/GP norms from model forecast (NODE files)

Usage:
	plot.sh FILE -o HTML [-h]

Options:
	FILE: text file containing lines made of NODE file names and their description
	HTML: output file containing HTML content (text and links to images)
	-h: print this help and exit normally

Details:
	Lines in the input text file must match 'NODEFILE: description...'. Within lines, \
blanks are ignored and colon (':') is the separator. Case is also ignored for NODE \
filename. These NODE file names must begin with 'NODE' (or 'node'). Description is what \
follows the separator ':'.
	Plots are produced in directories named after NODE files, prefix 'NODE'/'node' being \
deleted.
	A single HTML file is produced, containing text and images from these NODE files.
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

ls -L $fin > /dev/null
file -L $fin | grep -q text

temp=$(mktemp -d -t plotXXX)

trap 'rm -r $temp' 0

if grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $fin
then
	# make fic local (no path management in info.txt)
	cd $(dirname $fin) > /dev/null
	fin=$(basename $fin)
	echo "$fin: no description" > $temp/info.txt
elif ! grep -Ei '^\s*node\.?\w+ *: *\w' $fin > $temp/info.txt
then
	echo "Error: pattern 'NODExxx: ...' not found in $fin" >&2
	exit 1
fi

loc=$(dirname $fin)
xpdir=$(cd $loc > /dev/null && pwd)

cd $temp > /dev/null

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

while read -a ff
do
	echo "Line read: '${ff[*]}'"
	fic=$(echo ${ff[0]} | sed -re 's/^\s*(node\.?\w+) *:?/\1/i')
	desc=$(echo "${ff[*]}" | sed -re 's/.+: +//')
	ls -L $xpdir/$fic > /dev/null

	if ! { file -L $xpdir/$fic | grep -qE "(ASCII|UTF-8 Unicode) text" &&
		grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $xpdir/$fic; }
	then
		echo "--> file not text nor model forecast" >&2
		continue
	fi

	dd=$(echo $fic | sed -re 's:^node\.?(\w+):\1:i')
	mkdir -p $xpdir/$dd
	echo "--> output sent to $loc/$dd"

	grep -q 'END CNT0' $xpdir/$fic || echo "Warning: no 'END CNT0', program may crash" >&2

	rm -f sp.txt
	R --slave -f $diags/norms.R --args $xpdir/$fic lev=0 plot=sp \
		spre="spnorm .*t1(tr)?\>:spnorm .*t1si" png=png > sp.txt

	if grep -qE "gpnorm gmvt0" $xpdir/$fic
	then
		rm -f gpgmv
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpgmv \
			gpre="gpnorm gmvt1 sl:gpnorm gmvt1 lag" png=png > gpgmv.txt
	fi

	if grep -qE "gpnorm adiab" $xpdir/$fic
	then
		rm -f gpadiab.txt
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpadiab \
			gpref="gpnorm adiab" png=png > gpadiab.txt
	fi

	if grep -qE "gpnorm zb2" $xpdir/$fic
	then
		rm -f gpsi.txt
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=0 plot=gp type=gpsi \
			gpref="gpnorm zb2 cpg" gpre="gpnorm zb2 sl" png=png > gpsi.txt
	fi

	if ! grep -qE "NGDITS *= *-?[12]\>" $xpdir/$fic
	then
		rm -f gpgfl.txt
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=0 plot=gp png=png > gpgfl.txt
	fi

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

			if [ -s $pre.txt ]
			then
				echo -e "\t<td><pre>\n"
				grep -ivE "Read|no X11 device|(GP|Spectral) norms" $pre.txt
				echo "</pre></td>"
			fi

			echo "</tr>"
		} > $pre.html
	done

	[ $(wc -l $pre.html | awk '{print $1}') -le 2 ] && echo > $pre.html

	date=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
	res=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
	base=$(printf "%s %dh" $date $((res/3600)))
	sed -re "s:TAG NODE:$fic:" -e "s:TAG BASE:$base:" -e "s:TAG DESC:$desc:" \
		-e '/TAG SP/r sp.html' -e '/TAG GPGMV/r gpgmv.html' -e '/TAG GPGFL/r gpgfl.html' \
		-e '/TAG GPADIAB/r gpadiab.html' -e '/TAG GPSI/r gpsi.html' \
		$diags/img.html >> img.html
done < info.txt

sed -re '/TAG IMG/r img.html' $diags/norms.html > out.html

cd $OLDPWD > /dev/null
mv $temp/out.html $fout
