#!/bin/sh

diags=~/util/diags

set -e

usage()
{
	printf "
Description:
	Produce an HTML file related to plots of SP/GP norms from model forecast (NODE files)

Usage:
	plot.sh FILE -o HTML [-lev (LEV|I1:I2:...)] [-detail] [-opt OPTIONS] [-h]

Options:
	FILE: text file containing lines made of NODE file names and their description
	HTML: output file containing HTML content (text and links to images)
	-lev LEV: produce plots for level LEV, an index among the model levels. Default is 0 \
and stands for the summary level (ie level 'AVE')
	-lev I1:I2:...: produce plots for groups of levels
	-detail: activate 'detailed' plots for gpnorm (graphics contain 4 plots instead of 2)
	OPTIONS: options to pass to the R script, in the form key=value, separated by space \
and with keys among hmin, hmax (see Options)
	-h: print this help and exit normally

Details:
	Lines in the input text file must match 'NODEFILE: description...'. Within lines, \
blanks are ignored and colon (':') is the separator. Case is also ignored for NODE \
filename. These NODE file names must begin with 'NODE' (or 'node'). Description is what \
follows the separator ':'.
	Plots are produced in directories named after NODE files, prefix 'NODE'/'node' being \
deleted.
	A single HTML file is produced, containing text and images from these NODE files.
	When option '-lev' is passed, level norms must exist in NODE files. \
For option in form '-lev I1:I2:...', values are indicated by integer values, separated \
with colon. Values indicate 'breaks', where group i contains levels from 1 to i, but \
not belonging to preceding groups (ie groups do not overlap).

Options:
	Options passed with argument '-opt' consist in 1 character string (-> use quotes) \
to be passed to R script norms.R. Each option is a key/value couple just like \
'key1=value key2=value...'. Currently, keys can be hmin or hmax, setting min and max \
values (in hours) to plot on time-step diagrams. In other words, the full range of time \
found in norms from input files is focused between the given min and/or max hours.
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

norms="gmvgfladiabslb2"
fin=""
fout=""
suf=""
lev=0
ropt=""

while [ $# -ne 0 ]
do
	case $1 in
	-o)
		fout=$2
		shift
		;;
	-lev)
		lev=$2
		shift
		;;
	-detail) suf=detail;;
	-nogfl) norms=$(echo $norms | sed -re 's:gfl::');;
	-nogmv) norms=$(echo $norms | sed -re 's:gmv::');;
	-noadiab) norms=$(echo $norms | sed -re 's:adiab::');;
	-noslb2) norms=$(echo $norms | sed -re 's:slb2::');;
	-opt)
		ropt=$2
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
elif ! grep -Ei '^\s*node\.?\w+ *:? *\w' $fin > $temp/info.txt
then
	echo "Error: pattern 'NODExxx: ...' not found in $fin" >&2
	exit 1
fi

loc=$(dirname $fin)
xpdir=$(cd $loc > /dev/null && pwd)

cd $temp > /dev/null

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1
env | grep -qw R_LIBS || export R_LIBS=~/lib

while read fnode desc
do
	echo "Line read: '$fnode' '$desc'"
	fic=$(echo $fnode | sed -re 's/^\s*(node\.?\w+) *:?/\1/i')
	ls -L $xpdir/$fic > /dev/null

	if ! { file -L $xpdir/$fic | grep -qE "(ASCII|UTF-8 Unicode) text" &&
		grep -iqE '^ \w+:\w+:\w+ +STEP +[0-9]+' $xpdir/$fic; }
	then
		echo "--> file not text nor model forecast" >&2
		continue
	fi

	dd=$(echo $fic | sed -re 's:^node\.?(\w+):\1:i')
	mkdir -p $xpdir/$dd
	echo "--> output sent to $loc/$dd"

	grep -q 'END CNT0' $xpdir/$fic || echo "Warning: no 'END CNT0', program may crash" >&2

	echo "SP norms for SPEC"
	R --slave -f $diags/norms.R --args $xpdir/$fic lev=$lev type=spec \
		spre="spnorm t1 *tr(ansdir)?:spnorm t1si:spnorm t1 spcm" png=$xpdir/$dd $ropt

	if echo $norms | grep -q gfl
	then
		echo "GP norms for GFL"
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=$lev type=gpgfl$suf \
			gpre="gpnorm gflt1 (call_)?sl$:gpnorm gflt1 slmf:gpnorm gflt1 (cpg)?lag:gpnorm gfl tstep" png=$xpdir/$dd $ropt
	fi

	if grep -iqE "gpnorm gmvt0" $xpdir/$fic && echo $norms | grep -q gmv
	then
		echo "GP norms for GMV"
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=$lev type=gpgmv$suf \
			gpref="gpnorm gmvt0" gpre="gpnorm gmvt1 sl:gpnorm gmvt1 slmf:gpnorm gmvt1 cpglag" \
			png=$xpdir/$dd $ropt
	fi

	if grep -iqE "gpnorm adiab" $xpdir/$fic && echo $norms | grep -q adiab
	then
		echo "GP norms for ADIAB"
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=$lev type=gpadiab$suf \
			gpref="gpnorm adiab" png=$xpdir/$dd $ropt
	fi

	if grep -iqE "gpnorm zb2" $xpdir/$fic && echo $norms | grep -q slb2
	then
		echo "GP norms for ZB2"
		R --slave -f $diags/norms.R --args $xpdir/$fic lev=$lev type=gpsi$suf \
			gpref="gpnorm zb2 cpg$" gpre="gpnorm zb2 (call_)?sl$:gpnorm zb2 slmf" png=$xpdir/$dd \
			$ropt
	fi

	for pre in sp gpgmv$suf gpgfl$suf gpadiab$suf gpsi$suf
	do
		echo "HTML files for $pre"
		for s in "" v
		do
			echo "<tr>"

			n=0
			for ficp in $(ls -1 $xpdir/$dd | grep -E "${pre}norm$s[[:digit:]]+\.png")
			do
				if [ $n -eq 2 ]
				then
					echo -e "</tr>\n<tr>\n"
					n=0
				fi

				n=$((n+1))
				printf "\t<td><img src=\"%s\"/></td>\n" $loc/$dd/$ficp
			done

			if [ -s $xpdir/$dd/$pre$s.txt ]
			then
				echo -e "\t<td><pre>\n"
				cat $xpdir/$dd/$pre$s.txt
				echo "</pre></td>"
			fi

			echo "</tr>"
		done > $pre.html

		[ $(wc -l $pre.html | awk '{print $1}') -le 2 ] && echo > $pre.html
	done

	echo "Adding thematic HTML files"
	date=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
	res=$(grep -E 'NUDATE *=' $xpdir/$fic | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
	base=$(printf "%s %dh" $date $((res/3600)))
	sed -re "s:TAG NODE:$fic:" -e "s:TAG BASE:$base:" -e "s:TAG DESC:$desc:" \
		-e '/TAG SP/r sp.html' -e "/TAG GPGMV/r gpgmv$suf.html" \
		-e "/TAG GPGFL/r gpgfl$suf.html" -e "/TAG GPADIAB/r gpadiab$suf.html" \
		-e "/TAG GPSI/r gpsi$suf.html" $diags/img.html >> img.html
done < info.txt

sed -re '/TAG IMG/r img.html' $diags/norms.html > out.html

cd $OLDPWD > /dev/null
mv $temp/out.html $fout

