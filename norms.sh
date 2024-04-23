#!/bin/sh

node=~/util/node

set -e

usage()
{
	printf "
Description:
	Produce an HTML file showing plots of SP/GP norms from model forecast (NODE files)

Usage:
	norms.sh FILE1 FILE2... -o HTML [-lev (LEV|I1:I2:...)] [-detail] [-nogfl] [-nogmv] \
[-noadiab] [-noslb2] [-opt OPTIONS] [-h]

Options:
	FILE1, FILE2,...: NODE files containing SP and/or GP norms to plot (see Details)
	HTML: output file containing HTML content and links to PNG files
	-lev LEV: produce plots for level LEV, an index among the model levels. Default is 0 \
and stands for the summary level (ie level 'AVE'). If LEV is -1, plots are vertical \
profiles.
	-lev I1:I2:...: produce plots for groups of levels
	-detail: activate 'detailed' plots for gpnorm (graphics contain 4 plots instead of 2)
	-nogfl,nogmv,...: deactivate plots for this kind of norms. SP norms cannot!
	OPTIONS: options to pass to the R script, in the form key=value, separated by space \
and with keys among hmin, hmax (see Options)
	-h: print this help and exit normally

Details:
	NODE file names must begin with 'NODE' (or 'node'). One file at least is mandatory.
	Plots are produced in one single directory named after the 1st NODE file, its prefix \
'NODE'/'node' being deleted. If several NODE files are passed, plots show curves for \
all the files but for main norms only. Il one single file is passed, plots show curves \
for all the norms printed.
	A single HTML file is produced, containing text and images from all NODE files.
	When option '-lev' is passed, level norms must exist in NODE files. For option \
in form '-lev I1:I2:...', values are indicated by integer values, separated \
with colon. These values indicate 'breaks', where group i contains levels from 1 to i, \
but not belonging to preceding groups (ie groups do not overlap).

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
fin2=""
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
		[ -z "$fin" ] && fin=$1 || fin2="$fin2 $1"
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

if ! grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $fin
then
	echo "$fin: no 'STEP...' in $fin" >&2
	exit 1
fi

if [ -d $fout ]
then
	fout=$fout/norms.html
	echo "HTML output file: $fout"
fi

loc=$(dirname $fout)

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1
env | grep -qw R_LIBS || export R_LIBS=~/lib

fic=$(echo $fin | sed -re 's/^\s*(node\.?\w+[^ :]*) *:?/\1/i')
ls -L $fic > /dev/null

if ! { file -L $fic | grep -qE "(ASCII|UTF-8 Unicode) text" &&
	grep -iqE '^ \w+:\w+:\w+ +STEP +[0-9]+' $fic; }
then
	echo "--> file not text nor model forecast" >&2
	continue
fi

dd=$loc/$(echo $fic | sed -re 's:^(.+/)?node\.?(\w+):\2:i')
mkdir -p $dd
echo "--> output sent to $dd"

grep -q 'END CNT0' $fic || echo "Warning: no 'END CNT0', program may crash" >&2

echo "SP norms for SPEC"
R --slave -f $node/spnorms.R --args $fic $fin2 lev=$lev \
	spre="spnorm t1 *tr(ansdir)?:spnorm t1si:spnorm t1 spcm" png=$dd $ropt

if echo $norms | grep -q gfl
then
	echo "GP norms for GFL"
	R --slave -f $node/gpnorms.R --args $fic $fin2 lev=$lev type=gpgfl$suf \
		gpre="gpnorm gflt1 (call_)?sl$:gpnorm gflt1 slmf:gpnorm gflt1 (cpg)?lag:gpnorm gfl tstep" png=$dd $ropt
fi

if grep -iqE "gpnorm gmvt0" $fic && echo $norms | grep -q gmv
then
	echo "GP norms for GMV"
	R --slave -f $node/gpnorms.R --args $fic $fin2 lev=$lev type=gpgmv$suf \
		gpref="gpnorm gmvt0" gpre="gpnorm gmvt1 sl:gpnorm gmvt1 slmf:gpnorm gmvt1 cpglag" \
		png=$dd $ropt
fi

if grep -iqE "gpnorm adiab" $fic && echo $norms | grep -q adiab
then
	echo "GP norms for ADIAB"
	R --slave -f $node/gpnorms.R --args $fic $fin2 lev=$lev type=gpadiab$suf \
		gpref="gpnorm adiab" png=$dd $ropt
fi

if grep -iqE "gpnorm zb2" $fic && echo $norms | grep -q slb2
then
	echo "GP norms for ZB2"
	R --slave -f $node/gpnorms.R --args $fic $fin2 lev=$lev type=gpsi$suf \
		gpref="gpnorm zb2 cpg$" gpre="gpnorm zb2 (call_)?sl$:gpnorm zb2 slmf" png=$dd \
		$ropt
fi

for pre in sp gpgmv$suf gpgfl$suf gpadiab$suf gpsi$suf
do
	echo "HTML files for $pre"
	for s in "" v
	do
		echo "<tr>"

		n=0
		for ficp in $(ls -1 $dd | grep -E "${pre}norm$s[[:digit:]]+\.png")
		do
			if [ $n -eq 2 ]
			then
				echo -e "</tr>\n<tr>\n"
				n=0
			fi

			n=$((n+1))
			printf "\t<td><img src=\"%s\"/></td>\n" $dd/$ficp
		done

		if [ -s $dd/$pre$s.txt ]
		then
			echo -e "\t<td><pre>\n"
			cat $dd/$pre$s.txt
			echo "</pre></td>"
		fi

		echo "</tr>"
	done > $temp/$pre.html

	[ $(wc -l $temp/$pre.html | awk '{print $1}') -le 2 ] && echo > $temp/$pre.html
done

echo "Adding thematic HTML files from $temp"
date=$(grep -E 'NUDATE *=' $fic | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
res=$(grep -E 'NUDATE *=' $fic | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
base=$(printf "%s %dh" $date $((res/3600)))
sed -re "s:TAG NODE:$fic:" -e "s:TAG BASE:$base:" \
	-e "/TAG SP/r $temp/sp.html" -e "/TAG GPGMV/r $temp/gpgmv$suf.html" \
	-e "/TAG GPGFL/r $temp/gpgfl$suf.html" -e "/TAG GPADIAB/r $temp/gpadiab$suf.html" \
	-e "/TAG GPSI/r $temp/gpsi$suf.html" $node/norms.html > $fout
