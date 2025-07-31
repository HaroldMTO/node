#!/bin/sh

node=~/util/node

set -e

usage()
{
	printf "
Description:
	Produce an HTML file showing plots of SP/GP norms from model forecast (NODE files)

Usage:
	norms.sh NODE1 NODE2... -o HTML [-lev (LEV|I1:I2:...)] [-detail] [-nogfl] [-nogmv] \
[-noadiab] [-noslb2] [-nospt1] [-nogpt1] [-not1] [-ref FILE] [-opt OPTIONS] [-h]

Options:
	NODE1, NODE2,...: NODE files containing SP and/or GP norms to plot (see Details)
	HTML: output file containing HTML content and links to PNG files
	-lev LEV: produce plots for level LEV, an index among the model levels. Default is 0 \
and stands for the summary level (ie level 'AVE'). If LEV is -1, plots are vertical \
profiles.
	-lev I1:I2:...: produce plots for groups of levels
	-detail: activate 'detailed' plots for gpnorm (graphics contain 4 plots instead of 2)
	-nogfl,nogmv,...: deactivate plots for this kind of norms. SP norms cannot!
	-nospt1/nogpt1/not1: deactivate plots for SP, GP or SP+GP norms at t1 (resp.)
	-ref: use FILE (an R image of SP/GP norms) as a reference to plot in addition to NODE \
file(s)
	OPTIONS: options to pass to the R script, in the form key=value, separated by space \
and with keys among hmin, hmax (see Options)
	-h: print this help and exit normally

Details:
	NODE file names must begin with 'NODE' (or 'node') after any path part. \
One file at least is mandatory.
	A single HTML file is produced, containing text and images from all NODE files.
	Plots are produced in one single directory named after the 1st NODE file, its prefix \
'NODE'/'node' being deleted (an exception is made for classical NODE.001_01: the directory \
name is the one of the output HTML file without extension 'html'). If several NODE files \
are passed, plots show curves for all the files but for main norms only. If one single \
file is passed, plots show curves for all the norms printed.
	When option '-lev' is passed, level norms must exist in NODE files. For option \
in form '-lev I1:I2:...', values are indicated by integer values, separated \
with colon. These values indicate 'breaks', where group i contains levels from 1 to i, \
but not belonging to preceding groups (ie groups do not overlap).
	When option '-ref' is passed, FILE is an R image file to be loaded in R scripts. It \
must load objects ddt, sps and gps. See runtimes.sh for more info.

Options:
	Options passed with argument '-opt' consist in 1 character string (-> use quotes) \
to be passed to R script norms.R. Each option is a key/value couple just like \
'key1=value key2=value...'. Currently, keys can be hmin or hmax, setting min and max \
values (in hours) to plot on time-step diagrams. In other words, the full range of time \
found in norms from input files is focused between the given min and/or max hours.
"
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
then
	usage
	exit
fi

norms="spt1gpt1gmvgfladiabslb2fp"
fin=""
fin2=""
fout=""
suf=""
lev=0
ropt=""
spref=""
gpref=""
ref=""

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
	-spref)
		spref="spref=$2"
		shift
		;;
	-gpref)
		gpref="gpref=$2"
		shift
		;;
	-detail) suf=detail;;
	-nogfl) norms=$(echo $norms | sed -re 's:gfl::');;
	-nogmv) norms=$(echo $norms | sed -re 's:gmv::');;
	-noadiab) norms=$(echo $norms | sed -re 's:adiab::');;
	-noslb2) norms=$(echo $norms | sed -re 's:slb2::');;
	-nofp) norms=$(echo $norms | sed -re 's:fp::');;
	-nospt1) norms=$(echo $norms | sed -re 's:spt1::');;
	-nogpt1) norms=$(echo $norms | sed -re 's:gpt1::');;
	-not1) norms=$(echo $norms | sed -re 's:[sg]pt1::g');;
	-ref)
		ref=$2
		shift
		;;
	-opt)
		ropt=$2
		shift
		;;
	-*)
		echo "Error: unknown option '$1'" >&2
		exit 1
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

if ! file -L $fin | grep -qE "(ASCII|UTF-8 Unicode) text"
then
	echo "Error: $fin is not a text file" >&2
	exit 1
elif ! grep -qEi '^ \w+:\w+:\w+ +STEP +[0-9]+' $fin
then
	echo "Warning: no 'STEP...' in $fin" >&2
	#exit 1
fi

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1
if ! env | grep -qw R_LIBS
then
	export R_LIBS=~petithommeh/lib
	echo "--> setting R_LIBS: $R_LIBS"
else
	R_LIBS=$R_LIBS:~petithommeh/lib
fi

loc=$(dirname $fout)
dd=$loc/$(echo $fin | sed -re 's:^(.+/)?node\.?(\w+):\2:i')
echo $dd | grep -qE '(.+/)?001_01$' && dd=$loc/$(basename $fout .html)
dd=$(echo $dd | sed -re 's:^\./::')
mkdir -p $dd
echo "--> output sent to $dd"

grep -q 'END CNT0' $fin || echo "Warning: no 'END CNT0', program may crash" >&2

[ -n "$ref" ] && ropt="$ropt ref=$ref"

echo "SP norms for SPEC"
spre=""
st1="spnorm t1"
echo $norms | grep -q spt1 && spre="$st1 *tr(ansdir)?:$st1 *si:$st1 spcm"
R --slave -f $node/spnorms.R --args $fin $fin2 lev=$lev $spref spre="$spre" png=$dd $ropt

if echo $norms | grep -q gfl
then
	echo "GP norms for GFL"
	gpre="gpnorm gfl tstep"
	gt1="gpnorm gflt1"
	echo $norms | grep -q gpt1 && gpre="$gt1 (call_)?sl$:$gt1 slmf:$gt1 (cpg)?lag:$gpre"
	R --slave -f $node/gpnorms.R --args $fin $fin2 lev=$lev $gpref type=gpgfl$suf \
		gpre="$gpre" png=$dd $ropt

	if [ -z "$fin2" ]
	then
		echo "GP norms for GFL+SPEC"
		R --slave -f $node/gpnorms.R --args $fin lev=$lev type=gflspec$suf png=$dd $ropt
	fi
fi

if grep -iqE "gpnorm gmvt0" $fin && echo $norms | grep -q gmv
then
	echo "GP norms for GMV"
	gpre=""
	gt1="gpnorm gmvt1"
	echo $norms | grep -q gpt1 && gpre="$gt1 sl:$gt1 slmf:$gt1 cpglag"
	R --slave -f $node/gpnorms.R --args $fin $fin2 lev=$lev type=gpgmv$suf \
		gpref="gpnorm gmvt0" gpre="$gpre" png=$dd $ropt
fi

if grep -iqE "gpnorm adiab call_sl" $fin && echo $norms | grep -q adiab
then
	echo "GP norms for ADIAB"
	R --slave -f $node/gpnorms.R --args $fin $fin2 lev=$lev type=gpadiab$suf \
		gpref="gpnorm adiab call_sl" gpre="gpnorm adiab cpg" png=$dd $ropt
fi

if grep -iqE "gpnorm zb2 cpg$" $fin && echo $norms | grep -q slb2
then
	echo "GP norms for ZB2"
	R --slave -f $node/gpnorms.R --args $fin $fin2 lev=$lev type=gpsi$suf \
		gpref="gpnorm zb2 cpg$" gpre="gpnorm zb2 (call_)?sl$:gpnorm zb2 slmf" png=$dd \
		$ropt
fi

if grep -iqE "FULL-POS GPNORMS" $fin && echo $norms | grep -q fp
then
	if grep -iqE "gpnorm dynfpos z" $fin
	then
		echo "GP norms for FullPOS"
		fpre="gpnorm dynfpos di:gpnorm dynfposlag:gpnorm endvpos z:gpnorm endvpos fp:allfpos"
		R --slave -f $node/fpnorms.R --args $fin $fin2 lev=$lev type=fp$suf \
			fpref="gpnorm dynfpos z" fpre="$fpre" png=$dd $ropt
	fi

	if [ -z "$fin2" ] && grep -iq "allfpos" $fin
	then
		echo "GP norms for FullPOS+SPEC+GFL"
		R --slave -f $node/fpnorms.R --args $fin lev=$lev type=fpgp$suf fpref="allfpos" \
			png=$dd $ropt
	fi
fi

if grep -iqE "FULL-POS SPNORMS" $fin && echo $norms | grep -q fp
then
	echo "SP norms for FullPOS"
	R --slave -f $node/fpspnorms.R --args $fin $fin2 lev=$lev type=fpsp \
		fpref="full-pos spnorms" fpre="xxxyyyzzz" png=$dd $ropt

	if [ -z "$fin2" ]
	then
		echo "SP norms for FullPOS+SPEC"
		R --slave -f $node/fpspnorms.R --args $fin lev=$lev type=fpspec png=$dd $ropt
	fi
fi

for pre in sp gpgmv$suf gpgfl$suf gpadiab$suf gpsi$suf fp$suf fpsp \
	gflspec$suf fpgp$suf fpspec
do
	echo "HTML files for $pre"
	{
	echo "<table>"
	for s in "" v
	do
		echo "<tr>"

		nc=0
		# split lists in 2 so that 1...9 and 10... remain ordered
		for ficp in $(ls -1 $dd | grep -E "^${pre}norm$s[[:digit:]]\.png") \
			$(ls -1 $dd | grep -E "^${pre}norm$s[[:digit:]]{2,}\.png")
		do
			[ $nc -gt 0 -a $((nc%2)) -eq 0 ] && echo -e "</tr>\n<tr>"
			echo -e "\t<td><img src=\"$dd/$ficp\"/></td>"
			nc=$((nc+1))
		done

		echo "</tr>"
	done
	echo "</table>"
	} > $temp/$pre.html

	if [ -s $dd/$pre.txt ]
	then
		echo "<pre>"
		cat $dd/$pre.txt
		echo "</pre>"
	fi >> $temp/$pre.html

	[ $(wc -l $temp/$pre.html | awk '{print $1}') -le 4 ] && echo > $temp/$pre.html
done

echo "Adding thematic HTML files from $temp"
date=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
res=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
base=$(printf "%s %dh" $date $((res/3600)))
sed -re "s:TAG NODE:$fin:" -e "s:TAG BASE:$base:" \
	-e "/TAG SP\>/r $temp/sp.html" -e "/TAG GPGMV/r $temp/gpgmv$suf.html" \
	-e "/TAG GPGFL/r $temp/gpgfl$suf.html" -e "/TAG GPADIAB/r $temp/gpadiab$suf.html" \
	-e "/TAG GPSI/r $temp/gpsi$suf.html" -e "/TAG SPFPOS/r $temp/fpsp.html" \
	-e "/TAG FPOS/r $temp/fp$suf.html" -e "/TAG GFLSPEC/r $temp/gflspec$suf.html" \
	-e "/TAG SPECFPOS/r $temp/fpspec.html" -e "/TAG GPFPOS/r $temp/fpgp$suf.html" \
	$node/norms.html > $fout
