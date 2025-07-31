#!/bin/sh

node=~/util/node

set -e

usage()
{
	echo "
Description:
	Produce an HTML file showing time-series of some measure from NODE files

Usage:
	runtimes.sh PATH PATTERN -o PNG [HTML] [-start START] [-end END] [-res RES] \
[-save FILE] [-h]

Options:
	PATH: path to a root directory containing NODE files in its subdirectories
	PATTERN: pattern (glob) for accessing the NODE files from a MASTERODB job
	HTML: HTML output file containing text and links to images related to the NODE file \
(default: runtimes.html)
	-start: limit base date-times to those after START
	-end: limit base date-times to those before END
	-res: limit base times to RES, a colon (ie ':') separated list of times
	-save: save an R image of runtime data (SP/GP norms only)
	-h: print this help and exit normally

Details:
	NODE files are parsed, looking for the measures from the runtime part of the job.
	An HTML output file is produced, containing text and images built from the NODE file.
	Plots are produced in the directory PNG.
	Arguments START and END should be in POSIX time format, eg 'YYYY-MM-DD [HH[:MM:SS]]'. \
If a time part is passed, the whole sequence must be quoted."
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
then
	usage
	exit
fi

path=""
patt=""
png=""
fout=""
start=""
end=""
res=""
save=""

while [ $# -ne 0 ]
do
	case $1 in
	-o)
		png=$2
		shift
		;;
	-start)
		start=$2
		shift
		;;
	-end)
		end=$2
		shift
		;;
	-res)
		res=$2
		shift
		;;
	-save)
		save=$2
		shift
		;;
	-*)
		echo "Error: unknown option '$1'" >&2
		exit 1
		;;
	*)
		if [ -z "$path" ]
		then
			path=$1
		elif [ -z "$patt" ]
		then
			patt=$1
		elif [ -z "$fout" ]
		then
			fout=$1
		else
			echo "Error: output file already set as '$fout', unknown option '$1'" >&2
			exit 1
		fi
		;;
	esac

	shift
done

[ -n "$fout" ] || fout=runtimes.sh

if [ -z "$path" -o -z "$patt" -o -z "$png" -o -z "$fout" ]
then
	echo "Error: mandatory arguments missing
path: '$path'
patt: '$patt'
png: '$png'
fout: '$fout'" >&2
	exit 1
fi

set -e

png=$(echo $png | sed -re 's:^\./::')
mkdir -p $png
echo "Parse NODE files (PNG graphics in $png/)"

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1
if ! env | grep -qw R_LIBS
then
	export R_LIBS=~petithommeh/lib
	echo "--> setting R_LIBS: $R_LIBS"
else
	R_LIBS=$R_LIBS:~petithommeh/lib
fi

ropt=""
[ -n "$start" ] && ropt="start=$start"
[ -n "$end" ] && ropt="$ropt end=$end"
[ -n "$res" ] && ropt="$ropt res=$res"
[ -n "$save" ] && ropt="$ropt save=$save"

R --slave -f $node/runtimes.R --args path=$path patt=$patt png=$png $ropt

nc=0
{
	echo "<tr>"
	for fic in $png/*.png
	do
		[ $nc -gt 0 -a $((nc%2)) -eq 0 ] && echo -e "</tr>\n<tr>"
		echo "<td><img src='$fic'/></td>"
		nc=$((nc+1))
	done
	echo "</tr>"
} > $png/img.html

path_=$(echo $path | sed -e 's/:/_/g')
sed -re "s:TAG PATH:$path_:" -e "s:TAG PATT:$patt:" -e "s:TAG DIR:$png:g" \
	-e "/TAG IMG/r $png/img.html" $node/runtimes.html > $fout
