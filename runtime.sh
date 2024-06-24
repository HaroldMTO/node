#!/bin/sh

node=~/util/node

set -e

usage()
{
	echo "
Description:
	Produce an HTML file gathering information from the NODE file, 'run' part

Usage:
	runtime NODE [-o HTML] [-h]

Options:
	NODE: a NODE file from a MASTERODB job
	HTML: HTML output file containing text and links to images related to the NODE file \
(default: runtime.html)
	-h: print this help and exit normally

Details:
	NODE file is parsed, looking for information from the runtime part of the job.
	An HTML output file is produced, containing text and images built from the NODE file.
	Plots are produced in one directory named after the NODE file, its prefix \
'NODE'/'node' being deleted. This directory is created (if needed) aside from HTML."
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
then
	usage
	exit
fi

fin=""
fin2=""
fout=runtime.html
png=""

while [ $# -ne 0 ]
do
	case $1 in
	-o)
		fout=$2
		shift
		;;
	-*)
		echo "Error: unknown option '$1'" >&2
		exit 1
		;;
	*)
		if [ -z "$fin" ]
		then
			fin=$1
		else
			fin2="$fin2 $1"
		fi
		;;
	esac

	shift
done

if [ -z "$fin" -o -z "$fout" ]
then
	echo "Error: mandatory arguments missing
fin: '$fin'
fout: '$fout'" >&2
	exit 1
fi

set -e

for f in $fin $fin2
do
	ls -L $f > /dev/null
done

if ! file -L $fin | grep -q text
then
	ftmp=$(mktemp --tmpdir)
	cat $fin | tr -d '\0' > $ftmp

	if ! file -L $ftmp | grep -q text
	then
		echo "Error: $fin is not a text file" >&2
		exit 1
	fi

	fin=$ftmp
fi

loc=$(dirname $fout)

if echo $fin | grep -qEi '(.+/)?\<node\.?\w+'
then
	png=$loc/$(echo $fin | sed -re 's:(.+/)?\<node\.?(\w+):\2:i')
else
	png=$(mktemp -t $loc -d runtimeXXX)
fi

png=$(echo $png | sed -re 's:^\./::')
mkdir -p $png
echo "Parse NODE file (PNG graphics in $png/)"

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1
if ! env | grep -qw R_LIBS
then
	export R_LIBS=~petithommeh/lib
	echo "--> setting R_LIBS: $R_LIBS"
else
	R_LIBS=$R_LIBS:~petithommeh/lib
fi

R --slave -f $node/runtime.R --args $fin $fin2 png=$png

n=0
{
echo "<table><tr>"
for fic in $(ls $png | grep anaguess_.+.png)
do
	[ $((n%2)) -eq 1 ] && echo "</tr>\n<tr>"
	echo "<td><img src='$fic'/></td>"
	n=$((n+1))
done
echo "</tr></table>"
} > $png/anaguess.html

[ $n -eq 0 ] && echo > $png/anaguess.html

if [ -s $png/canari.txt ]
then
	echo "<br><pre>"
	cat $png/canari.txt
	echo "</pre>"
fi >> $png/anaguess.html

date=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
res=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
base=$(printf "%s %dh" $date $((res/3600)))
sed -re "s:TAG NODE:$fin:" -e "s:TAG BASE:$base:" -e "s:TAG DIR:$png:g" \
	-e "/TAG ANA/r $png/anaguess.html"  -e "/TAG JO/r $png/jo.txt" $node/runtime.html > \
	$fout
