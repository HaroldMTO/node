#!/bin/sh

diags=~/util/diags

set -e

usage()
{
	echo "
Description:
	Print HTML output related to a few plots from setup information of any model forecast

Usage:
	setup.sh FNODE [FHTML] [-o PNG] [-h]

Options:
	FNODE: input NODE file (as produced by MASTERODB)
	FHTML: HTML output file containing text and links to images related to FNODE \
(default: setup.html)
	PNG: local directory name (no '/') where output PNG files (+ text files) are produced
	-h: print this help and exit normally

Details:
	Input file is a 'NODE file' as produced by MASTERODB.
	An HTML output file is produced, containing text and images built from the NODE file. \
Its default name is 'setup.html'.
	Plots are placed into directory [PNG] (created if needed) aside (ie local) to the \
HTML file. If the name of this NODE file begins with 'node', the directory is named \
after 'nodexxx', prefix 'node' being removed. Otherwise, the directory is created by \
mktemp, with prefix 'setup'. Be aware that in this latter case, the directory is not \
removed."
}

if [ $# -eq 0 ] || echo " $*" | grep -qE '\-\w*h'
then
	usage
	exit
fi

fin=""
fout="setup.html"

while [ $# -ne 0 ]
do
	case $1 in
	-o)
		png=$2
		shift
		;;
	*)
		[ -z "$fin" ] && fin=$1 || fout=$1
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

echo "NODE file $fin"
ls -L $fin > /dev/null

if ! file -L $fin | grep -q text
then
	ftmp=$(mktemp --tmpdir)
	cat $fin | tr -d '\0' > $ftmp

	if ! file -L $ftmp | grep -q text
	then
		echo "$fin is not a text file" >&2
		exit 1
	fi

	fin=$ftmp
fi

if ! grep -qi 'END OF SETUPS' $fin
then
	echo "Warning: no text \"END OF SETUPS\" in $fin" >&2
fi

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

if [ -z "$png" ]
then
	if echo $fin | grep -qEi '(.+/)?\<node\.?\w+'
	then
		png=$(echo $fin | sed -re 's:(.+/)?\<node\.?(\w+):\2:i')
	else
		loc=$(dirname $fin)
		png=$(mktemp -d setupXXX)
	fi
fi

mkdir -p $png
echo "Parse NODE file (PNG graphics in $png)"

R --slave -f $diags/procmap.R --args ficin=$fin png=$png

echo "Write HTML file $fout"
{
	echo "<pre>"
	grep -A 1 -iw values $png/out.txt || true
	grep -w SL $png/out.txt || true
	echo "</pre>"
} > $png/map.txt

[ -s $png/jo.txt ] || echo "no Jo tables" > $png/jo.txt

date=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
res=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
base=$(printf "%s %dh" $date $((res/3600)))
sed -re "s:TAG NODE:$fin:" -e "s:TAG BASE:$base:" -e "s:TAG DIR:$png:g" \
	-e "/TAG MAP/r $png/map.txt" -e "/TAG JO/r $png/jo.txt" $diags/setup.html > $fout
