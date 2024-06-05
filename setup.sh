#!/bin/sh

node=~/util/node

set -e

usage()
{
	echo "
Description:
	Create an HTML output and a few PNG plots from setup information of any log from \
MASTERODB

Usage:
	setup.sh FNODE [-o HTML] [-h]

Options:
	FNODE: input NODE file (ie log file from MASTERODB)
	HTML: HTML output file containing text and links to images related to FNODE \
(default: setup.html)
	-h: print this help and exit normally

Details:
	The file name FNODE must begin with 'NODE' (or 'node') after any path part.
	An HTML output file is produced, containing text and images built from the NODE file.
	Plots are produced in one directory named after the NODE file, its prefix \
'NODE'/'node' being deleted. This directory is created (if needed) aside from HTML. \
If the name of this NODE file begins with 'node', the directory is named \
after 'nodexxx', prefix 'node' being removed. Otherwise, the directory is created by \
mktemp, with prefix 'setup'."
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
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
			echo "Error: input file already set as '$fin', unknown option '$1'" >&2
			exit 1
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

echo "NODE file $fin"
ls -L $fin > /dev/null

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

if ! grep -qiE '\-+ Set up' $fin
then
	echo "Error: $fin is not a NODE file (no 'Set up...' content)" >&2
	exit 1
fi

if ! grep -qi 'END OF SETUPS' $fin
then
	echo "Warning: no text \"END OF SETUPS\" in $fin" >&2
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

if echo $fin | grep -qEi '(.+/)?\<node\.?\w+'
then
	png=$loc/$(echo $fin | sed -re 's:(.+/)?\<node\.?(\w+):\2:i')
else
	png=$(mktemp -t $loc -d setupXXX)
fi

png=$(echo $png | sed -re 's:^\./::')
mkdir -p $png
echo "Parse NODE file (PNG graphics in $png/)"

R --slave -f $node/procmap.R --args ficin=$fin png=$png

echo "Write HTML file $fout"
if [ -s $png/out.txt ]
then
	echo "<pre>"
	grep -A 1 -iw values $png/out.txt || true
	grep -w SL $png/out.txt || true
	echo "</pre>"
fi > $png/map.txt

[ -s $png/jo.txt ] || echo "no Jo tables" > $png/jo.txt

date=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
res=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
base=$(printf "%s %dh" $date $((res/3600)))
sed -re "s:TAG NODE:$fin:" -e "s:TAG BASE:$base:" -e "s:TAG DIR:$png:g" \
	-e "/TAG MAP/r $png/map.txt" -e "/TAG JO/r $png/jo.txt" $node/setup.html > $fout
