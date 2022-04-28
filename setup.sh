#!/bin/sh

diag=~/util/diags

set -e

usage()
{
	printf "
Description:
	Print HTML output related to a few plots from setup information of any model forecast

Usage:
	setup.sh FILE [-o HTML] [-h]

Options:
	FILE: NODE file
	HTML: output file containing HTML content (text and links to images)
	-h: print this help and exit normally

Details:
	An HTML output is produced, containing text and images built from the NODE file. \
Plots are placed into a directory (created if needed) aside to the NODE file. \
If the name of this NODE file begins with 'node', the directory is named after \
'nodexxx', prefix 'node' being removed. Otherwise, the directory is created by mktemp, \
with prefix 'setup'. Be aware that in this latter case, the directory is not removed.
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

if [ -z "$fin" ]
then
	echo "Error: input option missing" >&2
	exit 1
fi

file $fin | grep -q text
grep -qi 'END OF SETUPS' $fin

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

if [ -z "$fout" ]
then
	R --slave -f $diag/procmap.R --args $fin
else
	loc=$(dirname $fin)
	xpdir=$(cd $loc > /dev/null && pwd)
	fic=$xpdir/$(basename $fin)

	if echo $fic | grep -qEi '\<node\w'
	then
		temp=$(echo $fic | sed -re 's:\<node::i')
		mkdir -p $temp
	else
		temp=$(mktemp -d -p $loc setupXXX)
	fi

	echo "--> output sent to $temp"
	cd $temp > /dev/null

	R --slave -f $diag/procmap.R --args $fic png > out.txt

	{
		echo "<pre>"
		grep -A 1 -iw values out.txt
		grep -w SL out.txt
		echo "</pre>"
	} > map.txt

	date=$(grep -E 'NUDATE *=' $fic | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
	res=$(grep -E 'NUDATE *=' $fic | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
	base=$(printf "%s %dh" $date $((res/3600)))
	sed -re "s:TAG NODE:$1:" -e "s:TAG BASE:$base:" -e "s:TAG DIR:$temp:g" \
		-e '/TAG MAP/r map.txt' $diag/setup.html > out.html

	cd $OLDPWD
	mv $temp/out.html $fout
fi
