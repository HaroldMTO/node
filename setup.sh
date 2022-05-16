#!/bin/sh

diags=~/util/diags

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

ls -L $fin > /dev/null
file -L $fin | grep -q text
grep -qi 'END OF SETUPS' $fin

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

if [ -z "$fout" ]
then
	echo "--> output sent to Rplots.pdf or PNG files"
	R --slave -f $diags/procmap.R --args $fin
else
	if echo $fin | grep -qEi '(.+/)?\<node\.?\w+'
	then
		dd=$(echo $fin | sed -re 's:(.+/)?\<node\.?(\w+):\1\2:i')
		mkdir -p $dd
	else
		loc=$(dirname $fin)
		dd=$(mktemp -d -p $loc setupXXX)
	fi

	cd $dd > /dev/null
	echo "--> output sent to $dd"

	fic=$(basename $fin)
	R --slave -f $diags/procmap.R --args ../$fic png > out.txt

	{
		echo "<pre>"
		grep -A 1 -iw values out.txt
		grep -w SL out.txt
		echo "</pre>"
	} > map.txt

	cd $OLDPWD > /dev/null

	date=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUDATE *= *([0-9]+) .+:\1:')
	res=$(grep -E 'NUDATE *=' $fin | sed -re 's:.*\<NUSSSS *= *([0-9]+).*:\1:')
	base=$(printf "%s %dh" $date $((res/3600)))
	sed -re "s:TAG NODE:$fin:" -e "s:TAG BASE:$base:" -e "s:TAG DIR:$dd:g" \
		-e "/TAG MAP/r $dd/map.txt" $diags/setup.html > $fout
fi
