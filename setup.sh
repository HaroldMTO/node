#!/bin/sh

diag=~/util/diag

set -e

usage()
{
	printf "
Description:
	Print HTML output related to a few plots from setup information of any model forecast

Usage:
	setup.sh FILE [-h]

Options:
	FILE: NODE file
	-h: print this help and exit normally

Details:
	An HTML output is produced, containing text and images built from the NODE file. \
Plots are placed into a directory (created if needed) aside to the NODE file. \
If the name of this NODE file begins with 'node', the directory is named after \
'nodexxx', prefix 'node' being removed. Otherwise, the directory is created by mktemp, \
with prefix 'setup'. Be aware that in this latter case, the directory is not removed.
"
}

if [ $# -eq 0 ] || echo $* | grep -qE '(^| )\-h\>'
then
	usage
	exit
fi

file $1 | grep -q text
grep -qi 'END OF SETUPS' $1

loc=$(dirname $1)
xpdir=$(cd $loc > /dev/null && pwd)
fic=$xpdir/$(basename $1)

if echo $fic | grep -qE '\<node\w'
then
	temp=$(echo $1 | sed -re 's:\<node::')
else
	temp=$(mktemp -d -p $loc setupXXX)
fi

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
	-e '/TAG MAP/r map.txt' $diag/setup.html
