#!/bin/sh

diags=~/util/diags

usage()
{
	printf "
Description:
	Produce a list of diagnostics on model forecasts

Usage:
	diag.sh PATH [-h]

Arguments:
	PATH: file or path to files ('.diff') of field differences
	-h: print this help and exits normally

Details:

Dependencies:
	R software
"
}

fic=""
help=0

[ $# -eq 0 ] && help=1

while [ $# -ne 0 ]
do
	case $1 in
		-h)
			help=1
			;;
		*)
			[ -z "$fic" ] && fic=$1
			;;
	esac

	shift
done

if [ $help -eq 1 ]
then
	usage
	exit
elif [ -z "$fic" ]
then
	echo "usage: diag.sh PATH [-h]" >&2
	exit 1
fi

set -e

if [ -d $fic ]
then
	path=$fic
	find -L $fic -maxdepth 1 -name ICMSHARPE+\* | grep -E '/ICMSHARPE\+[0-9]{4,}$' | \
		sort > list.txt
else
	path=$(dirname $fic)
	find -L -maxdepth 1 -name $fic\* | sort > list.txt
fi

date=$(grep '^ *NUDATE *=' $path/NODE.001_01 | \
	sed -re 's:^ *\w+ *= *([0-9]{8}) \w+ *= *([0-9]+):\1/\2:')

type R > /dev/null 2>&1 || module -s load intel R > /dev/null 2>&1

if true
then
	R --slave -f $diags/diag.R --args fic=list.txt params=SURFTEMPERATURE:SURFPRESSION
else
	R --slave -f $diags/diff.R --args fic=list.txt params=SURFTEMPERATURE:SURFPRESSION
		date=$date path=/scratch/work/petithommeh/oper/arpege/4dvarfr/OPER
fi
