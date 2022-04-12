#!/bin/sh

diags=~/util/diags

usage()
{
   printf "
Description:
	Produce vertical profiles of T/Qv in case of error found by checkmv

Usage:
	error.sh ERRFILE [-nlev NLEV|-node NODEFILE]

Arguments:
	ERRFILE: file where checkmv output can be found
	NODEFILE: path to (classical) file 'NODE.001_01'
	NLEV: number of vertical levels in checkmv output

Details:
	Number of vertical levels (aka NFLEVG) is provided either by way of NLEV or \
by way of NODEFILE, where it would be found in this very script.

Dependencies:
	R software
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

fic=""
fnode=NODE.001_01
nl=0

while [ $# -ne 0 ]
do
	case $1 in
		-nlev)
			nl=$2
			shift
			;;
		-node)
			fnode=$2
			shift
			;;
		*)
			fic=$1
			;;
	esac

	shift
done

if [ -z "$fic" ] || [ -z "$fnode" -a $nl -eq 0 ]
then
	echo "Error: arguments missing
ERRFILE: '$fic'
NODEFILE: '$fnode'
NLEV: '$nl'
" >&2
	exit 1
fi

set -e

ls $fic >/dev/null
[ $nl -eq 0 ] && ls $fnode >/dev/null

if ! grep -qE ' out of physical bounds|wind too strong at|velocity out at' $fic
then
	echo "--> no error found by checkmv/tq/gmv"
	exit
fi

nl=$(grep -E 'NFLEVG *=' $fnode | sed -re 's: *NFLEVG *= *([0-9]+).*:\1:')

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

R --slave -f $diags/error.R --args $fic $nl

