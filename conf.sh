#!/bin/sh

diags=~/util/diags

usage()
{
   printf "
Description:
	Print logical keys appearing in NODE file

Usage:
	conf.sh NODEFILE [-h]

Arguments:
	NODEFILE: path to (classical) file 'NODE.001_01'
	-h: print this message and exit normally

Details:
	/

Dependencies:
	None
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

fnode=NODE.001_01
opts="$*"

while [ $# -ne 0 ]
do
	case $1 in
		-node)
			fnode=$2
			shift
			;;
		-h)
			opts=""
			;;
		*)
			fnode=$1
			;;
	esac

	shift
done

if [ -z "$opts" ]
then
	usage
	exit
elif [ -z "$fnode" ]
then
	echo "Error: arguments missing
NODEFILE: '$fnode'
" >&2
	exit 1
fi

set -e

ls $fnode >/dev/null

# works pretty well..
grep -E ' *[[:alnum:]]+ *= *(T(RUE)?|F(ALSE)?)\>' $fnode | \
	sed -re 's:(\<([[:alnum:]]+%)*[[:alnum:]]+) *= *([TF])(RUE|ALSE)?\>:\n\1 = \3\n:g' | \
	grep -E '[[:alnum:]]+ = [TF]\>' | \
	sed -re 's:.*(\<([[:alnum:]]+%)*[[:alnum:]]+ = [TF])\>:\1:' | sort -u
