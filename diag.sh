#!/bin/sh

diags=~/util/diags

usage()
{
	printf "
Description:
	Produce a list of diagnostics on model forecasts

Usage:
	diag.sh (-conf CONFIG|-file FILE -param PAR -level LEV) [-o HTML] [-ref PATH] [-h]

Arguments:
	PATH: path to where to find reference files
	-h: print this help and exits normally

Details:
	PATH is the path to a hierarchy of files and directories following \
'HH/YYYYMMDD/file_ref'.

Dependencies:
	R software
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

conf="config"
fic=""
lev=0
par=""
html=""

while [ $# -ne 0 ]
do
	case $1 in
		-conf)
			conf=$2
			shift
			;;
		-level)
			lev=$2
			shift
			;;
		-param)
			par=$2
			shift
			;;
		-file)
			fic=$2
			shift
			;;
		-ref)
			ref=$2
			shift
			;;
		-o)
			html=$2
			shift
			;;
		-h)
			usage
			exit
			;;
		*)
			echo "Warning: unknown option '$1', ignored" >&2
			;;
	esac

	shift
done

if [ -n "$conf" ]
then
	fic=$conf/file.txt
	lev=$conf/level.txt
	par=$conf/param.txt
elif [ -z "$fic" -o -z "$lev" -o -z "$par" ]
then
	echo "Error: input options missing
conf: '$conf'
lev: '$lev'
par: '$par'
fic: '$fic'" >&2
	exit 1
fi

set -e

ls $par > /dev/null
if [ -n "$ref" ]
then
	ls -d $ref > /dev/null
	opt="ref=$ref"
fi

if [ -n "$html" ]
then
	loc=$html
else
	loc=$(mktemp -d diagXXX)
fi

echo "--> sending output to $loc"
mkdir -p $loc

if [ ! -s "$lev" ] && ! echo $lev | grep -qE '[0-9]+(:[0-9]+)*'
then
	echo "Error: option '-level' uncorrectly defined" >&2
	exit 1
	#echo $lev | tr ':' '\n' > $loc/levels.txt
fi

type R > /dev/null 2>&1 || module -s load intel R > /dev/null 2>&1

R --slave -f $diags/diag.R --args fic=$fic params=$par level=$lev png=$loc $opt

if [ -n "$html" ]
then
	for type in map hist err score spec prof
	do
	{
		echo "<tr>"

		n=0
		for ficp in $(ls -1 $loc | grep -E "$type.+.png")
		do
			if [ $n -eq 2 ]
			then
				echo -e "</tr>\n<tr>\n"
				n=0
			fi

			n=$((n+1))
			printf "\t<td><img src=\"%s\"/></td>\n" $loc/$ficp
		done

		echo "</tr>"
	} > $loc/$type.html
	done

	sed -re "/TAG MAP/r $loc/map.html" -e "/TAG HIST/r $loc/hist.html" \
		-e "/TAG SPEC/r $loc/spec.html" -e "/TAG PROF/r $loc/prof.html" \
		-e "/TAG ERR/r $loc/err.html" -e "/TAG SCORE/r $loc/score.html" $diags/diag.html > $html/diag.html
fi
