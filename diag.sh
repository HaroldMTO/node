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
	HTML: path to output files (graphics and HTML files)
	-h: print this help and exits normally

Details:
	PATH is the path to a hierarchy of files and directories following \
'HH/YYYYMMDD/file_ref'.

Dependencies:
	R software
"
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
then
	usage
	exit
fi

conf="config"
fic=""
lev=""
par=""
graph=1
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
		-nograph)
			graph=0
			;;
		*)
			echo "Warning: unknown option '$1', ignored" >&2
			;;
	esac

	shift
done

if [ -n "$conf" ]
then
	[ -z "$fic" ] && fic=$conf/file.txt
	[ -z "$lev" -a -s $conf/level.txt ] && lev=$conf/level.txt || lev=0
	[ -z "$par" ] && par=$conf/param.txt
fi

if [ -z "$fic" -o -z "$lev" -o -z "$par" ]
then
	echo "Error: input options missing
conf: '$conf'
lev: '$lev'
par: '$par'
fic: '$fic'" >&2
	exit 1
fi

if [ ! -s $lev ] && ! echo $lev | grep -qE '[0-9]+(:[0-9]+)*'
then
	echo "Error: option '-level' uncorrectly defined" >&2
	exit 1
	#echo $lev | tr ':' '\n' > $loc/levels.txt
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

if [ $graph -eq 1 ]
then
	type R > /dev/null 2>&1 || module -s load intel R > /dev/null 2>&1
	R --slave -f $diags/diag.R --args fic=$fic params=$par level=$lev png=$loc $opt > \
		$loc/diag.log
fi

[ -n "$html" ] || exit

ficdom=$conf/domain.txt
if [ ! -s $ficdom ]
then
	echo "Error: no domains file '$ficdom'" >&2
	exit 1
fi

doms=$(grep -E '^ *\w+' $ficdom | awk 'NR > 1 {print $1}' | xargs)
params=$(ls -1 $loc | grep -E 'map[1-9].+_\w+\.png$' | sed -re 's:.+_(\w+)\.png:\1:' | \
	sort -u | xargs)
echo "HTML files: $(cat $loc/steps.txt | wc -l) forecast steps
domains: $doms"

for par in $params
do
	ficpar=$(dirname $html)/$par.html
	echo ". creating $ficpar"

	{
	it=0
	while read -a tt
	do
		it=$((it+1))
		echo ${tt[*]} | grep -qE 'graph: TRUE' || continue

		echo "<table>"
		echo "<tr><th>Maps, selection of levels</th><th>Cross-sections and diagrams</th>"
		att="colspan='2'"
		if [ -n "$ref" ]
		then
			echo "<th>Maps for reference</th><th>Cross-sections and diagrams for ref</th>"
			att="colspan='4'"
		fi

		echo "</tr>"

		for dom in $doms
		do
			echo "<tr><th name='title' $att>Param '$par' - ${tt[*]} - Domain $dom</th></tr>"
			echo "<tr>"

			for typ in map hist mapdiff histdiff
			do
				fic=$loc/$typ$it${dom}_$par.png
				[ -s $fic ] || continue

				printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
			done

			echo "</tr>"

			cat $diags/step.html

			for stat in bias rmse errx
			do
				echo "<tr><th colspan='2'>Param '$par' $stat - Domain $dom</th></tr>"
				echo "<tr>"

				for typ in map$stat hist$stat
				do
					fic=$loc/$typ$it${dom}_$par.png
					[ -s $fic ] || continue

					printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
				done

				echo "</tr>"
			done
		done > $loc/idom.html

		if grep -qE '<img .+ src=' $loc/idom.html
		then
			cat $loc/idom.html
			echo "</table>"
			break
		fi

		echo "</table>"
	done < $loc/steps.txt

	sed -re 's:TAG NAME:step:' -e "/TAG OPT/r $loc/steps.html" $diags/select.html
	for dom in $doms
	do
		for stat in "" diff bias rmse errx
		do
			for typ in map$stat hist$stat
			do
				fic=$loc/$typ${dom}_$par.html
				[ -s $fic ] || continue

				sed -re 's:TAG NAME:map:' -e "/TAG OPT/r $fic" $diags/select.html
			done
		done
	done

	typd=("stat" "statv" "err" "errv" "score" "scorev" "scoret" "rmsevt")
	titd=("Statistics of forecast" "Statistics of profile" "Statistics of forecast error" \
		"Statistics of profile error" "Scores of forecasts" "Scores of profile" \
		"Daily scores" "Daily RMSE")
	i=0
	while [ $i -lt ${#typd[*]} ]
	do
		echo "<h2>${titd[i]} on domains</h2>
<table>
<tr>"

		for ficp in $(ls -1 $loc | grep -E "${typd[i]}[0-9]+_$par.png" | sort)
		do
			printf "\t<td><img src=\"%s\" alt=\"missing image\"/></td>\n" $loc/$ficp
		done

		echo "</tr>
</table>"
		i=$((i+1))
	done
	} > $loc/$par.html

	sed -re "s:TAG PAR:$par:" -e "/TAG BODY/r $loc/$par.html" $diags/par.html > $ficpar
done
