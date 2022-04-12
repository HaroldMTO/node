#!/bin/sh

diags=~/util/diags

usage()
{
   printf "
Description:
   Produce plots of spectral norms

Usage:
   (1) spnorm NODEFILE1 [NODEFILE2] [NLEV] [NMAX]
   (2) spnorm CMDFILE [NLEV] [NMAX]

Arguments:
   NODEFILE[12]: 1 or 2 paths to (classical) files 'NODE.001_01'
   CMDFILE: path to a command file listing node files, colors and legends
	NLEV: index (1...n) of vertical level for norms plot
	NMAX: maximum number or time-steps to plot norms for

Details:
	Plots show spectral norms along with time-steps. Plots are gathered on 1 or 2 \
pages, depending on the number of spectral fields (NH model has more than H). \
Each plot shows one curve for each node file.
	Node files are passed as files (option 1) or in a command file (option 2). \
In option 2, the command file is a tabular text file with txt as extension. It \
should contain file paths, colors and legend names (for plot legends). Values \
are those of R language.
	When option nlev is not used, the norms averages are plotted.  When option \
nlev is used, plots show norms for this level only. In this case, printing of \
levels must have been activated in the model run. As it is a 2D spectral \
field, SP norms are always shown, whatever the level (if present).
	Option nmax lets the user limit the number of time-steps, thus limiting \
the number of points in plots.

Dependencies:
   R software
"
}

if [ $# -eq 0 ]
then
	usage
	exit
fi

set -e

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

nf=$(echo $args | tr ' ' '\n' | grep -v '=' | wc -l | awk '{print $1}')
if [ $nf -eq 1 ] && ! echo $args | grep -qE '\<(time|nmax)='
	R --slave -f $diags/spnorm1.R --args $*
else
	R --slave -f $diags/spnorm.R --args $*
fi
