#!/bin/sh

diags=~/util/diags

usage()
{
   printf "
Description:
   Produce a map of MPI tasks and their grid-points from a NODE file

Usage:
   procmap.sh NODEFILE

Arguments:
   NODEFILE: path to a file 'NODE.001_01'
   -h: print this help message and exit normally

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

R --slave -f $diags/procmap.R --args $*
