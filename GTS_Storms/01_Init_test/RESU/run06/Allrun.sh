#!/bin/bash
source /home/finzeo/.bashrc
cs800
## $1: case name
if [ -z "$1" -o ! \( -z "$2" \) ]; then
	echo "Ingrese solo nombre de caso. Vuelva a ejecutar.\n"
	exit 1
fi
code_saturne8 run --id=$1 --param=setup.xml --initialize
cd ../RESU/$1
sbatch runjob.sh
