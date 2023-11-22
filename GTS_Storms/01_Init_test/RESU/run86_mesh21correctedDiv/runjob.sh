#!/bin/bash
#SBATCH --job-name=SS8621 # nombre para identificar el trabajo. Por defecto es el nombre del script # Autor-SimMall-NumCaso-NumMalla
#SBATCH --ntasks=96 # cantidad de cores pedidos
##SBATCH --nodes=1 # cantidad de nodos pedidos
#SBATCH --ntasks-per-node=32 # cantidad de cores por nodo, para que agrupe o distribuya procesos
#SBATCH --cpus-per-task=1
# la linea siguiente es ignorada por Slurm porque empieza con ##
##SBATCH --mem-per-cpu=4G # cantidad de memoria por core
#SBATCH --output=trabajo-%j-salida.txt # la salida y error estandar van a este archivo. Si no es especifca es slurm-%j.out (donde %j es el Job ID)
#SBATCH --error=trabajo-%j-error.txt # si se especifica, la salida de error va por separado a este archivo
#SBATCH --time=13-23 # tiempo máximo de ejecución, el formato es: dias-horas / dias-horas:minutos / horas:minutos:segundos
##SBATCH --nodelist=n-8,n-9,n-10,n-11
#SBATCH --partition=secondary

export OMP_NUM_THREADS=1
export LD_PRELOAD=/share/apps/easybuild/software/GCCcore/10.3.0/lib64/libgfortran.so
mpiexec ./cs_solver --mpi $@
