#!/bin/bash

## Declare a name for this job (replace jobname with something more descriptive)
#PBS -N AlSm
## Request the queue for this job (replace queue with morgan1, izabela1, morgan2, izabela2, morgan3, or izabela3)
#PBS -q morgan2
## Request computational resources for this job as follows
##  nodes - specifies how many nodes to request
##  ppn   - specifies how many processors per node to request. 
##          Set it as follows:
##          Replace <num> below with 8 for morgan1 and izabela1
##          Replace <num> below with 12 for morgan2, izabela2 and izabela3
##          Replace <num> below with 32 for morgan3
## Request up to 2GB virtual memory per job, or 1GB for morgan3
#PBS -l nodes=12:ppn=12,pvmem=2GB
## Request walltime, with limits as follows:
##      72:00:00 for morgan1
##      96:00:00 for morgan2
##     168:00:00 for morgan3
#PBS -l walltime=96:00:00
## Combine PBS standard output and error files
##PBS -j oe
##PBS -k eo
## These are PBS standard output and error files.  Uncomment only if you don't want the defaults.
##PBS -o output.$PBS_JOBID
##PBS -e error.$PBS_JOBID

## How many procs do I have?
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"

## cd into the directory where I typed qsub
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE

## If you are not using the general purpose mpiexec, make sure your mpi 
## environment is properly set up such that the correct mpirun is found 
## (you should use the mpirun provided with the compiler
## used to compile the program you are running).

/opt/mpiexec/bin/mpiexec lmp_bardeeneth < system.in > out.txt
