#!/bin/bash
## RENAME FOR YOUR JOB
#PBS -N NC1-SO4-OPT

## EDIT FOR YOUR JOB
## Request 8 CPUs (cores) on 4 nodes, 16 total cores
#PBS -l nodes=1:ppn=8,mem=24gb,feature=8core

## WALLTIME DEFAULTS TO ONE HOUR - ALWAYS SPECIFY FOR LONGER JOBS
## If the job doesn't finish in 10 minutes, cancel it
#PBS -l walltime=80:30:00

## EDIT FOR YOUR JOB
## Put the output from jobs into the below directory
#PBS -o /gscratch/esci/qshao/gaussian/zw/NC1-SO4-OPT
## Put both the stderr and stdout into a single file
#PBS -j oe

## EDIT FOR YOUR JOB
## Sepcify the working directory for this job
#PBS -d /gscratch/esci/qshao/gaussian/zw/NC1-SO4-OPT

## Some applications, particularly FORTRAN applications require
##  a larger than usual data stack size. Uncomment if your
##  application is exiting unexpectedly.
#ulimit -s unlimited

 
### Specify the app to run here                           ###
###                                                       ###
# EDIT FOR YOUR JOB
# ALWAYS use mpiexec.hydra NOT mpiexec
# ALWAYS include the "-rmk pbs" directive
# run YOURJOB, it's in the working directory specified above, /gscratch/GROUPNAME/USERNAME/JOB_DIR
# this script is to generate the necessary tpr file and make the mdrun file
# variables
# the gromacs binary file
export GAUSS_EXEDIR="/gscratch/esci/qshao/g09"
FILE=*.com
for file in $FILE
do
echo $file
file1=${file%.*}
echo $file1
mpiexec.hydra -n 8 g09 <$file> $file1.out 
done
### include any post processing here                      ###
###                                                       ###
 
