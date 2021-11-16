#!/bin/bash
#
#$ -cwd
##$ -j y
#$ -S /bin/bash
#$ -N WF7.85
#$ -pe orte 16
#$ -o WF.out
#$ -e WF.err
#
 
#
# Use modules to setup the runtime environment
#
. /etc/profile

ulimit -a
 
#
# Execute the run
#
mpirun -np $NSLOTS ./WF_measure -l conflist -i input_file -o Sp4WFB7.85.out
