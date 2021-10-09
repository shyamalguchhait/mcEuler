#!/bin/bash
# Parallel job submission script:
# Usage: qsub  <this_script>
#
# The shell used to run the job
#$ -S /bin/bash
#
# The name of the parallel queue to submit the job to
#
# Define the parallel runtime environment and number of nodes
# NB: number of nodes is one more than needed as one copy resideson the master node
#$ -pe mpi 32
# Use location that job was submitted as working directory
#$ -cwd
#
# Export all environment variables to the slave jobs
#$ -V
#
# Put stdout & stderr into the same file
#$ -j y
#
# The name of the SGE logfile for this job
#$ -o output_p.log
#mpiexec -n $NSLOTS ./heat
gcc *.c -lm -o ./a_p.out
./a_p.out
