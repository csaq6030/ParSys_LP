#!/bin/bash

# Execute job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# The batch system should use the current directory as working directory.
#$ -cwd

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N name

# Redirect output stream to this file.
#$ -o outNP32.dat

# Join the error stream to the output stream.
#$ -j yes

# Send status information to this email address. 
#$ -M Chhong.Lee@uibk.ac.at

# Send me an e-mail when the job has finished. 
#$ -m e

# Specify the amount of virtual memory given to each MPI process
# in the job.
#$ -l h_vmem=1G

# Use the parallel environment "openmpi-fillup", which assigns as many processes
# as available on each host. Start 16 MPI processes across an arbitrary number of
# hosts. For each process, SGE will reserve one CPU-core.
#$ -pe openmpi-fillup 2

## ALTERNATIVE
# Use the parallel environment "openmpi-fillup", which assigns as many processes
# as available on each host. If there are not enough machines to run the MPI job
# on the maximum of 16 cores, the batch system will gradually try to run the job
# on fewer cores, but not less than 8. 
##  #$ -pe openmpi-fillup 8-16

module load gcc/5.1.0
module load openmpi/2.1.0
./jsHelper.sh
