#!/bin/bash

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N monteCarloPi

# Execute job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# The batch system should use the current directory as working directory.
#$ -cwd

# Redirect output stream to this file.
#$ -o output.txt

# Join the error stream to the output stream.
#$ -j yes

# Send status information to this email address. 
#$ -M Chhong.Lee@student.uibk.ac.at

# Send me an e-mail when the job has finished. 
#$ -m e

# Use the parallel environment "openmp" with 8 job slots. Each requested
# job slot will be assigned one core on the execution host.
#$ -pe openmpi-fillup 16

# Allocate 2 Gigabytes per job slot.
# The total memory available to your program
# (i.e. the UNIX "ulimit -v" value) will be the
# product of job slots from the -pe directive
# times the h_vmem requirement. For the present
# example, the job will get 16GB of virtual memory.
#$ -l h_vmem=2G

# tell OpenMP how many threads to start

export OMP_NUM_THREADS=8
module load gcc/5.1.0
module load openmpi/2.1.1

./job_script_helper.sh
