#!/bin/bash
##############################################
# An example jobarray sbatch file
##############################################

#SBATCH -J ABCDtask		# Job name

#SBATCH --partition=small	# Name of partition to assign jobs to

#SBATCH --array=1-24404%100	# Read in the first 24,404 lines in script_list.sh by saying array=1-24,404
				# Assign to nodes for them to run
				# %100 tells SLURM to only have 100 jobs running at any given time

#SBATCH --time=03:00:00		# Walltime/time limit for jobs = 3 hours

#SBATCH -n 1			# 1 CPU per job

#SBATCH --mem=6G		# 6 Gb of RAM per job

				# Use the job number as the prefix for .err/.out files

# Outputs ----------------------------------
#SBATCH -o %A-%a.out
#SBATCH -e %A-%a.err
# ------------------------------------------

				# Pull each line from script_list, create array, and assign a task number

eval $( sed "${SLURM_ARRAY_TASK_ID}q;d" script_list.sh )

