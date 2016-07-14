#!/bin/bash
#SBATCH -J CRVE		 		# Job name
#SBATCH -o CRVE.o%j 		# Name of stdout output file (%j expands to jobId)
#SBATCH -e CRVE.o%j 		# Name of stderr output file(%j expands to jobId)
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -N 1				# Total number of nodes
#SBATCH -n 16 		   		# Total number of mpi tasks requested
#SBATCH -t 0:15:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < panel_simulation_for_TACC.R
