#!/bin/bash
#SBATCH -J CRVE		 		# Job name
#SBATCH -o CRT-IV.o%j 		# Name of stdout output file (%j expands to jobId)
#SBATCH -e CRT-IV.o%j 		# Name of stderr output file(%j expands to jobId)
#SBATCH -p development	 	# Submit to the 'normal' or 'development' queue
#SBATCH -N 3				# Total number of nodes
#SBATCH -n 180 		   		# Total number of mpi tasks requested
#SBATCH -t 1:00:00 			# Run time (hh:mm:ss)
#SBATCH --mail-user=jepusto@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < "./iv simulation - cluster-level non-compliance.R"
