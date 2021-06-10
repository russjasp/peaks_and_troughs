
rm(list=ls(all=TRUE))
options("scipen"=100, "digits"=6) 

setwd()

temp <- read.table ("init_template.txt", sep = "&")

basename <- "example"


# enter variables here #
m <- c (0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5) ## this migration is island migration ## stepping stone dispersal matrix is at bottom if desired instead ##
N <- c (1000) ## population size per patch
mut <- c (0.00001) ## mutation rate of quanti loci
Vs <- c (2,5,10) ## variance of fitness function, ie "strength of selection"


for (i in 1:length(m)){
  
  migration_rate <- m[i]
  
  ### ### ### Comment in for dispersal matrix, change "dispersal_rate" to "dispersal_matrix" lines 58-62
  # init_migration_rate <- dispersal_matrix_list[i] ### Make dispersal matrices at the bottom of script ###
  
  for (j in 1:length(N)){
    
    N_value <- N[j]
    
    for (k in 1:length(mut)){
      
      mutation_rate <- mut[k]
      
      for (l in 1:length(Vs)){
        
        strength_of_selection <- Vs [l]
        
        
        
        # init filename
        fname <- paste (basename, "_m", as.character(migration_rate), "_Vs", as.character(strength_of_selection), "_u", as.character(mutation_rate), "_N", as.character(N_value), ".init", sep = "")
        # corresponding sh file
        fname3 <- paste (basename, "_m", as.character(migration_rate), "_Vs", as.character(strength_of_selection), "_u", as.character(mutation_rate), "_N", as.character(N_value), ".slurm", sep = "")
        # output filename for simulations
        fname2 <- paste (basename,"_m", as.character(migration_rate), "_Vs", as.character(strength_of_selection), "_u", as.character(mutation_rate), "_N", as.character(N_value), sep = "")
        
        write.table (temp, fname, row.names =  FALSE, col.names = FALSE, quote = FALSE )
        
        out1 <- paste ("filename ", fname2, sep = "")
        write (out1, fname, append = TRUE)
        
        out1 <- paste ("patch_nbfem ", N_value, sep = "")
        write (out1, fname, append = TRUE)
        
        out1 <- paste ("patch_nbmal ", 0, sep = "")
        write (out1, fname, append = TRUE)
        
        out1 <- paste ("dispersal_rate ", migration_rate, sep = "") ## island migration
        write (out1, fname, append = TRUE)                          ## island migration
        
        # out1 <- paste ("dispersal_matrix ", init_migration_rate, sep = "") ## stepping stone migration
        # write (out1, fname, append = TRUE)                                 ## stepping stone migration
        
        out1 <- paste ("quanti_mutation_rate ", mutation_rate, sep = "")
        write (out1, fname, append = TRUE)
        
        out1 <- paste ("selection_matrix {{", strength_of_selection,"}}", sep = "")
        write (out1, fname, append = TRUE)
        
        out1 <- paste ("random_seed", floor(runif(1)*10^8))
        write (out1, fname, append = TRUE)
        
        
        
        
        write ("#!/bin/bash", fname3)                       ## adjust slurm parameters here
        write ("#SBATCH --nodes=1", fname3, append = TRUE)
        write ("#SBATCH --mem=2000M", fname3, append = TRUE)
        write ("#SBATCH --ntasks=1", fname3, append = TRUE)
        write ("#SBATCH --time=24:00:00", fname3, append = TRUE)
        write ("#SBATCH --cpus-per-task=1", fname3, append = TRUE)
        write ("#SBATCH --account=PLACEHOLDER", fname3, append = TRUE)
        write ("#SBATCH --mail-user=<PLACEHOLDER>", fname3, append = TRUE)
        write ("#SBATCH --mail-type=FAIL", fname3, append = TRUE)
        
        write ("NEMOBIN=$HOME/software/nemo/Nemo-2.3.46/bin/", fname3, append = TRUE) ## change to directory pointing to Nemo program
        write ("WD=$SLURM_TMPDIR", fname3, append = TRUE)
        write (paste0("INPUT=",fname), fname3, append = TRUE)
        write ("ORIG=$PWD", fname3, append = TRUE)
        write ("cp $INPUT $WD/", fname3, append = TRUE)
        write ("cd $WD", fname3, append = TRUE)
        write ("$NEMOBIN/nemo2.3.46 $INPUT", fname3, append = TRUE)
        write ("cp -r * $ORIG/", fname3, append = TRUE)
        write ("rm -rf *", fname3, append = TRUE)
        write ("cd $ORIG", fname3, append = TRUE)
        
        write ("", fname3, append = TRUE)
        
        
      }
    }
  }
}







### ### ### ### ### ###

# Use for stepping stone migration #

# Dispersal matrices #

number_of_patches <- 10

dispersal_matrix_list <- NA

dispersal_matrix_var <- array(0,dim = c(number_of_patches,number_of_patches))

for(mm in 1:length(m)){
  
  choose_m <- m[mm]
  
  diag(dispersal_matrix_var) <- (1-2*choose_m)
  
  dispersal_matrix_var[row(dispersal_matrix_var) == col(dispersal_matrix_var)-1] <- choose_m
  dispersal_matrix_var[row(dispersal_matrix_var) == col(dispersal_matrix_var)+1] <- choose_m
  
  dispersal_matrix_var[1,1] <- 1-choose_m
  dispersal_matrix_var[nrow(dispersal_matrix_var),ncol(dispersal_matrix_var)] <- 1-choose_m
  
  
  empty_var <- vector()
  
  for(i in 1:nrow(dispersal_matrix_var)){
    
    empty_var[i] <- paste0("{",paste0(dispersal_matrix_var[i,], collapse = ","),"}")
    
  }
  
  dispersal_matrix_list[[mm]] <- paste0("{",paste0(empty_var, collapse = ""),"}")
  
}

