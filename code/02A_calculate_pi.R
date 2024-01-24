rm(list=ls(all=TRUE))
options(scipen = 999) 

Total_list_of_all_nemo_raw_results_files <- list.files(pattern = "*.txt") ## The .txt Nemo output files in ntrl_output_dir will be used to calculate the Pi

DATE <- "jan24" ## input date for naming results files



neutral_genetic_map <- c(-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) 
neutral_genetic_map_no_sign <- abs(neutral_genetic_map)
transformed_neutral_genetic_map <- log10(neutral_genetic_map_no_sign)


Master_results <- array(data = NA, dim = c(0,12), dimnames = list(c(),c("slope","R.sq","log.slope","log.R.sq","mean.pi","patch","qtl","log.migration.rate","log.mutation.rate","Vs","pop.size","time")))

Temp_results <- array(data = NA, dim = c(1,12), dimnames = list(c(),c("slope","R.sq","log.slope","log.R.sq","mean.pi","patch","qtl","log.migration.rate","log.mutation.rate","Vs","pop.size","time")))


different_pop_sizes_1 <- gsub(".*_N","",Total_list_of_all_nemo_raw_results_files)
different_pop_sizes <- unique(as.numeric(gsub("_.*","",different_pop_sizes_1)))

number_of_different_pop_sizes <- length(different_pop_sizes)


for(n in 1:length(number_of_different_pop_sizes)){
  
  choose_pop <- different_pop_sizes[n]

  pop_size_probe <- paste0("_N",choose_pop,"_")
  
  all_files_this_pop_size <- Total_list_of_all_nemo_raw_results_files[grepl(pop_size_probe, Total_list_of_all_nemo_raw_results_files)]

  different_times1 <- gsub(paste0(".*",pop_size_probe),"",all_files_this_pop_size)
  different_times <- unique(gsub("_.*","",different_times1))
  
  number_of_different_times <- length(different_times)
  
  for(t in 1:number_of_different_times){
    
    choose_time <- different_times[t]
    
    time_probe <- paste0(pop_size_probe,choose_time,"_")
    
    all_files_this_time <- all_files_this_pop_size[grepl(time_probe, all_files_this_pop_size)]
    
    different_migrations <- unique(gsub(".*_m|_Vs.*","",all_files_this_time))
    
    number_of_different_migrations <- length(different_migrations)
    
    for(m in 1:number_of_different_migrations){
      
      choose_migration <- different_migrations[m]
      
      migration_probe <- paste0("_m",choose_migration,"_Vs")
      
      all_files_this_migration <- all_files_this_time[grepl(migration_probe, all_files_this_time)]
      
      different_Vs <- unique(gsub(".*_Vs|_u.*","",all_files_this_migration))
      
      number_of_different_Vs <- length(different_Vs)
      
      for(v in 1:number_of_different_Vs){
      
        choose_Vs <- different_Vs[v]
        
        Vs_probe <- paste0("_Vs",choose_Vs,"_")
        
        all_files_this_Vs <- all_files_this_migration[grepl(Vs_probe, all_files_this_migration)]
        
        different_mutation <- unique(gsub(".*_u|_N.*","",all_files_this_Vs))
        
        number_of_different_mutation <- length(different_mutation)
      
        for(u in 1:number_of_different_mutation){
      
          choose_mutation <- different_mutation[u]
          
          mutation_probe <- paste0("_u",choose_mutation,"_N")
          
          all_files_this_mutation <- all_files_this_Vs[grepl(mutation_probe, all_files_this_Vs)]
          
          Different_qtl <- 1 ## adjust when using multilocus
          
          Numbers_of_different_qtl <- length(Different_qtl)
          
          for(q in 1:Numbers_of_different_qtl){
          
            choose_qtl <- Different_qtl[q]
            Number_of_ntrl_loci <- choose_qtl*74
            
            # qtl_probe <- paste0("") ## adjust when using multilocus
            
            Subset_list_of_all_nemo_raw_results_files <- all_files_this_mutation
            
            Number_of_different_nemo_files <- length(Subset_list_of_all_nemo_raw_results_files)
            
            
            ## Let there be pi
            for(r in 1:Number_of_different_nemo_files){
              
              master_results_variable <- tryCatch(read.table(Subset_list_of_all_nemo_raw_results_files[r], header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
              
              if( is.null(master_results_variable) == TRUE ) { print(paste("Results text file is completely empty",Subset_list_of_all_nemo_raw_results_files[r])) } else {
                
                master_results_variable <- as.matrix(master_results_variable)
                
                
                if(dim(master_results_variable)[2] < (c(2+2*Number_of_ntrl_loci))){ print("Results test file has an incomplete header line") } else {
                  
                  
                  master_results_only_adults <- master_results_variable[master_results_variable[,c(2+2*Number_of_ntrl_loci)] == 2,]
                  
                    nonNA_patches <- unique(master_results_only_adults[,1])[!is.na(unique(master_results_only_adults[,1]))]
                    
                    Number_of_patches <- length(unique(nonNA_patches))
                  
                  
                  if(r == 1){
                    
                    
                    Nucleotide_diversity_results <- array(NA, 
                                                          dim = c(Number_of_ntrl_loci,Number_of_patches,Number_of_different_nemo_files),
                                                          dimnames = list(c(paste0("l",1:Number_of_ntrl_loci)),
                                                                          c(paste0("p",1:Number_of_patches)),
                                                                          c(paste0("r",1:Number_of_different_nemo_files)))
                                                          
                    )
                  }
                  
                  
                  
                  for(p in 1:Number_of_patches){
                    
                    
                    master_results_only_adults_by_patch <- master_results_only_adults[!is.na(master_results_only_adults[,1]) & master_results_only_adults[,1] == p,]
                    
                    Allele_population_size <- dim(master_results_only_adults_by_patch)[1] * 2
                    
                    
                    for(l in 1:Number_of_ntrl_loci){
                      
                      
                      allele_counts_at_one_locus <- c(master_results_only_adults_by_patch[,2*l],
                                                      master_results_only_adults_by_patch[,2*l+1])
                      
                      
                      allele_table <- table(allele_counts_at_one_locus)
                      
                      
                      if( dim(allele_table) == 0 ) { temp_Nucleotide_diversity_results[l,p] <- NA } else {
                        
                        
                        An_allele_frequency <- allele_table[[1]] / Allele_population_size
                        
                        pairwise_differences_new_method <- (An_allele_frequency) * (1 - An_allele_frequency) * (Allele_population_size)^2
                        
                        
                        Nucleotide_diversity_results[l,p,r] <- (pairwise_differences_new_method /
                                                                  ( choose(Allele_population_size, 2) ) )
                        
                        
                        
                        
                      } # end of else (missing locus information)
                      
                    } # end of loci for loop
                    
                  } # end of patch for loop
                  
                } # end of else (incomplete header check)
                
              } # end of else (empty results file check)
              
            } # end of replicates for loop
            
            ## Write results for this particular parameter set
            filename_var <- paste0("Pi_mig_", choose_migration, "_mut_", choose_mutation, "_Vs_", choose_Vs, "_qtl_", choose_qtl, "_N_", choose_pop, "_t_",choose_time, "_", DATE, ".csv")
            # setwd()
            write.csv(Nucleotide_diversity_results, filename_var, row.names = FALSE)

          
          } ## qtl
        } ## mutation
      } ## Vs
    } ## migration
  } ## time
} ## pop size

## Write total results for all parameter sets
## setwd()
write.csv(Master_results, paste0("PI_master_results_",DATE,".csv"), row.names = FALSE)
