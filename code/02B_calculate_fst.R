rm(list=ls(all=TRUE))
options(scipen = 999) 

setwd()
Total_list_of_all_nemo_raw_results_files <- list.files(pattern = "*freq") ## The .freq Nemo output files in ntrl_output_dir will be used to calculate the Fst

DATE <- "" ## input date for naming results files



neutral_genetic_map <- c(-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) 
neutral_genetic_map_no_sign <- abs(neutral_genetic_map)
transformed_neutral_genetic_map <- log10(neutral_genetic_map_no_sign)


Master_results <- array(data = NA, dim = c(0,11), dimnames = list(c(),c("slope","R.sq","log.slope","log.R.sq","mean.fst","qtl","log.migration.rate","log.mutation.rate","Vs","pop.size","time")))

Temp_results <- array(data = NA, dim = c(1,11), dimnames = list(c(),c("slope","R.sq","log.slope","log.R.sq","mean.fst","qtl","log.migration.rate","log.mutation.rate","Vs","pop.size","time")))


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
            
            
            ## Let there be fst
            for(r in 1:Number_of_different_nemo_files){
              
              
              if(r == 1){
                
                
                Fst_results <- array(NA,dim = c(Number_of_ntrl_loci,Number_of_different_nemo_files),
                                     dimnames = list(c(paste0("l",1:Number_of_ntrl_loci)),
                                                     c(paste0("r",1:Number_of_different_nemo_files)))
                                     
                )
              }
              
              
              
              
              master_results_variable <- tryCatch(read.table(Subset_list_of_all_nemo_raw_results_files[r], header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
              
              if( is.null(master_results_variable) == TRUE ) { print(paste("Results text file is completely empty",Subset_list_of_all_nemo_raw_results_files[r])) } else {
                
                master_results_variable <- as.matrix(master_results_variable)
                
                
                if(dim(master_results_variable)[2] < 13){ print("Results test file has an incomplete header line") } else {
                  
                  
                  Fst_var <- master_results_variable[,8]
                  
                  Fst_var[is.nan(Fst_var)] <- 0
                  
                  if(length(Fst_var) != Number_of_ntrl_loci){ print("Incomplete txt file (missing some loci)") }
                  
                  Fst_results[,r] <- Fst_var
                  
                  
                  
                  
                  
                  
                } # end of else (incomplete header check)
                
              } # end of else (empty results file check)
              
            } # end of replicates for loop
            
            ## Write results for this particular parameter set
            filename_var <- paste0("Fst_mig_", choose_migration, "_mut_", choose_mutation, "_Vs_", choose_Vs, "_qtl_", choose_qtl, "_N_", choose_pop, "_t_", choose_time, "_", DATE, ".csv")
            # setwd()
            write.csv(Fst_results, filename_var, row.names = FALSE)
            
            
            
            genetic_map_qtl_specific <- rep(transformed_neutral_genetic_map, times = choose_qtl)
            
            
            
            
            
            Fst_means <- rowMeans(Fst_results, na.rm = TRUE)
            
            
            lm_mod_variable <- lm(Fst_means~genetic_map_qtl_specific)
            
            
            slope_variable <- summary(lm_mod_variable)[[4]][2]
            
            
            R_sq_variable <- summary(lm_mod_variable)$r.squared
            
            
            Temp_results[,1] <- slope_variable
            
            Temp_results[,2] <- R_sq_variable
            
            
            if(sum(Fst_means == 0) > 0) { 
              
              Temp_results[,3:4] <- NA
              
            } else {          
              
              
              log_lm_mod_variable <- lm(log10(Fst_means)~genetic_map_qtl_specific)
              
              
              
              log_slope_variable <- summary(log_lm_mod_variable)[[4]][2]
              
              log_R_sq_variable <- summary(log_lm_mod_variable)$r.squared
              
              Temp_results[,3] <- log_slope_variable
              
              Temp_results[,4] <- log_R_sq_variable
              
            }
            
            
            average_Fst <- mean(Fst_means, na.rm = TRUE)
            
            Temp_results[,5] <- average_Fst
            
            
            Temp_results[,6] <- as.numeric(choose_qtl)
            Temp_results[,7] <- log10(as.numeric(choose_migration))
            Temp_results[,8] <- log10(as.numeric(choose_mutation))
            Temp_results[,9] <- as.numeric(choose_Vs)
            Temp_results[,10] <- as.numeric(choose_pop)
            Temp_results[,11] <- as.numeric(choose_time)
            
            
            Master_results <- rbind(Master_results, Temp_results)
              

            
          
          } ## qtl
        } ## mutation
      } ## Vs
    } ## migration
  } ## time
} ## pop size

## Write total results for all parameter sets
setwd()
write.csv(Master_results, paste0("FST_master_results_",DATE,".csv"), row.names = FALSE)



