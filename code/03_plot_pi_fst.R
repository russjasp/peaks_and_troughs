

## I generate plots using the parameter-specific results files (Objects: Nucleotide_diversity_results/Fst_results in 02A/02B), so that I can continuously add data as I go over time
## One could simply use the "Master_results" files generated in 02A/02B to plot if they were to generate ALL there data in one fell swoop


rm(list=ls(all=TRUE))
library(ggplot2)
options(scipen = 999) 

setwd()
Total_list_of_all_nemo_raw_results_files <- list.files(pattern = "Pi*") ## Choose Pi or Fst here

DATE <- "" ## input date for naming results files



template_genetic_map <- c(90000, 91000, 92000, 93000, 94000, 95000, 96000, 97000, 98000, 99000, 99100, 99200, 99300, 99400,  99500, 99600, 99700, 99800, 99900, 99910, 99920, 99930, 99940, 99950, 99960, 99970, 99980, 99990, 99991, 99992, 99993, 99994, 99995, 99996, 99997, 99998, 99999, 100001, 100002, 100003, 100004, 100005, 100006, 100007, 100008, 100009, 100010, 100020, 100030, 100040, 100050, 100060, 100070, 100080, 100090, 100100, 100200, 100300, 100400, 100500, 100600, 100700, 100800, 100900, 101000, 102000, 103000, 104000, 105000, 106000, 107000, 108000, 109000, 110000)
transformed_neutral_genetic_map <- log10(abs(template_genetic_map - 100000))


Total_Plot_Results <- array(NA, dim = c(0,12), dimnames = list(c(),c("slope","log_slope","mean_pi_fst","close_pi_fst","far_pi_fst","log_migration","log_mutation","number_of_QTL","Vs","pop","time","reps")))


different_pop_sizes_1 <- gsub(".*_N_","",Total_list_of_all_nemo_raw_results_files)
different_pop_sizes <- unique(as.numeric(gsub("_.*","",different_pop_sizes_1)))

number_of_different_pop_sizes <- length(different_pop_sizes)

for(n in 1:length(number_of_different_pop_sizes)){
  
  choose_pop <- different_pop_sizes[n]
  
  pop_size_probe <- paste0("_N_",choose_pop,"_")
  
  all_files_this_pop_size <- Total_list_of_all_nemo_raw_results_files[grepl(pop_size_probe, Total_list_of_all_nemo_raw_results_files)]
  
  different_times1 <- gsub(".*_t_","",all_files_this_pop_size)
  different_times <- unique(gsub("_.*","",different_times1))
  
  number_of_different_times <- length(different_times)
  
  for(t in 1:number_of_different_times){
    
    choose_time <- different_times[t]
    
    time_probe <- paste0("t_",choose_time,"_")
    
    all_files_this_time <- all_files_this_pop_size[grepl(time_probe, all_files_this_pop_size)]
    
    different_migrations <- unique(gsub(".*_mig_|_mut.*","",all_files_this_time))
    
    number_of_different_migrations <- length(different_migrations)
    
    for(m in 1:number_of_different_migrations){
      
      choose_migration <- different_migrations[m]
      
      migration_probe <- paste0("_mig_",choose_migration,"_mut")
      
      all_files_this_migration <- all_files_this_time[grepl(migration_probe, all_files_this_time)]
      
      different_Vs <- unique(gsub(".*_Vs_|_qtl.*","",all_files_this_migration))
      
      number_of_different_Vs <- length(different_Vs)
      
      for(v in 1:number_of_different_Vs){
        
        choose_Vs <- different_Vs[v]
        
        Vs_probe <- paste0("_Vs_",choose_Vs,"_")
        
        all_files_this_Vs <- all_files_this_migration[grepl(Vs_probe, all_files_this_migration)]
        
        different_mutation <- unique(gsub(".*_mut_|_Vs.*","",all_files_this_Vs))
        
        number_of_different_mutation <- length(different_mutation)
        
        for(u in 1:number_of_different_mutation){
          
          choose_mutation <- different_mutation[u]
          
          mutation_probe <- paste0("_mut_",choose_mutation,"_Vs")
          
          all_files_this_mutation <- all_files_this_Vs[grepl(mutation_probe, all_files_this_Vs)]
          
          Different_qtl <- 1 ## adjust when using multilocus
          
          Numbers_of_different_qtl <- length(Different_qtl)
          
          for(q in 1:Numbers_of_different_qtl){
            
            choose_qtl <- Different_qtl[q]
            Number_of_ntrl_loci <- choose_qtl*74
            Temporary_Array_to_Collect_Results <- array(NA, dim = c(Number_of_ntrl_loci,2), dimnames = list(c(),c("mean Pi Fst","distance")))
            
            # qtl_probe <- paste0("") ## adjust when using multilocus
            
            Subset_list_of_all_nemo_raw_results_files <- all_files_this_mutation
            
            Number_of_different_nemo_files <- length(Subset_list_of_all_nemo_raw_results_files)
            
            
          
            csv_results_list <- vector("list", length = Number_of_different_nemo_files)
            
            
            for(r in 1:Number_of_different_nemo_files){
              
              csv_results_var <- read.csv(Subset_list_of_all_nemo_raw_results_files[r], header = TRUE, stringsAsFactors = FALSE)
              
              colnames(csv_results_var) <- paste0(colnames(csv_results_var),"_csv",r)
              
              csv_results_list[[r]] <- csv_results_var
              
            }
            
            csv_results <- do.call("cbind", csv_results_list)
            
            remove_NA_columns_dataframe_variable <- Filter(function(x)!all(is.na(x)), csv_results)
            
            remove_NA_columns_dataframe_variable <- as.matrix(remove_NA_columns_dataframe_variable)
            
            if("mean" %in% colnames(remove_NA_columns_dataframe_variable)){ remove_NA_columns_dataframe_variable <- remove_NA_columns_dataframe_variable[ , -which(colnames(remove_NA_columns_dataframe_variable) %in% "mean")] }
            
            number_of_reps_variable <- dim(remove_NA_columns_dataframe_variable)[2]
            
            
            TEMP_Plot_Results <- array(NA, dim = c(1,12), dimnames = list(c(),c("slope","log_slope","mean_pi_fst","close_pi_fst","far_pi_fst","log_migration","log_mutation","number_of_QTL","Vs","pop","time","reps")))
            
            ## SLOPE METHODS 1 ##
            
            ### Calculate Mean Pi or Fst ### POPULATIONS ARE COMBINED ###
            
            Temporary_Array_to_Collect_Results[,1] <- rowMeans(remove_NA_columns_dataframe_variable,na.rm = TRUE)
            
            genetic_map_variable <- rep(template_genetic_map, times = choose_qtl)
            
            transformed_map_variable <- log10(abs(genetic_map_variable - 100000))
            

              
            Temporary_Array_to_Collect_Results[,2] <- transformed_map_variable
              

            Mod_variable11 <- lm(Temporary_Array_to_Collect_Results[,1]~Temporary_Array_to_Collect_Results[,2])
              
            Slope_variable11 <- summary(Mod_variable11)[[4]][2]
              
            TEMP_Plot_Results[1,1] <- Slope_variable11
              
            if(sum(Temporary_Array_to_Collect_Results[,1]<=0) > 0 ){print("Trying to log(0)")} else {
            
              Log10_Mod_variable11 <- lm(log10(Temporary_Array_to_Collect_Results[,1])~Temporary_Array_to_Collect_Results[,2])
              Log10_Slope_variable11 <- summary(Log10_Mod_variable11)[[4]][2]
              TEMP_Plot_Results[1,2] <- Log10_Slope_variable11
            
              }
            
            Average_Pi_Fst_variable11 <- mean(Temporary_Array_to_Collect_Results[,1])
              
            TEMP_Plot_Results[1,3] <- Average_Pi_Fst_variable11
              
            Close_Pi_Fst <- (Temporary_Array_to_Collect_Results[c(Number_of_ntrl_loci/2),1] + Temporary_Array_to_Collect_Results[c(Number_of_ntrl_loci/2+1),1]) / 2
              
            TEMP_Plot_Results[1,4] <- Close_Pi_Fst
              
            Far_Pi_Fst <- (Temporary_Array_to_Collect_Results[1,1] + Temporary_Array_to_Collect_Results[Number_of_ntrl_loci,1]) / 2
              
            TEMP_Plot_Results[1,5] <- Far_Pi_Fst
              
              
            TEMP_Plot_Results[,6] <- log10(as.numeric(choose_migration))
            TEMP_Plot_Results[,7] <- log10(as.numeric(choose_mutation))
            TEMP_Plot_Results[,8] <- as.numeric(choose_qtl)
            TEMP_Plot_Results[,9] <- as.numeric(choose_Vs)
            TEMP_Plot_Results[,10] <- as.numeric(choose_pop)
            TEMP_Plot_Results[,11] <- as.numeric(choose_time)
            TEMP_Plot_Results[,12] <- as.numeric(number_of_reps_variable)
              
            Total_Plot_Results <- rbind(Total_Plot_Results,TEMP_Plot_Results)
              
              
            
          } ## qtl
        } ## mutation
      } ## Vs
    } ## migration
  } ## time
} ## pop size

setwd()
write.csv(Total_Plot_Results, paste0("PI_or_FST_plotting_results_",DATE,".csv"), row.names = FALSE)


## massage data
Total_Plot_Results <- data.frame(Total_Plot_Results)
Total_Plot_Results$log_mutation <- as.factor(Total_Plot_Results$log_mutation)
Total_Plot_Results$Vs <- as.factor(Total_Plot_Results$Vs)
Total_Plot_Results$number_of_QTL <- as.factor(Total_Plot_Results$number_of_QTL)
Total_Plot_Results$pop <- as.factor(Total_Plot_Results$pop)
Total_Plot_Results$time <- as.factor(Total_Plot_Results$time)

Close_DF <- Total_Plot_Results[!grepl("far_pi_fst", names(Total_Plot_Results))]
names(Close_DF)[which(names(Close_DF)=="close_pi_fst")] <- "pi_fst";Close_DF$type <- "Nearby"
Far_DF <- Total_Plot_Results[!grepl("close_pi_fst", names(Total_Plot_Results))]
names(Far_DF)[which(names(Far_DF)=="far_pi_fst")] <- "pi_fst";Far_DF$type <- "Distant"
Close_Far_DF <- rbind(Close_DF,Far_DF)
Close_Far_DF$type <- as.factor(Close_Far_DF$type);Close_Far_DF$type <- factor(Close_Far_DF$type, levels = c("Nearby", "Distant"))



## PLOTS ##

# Figure 1
# Slope of metric against migration rate (log10) and Vs
ggplot(Total_Plot_Results, aes(x=log_migration, y=slope, colour = Vs)) + geom_point(aes(), size=2.25) + geom_line(aes(), linetype=2, size=0.875) +
  ylab(expression("Slope of Nucleotide Diversity versus Distance ("*Log[10]*" cM)")) +
  # ylab(expression("Slope of "*F[ST]*" versus Distance ("*Log[10]*" cM)")) +
  xlab(expression("Migration Rate ("*Log[10]*")")) +
  labs(colour = bquote("Strength of\nSelection ("*italic(V)[S]*")")) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.5)

# Figure 2
# Metric at loci 1/1000 cM and 10 cM away from adaptive locus against migration rate (log10) and Vs
ggplot(Close_Far_DF, aes(x=log_migration, y=pi_fst, colour = Vs)) +
  geom_line(aes(linetype=type), size=0.875) + geom_point(aes(shape = type), size=2.25) +
  scale_linetype_manual(name="Locus Position", values = c("Distant" = 4, "Nearby" = 2), labels = c(bquote(10^-3~"cM"),"10 cM")) +
  scale_shape_manual(name="Locus Position", values = c("Distant" = 1, "Nearby" = 16), labels = c(bquote(10^-3~"cM"),"10 cM")) +
  ylab("Nucleotide Diversity") +
  # ylab(expression(F[ST])) +
  xlab(expression("Migration Rate ("*Log[10]*")")) +
  labs(colour = bquote("Strength of\nSelection ("*italic(V)[S]*")"))

# Figure 3
# Mean metric per locus against migration rate (log10) and Vs
ggplot(Total_Plot_Results, aes(x=log_migration, y=mean_pi_fst, colour = Vs)) + geom_point(aes(), size=2.25) + geom_line(aes(), linetype=2, size=0.875) +
  ylab("Per Locus Mean Nucleotide Diversity") +
  # ylab(expression("Per Locus Mean "*F[ST])) +
  xlab(expression("Migration Rate ("*Log[10]*")")) +
  labs(colour = bquote("Strength of\nSelection ("*italic(V)[S]*")"))





