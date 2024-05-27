


# Define function to generate transition matrices
generate_transition_matrices <- function(input_parameters) {
  
  n_samples <- dim(input_parameters)[1]
  
  
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  # Generate empty 'transition matrices' array. 
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),  
                               dimnames = list(NULL, NULL, state_names, state_names))
  
  
  
  # Osteoporosis probabilities 
  # Generate probabilities based on these rates and ratios
  osteoporosis_probability_GFD <- osteoporosis_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    osteoporosis_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10*(i_age_category - 1))] + input_parameters$log_or_osteoporosis_GFD))
    osteoporosis_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_or_osteoporosis_noGFD))
  }
  osteoporosis_probability_GFD_all <- as.data.frame(osteoporosis_probability_GFD)
  colnames(osteoporosis_probability_GFD_all) <- c("osteoporosis_probability_GFD_0","osteoporosis_probability_GFD_10","osteoporosis_probability_GFD_20","osteoporosis_probability_GFD_30",
                                                  "osteoporosis_probability_GFD_40","osteoporosis_probability_GFD_50","osteoporosis_probability_GFD_60",
                                                  "osteoporosis_probability_GFD_70","osteoporosis_probability_GFD_80","osteoporosis_probability_GFD_90")
  osteoporosis_probability_noGFD_all <- as.data.frame(osteoporosis_probability_noGFD)
  colnames(osteoporosis_probability_noGFD_all) <- c("osteoporosis_probability_noGFD_0","osteoporosis_probability_noGFD_10","osteoporosis_probability_noGFD_20","osteoporosis_probability_noGFD_30",
                                                    "osteoporosis_probability_noGFD_40","osteoporosis_probability_noGFD_50","osteoporosis_probability_noGFD_60",
                                                    "osteoporosis_probability_noGFD_70","osteoporosis_probability_noGFD_80","osteoporosis_probability_noGFD_90")
  
  
  # Calculate NHL probabilities
  NHL_probability_GFD <- NHL_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    NHL_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("NHL_lograte_", 10*(i_age_category - 1))] + input_parameters$log_rr_NHL_GFD))
    NHL_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("NHL_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_rr_NHL_noGFD))
  }
  NHL_probability_GFD_all <- as.data.frame(NHL_probability_GFD)
  colnames(NHL_probability_GFD_all) <- c("NHL_probability_GFD_0","NHL_probability_GFD_10","NHL_probability_GFD_20","NHL_probability_GFD_30",
                                         "NHL_probability_GFD_40","NHL_probability_GFD_50","NHL_probability_GFD_60",
                                         "NHL_probability_GFD_70","NHL_probability_GFD_80","NHL_probability_GFD_90")
  NHL_probability_noGFD_all <- as.data.frame(NHL_probability_noGFD)
  colnames(NHL_probability_noGFD_all) <- c("NHL_probability_noGFD_0","NHL_probability_noGFD_10","NHL_probability_noGFD_20","NHL_probability_noGFD_30",
                                           "NHL_probability_noGFD_40","NHL_probability_noGFD_50","NHL_probability_noGFD_60",
                                           "NHL_probability_noGFD_70","NHL_probability_noGFD_80","NHL_probability_noGFD_90")
  
  #calulate GIC probabilities
  GIC_probability_GFD <- GIC_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    GIC_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("GIC_lograte_", 10*(i_age_category - 1))] + input_parameters$log_rr_GIC_GFD))
    GIC_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("GIC_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_rr_GIC_noGFD))
  }
  GIC_probability_GFD_all <- as.data.frame(GIC_probability_GFD)
  colnames(GIC_probability_GFD_all) <- c("GIC_probability_GFD_0","GIC_probability_GFD_10","GIC_probability_GFD_20","GIC_probability_GFD_30",
                                         "GIC_probability_GFD_40","GIC_probability_GFD_50","GIC_probability_GFD_60",
                                         "GIC_probability_GFD_70","GIC_probability_GFD_80","GIC_probability_GFD_90")
  GIC_probability_noGFD_all <- as.data.frame(GIC_probability_noGFD)
  colnames(GIC_probability_noGFD_all) <- c("GIC_probability_noGFD_0","GIC_probability_noGFD_10","GIC_probability_noGFD_20","GIC_probability_noGFD_30",
                                           "GIC_probability_noGFD_40","GIC_probability_noGFD_50","GIC_probability_noGFD_60",
                                           "GIC_probability_noGFD_70","GIC_probability_noGFD_80","GIC_probability_noGFD_90")
  

  
  lifetables <- read_xlsx("lifetable/genpop_mort_rates.xlsx", sheet = "input")
  lifetables$lograte <- log(lifetables$mx)
  # Lifetables are rates so must be converted to probabilities
  lifetables$Overall <- 1 - exp(-lifetables$mx)
  death_probability_nocomplications	<- data.frame(lifetables$age, lifetables$Overall)

  
  

  # Combine on log scale
  # Rows are for samples, columns are for ages
  death_probability_osteoporosis <- list()
  death_probability_osteoporosis$Overall <- 1 - exp(-exp(matrix(rep(input_parameters$death_log_hr_osteoporosis_mixedj, length(lifetables$lograte)) +   rep(lifetables$lograte, each = n_samples), nrow = n_samples)))
  colnames(death_probability_osteoporosis$Overall) <- lifetables$age

  
  
  #Fill in the array with the correct transition probabilities over the cycles
  
  #transition probabilities over cycles
  for (i_cycle_tmp in 1:n_cycles) {
    i_age_tmp <- sprintf("%02d", i_cycle_tmp + (starting_age-1) )
    col_ref_tmp <- as.integer(substr(i_age_tmp, 1,1))
    if (as.numeric(i_age_tmp)<18) {
      probability_late_diag_tmp <- input_parameters$probability_late_diagnosis_children
    } else {
      probability_late_diag_tmp <- input_parameters$probability_late_diagnosis_adults
    }
    #from state to GIC
    transition_matrices[,i_cycle_tmp, "diagnosed no complications", "diagnosed GIC"] <- 
      transition_matrices[, i_cycle_tmp, "diagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, i_cycle_tmp, "diagnosed NHL", "diagnosed GIC"] <- 
      GIC_probability_GFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "Undiagnosed GIC"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed osteoporosis", "Undiagnosed GIC"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed NHL", "Undiagnosed GIC"] <-
      (1 - probability_late_diag_tmp) * GIC_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed GIC"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed NHL", "diagnosed GIC"] <-
      probability_late_diag_tmp * GIC_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp,"Undiagnosed GIC", "diagnosed GIC"] <- probability_late_diag_tmp 
    
    #from 'state' to 'NHL'
    transition_matrices[, i_cycle_tmp, "diagnosed no complications", "diagnosed NHL"] <- 
      transition_matrices[, i_cycle_tmp, "diagnosed osteoporosis", "diagnosed NHL"] <-
      NHL_probability_GFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "Undiagnosed NHL"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed osteoporosis", "Undiagnosed NHL"] <-
      (1 - probability_late_diag_tmp) * NHL_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed NHL"] <-
      transition_matrices[, i_cycle_tmp, "Undiagnosed osteoporosis", "diagnosed NHL"] <-
      probability_late_diag_tmp * NHL_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp,"Undiagnosed NHL", "diagnosed NHL"] <- probability_late_diag_tmp - 
      (transition_matrices[, i_cycle_tmp, "Undiagnosed NHL","diagnosed GIC"])
    #from state to Osteoporosis
    transition_matrices[, i_cycle_tmp, "diagnosed no complications", "diagnosed osteoporosis"] <- osteoporosis_probability_GFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "Undiagnosed osteoporosis"] <- (1 - probability_late_diag_tmp) * osteoporosis_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed osteoporosis"] <-  probability_late_diag_tmp * osteoporosis_probability_noGFD_all[,col_ref_tmp+1]
    transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed no complications"] <- probability_late_diag_tmp - (
      (transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed osteoporosis"]) +
        (transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed NHL"]) +
        (transition_matrices[, i_cycle_tmp, "Undiagnosed no complications", "diagnosed GIC"])
    ) 
    transition_matrices[,i_cycle_tmp ,"Undiagnosed osteoporosis", "diagnosed osteoporosis"] <- probability_late_diag_tmp - (
      (transition_matrices[, i_cycle_tmp,"Undiagnosed osteoporosis", "diagnosed NHL"]) +
        (transition_matrices[, i_cycle_tmp,"Undiagnosed osteoporosis", "diagnosed GIC"])
    )
    
    
  }
 
  
  
  n_ages <- 100 - starting_age
  

  death_state_index <- which(state_names == "Death")
  for (i_age in 1:n_ages){
    transition_matrices[, i_age, "diagnosed no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "diagnosed osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "diagnosed osteoporosis", -death_state_index] <-  transition_matrices[, i_age, "diagnosed osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
    transition_matrices[, i_age, "Undiagnosed no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "Undiagnosed osteoporosis", -death_state_index] <- transition_matrices[, i_age, "Undiagnosed osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
  }
  
  # Probabilities of death form NHL
  transition_matrices[, , "diagnosed NHL", "Death"] <- 
    transition_matrices[, , "Undiagnosed NHL", "Death"] <- input_parameters$death_probability_NHL
  
  # Probabilities of death form GIC
  transition_matrices[, , "diagnosed GIC", "Death"] <- 
    transition_matrices[, , "Undiagnosed GIC", "Death"] <- input_parameters$death_probability_GIC
  
  # Multiply other transitions to prevent probabilities exceeding 1  
  for(i_sample in 1:n_samples) {
    transition_matrices[i_sample, , "Undiagnosed NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "Undiagnosed NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
    transition_matrices[i_sample, , "diagnosed NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "diagnosed NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
  }
  
  for(i_sample in 1:n_samples) {
    transition_matrices[i_sample, , "Undiagnosed GIC", -death_state_index] <- 
      transition_matrices[i_sample, , "Undiagnosed GIC", -death_state_index] * (1 - input_parameters$death_probability_GIC[i_sample])
    transition_matrices[i_sample, , "diagnosed GIC", -death_state_index] <- 
      transition_matrices[i_sample, , "diagnosed GIC", -death_state_index] * (1 - input_parameters$death_probability_GIC[i_sample])
  }
  
  
  #Complete the matrix by adding the complement of all probabilities  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , i_state, -i_state], c(1,2), sum, na.rm=TRUE)
  }
  
  
  transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
  
  # Check that rows sum to 1
  rowSums (transition_matrices[1, 4, , ], na.rm = FALSE , dims = 1)
  
  return(transition_matrices[, , , ])
}

