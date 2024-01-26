


generate_state_qalys <- function(input_parameters) {
  n_samples <- dim(input_parameters)[1]
  
  
  # define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  
  state_qalys <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  eq5d_norms <- read.csv("eq5d_norms_nl.csv")
  eq5d_norms$age <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)

  
  
  
  #costs per state over cycles
  
  #call correct values to use in cost calculation for state
  for (i_cycle_tmp in 1:n_cycles) {
    i_age_tmp <- sprintf("%02d", i_cycle_tmp + (starting_age-1) )
    col_ref_tmp <- as.integer(substr(i_age_tmp, 1,1))
    if (as.numeric(i_age_tmp)<18) {
      probability_late_diagnosis <- input_parameters$probability_late_diagnosis_children
      utility_GFD <- input_parameters$utility_GFD_children
      utility_undiagnosedCD <- input_parameters$utility_undiagnosedCD_children
      probability_biopsy <- input_parameters$probability_biopsy_children 
      disutility_biopsy <- input_parameters$disutility_biopsy_children
    } else {
      probability_late_diagnosis <- input_parameters$probability_late_diagnosis_adults
      utility_GFD <- input_parameters$utility_GFD_adults
      utility_undiagnosedCD <- input_parameters$utility_undiagnosedCD_adults
      probability_biopsy <- input_parameters$probability_biopsy_adults
      disutility_biopsy <- input_parameters$disutility_biopsy_adults
      
    }
  
    state_qalys[, i_cycle_tmp, "diagnosed no complications"] <- eq5d_norms[col_ref_tmp + 1, 2] * utility_GFD
    state_qalys[, i_cycle_tmp, "Undiagnosed no complications"] <- eq5d_norms[col_ref_tmp + 1, 2] * utility_undiagnosedCD - (probability_late_diagnosis * probability_biopsy * disutility_biopsy) 
    
    
  }
  
  
  
  state_qalys[, , "diagnosed osteoporosis"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_osteoporosis 
  state_qalys[, , "diagnosed NHL"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "diagnosed GIC"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_GIC
  state_qalys[, , "Undiagnosed osteoporosis"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_osteoporosis
  state_qalys[, , "Undiagnosed NHL"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Undiagnosed GIC"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_GIC
  state_qalys[, , "Death"] <- 0
  return(state_qalys[, , ])
}