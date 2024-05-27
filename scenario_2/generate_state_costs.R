


generate_state_costs <- function(input_parameters) {
  
  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  #copy probability_IDA for inclusion in state costs
  probability_IDA <- data.frame(input_parameters$probability_IDA_0, input_parameters$probability_IDA_10, input_parameters$probability_IDA_20, input_parameters$probability_IDA_30, input_parameters$probability_IDA_40, 
                                input_parameters$probability_IDA_50, input_parameters$probability_IDA_60, input_parameters$probability_IDA_70, input_parameters$probability_IDA_80, input_parameters$probability_IDA_90)
  
  
  # costs per state over cycles
  
  #call correct values to use in cost calculation for state
  for (i_cycle_tmp in 1:n_cycles) {
    i_age_tmp <- sprintf("%02d", i_cycle_tmp + (starting_age-1) )
    col_ref_tmp <- as.integer(substr(i_age_tmp, 1,1))
    if (as.numeric(i_age_tmp)<18) {
      probability_late_diagnosis <- input_parameters$probability_late_diagnosis_children
      cost_CDGFD <- input_parameters$cost_CDGFD_children
      cost_undiagnosedCD  <- input_parameters$cost_undiagnosedCD_children
      cost_diagnosis_ns <- input_parameters$cost_diagnosis_children_ns
    } else {
      probability_late_diagnosis <- input_parameters$probability_late_diagnosis_adults
      cost_CDGFD <- input_parameters$cost_CDGFD_adults
      cost_undiagnosedCD  <- input_parameters$cost_undiagnosedCD_adults
      cost_diagnosis_ns <- input_parameters$cost_diagnosis_adults_ns
    }
    #account for retirement at age 67
    if (as.numeric(i_age_tmp)>20 & as.numeric(i_age_tmp)<67) {
      cost_abs_diagnosed_fc <- input_parameters$cost_abs_diagnosed_fc
      cost_abs_undiagnosed_fc <- input_parameters$cost_abs_undiagnosed_fc
    }else{
      cost_abs_diagnosed_fc <- 0
      cost_abs_undiagnosed_fc <- 0 
    }
    
    if (perspective == "healthcare") {
      cost_abs_diagnosed_fc <- 0
      cost_abs_undiagnosed_fc <- 0 
      input_parameters$cost_gfp <- 0
    }
    
    state_costs[, i_cycle_tmp, "diagnosed no complications"] <- cost_CDGFD +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      input_parameters$cost_gfp + cost_abs_diagnosed_fc

    state_costs[, i_cycle_tmp, "diagnosed osteoporosis"]  <- input_parameters$cost_osteoporosis +
      cost_CDGFD +
      input_parameters$cost_gfp +
      (probability_IDA[,col_ref_tmp+1] * input_parameters$cost_IDA) +
      cost_abs_diagnosed_fc

    state_costs[, i_cycle_tmp, "diagnosed GIC"]  <- input_parameters$cost_GIC +
      cost_CDGFD +
      input_parameters$cost_gfp +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_diagnosed_fc

    state_costs[, i_cycle_tmp, "diagnosed NHL"] <-  input_parameters$cost_NHL +
      cost_CDGFD +
      input_parameters$cost_gfp +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_diagnosed_fc
    
    
    
    state_costs[, i_cycle_tmp, "Undiagnosed no complications"] <- cost_undiagnosedCD +
      (probability_late_diagnosis * cost_diagnosis_ns) +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_undiagnosed_fc

    state_costs[, i_cycle_tmp, "Undiagnosed osteoporosis"] <- cost_undiagnosedCD +
      input_parameters$cost_osteoporosis +
      (probability_late_diagnosis * cost_diagnosis_ns) +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_undiagnosed_fc

    state_costs[, i_cycle_tmp, "Undiagnosed GIC"] <- input_parameters$cost_GIC +
      cost_undiagnosedCD +
      (probability_late_diagnosis * cost_diagnosis_ns) +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_undiagnosed_fc

    state_costs[, i_cycle_tmp, "Undiagnosed NHL"] <- input_parameters$cost_NHL +
      cost_undiagnosedCD +
      (probability_late_diagnosis * cost_diagnosis_ns) +
      (probability_IDA[, col_ref_tmp + 1] * input_parameters$cost_IDA) +
      cost_abs_undiagnosed_fc
    
    

  }
  
  
  
  
  state_costs[, , "Death"] <- 0
  return(state_costs[, , ])
}

