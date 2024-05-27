

require(Rcpp)

# Compiles the C++ file for the Markov loop
Rcpp::sourceCpp(file = "rcpp_loop_full.cpp")


generate_net_benefit <- function(input_parameters, 
                                 strategies_of_interest, 
                                 transition_matrices,
                                 combinations,
                                 lambda)
{
  
  # Input costs and utilities
  
  state_qalys <- generate_state_qalys(input_parameters) 
  
  state_costs <- generate_state_costs(input_parameters) 
  
  
  # Sensitivity and specificity parameters
  spec_POCT <- input_parameters$spec_POCT
  sens_POCT <- input_parameters$sens_POCT
  
  
  pre_test_probability_overall <- input_parameters$pre_test_probability_overall 
  
  
  ## Define the true positives, false negatives etc for the symptom questionnaire
  
  tp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tp_riskfactor")))
  fn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fn_riskfactor")))
  fp_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "fp_riskfactor")))
  tn_riskfactor <- array(dim=c(n_samples, n_combinations), dimnames = list(NULL, paste(combinations_names, "tn_riskfactor")))
  
  for (i in 1:n_combinations) {
    tp_riskfactor[,i] <- pre_test_probability_overall * combinations$sens_riskfactor[i]
    fn_riskfactor[,i] <- pre_test_probability_overall - tp_riskfactor[,i]  
    tn_riskfactor[,i] <- (1 - pre_test_probability_overall) * combinations$spec_riskfactor[i]
    fp_riskfactor[,i] <- (1 - pre_test_probability_overall) - tn_riskfactor[,i]
  }
  
  
  tp <- fn <- fp <- tn <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  
  
  for (i in 1:n_combinations) {
    
    #Probabilities for no screening
    tp[,1] <- 0 
    fn[,1] <- pre_test_probability_overall  #in no screening all False Negatives
    tn[,1] <- 1 - pre_test_probability_overall
    fp[,1] <- 0
    
    # Probabilities for POCT
    tp[,i+1] <- (tp_riskfactor[,i] * sens_POCT)
    fn[,i+1] <- tp_riskfactor[,i] - tp[,i+1]  
    tn[,i+1] <- ((fp_riskfactor[,i]) * spec_POCT)
    fp[,i+1] <- (fp_riskfactor[,i]) - tn[,i+1]
  }
  
  
  #Adding costs   

  diagnosis_costs <- array(dim = c(n_samples, n_tests),
                           dimnames = list(NULL, t_names))
  
  diagnosis_costs[,] <- ((tp[,] + fp[,]) * input_parameters$cost_diagnosis_s)/pre_test_probability_overall  #costs of post-POCT confirmatory diagnostics (the input parameter includes cost of biopsy*probability of biopsy in this group)
  
  
  test_costs <- array(dim=c(n_samples, n_tests), dimnames=list(NULL, t_names))
  for (i in 1:n_combinations) {
    test_costs[, 1] <- 0
    
    # POCT
    test_costs[, i+1] <-  input_parameters$cost_POCT
    
  }
  
  
  #Adding disutilities
  
  disutility_biopsy_screen <- array(dim=c(n_samples,n_tests),dimnames=list(NULL, t_names))
  for (i in 1:n_combinations) {
    disutility_biopsy_screen[, 1] <- 0
    disutility_biopsy_screen[, i+1] <-  input_parameters$disutility_biopsy_children 
  }
  
  
  # Convert to tables for convenience
  tp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  fp_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  fn_riskfactor_table <- array(dim=c(n_samples, n_tests), dimnames = list(NULL, t_names))
  for (i in 1:n_combinations) {
    tp_riskfactor_table[,1] <- 1
    fp_riskfactor_table[,1] <- 1
    fn_riskfactor_table[,1] <- 1
    tp_riskfactor_table[,i + 1] <- tp_riskfactor[, i]
    fp_riskfactor_table[,i + 1] <- fp_riskfactor[, i]
    fn_riskfactor_table[,i + 1] <- fn_riskfactor[, i]
  }
  
  test_costs_applied <- array(dim = c(n_samples, n_tests),
                              dimnames = list(NULL, t_names))
  
  test_costs_applied <- (test_costs * (tp_riskfactor_table[,] + fp_riskfactor_table[,]))/pre_test_probability_overall #anyone meeting the symptom criteria (i.e. tp or fp on questionnaire) gets a POCT
  false_positive_costs_applied <- (test_costs * (fp_riskfactor_table[,]))/pre_test_probability_overall #these Costs represent costs of testing in those without CD after risk factor test [in our case the questionnaire] 
  biopsy_disutility_applied <- (disutility_biopsy_screen * ((tp[,] + fp[,])*5/61))/pre_test_probability_overall # only about 8% of those who tested positive on the POCT will receive a biopsy

  
  fn_riskfactor_table <- fn_riskfactor_table 
  fn_all <- fn + fn_riskfactor_table
  
  
  #Scale up true positives and false negatives to obtain only the cohort of CD cases
  tp <- tp/(tp + fn_all)
  fn <- 1 - tp
  
  
  # Get transition matrices as a data frame
  transition_matrices_df <- convert_transition_matrices_to_df(transition_matrices)

  
  # Build an array to store the cohort vector at each cycle
  # Each cohort vector has n_states elements: probability of being in each state,
  # There is one cohort vector for each test, for each PSA sample, for each cycle.
  # Faster to first create cohort vectors as an array and then convert to data frame
  cohort_vectors<-array(dim=c(n_tests,n_samples,n_cycles,n_states),  
                        dimnames=list(t_names,NULL,NULL, state_names))
  
  
  for (i_test in c(1:n_tests)) { 
    cohort_vectors[i_test, , 1, "diagnosed no complications"] <- tp[, i_test] * input_parameters$probability_nocomplications 
    cohort_vectors[i_test, , 1, "diagnosed osteoporosis"] <- tp[, i_test] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_test, , 1, "diagnosed NHL"] <- tp[, i_test] * input_parameters$probability_NHL
    cohort_vectors[i_test, , 1, "diagnosed GIC"] <- tp[, i_test] * input_parameters$probability_GIC
    
    cohort_vectors[i_test, , 1, "Undiagnosed no complications"] <- fn[, i_test] * input_parameters$probability_nocomplications 
    cohort_vectors[i_test, , 1, "Undiagnosed osteoporosis"] <- fn[, i_test] * input_parameters$probability_osteoporosis 
    cohort_vectors[i_test, , 1, "Undiagnosed NHL"] <- fn[, i_test] * input_parameters$probability_NHL 
    cohort_vectors[i_test, , 1, "Undiagnosed GIC"] <- fn[, i_test] * input_parameters$probability_GIC 
    
    
    cohort_vectors[i_test, , 1, "Death"] <- 0
  }
  
  
  cohort_vectors[, , , ] [is.na(cohort_vectors[, , , ] )] <- 0
  cohort_vectors_df <- convert_cohort_vectors_to_df(cohort_vectors)
  
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each test, for each PSA sample, for each cycle
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_tests, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  cycle_qalys <- array(dim = c(n_tests, n_samples, n_cycles),
                       dimnames = list(t_names, NULL, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each test and each PSA sample
  # These are filled in below using cycle_costs, 
  # test_costs, and cycle_qalys
  total_costs <- array(dim = c(n_tests, n_samples),
                       dimnames = list(t_names, NULL))
  total_qalys <- array(dim = c(n_tests, n_samples),
                       dimnames = list(t_names, NULL))
  
  
  disc_vec_costs <- (1/1.04) ^ rep(c(0:(n_cycles-1)))
  disc_vec_outcomes <- (1/1.015) ^ rep(c(0:(n_cycles-1)))
  
  
  
  # Function currently requires conversion to matrices and returns a matrix, not a dataframe
  cohort_vectors_df <- 
    (rcpp_loop_full(cohort_vectors_in = as.matrix(cohort_vectors_df), 
                    transition_matrices = as.matrix(transition_matrices_df),
                    n_cycles = n_cycles, n_tests = n_tests, n_samples = n_samples,
                    n_states = n_states))
  
 
  lapply(c(1:n_tests), function(i_test){
    # Pre-index to reduce runtime
    cycle_costs_tr <- cycle_costs[i_test, , ]
    cycle_qalys_tr <- cycle_qalys[i_test, , ]
    total_costs_tr <- total_costs[i_test, ]
    total_qalys_tr <- total_qalys[i_test, ]
    test_costs_applied_tr <- test_costs_applied[ , i_test]
    diagnosis_costs_tr <- diagnosis_costs[, i_test]
    biopsy_disutility_applied_tr <- biopsy_disutility_applied[, i_test]

    # In this case state qalys and costs are the same for all strategies 
    state_costs_tr <- state_costs[, , ]
    state_qalys_tr <- state_qalys[, , ]
    
    for(i_sample in 1:n_samples) {
      cohort_vectors_tr_sample <- cohort_vectors_df[c(0:(n_cycles-1)) * (n_tests * n_samples) +
                                                      (i_test - 1) * n_samples + i_sample, 3 + c(1:n_states) ]
      
      # Use cohort vectors to calculate cycle costs and qalys
      cycle_costs_tr[i_sample, ] <- rowSums(cohort_vectors_tr_sample[,  ] * state_costs_tr[i_sample, , ])
      cycle_qalys_tr[i_sample, ] <- rowSums(cohort_vectors_tr_sample[,  ] * state_qalys_tr[i_sample, , ])
      
      # Sum  and discount to get total costs and qalys
      total_costs_tr[i_sample] <- cycle_costs_tr[i_sample, ] %*% disc_vec_costs + 
        test_costs_applied_tr[i_sample] + 
        diagnosis_costs_tr[i_sample] 
      
      total_qalys_tr[i_sample] <- cycle_qalys_tr[i_sample, ] %*% disc_vec_outcomes -
        biopsy_disutility_applied[i_sample, i_test] 
      
    } # End loop over samples
    
    return(list(total_qalys = total_qalys_tr, total_costs = total_costs_tr, cycle_qalys = cycle_qalys_tr, cycle_costs = cycle_costs_tr))
    
  }) -> output_list # End lapply   
  
  names(output_list) <- t_names
  
  # Reconvert result to a matrix
  total_costs <- sapply(t_names, function(test) {total_costs[test, ] <- output_list[[test]]$total_costs})
  total_qalys <- sapply(t_names, function(test) {total_qalys[test, ] <- output_list[[test]]$total_qalys})
  cycle_costs <- sapply(t_names, function(test) {cycle_costs[test,, ] <- output_list[[test]]$cycle_costs})
  cycle_qalys <- sapply(t_names, function(test) {cycle_qalys[test,, ] <- output_list[[test]]$cycle_qalys})
  # sapply inverts the matrices to uninvert
  total_costs <- t(total_costs)  #add Questionnaire cost per CD case here (below)
  total_costs[2:n_tests,] <- total_costs[2:n_tests,]+(input_parameters$cost_quest[1]/mean(pre_test_probability_overall)) 
  total_costs['0.9999 0 POCT',] <- total_costs['0.9999 0 POCT',]-(input_parameters$cost_quest[1]/mean(pre_test_probability_overall)) 
  total_qalys <- t(total_qalys)
  
  
  
  
  #############################################################################
  ## Analysis of results ######################################################
  #############################################################################
  output <- list()
  
  
  #Costs and QALYs
  output$total_costs <- total_costs  
  output$total_qalys <- total_qalys
  # Average costs
  output$average_costs <- rowMeans(total_costs)
  # Average effects (in QALY units)
  output$average_effects <- rowMeans(total_qalys)
  
  #for all TPs/FNs comparison
  quantile(output$total_costs, 0.025)
  quantile(output$total_costs, 0.975)
  quantile(output$total_qalys, 0.025)
  quantile(output$total_qalys, 0.975)
  
  output$incremental_costs <-  output$average_costs - output$average_costs["No screening"]
  output$incremental_effects <-  output$average_effects -  output$average_effects["No screening"]
  
  output$all_incremental_costs <-   output$all_incremental_effects <- array(dim=c(n_tests,n_samples),  
                                                                            dimnames=list(t_names,NULL))
  
  for (i_sample in 1:n_samples) {
    for (i_test in 1:n_tests) {
      output$all_incremental_costs[i_test, ] <-  output$total_costs[i_test,] - output$total_costs["No screening",]
      output$all_incremental_effects[i_test, ] <-  output$total_qalys[i_test,] -  output$total_qalys["No screening",]
    }}
  
  
  output$ICER <- output$incremental_costs/output$incremental_effects
  
  # Incremental net benefit at the 20,000 willingness-to-pay
  
  output$incremental_net_benefit <- lambda*output$incremental_effects - output$incremental_costs
  
  output$net_benefit<- lambda*output$average_effects - output$average_costs
  output$all_net_benefit <- lambda*total_qalys - total_costs
  
  quantile(output$all_net_benefit, 0.025)
  quantile(output$all_net_benefit, 0.975)
  
  output$ceac_calculation <- array(dim=c(n_tests,n_samples),  
                                   dimnames=list(t_names,NULL))
  
  for (i_sample in 1:n_samples) {
    for (i_test in 1:n_tests) {
      output$ceac_calculation[i_test,i_sample] <- ifelse(output$all_net_benefit[i_test,i_sample] == max(output$all_net_benefit[,i_sample]), 1, 0)
    }}
  
  output$probability_best <- rowMeans(output$ceac_calculation)
  
  #time spent in states
  time_in_states <- matrix(nrow = length(t_names), ncol = n_states)
  rownames(time_in_states) <- t_names
  colnames(time_in_states) <- state_names
  for(test_name in t_names) {
    time_in_states[test_name, ] <- colMeans(cohort_vectors_df[which(cohort_vectors_df[, "test"] == which(t_names == test_name)), 3 + c(1:n_states)])
  }      
  output$time_in_states <- time_in_states[t_names, ]
  
  
  
  
  #cost breakdown
  output$test_costs_applied <- test_costs_applied
  output$test_costs <- colMeans(test_costs_applied) #costs of test 
  output$false_positive_costs_applied <- false_positive_costs_applied
  output$diagnosis_costs <- diagnosis_costs #cost if diagnostic procedure (incl. biopsies) after positive POCT
  output$cycle_costs <- colMeans(cycle_costs)
  
  
  #utility breakdown
  output$cycle_qalys <- colMeans(cycle_qalys)
  output$disutility_biopsy <- biopsy_disutility_applied

  
  
  # Probability cost-effective
  # This is the proportion of samples for which the incremental net benefit is positive
  output$all_incremental_net_benefit <- lambda*output$all_incremental_effects - output$all_incremental_costs
  output$incremental_net_benefit <- rowMeans(output$all_incremental_net_benefit)
  
  
  output$inb_lci <-   output$inb_uci <- array(dim=c(n_tests,n_samples),  
                                              dimnames=list(t_names,NULL))
  
  for (i_sample in 1:n_samples) {
    for (i_test in 1:n_tests) {
      output$inb_lci[i_test, i_sample] <- quantile(output$all_incremental_net_benefit[i_test,], 0.025)
      output$inb_uci[i_test, i_sample] <- quantile(output$all_incremental_net_benefit[i_test,], 0.975)
      
    }
  }
  
  output$inb_lci <- rowMeans(output$inb_lci)
  
  output$inb_uci <- rowMeans(output$inb_uci)
  
  output$pce_calculation <- array(dim=c(n_tests,n_samples),  
                                  dimnames=list(t_names,NULL))
  
  for (i_sample in 1:n_samples) {
    for (i_test in 1:n_tests) {
      output$pce_calculation[i_test,i_sample] <- sum(output$all_incremental_net_benefit[i_test,] > 0)/n_samples
    }}
  
  output$probability_cost_effective <- rowMeans(output$pce_calculation)
  output$pce_lci <-   output$pce_uci <- array(dim=c(n_tests,n_samples),  
                                              dimnames=list(t_names,NULL))
  
  for (i_sample in 1:n_samples) {
    for (i_test in 1:n_tests) {
      output$pce_lci[i_test, i_sample] <- quantile(output$pce_calculation[i_test,], 0.025)
      output$pce_uci[i_test, i_sample] <- quantile(output$pce_calculation[i_test,], 0.975)
      
    }
  }

  
  return(output)
}






