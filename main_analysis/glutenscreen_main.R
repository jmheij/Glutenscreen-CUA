
#######################################################################################
#This script calls all the necessary functions and inputs to conduct the main analysis. 
#Run this script to recreate results from the main analysis.
#######################################################################################
rm(list = ls())
#load required packages (install if necessary)
library(Rcpp, lib.loc = "C:/Program Files/R/R-4.3.1/library")
library(readxl)
library(ggplot2)
library(reshape)
library(Rmisc)
library(dplyr)
library(tictoc)
library(SimDesign)
library(BCEA)
library(EnvStats)
library(car)


set.seed(14143)

#load necessary functions
source('generate_input_parameters.R')
source('generate_transition_matrices.R')
source('generate_state_qalys.R')
source('generate_state_costs.R')
source('convert_transition_matrices_to_df.R')
source('convert_cohort_vectors_to_df.R')
source('generate_net_benefit.R')


# Define global simulation parameters
perspective <- "societal" #can also set to "healthcare"

test <- "POCT" #POCT = point of care test

n_samples <- 1000

n_states <- 9
state_names <- c("diagnosed no complications",
                 "diagnosed osteoporosis",
                 "diagnosed NHL",
                 "diagnosed GIC",
                 "Undiagnosed no complications",
                 "Undiagnosed osteoporosis",
                 "Undiagnosed NHL",
                 "Undiagnosed GIC",
                 "Death")


starting_age <- 3 

n_cycles <- 100-starting_age # 83 = average LE in the Netherlands


# Define strategies based on sensitivity and specificity of symptom questionnaire
# Note: the main analysis compares two testing strategies out of the table below: "0.9999 0" and "0.58 0.63" are mass-screening and case-finding respectively
# rows 3-27 in the 'combinations' data frame below are further (two-way sensitivity) analyses.
combinations <- as.data.frame(read_xlsx("strategies.xlsx"))
combinations_names <- combinations$x
n_combinations <- length(combinations$x)


# Define the number and names of tests
n_tests <- (1 * n_combinations) + 1 # +1 for no screening
t_names <-  c("No screening", outer(combinations_names, test, FUN = "paste")[1:n_tests-1])


# Generate the input parameters
input_parameters <- generate_input_parameters(n_samples = n_samples,
                                                starting_age = starting_age)

#generate transition matrices
transition_matrices <- generate_transition_matrices(input_parameters)

# Run the Markov model to get the model outputs (note 'lambda' = WTP threshold)
output <- generate_net_benefit(input_parameters = input_parameters, 
                                 strategies_of_interest = strategies_of_interest, 
                                 transition_matrices = transition_matrices,
                                 combinations = combinations,
                                 lambda=20000)


#to compare only the strategies that form part of the main analysis (i.e. mass-screening and case finding)
holder <- t_names[t_names %in% c("No screening"  ,   "0.9999 0 POCT"  ,    "0.58 0.63 POCT")] 

m <- bcea(e = t(output$total_qalys[holder,])  , 
          c = t(output$total_costs[holder,]) , 
          ref = 1, 
          interventions = holder,
          Kmax = 30000)
summary(m, wtp=30000)

#separate ceacs compared to no screening
ceac.plot(m, pos = "topright",
          theme = theme_minimal(),
          line = list(color=c("red", "blue")),
          comparison = c(2,3))


#multiple comparisons
mce <- multi.ce(m)
mce$interventions <- c("No screening", "Mass screening", "Case-finding")
ceac.plot(mce, pos = "topright", 
               theme = theme_minimal(),
               line = list(color = c(1:length(holder)), types=c(1,4,2))) 

#redone manually for better display
mcomps <- as.data.frame(mce$p_best_interv) 
colnames(mcomps) <- c("Mass_screening", "Case_finding", "No_screening")
mcomps$xax <- mce$k

personal_theme = theme(plot.title = element_text(hjust = 0.5, size=20), axis.text = element_text(size=14), axis.title = element_text(size=18), legend.text = element_text(size=14))
jpeg(file = "results/CEAC.jpeg", quality = 100, height = 15, width = 20, units = "cm", res = 300)

CEAC <- ggplot(mcomps, aes(x = xax)) +
  geom_line(aes(y = Case_finding, color = "Case finding", linetype = "Case finding")) +
  geom_line(aes(y = No_screening, color = "No screening", linetype = "No screening")) +
  geom_line(aes(y = Mass_screening, color = "Mass screening", linetype = "Mass screening")) +
  labs(title = "Cost-effectiveness acceptability curves for main analysis",
       x = "Willingness to pay",
       y = "Probability of cost effectiveness") +
  theme_minimal() +
  scale_color_manual(name = "", values = c("Case finding" = "green", "No screening" = "black", "Mass screening" = "red")) +
  scale_linetype_manual(name = "", values = c("Case finding" = "dashed", "No screening" = "solid", "Mass screening" = "dotted"))

CEAC + personal_theme + theme(legend.position = "bottom")

dev.off()
#old equivalent of the above code. The below code stopped working after updating to latest versions of R and ggplot2 (should work in older versions)
# jpeg(file = "results/CEAC.jpeg", quality = 100, height = 15, width = 20, units = "cm", res = 300)
# CEAC <- ceac.plot(mce, graph = "ggplot2", pos = "bottom",
#                   theme = theme_minimal(),
#                   line = list(colors = c(1:length(holder)), types=c(1,4,2))) 
# CEAC + personal_theme 
# dev.off()


# Create a table (i.e. table of results) comparing costs, QALYs, IBN, and ICERs for strategies
format_results<-function(x, n_digits = 2) {
  paste(round(mean(x, na.rm = TRUE),digits = n_digits), " (",
        round(quantile(x, probs = 0.025, na.rm = TRUE), digits = n_digits), ", ", 
        round(quantile(x, probs = 0.975, na.rm = TRUE), digits = n_digits),")",sep="")
}


t_names_qalys <- output$total_qalys[t_names, ]
t_names_costs <- output$total_costs[t_names, ]
t_names_inetbenefit <- output$all_incremental_net_benefit[t_names, ]
t_names_qalys_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_costs_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_inetbenefit_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))

for (i in 1:length(t_names)) { 
  t_names_qalys_table[, i] <- format_results(t_names_qalys[i, ], n_digits = 2)
  t_names_costs_table[, i] <- format_results(t_names_costs[i, ], n_digits = 0)
  t_names_inetbenefit_table[, i] <- format_results(t_names_inetbenefit[i, ], n_digits = 0)
}

results_table <- data.frame(
  c("", rep(combinations$sens_riskfactor, length(n_tests))),
  c("", rep(combinations$spec_riskfactor, length(n_tests))),
  t(t_names_costs_table), t(t_names_qalys_table), t(t_names_inetbenefit_table))
colnames(results_table) <- c("Sensitivity", "Specificity", "Costs", "QALYs", "Incremental net benefit v no screening")

# Change the sensitivity and specificity of 0.9999 back to 1 for ease of presentation
results_table$Sensitivity[results_table$Sensitivity == 0.9999] <- 1
results_table$Specificity[results_table$Specificity == 0.9999] <- 1
rownames(results_table) <- t_names
results_table$ICER <- output$ICER
writexl::write_xlsx(results_table, "results/table of results.xlsx")



#Cost breakdown table (only breaks down test costs, diagnosis costs, and cycle costs)
  #Note: for the cost breakdown presented in the manuscript (i.e. Tables 1 and 2), model cost parameters need to be changed.
t_names_testcosts <- output$test_costs_applied[, t_names]
t_names_diagnosiscosts <- output$diagnosis_costs[, t_names]

t_names_cyclecosts <- output$cycle_costs[t_names]
t_names_testcosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_diagnosiscosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_cyclecosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
for (i in 1:length(t_names)) { 
  t_names_testcosts_table[, i] <- format_results(t_names_testcosts[, i], n_digits = 0)
  t_names_diagnosiscosts_table[, i] <- format_results(t_names_diagnosiscosts[, i], n_digits = 0)
  t_names_cyclecosts_table[i] <- format_results(t_names_cyclecosts[i], n_digits = 0)
}

cost_breakdown_table <- data.frame(t(t_names_testcosts_table), t(t_names_diagnosiscosts_table), t(t_names_cyclecosts_table))
colnames(cost_breakdown_table) <- c("Test costs", "Diagnosis costs", "Cycle costs")
rownames(cost_breakdown_table) <- t_names

writexl::write_xlsx(cost_breakdown_table, "results/cost breakdown table.xlsx")



#Utility breakdown table
t_names_disutilitybiopsy <- output$disutility_biopsy[, t_names]
t_names_cycleqalys <- output$cycle_qalys[t_names]
t_names_disutilitybiopsy_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_cycleqalys_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
for (i in 1:length(t_names)) { 
  t_names_disutilitybiopsy_table[, i] <- format_results(t_names_disutilitybiopsy[, i], n_digits = 4)
  t_names_cycleqalys_table[i] <- format_results(t_names_cycleqalys[i], n_digits = 4)
}

qaly_breakdown_table <- data.frame(t(t_names_cycleqalys_table), t(t_names_disutilitybiopsy_table))
colnames(qaly_breakdown_table) <- c("Cycle QALYs", "Disutility biopsy")
rownames(qaly_breakdown_table) <- t_names

writexl::write_xlsx(qaly_breakdown_table, "results/qaly breakdown table.xlsx")



#Time in states table
time_in_states <- output$time_in_states
writexl::write_xlsx(as.data.frame(time_in_states), "results/time in states.xlsx")


################################################################################################
# The below code is relevant only for the further two-way analyses looking at strategies 4-28 in t_names
# i.e. the strategies that were not part of the main analysis but represent combinations of
# implied sensitivities and specificities of a symptom questionnaire depending on the criteria set for the 
# numbers of symptoms.
################################################################################################

# Need a simple function to pick out sensitivity, specificity and POCT strategy
inb_strategies <- function(sensitivity_vector, specificity_vector, sero_test) {
  # Identify results for this POCT strategy
  temp <- grepl(sero_test, names(output$incremental_net_benefit))
  # And remove serological tests that might have strings in common
  for(wrong_test in tests[tests != sero_test]) {
    # Only remove if wrong_test is not a substring of sero_test
    if(!grepl(wrong_test, sero_test)) {
      temp <- temp & !grepl(wrong_test, names(output$incremental_net_benefit))  
    }
  }
  # And allow user to pass vector of sensitivity and specificity
  sens_spec_temp <- rep(FALSE, length(output$incremental_net_benefit))
  for(sensitivity in sensitivity_vector) {
    for(specificity in specificity_vector) {
      sens_spec_temp <- sens_spec_temp |
        grepl(paste0(sensitivity, " ", specificity, " "), names(output$incremental_net_benefit))
    }
  }
  # Right serological test and sens/spec in user provided vectors
  temp <- temp & sens_spec_temp
  return(temp)
}


jpeg(file = "results/ibnplot.jpeg", quality = 100, height = 15, width = 20, units = "cm", res = 300)
par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(5, 1, 0))
# All tests
formatted_test_names <- c("POCT")
tests <- "POCT"
names(formatted_test_names) <- tests
# Each of the lines includes six sensitivities
# One line for each specificity
for(sero_test in tests) {
  plot(output$incremental_net_benefit[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999), 
                                                     specificity_vector = 0.2, sero_test = sero_test)],
       pch=19, ylim=c(-10000,80000), ylab = "Incremental net benefit", xlab = "Sensitivity", xaxt="n" , col = 0, font.main = 1, cex.lab=1.5, cex.axis=1.2 , las = 1 )
  abline(h=0)
  specificity_vector <- c(0.2, 0.4, 0.6, 0.8)
  # All specificities
  for(i_specificity in 1:length(specificity_vector)) {
    lines(output$incremental_net_benefit[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999), 
                                                        specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
          col = i_specificity, lwd=2, lty= i_specificity)  
    
    # # Lower credible limit
    lines(output$inb_lci[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999), 
                                        specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
          col = i_specificity, lwd=1, lty= i_specificity) 
    # Upper credible limit
    lines(output$inb_uci[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999), 
                                        specificity_vector = specificity_vector[i_specificity], sero_test = sero_test)],
          col = i_specificity , lwd=1, lty= i_specificity) 
  } # End loop over specificity
}
# Add a legend
legend("right", legend=c("0.2", "0.4", "0.6", "0.8"),
       col= 1:length(specificity_vector), lwd=2, lty=1:length(specificity_vector), title = "Specificity", box.col = "transparent", cex=1.2)
axis(1,                         # Define x-axis manually
     at = 1:5,
     labels = c(0.2, 0.4, 0.6, 0.8, 1), font.axis = 1, cex.axis=1.2)

dev.off()



#plots of probability CE
pdf(paste0("results/probce.pdf"))

for(sero_test in tests) {
  plot(output$probability_cost_effective[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999),
                                                        specificity_vector = 0.2,
                                                        sero_test = sero_test)],pch=19, ylim=c(0,1), ylab = "Probability CE", xlab = "Sensitivity", xaxt="n" , cex.lab=1.5, cex.axis=1.2 )
  abline(h=0)
  
  specificity_vector <- c(0.2, 0.4, 0.6, 0.8)
  # All specificities
  for(i_specificity in 1:length(specificity_vector)) {
    lines(output$probability_cost_effective[inb_strategies(sensitivity_vector = c(0.2, 0.4, 0.6, 0.8, 0.9999),
                                                           specificity_vector = specificity_vector[i_specificity],
                                                           sero_test = sero_test)],
          col = i_specificity, lwd = 2, lty= i_specificity)
  }
  
}

# Add a legend
legend("right", legend=c("0.2", "0.4", "0.6", "0.8"),
       col= 1:length(specificity_vector), lwd=2, lty=1:length(specificity_vector), title = "Specificity", box.col = "transparent", cex=1.2)
axis(1,                         # Define x-axis manually
     at = 1:5,
     labels = c(0.2, 0.4, 0.6, 0.8, 1), font.axis = 1, cex.axis=1.2)
dev.off()


