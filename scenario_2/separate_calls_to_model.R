


# NOTE: ----------  Main script for Scenario Analysis 2  ----------------

# The aim and assumptions for this scenario analysis are described in detail in the Methods and Supplementary materials. 
# This script is based on the main analysis code but instead calls the procedure 3 times to generate the correct output
# information to compare the three strategies (i.e. no screening, mass-screening, and case-finding based on 1 symptom).

# Since several assumptions are made regarding the asymptomatic population in this scenario analysis (see Methods and Appendix), 
# and since the asymptomatic population in the diagnosed and undiagnosed CD states is present in different proportions depending 
# on the specific strategy, the input parameters need to be adapted accordingly for each strategy (thus requiring separate calls to the model).
# That is why each strategy has a different 'generate_input_parameters' function.

# Also, the last time the 'generate_net_benefit' function is called, an adapted version is called which collects the right 
# output values from the previous two calls to the function (see 'generate_net_benefit_B').




#load required libraries
library(Rcpp, lib.loc = "C:/Program Files/R/R-4.3.1/library")

library(BCEA)
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


# Define global simulation parameters

test <- "POCT"

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

n_cycles <- 100-starting_age 

perspective <- "societal"

#define strategies based on sensitivity and specificity of symptom questionnaire
combinations <- as.data.frame(read_xlsx("strategies.xlsx"))
combinations <- combinations[1:2,] # in this scenario analysis, only mass screening and case-finding on 1 symptom are considered.
combinations_names <- combinations$x
n_combinations <- length(combinations$x)


# Define the number and names of tests

n_tests <- n_combinations + 1 # +1 for no screening

t_names <-  c("No screening", outer(combinations_names, test, FUN = "paste")[1:n_tests-1])

set.seed(14143)
pre_test_probability_overall <- rtri( n = n_samples,
                                      min = 0.007,
                                      max = 0.015,
                                      mode = 0.01)

#call functions that will be used in all calls to model
source('generate_transition_matrices.R')
source('generate_state_qalys.R')
source('generate_state_costs.R')
source('convert_transition_matrices_to_df.R')
source('convert_cohort_vectors_to_df.R')
source('generate_net_benefit.R')

##############################################################################################
#call modified input parameters for no screening scenario
##############################################################################################
source('generate_input_parameters_no_screen.R')


# Generate the input parameters
input_parameters <- generate_input_parameters_no_screen(n_samples = n_samples,
                                                        starting_age = starting_age)

input_parameters$pre_test_probability_overall <- pre_test_probability_overall 

transition_matrices <- generate_transition_matrices(input_parameters)

# Run the Markov model to get the model outputs
output <- generate_net_benefit(input_parameters = input_parameters, 
                               strategies_of_interest = strategies_of_interest, 
                               transition_matrices = transition_matrices,
                               combinations = combinations,
                               lambda=20000)

#create separate output lists for this strategy (to combine later)
output_no_screening <- list()
output_no_screening$total_costs <- output$total_costs['No screening',]
output_no_screening$total_qalys <- output$total_qalys['No screening',]


output_no_screening$time_in_states <- output$time_in_states['No screening',]
output_no_screening$test_costs_applied <- output$test_costs_applied[,'No screening']
output_no_screening$test_costs <- output$test_costs['No screening']
output_no_screening$false_positive_costs_applied <- output$false_positive_costs_applied[,'No screening']
output_no_screening$diagnosis_costs <- output$diagnosis_costs[,'No screening']
output_no_screening$cycle_costs <- output$cycle_costs['No screening']
output_no_screening$cycle_qalys <- output$cycle_qalys['No screening']
output_no_screening$disutility_biopsy <- output$disutility_biopsy[,'No screening']
output_no_screening$unscaled_table <- output$unscaled_table['No screening',]



##############################################################################################
##### call for '0.9999 0 POCT' , i.e. mass-screening
##############################################################################################

source('generate_input_parameters_mass.R')

# Generate the (adapted) input parameters for this strategy and scenario
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters_mass(n_samples = n_samples,
                                                   starting_age = starting_age)

input_parameters$pre_test_probability_overall <- pre_test_probability_overall 

transition_matrices <- generate_transition_matrices(input_parameters)

# Run the Markov model to get the model outputs
output <- generate_net_benefit(input_parameters = input_parameters, 
                                 strategies_of_interest = strategies_of_interest, 
                                 transition_matrices = transition_matrices,
                                 combinations = combinations,
                                 lambda=20000)


output_0.9999_0_POCT <- list()
output_0.9999_0_POCT$total_costs <- output$total_costs['0.9999 0 POCT',]
output_0.9999_0_POCT$total_qalys <- output$total_qalys['0.9999 0 POCT',]


output_0.9999_0_POCT$time_in_states <- output$time_in_states['0.9999 0 POCT',]
output_0.9999_0_POCT$test_costs_applied <- output$test_costs_applied[,'0.9999 0 POCT']
output_0.9999_0_POCT$test_costs <- output$test_costs['0.9999 0 POCT']
output_0.9999_0_POCT$false_positive_costs_applied <- output$false_positive_costs_applied[,'0.9999 0 POCT']
output_0.9999_0_POCT$diagnosis_costs <- output$diagnosis_costs[,'0.9999 0 POCT']
output_0.9999_0_POCT$cycle_costs <- output$cycle_costs['0.9999 0 POCT']
output_0.9999_0_POCT$cycle_qalys <- output$cycle_qalys['0.9999 0 POCT']
output_0.9999_0_POCT$disutility_biopsy <- output$disutility_biopsy[,'0.9999 0 POCT']
output_0.9999_0_POCT$unscaled_table <- output$unscaled_table['0.9999 0 POCT',]


##############################################################################################
#####  call for '0.58 0.63 POCT' , i.e. Glutenscreen (case finding) strategy
##############################################################################################

source('generate_input_parameters_GLUTENSCREEN.R')
source('generate_net_benefit_B.R') #note difference here, this adapted net_benefit function will call the right costs, outcomes for each of the previous strategies. 

# Generate the input parameters
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters_GLUTENSCREEN(n_samples = n_samples,
                                                           starting_age = starting_age)

input_parameters$pre_test_probability_overall <- pre_test_probability_overall 

transition_matrices <- generate_transition_matrices(input_parameters)

# Run the Markov model to get the model outputs
output <- generate_net_benefit_B(input_parameters = input_parameters, 
                                  strategies_of_interest = strategies_of_interest, 
                                  transition_matrices = transition_matrices,
                                  combinations = combinations,
                                  lambda=20000)


##############################################################################################
##############################################################################################

holder <- t_names[t_names %in% c("No screening"  ,   "0.9999 0 POCT"  ,    "0.58 0.63 POCT")] 


m <- bcea(e = t(output$total_qalys[holder,])  , 
          c = t(output$total_costs[holder,]) , 
          ref = 1, 
          interventions = holder,
          Kmax = 30000)
summary(m, wtp = 30000)

m$interventions <- c("No screening", "Mass screening", "Case-finding")

ceac.plot(m, pos = "bottomright",
                cex.main = 0.2,
                theme = theme_minimal(),
                line = list(colors = c(1:length(holder )))) 


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
  labs(title = "Cost-effectiveness acceptability curves for scenario analysis",
       x = "Willingness to pay",
       y = "Probability of cost effectiveness") +
  theme_minimal() +
  scale_color_manual(name = "", values = c("Case finding" = "green", "No screening" = "black", "Mass screening" = "red")) +
  scale_linetype_manual(name = "", values = c("Case finding" = "dashed", "No screening" = "solid", "Mass screening" = "dotted"))

CEAC + personal_theme + theme(legend.position = "bottom")

dev.off()








#get all correct values
output$time_in_states['No screening',] <- output_no_screening$time_in_states
output$time_in_states['0.9999 0 POCT',] <- output_0.9999_0_POCT$time_in_states

output$test_costs_applied[,'No screening'] <- output_no_screening$test_costs_applied
output$test_costs_applied[,'0.9999 0 POCT'] <- output_0.9999_0_POCT$test_costs_applied

output$test_costs['No screening'] <- output_no_screening$test_costs
output$test_costs['0.9999 0 POCT'] <- output_0.9999_0_POCT$test_costs

output$diagnosis_costs[,'No screening'] <- output_no_screening$diagnosis_costs
output$diagnosis_costs[,'0.9999 0 POCT'] <- output_0.9999_0_POCT$diagnosis_costs

output$false_positive_costs_applied[, 'No screening'] <- output_no_screening$false_positive_costs_applied
output$false_positive_costs_applied[, '0.9999 0 POCT'] <- output_0.9999_0_POCT$false_positive_costs_applied

output$unscaled_table['No screening',] <- output_no_screening$unscaled_table
output$unscaled_table['0.9999 0 POCT',] <- output_0.9999_0_POCT$unscaled_table



strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB greater than no screening




#table of costs and QALYs
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

# Create a table comparing costs, QALYs, and incrmental net benefit for strategies of interest
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

#cost breakdown
t_names_testcosts <- output$test_costs_applied[, t_names]
t_names_diagnosiscosts <- output$diagnosis_costs[, t_names]

t_names_cyclecosts <- output$cycle_costs[t_names]
t_names_testcosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_diagnosiscosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_cyclecosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
for (i in 1:length(t_names)) { 
  t_names_testcosts_table[, i] <- format_results(t_names_testcosts[, i])
  t_names_diagnosiscosts_table[, i] <- format_results(t_names_diagnosiscosts[, i])
  t_names_cyclecosts_table[i] <- format_results(t_names_cyclecosts[i])
}

cost_breakdown_table <- data.frame(t(t_names_testcosts_table), t(t_names_diagnosiscosts_table), t(t_names_cyclecosts_table))
colnames(cost_breakdown_table) <- c("Test costs", "Diagnosis costs", "Cycle costs")
rownames(cost_breakdown_table) <- t_names

writexl::write_xlsx(cost_breakdown_table, "results/cost breakdown table.xlsx")


#utility breakdown
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











