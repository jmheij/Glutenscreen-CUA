
# Note: search for **ADAPTED** below to see what adaptations are made in this script in comparison to the 
# generate_input_parameters function used in the main analysis.



# Utility functions
expit <- function(logO) {
  return(exp(logO)/(1 + exp(logO)))
}
logit <- function(p) {
  return(log(p/(1-p)))
}



generate_input_parameters_GLUTENSCREEN <- function(n_samples,
                                        starting_age = starting_age) {
  
  set.seed(14143)
sens_POCT <- logit(0.94)
sens_POCT_lci <- logit(0.899)
sens_POCT_uci <- logit(0.965)
sd = (sens_POCT_uci-sens_POCT_lci)/3.92
sens_POCT <- rnorm(n=n_samples, sens_POCT, sd=sd)
sens_POCT <- expit(sens_POCT)

spec_POCT <- logit(0.944)
spec_POCT_lci <- logit(0.909)
spec_POCT_uci <- logit(0.965)
sd = (spec_POCT_uci-spec_POCT_lci)/3.92
spec_POCT <- rnorm(n=n_samples, spec_POCT, sd=sd)
spec_POCT <- expit(spec_POCT)

  
  # **ADAPTED**
  symp_prop <- 0.58
  asymp_prop <- 0.42
  tp_symp <- symp_prop*mean(sens_POCT) # symptomatics that immediately go to tp in markov
  symp_undiag <- symp_prop-tp_symp # symptomatics missed due to POCT sensitivity 
  undiag_all <- asymp_prop + symp_undiag
  

  prop_children_diagnosis <- 0.40 # based on: Bergman et al 2022. 

  
  # **ADAPTED**
  p_diag_ever <- (1/3)*symp_undiag  #lifetime probability of getting a late diagnosis
  p_diag_children <- prop_children_diagnosis*p_diag_ever  #probability of late diagnosis occurring in childhood.
   
  probability_late_diagnosis_children <- 1 - (1-p_diag_children)^(1/17) 
  probability_late_diagnosis_adults <- 1 - ((1-p_diag_ever)/(1-p_diag_children))^(1/35)  #the 35 is the result of calibrating in order to achieve the assumption that 1 in 3 people reach the diagnosed state throughout the model.  
  
  
  # Use the appropriate prevalence
  prevalence = as.data.frame(readxl::read_excel(path = "prevalence/adapted from CPRD prevalence.xlsx", sheet = "mixed")) 
  
  #Initial cohort at diagnosis - depends on age at diagnosis 
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_category <- 1
  # prevalence$Age.categories == starting_age
  probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "Osteoporosis_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "Osteoporosis_r"])
  probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "NHL_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "NHL_r"])
  probability_GIC <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "GIC_r"], shape2 = prevalence[starting_age_category, "GIC_N"] - prevalence[starting_age_category, "GIC_r"])
  probability_nocomplications <- 1 - probability_osteoporosis - probability_NHL - probability_GIC
  
  # IDA prevalence changes with age of cohort so is age stratified
  # Prevalence of osteoporosis and NHL are combined with incidence rates to model prevalence changing with age
  probability_IDA_0 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[1], shape2 = (prevalence$N[1] - prevalence$IDA_r[1]))
  probability_IDA_10 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[2], shape2 = (prevalence$N[2] - prevalence$IDA_r[2]))
  probability_IDA_20 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[3], shape2 = (prevalence$N[3] - prevalence$IDA_r[3]))
  probability_IDA_30 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[4], shape2 = (prevalence$N[4] - prevalence$IDA_r[4]))
  probability_IDA_40 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[5], shape2 = (prevalence$N[5] - prevalence$IDA_r[5]))
  probability_IDA_50 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[6], shape2 = (prevalence$N[6] - prevalence$IDA_r[6]))
  probability_IDA_60 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[7], shape2 = (prevalence$N[7] - prevalence$IDA_r[7]))
  probability_IDA_70 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[8], shape2 = (prevalence$N[8] - prevalence$IDA_r[8]))
  probability_IDA_80 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[9], shape2 = (prevalence$N[9] - prevalence$IDA_r[9]))
  probability_IDA_90 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[10], shape2 = (prevalence$N[10] - prevalence$IDA_r[10]))
  probability_IDA <- data.frame(probability_IDA_0 , probability_IDA_10 , probability_IDA_20 , probability_IDA_30 , probability_IDA_40
                                ,probability_IDA_50 , probability_IDA_60 , probability_IDA_70 , probability_IDA_80 , probability_IDA_90)
  
  
  #####Developing consequence
  
 ## Developing Osteoporosis (Osteoporosis probabilities On GFD)
  
  Osteoporosis_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/ostp_rates_NL.xlsx", sheet = "Sheet1")) 
  Osteoporosis_probability$Osteoporosis_rate_low <- Osteoporosis_probability$Osteoporosis_rate - (1.96*(Osteoporosis_probability$Osteoporosis_rate/10))
  Osteoporosis_probability$Osteoporosis_rate_high <- Osteoporosis_probability$Osteoporosis_rate + (1.96*(Osteoporosis_probability$Osteoporosis_rate/10))
 
  
  # Log odds ratio for osteoporosis
  log_or_osteoporosis_GFD <- rnorm(n_samples, mean = log(1.43), 
    sd = ((log(1.78) - log(1.15))/(2*1.96))) #Olmos 2008
  # **ADAPTED**
  log_or_osteoporosis_noGFD <- rnorm(n_samples, mean = log(
    ((1*(asymp_prop/undiag_all)+1.77*(symp_undiag/undiag_all))+
                                    ((asymp_prop/(undiag_all-p_diag_ever))*1+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*1.77)))/2), 
    sd = (((symp_undiag/undiag_all)*((log(2.65) - log(1.18))/(2*1.96))) + 
                                        ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*((log(2.65) - log(1.18))/(2*1.96)))))/2)  #Hujoel 2018
  
  
   # Log rates in general population
  osteoporosis_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(osteoporosis_lograte) <- paste0("osteoporosis_lograte_", seq(0, 90, 10))
  
  
  for(i_age_category in 1:10) {
    osteoporosis_lograte[, i_age_category] <- rnorm(n_samples, mean = log(Osteoporosis_probability$Osteoporosis_rate[i_age_category]),
                                                    sd = (log(Osteoporosis_probability$Osteoporosis_rate_high[i_age_category]) - log(Osteoporosis_probability$Osteoporosis_rate_low[i_age_category]))/(2*1.96))
  }
  
   ##developing NHL
  
  NHL_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/NHL_rates.xlsx", sheet = "input"))
  NHL_probability$NHL_rate_low <- NHL_probability$NHL_rate - (1.96*(NHL_probability$NHL_rate/10))
  NHL_probability$NHL_rate_high <- NHL_probability$NHL_rate + (1.96*(NHL_probability$NHL_rate/10))
  
  # Log incidence ratios for NHL
  log_rr_NHL_GFD <- rnorm(n_samples, mean = log(3.28), 
    sd =  ((log(6.28) - log(1.49))/(2*1.96)))
  # **ADAPTED**
  log_rr_NHL_noGFD <- rnorm(n_samples, mean = log(
    ((1*(asymp_prop/undiag_all) + 4.7*(symp_undiag/undiag_all))+
                                      ((asymp_prop/(undiag_all-p_diag_ever))*1+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*4.7)))/2),
    sd = (((symp_undiag/undiag_all)*((log(7.3) - log(2.9))/(2*1.96)))+
                                      ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*((log(7.3) - log(2.9))/(2*1.96)))))/2)
  
  
  
  # Log rates in general population
  NHL_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(NHL_lograte) <- paste0("NHL_lograte_", seq(0, 90, 10))

  for(i_age_category in 1:10) {
    NHL_lograte[, i_age_category] <- rnorm(n_samples, mean = log(NHL_probability$NHL_rate[i_age_category]),
                                                    sd = (log(NHL_probability$NHL_rate_high[i_age_category]) - log(NHL_probability$NHL_rate_low[i_age_category]))/(2*1.96))
  }
  
  
  ##developing GIC 
  GIC_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/GIC_rates.xlsx", sheet = "input"))
  GIC_probability$GIC_rate_low <- GIC_probability$GIC_rate - (1.96*(GIC_probability$GIC_rate/10))
  GIC_probability$GIC_rate_high <- GIC_probability$GIC_rate + (1.96*(GIC_probability$GIC_rate/10))

  
  
  # Log incidence ratios for GIC from studies
  log_rr_GIC_GFD <- rnorm(n_samples, mean = log(1.32), 
    sd =  ((log(1.43) - log(1.22))/(2*1.96)))
  # **ADAPTED**
  log_rr_GIC_noGFD <- rnorm(n_samples, mean = log(
    ((1*(asymp_prop/undiag_all)+2.33*(symp_undiag/undiag_all))+
                                   ((asymp_prop/(undiag_all-p_diag_ever))*1+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*2.33)))/2),
    sd = (((symp_undiag/undiag_all)*((log(4.04) - log(1.35))/(2*1.96)))+
                                    ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*((log(4.04) - log(1.35))/(2*1.96)))))/2)
  
  
  
  
 # Log rates in general population
  GIC_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(GIC_lograte) <- paste0("GIC_lograte_", seq(0, 90, 10))
  for(i_age_category in 1:10) {
    GIC_lograte[, i_age_category] <- rnorm(n_samples, mean = log(GIC_probability$GIC_rate[i_age_category]),
                                           sd = (log(GIC_probability$GIC_rate_high[i_age_category]) - log(GIC_probability$GIC_rate_low[i_age_category]))/(2*1.96))
  }
  
  
  ##Dying from consequence
  death_log_hazard_NHL <- exp(-2.6026)  # 2.6026 based on IKNL (indolent and aggressive NHL)
  death_probability_NHL <- 1-exp(-exp(death_log_hazard_NHL))
  
  
  # Death OSTP - from Klop 2014
  death_log_hr_osteoporosis_mixedj <- rnorm(n = n_samples, mean = log(2.80),
                                            sd = (log(3.04) - log(2.58))/(2*1.96))
  
  # Death GIC - from IKNL (spreadsheet "GIC rates death")
  death_log_hazard_GIC <- exp(-0.9729)
  death_probability_GIC <- 1-exp(-exp(death_log_hazard_GIC))
  
  
 
  #############################################################################
  ## Utilities ################################################################
  #############################################################################

  
  utility_GFD_adults <-  0.84 # from NCV data  
  utility_GFdse_adults <-  (0.94-0.74)/(2*1.96) # 0.001 
  utility_GFDalpha_adults <- (utility_GFD_adults ^ 2 * (1 - utility_GFD_adults)/utility_GFdse_adults ^ 2) - utility_GFD_adults
  utility_GFDbeta_adults <- (utility_GFDalpha_adults / utility_GFD_adults) - utility_GFDalpha_adults
  utility_GFD_adults <- rbeta(n = n_samples, shape1 = utility_GFDalpha_adults, shape2 = utility_GFDbeta_adults)
  
  utility_GFD_children <- 0.9 #from NCV data
  utility_GFdse_children <- (1-0.8)/(2*1.96) # 0.01
  utility_GFDalpha_children <- (utility_GFD_children ^ 2 * (1 - utility_GFD_children)/utility_GFdse_children ^ 2) - utility_GFD_children
  utility_GFDbeta_children <- (utility_GFDalpha_children / utility_GFD_children) - utility_GFDalpha_children
  utility_GFD_children <- rbeta(n = n_samples, shape1 = utility_GFDalpha_children, shape2 = utility_GFDbeta_children)
  
  # **ADAPTED**
  utility_undiagnosedCD_adults <- (((asymp_prop/undiag_all)*0.866+(symp_undiag/undiag_all)*0.65)+
                                                                ((asymp_prop/(undiag_all-p_diag_ever))*0.866+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*0.65)))/2 # 0.65 is from NCV data
  utility_undiagnosedCD_se_adults <- (((asymp_prop/undiag_all)*0+(symp_undiag/undiag_all)*((0.75-0.55)/(2*1.96)))+
                                                                 ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*((0.75-0.55)/(2*1.96)))))/2 
  utility_undiagnosedCD_alpha_adults <- (utility_undiagnosedCD_adults ^ 2 * (1 - utility_undiagnosedCD_adults)/utility_undiagnosedCD_se_adults ^ 2) - utility_undiagnosedCD_adults
  utility_undiagnosedCD_beta_adults <- (utility_undiagnosedCD_alpha_adults/utility_undiagnosedCD_adults) - utility_undiagnosedCD_alpha_adults
  utility_undiagnosedCD_adults <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_adults, shape2 = utility_undiagnosedCD_beta_adults)
  
  # **ADAPTED**
  utility_undiagnosedCD_children <- (((asymp_prop/undiag_all)*0.958+(symp_undiag/undiag_all)*0.59)+
                                                                ((asymp_prop/(undiag_all-p_diag_ever))*0.958+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*0.59)))/2 # 0.59 is from NCV data
  utility_undiagnosedCD_se_children <- (((asymp_prop/undiag_all)*0+(symp_undiag/undiag_all)*((0.69-0.49)/(2*1.96))) +
                                                                   ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*((0.69-0.49)/(2*1.96)))))/2 
  utility_undiagnosedCD_alpha_children <- (utility_undiagnosedCD_children ^ 2 * (1 - utility_undiagnosedCD_children)/utility_undiagnosedCD_se_children ^ 2) - utility_undiagnosedCD_children
  utility_undiagnosedCD_beta_children <- (utility_undiagnosedCD_alpha_children/utility_undiagnosedCD_children) - utility_undiagnosedCD_alpha_children
  utility_undiagnosedCD_children <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_children, shape2 = utility_undiagnosedCD_beta_children)
  
  
  #####
  

  probability_hipfracture <- 151/100000 # in NL, from Lotters et al 2016
  probability_vertebralfracture <- 53/100000 # in NL, from Lotters et al 2016
  probability_wristfracture <- 146/100000 # in NL, from Lotters et al 2016
  disutility_hipfracture <-   0.839 - 0.59  #0.839 corresponds to Dutch eq5d5l norm from Versteegh ages 60-70
  disutility_hipfractureSE <- (((0.839 - 0.65) - (0.839 - 0.54))/3.92)
  disutility_hipfracture_alpha <- (disutility_hipfracture ^ 2 * (1 - disutility_hipfracture)/disutility_hipfractureSE ^ 2) - disutility_hipfracture
  disutility_hipfracture_beta <- (disutility_hipfracture_alpha/disutility_hipfracture) - disutility_hipfracture_alpha
  disutility_hipfracture <- rbeta(n = n_samples, shape1 = disutility_hipfracture_alpha, shape2 = disutility_hipfracture_beta)
  disutility_wristfracture <-   0.839 - 0.78
  disutility_wristfractureSE <-  (((0.839 - 0.84) - (0.839 - 0.72))/3.92)
  disutility_wristfracture_alpha <- (disutility_wristfracture ^ 2 * (1 - disutility_wristfracture)/disutility_wristfractureSE ^ 2) - disutility_wristfracture
  disutility_wristfracture_beta <- (disutility_wristfracture_alpha/disutility_wristfracture) - disutility_wristfracture_alpha
  disutility_wristfracture <- rbeta(n = n_samples, shape1 = disutility_wristfracture_alpha, shape2 = disutility_wristfracture_beta)
  disutility_vertebralfracture <-   0.839 - 0.55
  disutility_vertebralfractureSE <- (((0.839 - 0.60) - (0.839 - 0.50))/3.92)
  disutility_vertebralfracture_alpha <- (disutility_vertebralfracture ^ 2 * (1 - disutility_vertebralfracture)/disutility_vertebralfractureSE ^ 2) - disutility_vertebralfracture
  disutility_vertebralfracture_beta <- (disutility_vertebralfracture_alpha/disutility_vertebralfracture) - disutility_vertebralfracture_alpha
  disutility_vertebralfracture <- rbeta(n = n_samples, shape1 = disutility_vertebralfracture_alpha, shape2 = disutility_vertebralfracture_beta)
  disutility_osteoporosis <- (probability_hipfracture * disutility_hipfracture) + (probability_wristfracture * disutility_wristfracture) + (probability_vertebralfracture * disutility_vertebralfracture)

  disutility_NHL <-   0.839 - 0.735 # Verstegh 2012 https://ars.els-cdn.com/content/image/1-s2.0-S1098301520320568-mmc1.pdf
  disutility_NHLSE <- (((0.839 - 0.784) - (0.839 - 0.686))/3.92)
  disutility_NHL_alpha <- (disutility_NHL ^ 2 * (1 - disutility_NHL)/disutility_NHLSE ^ 2) - disutility_NHL
  disutility_NHL_beta <- (disutility_NHL_alpha/disutility_NHL) - disutility_NHL_alpha
  disutility_NHL <- rbeta(n = n_samples, shape1 = disutility_NHL_alpha, shape2 = disutility_NHL_beta)
  
  
  disutility_GIC <- runif(n = n_samples, min =  0.839-0.74, max = 0.839-0.68) # min from Djalalov et al 2014 min = stage 1-3, max=stage 4 
  
  
  disutility_biopsy_adults <- rtri(n = n_samples, min = 0, max = 0.005, mode = 0.003)
  disutility_biopsy_children <- rtri(n = n_samples, min = 0, max = 0.010, mode = 0.006)
  
  
  
  #############################################################################
  ## Costs ####################################################################
  #############################################################################
  
  
  cost_osteoporosis <- (110700000/(29.2*17340))*1.0398 # 110,700,000 euros in 2019, 506328 (i.e. 29.2 per 1000 * 17340) people with OSTP in 2019 ((https://WWW.vzinfo.nl/osteoporose/zorguitgaven))
  cost_osteoporosis_se <- cost_osteoporosis*0.1
  cost_osteoporosis_alpha <- (cost_osteoporosis/cost_osteoporosis_se)^2
  cost_osteoporosis_beta <- (cost_osteoporosis_se^2)/cost_osteoporosis
  cost_osteoporosis <- rgamma(n=n_samples, shape = cost_osteoporosis_alpha, scale = cost_osteoporosis_beta)
  
  
  
  cost_IDA <- (52*(2*0.02))+(4*6) # 0.02 per tab, 2 times per week, 4 times the afleverkosten which = 6 eur
  
  cost_gfp <- 1506.02  # from NCV data, assuming year 2021
  cost_gfp_se <- 16.93  # from NCV data
  cost_gfp_alpha <- (cost_gfp/cost_gfp_se)^2
  cost_gfp_beta <- (cost_gfp_se^2)/cost_gfp
  cost_gfp <- rgamma(n = n_samples, shape = cost_gfp_alpha, scale = cost_gfp_beta) 
  
  ##
  
  ferritine <- 		6.48 #all below are 2021 euros
  ijzer <-		2.35
  transferrine <-		4.41
  tpo_antistoffen <-		31.71
  vit_b12 <-		6.45
  foliumzuur <-		5.86
  endo_iga <- 		34.96
  iga_anti_ttg <-		34.96
  hla <- 		263.44
  order_tarif <-  12.13
  biopsy <- 376.65
  poli_pediatric <- 112
  dietician <- 40
  phone_doc <- 19
  
  test_battery <- ferritine + ijzer + transferrine +tpo_antistoffen +  vit_b12 +foliumzuur +endo_iga +iga_anti_ttg +hla + order_tarif 
  cost_diagnosis_s <- test_battery + (biopsy*(5/61)) + poli_pediatric + dietician + phone_doc # 5/61 = prop receiving biopsy in POCT
  
  poli_general <- 101 #for use later, from costing manual, gewogen gemidelde
  
  probability_biopsy_adults <- rbeta(n=n_samples, shape1 = 2007, shape2 = (2324 - 2007)) # 2324=total with no errors and +18 in NCV, 2007=n who received biopsy
  probability_biopsy_children <- rbeta(n=n_samples, shape1 = 151, shape2 = (480 - 151)) # 480=total in NSCK, 151=n who received biopsy
  
  # **ADAPTED**
  cost_diag_asymp <-  test_battery + poli_general
  
  # **ADAPTED**
  cost_diagnosis_children_ns_tmp <- test_battery + poli_pediatric + (probability_biopsy_children*biopsy) #costs of clinical diagnosis if delayed diagnosis
  cost_diagnosis_adults_ns_tmp <- test_battery + poli_general + (probability_biopsy_adults*biopsy)  #costs of clinical diagnosis if delayed diagnosis
  
  # **ADAPTED**
  cost_diagnosis_children_ns <- (((((asymp_prop/undiag_all)*(cost_diag_asymp+biopsy*mean(probability_biopsy_children))) +((symp_undiag/undiag_all)*cost_diagnosis_children_ns_tmp)))+
                                   ((asymp_prop/(undiag_all-p_diag_ever))*(cost_diag_asymp+biopsy*mean(probability_biopsy_children))+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*cost_diagnosis_children_ns_tmp)))/2 
  
  # **ADAPTED**
  cost_diagnosis_adults_ns <- (((((asymp_prop/undiag_all)*(cost_diag_asymp+biopsy*mean(probability_biopsy_adults))) +((symp_undiag/undiag_all)*cost_diagnosis_adults_ns_tmp)))+
                                 ((asymp_prop/(undiag_all-p_diag_ever))*(cost_diag_asymp+biopsy*mean(probability_biopsy_adults))+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*cost_diagnosis_adults_ns_tmp)))/2  
  
  # **ADAPTED**
  cost_undiagnosedCD_adults_tmp <- ((649.96 - poli_general)/10.72)  #649 = pre_costs from NCV. 10.72 = mean duration of symptoms in NCV. 
  cost_undiagnosedCD_adults <- (((((asymp_prop/undiag_all)*0) +((symp_undiag/undiag_all)*cost_undiagnosedCD_adults_tmp)))+
                                        ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*cost_undiagnosedCD_adults_tmp)))/2 
  cost_undiagnosedCD_se_adults <- cost_undiagnosedCD_adults*0.1 # 1/10 of mean
  cost_undiagnosedCD_alpha_adults <- (cost_undiagnosedCD_adults/cost_undiagnosedCD_se_adults)^2
  cost_undiagnosedCD_beta_adults <- (cost_undiagnosedCD_se_adults^2)/cost_undiagnosedCD_adults
  cost_undiagnosedCD_adults <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_adults, scale = cost_undiagnosedCD_beta_adults)
  
  # **ADAPTED**
  cost_undiagnosedCD_children_tmp <- ((680.8 - cost_diagnosis_children_ns)/1.96) # 680.8 = costs from NSCK. 1.96 = mean duration of symptoms in NCV. 
  cost_undiagnosedCD_children <- (((((asymp_prop/undiag_all)*0) +((symp_undiag/undiag_all)*cost_undiagnosedCD_children_tmp)))+
                                          ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*cost_undiagnosedCD_children_tmp)))/2 
  cost_undiagnosedCD_se_children <- cost_undiagnosedCD_children*0.1 # 1/10 of mean
  cost_undiagnosedCD_alpha_children <- (cost_undiagnosedCD_children/cost_undiagnosedCD_se_children)^2
  cost_undiagnosedCD_beta_children <- (cost_undiagnosedCD_se_children^2)/cost_undiagnosedCD_children
  cost_undiagnosedCD_children <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_children, scale = cost_undiagnosedCD_beta_children)
  
  
  cost_CDGFD_adults <- 172.85  # 172.85 is from NCV data (condition >= 1 year diagnosed, 2021 eur).
  cost_CDGFD_se_adults <- cost_CDGFD_adults*0.1 # 1/10 of mean
  cost_CDGFD_alpha_adults <- (cost_CDGFD_adults/cost_CDGFD_se_adults)^2
  cost_CDGFD_beta_adults <- (cost_CDGFD_se_adults^2)/cost_CDGFD_adults
  cost_CDGFD_adults <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_adults, scale = cost_CDGFD_beta_adults)
  
  cost_CDGFD_children <- 205.65  # 205.65 from NCV data (condition >= 1 year diagnosed, 2021 eur) 
  cost_CDGFD_se_children <- cost_CDGFD_children*0.1 # 1/10 of mean
  cost_CDGFD_alpha_children <- (cost_CDGFD_children/cost_CDGFD_se_children)^2
  cost_CDGFD_beta_children <- (cost_CDGFD_se_children^2)/cost_CDGFD_children
  cost_CDGFD_children <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_children, scale = cost_CDGFD_beta_children)

  
  #######################################################################################
  
  
  
  cost_NHL <- ((207000000)/(13248 +	15131))*1.0398 # updated based on 20-year prevalence in 2019 and KvZ total in 2019 
  cost_NHL_se <- cost_NHL*0.1
  cost_NHL_alpha <- (cost_NHL/cost_NHL_se)^2
  cost_NHL_beta <- (cost_NHL_se^2)/cost_NHL
  cost_NHL <- rgamma(n=n_samples, shape = cost_NHL_alpha, scale = cost_NHL_beta)
  
  
  cost_GIC <- (558000000/108477)*1.0398 #based on 20-year prevalence for colon cancer in 2019 and total KvZ for colon cancer in 2019. Links: https://www.vzinfo.nl/dikkedarmkanker/zorguitgaven & https://iknl.nl/nkr-cijfers?fs%7Cepidemiologie_id=528&fs%7Ctumor_id=138&fs%7Cprevalentie_id=554&fs%7Cperiode_id=563&fs%7Cgeslacht_id=644&fs%7Cleeftijdsgroep_id=677&fs%7Cjaren_na_diagnose_id=687&cs%7Ctype=false&cs%7CxAxis=false&cs%7Cseries=epidemiologie_id&ts%7CrowDimensions=&ts%7CcolumnDimensions=&lang%7Clanguage=nl
  cost_GIC_se <- cost_GIC*0.1
  cost_GIC_alpha <- (cost_GIC/cost_GIC_se)^2
  cost_GIC_beta <- (cost_GIC_se^2)/cost_GIC
  cost_GIC <- rgamma(n=n_samples, shape = cost_GIC_alpha, scale = cost_GIC_beta)
  

  
  
  cost_POCT <- 11.51 + 13.87 # the 13.87 is for the POCT itself excl 6% VAT, the rest is for related costs. 2021 prices
  
  
  cost_quest <- 2 # cost of questionnaire per child
  
  
  #productivity losses (friction cost method)
  # **ADAPTED**
  cost_abs_undiagnosed_fc <-  ((((asymp_prop/undiag_all)*0+(symp_undiag/undiag_all)*1674.54))+
                                      ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))*1674.54)))/2  #WAS: 1674.54 #2021 eur 
  cost_abs_undiagnosed_fc_se <- ((((asymp_prop/undiag_all)*0+(symp_undiag/undiag_all)*76.83))+
                                         ((asymp_prop/(undiag_all-p_diag_ever))*0+(((symp_undiag-p_diag_ever)/(undiag_all-p_diag_ever))* 76.83)))/2  
  cost_abs_undiagnosed_fc_alpha <- (cost_abs_undiagnosed_fc/cost_abs_undiagnosed_fc_se)^2
  cost_abs_undiagnosed_fc_beta <- (cost_abs_undiagnosed_fc_se^2)/cost_abs_undiagnosed_fc
  cost_abs_undiagnosed_fc <- rgamma(n=n_samples, shape = cost_abs_undiagnosed_fc_alpha, scale = cost_abs_undiagnosed_fc_beta)
  


cost_abs_diagnosed_fc <- 332.12 #2021 eur . NCV data
cost_abs_diagnosed_fc_se <- 33.09
cost_abs_diagnosed_fc_alpha <- (cost_abs_diagnosed_fc/cost_abs_diagnosed_fc_se)^2
cost_abs_diagnosed_fc_beta <- (cost_abs_diagnosed_fc_se^2)/cost_abs_diagnosed_fc
cost_abs_diagnosed_fc <- rgamma(n=n_samples, shape = cost_abs_diagnosed_fc_alpha, scale = cost_abs_diagnosed_fc_beta)



  
  #######################
  input_parameters <- data.frame(probability_late_diagnosis_children, probability_late_diagnosis_adults, 
                                 probability_osteoporosis, probability_NHL, probability_GIC, probability_nocomplications, probability_IDA,
                                 osteoporosis_lograte, log_or_osteoporosis_GFD, log_or_osteoporosis_noGFD,
                                 NHL_lograte, log_rr_NHL_GFD, log_rr_NHL_noGFD,
                                 GIC_lograte, log_rr_GIC_GFD, log_rr_GIC_noGFD,
                                 death_probability_NHL, death_log_hr_osteoporosis_mixedj, death_probability_GIC, 
                                 utility_GFD_children, utility_GFD_adults, utility_undiagnosedCD_children, utility_undiagnosedCD_adults,  disutility_osteoporosis, disutility_NHL, disutility_GIC, #disutility_fp,
                                 probability_biopsy_children, probability_biopsy_adults, disutility_biopsy_adults, disutility_biopsy_children, 
                                 cost_CDGFD_children, cost_CDGFD_adults, cost_osteoporosis, cost_undiagnosedCD_children, cost_undiagnosedCD_adults, cost_IDA, cost_GIC, cost_NHL,
                                 cost_diagnosis_children_ns, cost_diagnosis_adults_ns, cost_diagnosis_s, #ns= no screening, s=screening  
                                 cost_abs_undiagnosed_fc, cost_abs_diagnosed_fc,
                                 sens_POCT, spec_POCT, cost_POCT, cost_gfp, cost_quest)
  
  return(input_parameters)
} #END FUNCTION

