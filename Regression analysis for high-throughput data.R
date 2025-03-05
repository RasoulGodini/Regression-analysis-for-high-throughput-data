library(truncnorm)
library(car)
library(tidyverse)

#Random data generation
set.seed(123)  # For reproducibility

n_participants <- 400
n_genes <- 500

sample_id <- sample(paste0("Sample_", sprintf("%03d", 1:n_participants)), n_participants, replace = FALSE)
sex <- sample(c("Male", "Female"), n_participants, replace = TRUE)
Age <- sample(c(40:60), n_participants, replace = TRUE)
Weight <- sample(c(50:90), n_participants, replace = TRUE)
Ethnicity <- sample(c("A", "B", "C", "D", "E"), n_participants, replace = TRUE)
blood_sugar <- sample(c(5:12), n_participants, replace = TRUE)


# Initialize matrix for gene expression
gene_data <- matrix(nrow = n_participants, ncol = n_genes)

# Generate gene expression values based on Sex
for (i in 1:n_participants) {
  if (sex[i] == "Male") {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 60, sd = 30)
  } else {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 40, sd = 20)
  }
}

# Convert to dataframe and add Sex column
colnames(gene_data) <- paste0("Gene_", 1:n_genes)  # Name genes
df <- data.frame(Sampe_ID = sample_id, Sex = sex, Age = Age, 
                 Weight = Weight, Ethnicity = Ethnicity, 
                 Blood_sugar = blood_sugar, gene_data)  # Make the datframe

genes_log10 <- log10(df[,7:ncol(df)])
df_log10 <- cbind(df[1:6], genes_log10)



# Unadjusted modelling
var_columns <- c(colnames(df_log10[7:ncol(df_log10)]))

lm_model_unadjusted <- function(var_columns){
  formula <- as.formula(paste0(var_columns," ~ Age"))
  model <- lm(formula, data = df_log10)
  model_summary <- summary(model)
  
  #Obtain the statistical values
  coeffs <- (model_summary$coefficients)[2, , drop = FALSE]
  Condition <- rownames(coeffs)
  estimates <- coeffs[, "Estimate"]
  p_values <- coeffs[, "Pr(>|t|)"]
  conf_int <- (confint(model))[2, , drop = FALSE]
  shapiro.test_pvalue <- shapiro.test(resid(model))$p.value
  
  #Prepare the results as column
  result <- data.frame(
    Comparison = Condition,
    Estimate = estimates,
    PValue = p_values,
    Lower_CI = conf_int[, 1],
    Upper_CI = conf_int[, 2],
    percent_change = (10^estimates - 1) * 100, #Convert the estimate to percent change
    Lower_CI_perc_change = (10^conf_int[, 1] - 1) * 100, #Calculate the CIs
    Upper_CI_perc_change = (10^conf_int[, 2] - 1) * 100, #Calculate the CIs
    shapiro.test_pvalue = shapiro.test_pvalue
  )
  
  return(result)
}
  
results_unadjusted <- sapply(var_columns, lm_model_unadjusted, simplify = FALSE)
result_unadjusted_2 <- data.frame(var_name = names(results_unadjusted), 
                                  do.call(rbind, results_unadjusted, ))

# Calculate FDR and add to the regression results
fdr=data.frame(adjusted_pvalue_BH = p.adjust(result_unadjusted_2$PValue, method = "BH")) 
final_result_unadjusted <- result_unadjusted_2 %>% mutate(FDR = fdr$adjusted_pvalue_BH, .after = PValue)

#write.csv(final_result, "final_result_unadjusted.csv")




# Adjusted regression model
# Make the ethnicity a factor as a dummy variable
df_log10$Ethnicity <- factor(df_log10$Ethnicity)


lm_model_adjusted <- function(var_columns){
  formula <- as.formula(paste0(var_columns," ~ Age + Weight + Sex + 
                               Ethnicity + Blood_sugar"))
  model <- lm(formula, data = df_log10)
  model_summary <- summary(model)
  
  
  #Obtain the statistical values
  coeffs <- (model_summary$coefficients)[2, , drop = FALSE]
  Condition <- rownames(coeffs)
  estimates <- coeffs[, "Estimate"]
  p_values <- coeffs[, "Pr(>|t|)"]
  conf_int <- (confint(model))[2, , drop = FALSE]
  shapiro.test_pvalue <- shapiro.test(resid(model))$p.value
  
  # Calculate Multicollinearity
  df_vif <- as.data.frame(vif(model))
  df_vif$names <- rownames(df_vif)
  names(df_vif)[names(df_vif) == "vif(model)"] <- "vif(model)"
  names(df_vif)[names(df_vif) == "GVIF"] <- "vif(model)"
  df_vif_2 <- pivot_wider(df_vif[, c("names","vif(model)")], 
                          names_from = "names", values_from = "vif(model)")
  
  #Prepare the results as column
  result <- data.frame(
    Comparison = Condition,
    Estimate = estimates,
    PValue = p_values,
    Lower_CI = conf_int[, 1],
    Upper_CI = conf_int[, 2],
    percent_change = (10^estimates - 1) * 100, #Convert the estimate to percent change
    Lower_CI_perc_change = (10^conf_int[, 1] - 1) * 100, #Calculate the CIs
    Upper_CI_perc_change = (10^conf_int[, 2] - 1) * 100, #Calculate the CIs
    shapiro.test_pvalue = shapiro.test_pvalue
  )
  
  return(cbind(result,df_vif_2))
}

results_adjusted <- sapply(var_columns, lm_model_adjusted, simplify = FALSE)
results_adjusted_2 <- data.frame(var_name = names(results_adjusted), 
                               do.call(rbind, results_adjusted, ))

# Calculate FDR and add to the regression results
fdr=data.frame(adjusted_pvalue_BH = p.adjust(results_adjusted_2$PValue, method = "BH"))
final_result_adjusted <- results_adjusted_2 %>% mutate(FDR = fdr$adjusted_pvalue_BH, .after = PValue)



