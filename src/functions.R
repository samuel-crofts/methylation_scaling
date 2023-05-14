filter_results <- function(results_dfs, cpg_list) {
  
  i <- 0
  
  for (df in results_dfs) {
    
    df <- dplyr::filter(df,cpg %in% cpg_list)
    
    if (nrow(df)!=0) {
      
      if (i==0) {
        result_df <- df
      } else {
        result_df <- rbind(result_df, df)
      }
    }
    
    i <- 1
    
  }
  
  return(result_df)
  
}

summarise_cpgs_comp_baseline <- function(df, baseline_species, comparison_species) {
  
  ratio_df <- data.frame(cpg=character(), ratio=numeric(), species=character())
  
  for (cpg_i in unique(df$cpg)) {
    
    #slope in baseline species
    baseline_est <- filter(df, cpg==cpg_i, organism==baseline_species) %>%
      dplyr::select(beta_est)
    
    #slope in comparison species
    comp_est <- filter(df, cpg==cpg_i, organism==comparison_species) %>%
      dplyr::select(beta_est)
    
    #absolute value of the ratio of the two
    ratio_i <- as.numeric(abs(mean_comp_est)/abs(mean_baseline_est))
    temp_df <- data.frame(cpg=cpg_i, ratio=ratio_i, species=comparison_species) 
    ratio_df <- rbind(ratio_df, temp_df)
    
  }
  return(ratio_df)
}

calc_slopes <- function(meth_df, meta_df, tissue_to_filt, cpgs) {
  
  results <- data.frame(cpg=character(), r2_age = numeric(), beta_est=numeric(), beta_low=numeric(), beta_high=numeric(), beta_p=numeric())
  meta_df_filt <- dplyr::filter(meta_df, tissue==tissue_to_filt)
  
  #make sure ordered correctly
  meth_df_filt <- meth_df[,meta_df_filt$sample_id]
  
  #for each cpg
  for (i in cpgs) {
    
    #if there are multiple measurements for some indivuduals, use a mixed-effects model with indivual_id as random effects. 
    #Otherwise, simple linear regression:
    if (any(names(meta_df_filt)=="individual_id")) {
      age_model <- lmer(unlist(meth_df_filt[i,]) ~ as.numeric(meta_df_filt$age) + (1|meta_df_filt$individual_id), data = NULL)  
      r2_age <- r.squaredGLMM(age_model)[1]
      sum_model <- summary(age_model)
      beta_est <- sum_model$coefficients[2,1]
      beta_p <- sum_model$coefficients[2,5]
      beta_low <- NA
      beta_high <- NA
    } else {
      age_model <- lm(unlist(meth_df_filt[i,]) ~ as.numeric(meta_df_filt$age), data = NULL)   
      sum_model <- summary(age_model)
      r2_age <- sum_model$adj.r.squared
      beta_est <- sum_model$coefficients[2,1]
      beta_p <- sum_model$coefficients[2,4]
      beta_low <- confint(age_model)[2,1]
      beta_high <- confint(age_model)[2,2]
    }
    
    #let's populate a dataframe with this info
    new_row <- data.frame(cpg=i, r2_age = r2_age, beta_est=beta_est, beta_low=beta_low, beta_high=beta_high, beta_p=beta_p)
    results <- rbind(results, new_row)
  }
  return(results)
}

make_plot_data <- function(meth_df, meta_df, cpg, organism) {
  
  meth_df <- meth_df[,meta_df$sample_id]
  
  x <- as.numeric(meta_df$age)
  y <- as.numeric(meth_df[cpg,])
  df <- data.frame(x,y)
  df$organism <- organism 
  
  return(df)
  
}
