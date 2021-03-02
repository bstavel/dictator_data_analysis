
load_interaction_results_and_combine <- function(epoch, freq, pred) {
  
  ### any subject general vars ###
  spec_vars <- c("self_var_payoff", "more_or_less")
  regions_to_combine <- c("All")
  
  ## IR57 ##
  choice_inter_reg_df_57 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR57",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR57", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_57 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR57",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR57", "post_eda", "interactions"))
  
  
  ## IR35 ##
  choice_inter_reg_df_35 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR35",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR35", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_35 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR35",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR35", "post_eda", "interactions"))
  
  
  ## IR28 ##
  choice_inter_reg_df_28 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR28",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR28", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_28 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR28",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR28", "post_eda", "interactions"))
  
  ## IR26 ##
  choice_inter_reg_df_26 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR26",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR26", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_26 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR26",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR26", "post_eda", "interactions"))
  
  
  ## IR16 ##
  choice_inter_reg_df_16 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR16",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR16", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_16 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR16",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR16", "post_eda", "interactions"))
  
  
  ## IR10 ##
  choice_inter_reg_df_10 <- compile_mult_reg_results(regions = regions_to_combine,
                                                     type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                     sub = "IR10",
                                                     spec_vars = spec_vars,
                                                     path = path(here(), "results", "IR10", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_10 <- compile_mult_reg_results(regions = regions_to_combine,
                                                           type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                           sub = "IR10",
                                                           spec_vars = spec_vars,
                                                           path = path(here(), "results", "IR10", "post_eda", "interactions"))
  
  
  ## IR9 ##
  choice_inter_reg_df_9 <- compile_mult_reg_results(regions = regions_to_combine,
                                                    type = paste0("interaction-", freq, "-choice-locked-hilbertRS"),
                                                    sub = "IR9",
                                                    spec_vars = spec_vars,
                                                    path = path(here(), "results", "IR9", "post_eda", "interactions"))
  
  
  presentation_inter_reg_df_9 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = paste0("interaction-", freq, "-pres-locked-hilbertRS"),
                                                          sub = "IR9",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR9", "post_eda", "interactions"))
  
  
  
  
  # presentation dfs #
  presentation_inter_reg_df_57 <- presentation_inter_reg_df_57 %>% mutate(sub = "IR57") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_35 <- presentation_inter_reg_df_35 %>% mutate(sub = "IR35") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_28 <- presentation_inter_reg_df_28 %>% mutate(sub = "IR28") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_26 <- presentation_inter_reg_df_26 %>% mutate(sub = "IR26") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_16 <- presentation_inter_reg_df_16 %>% mutate(sub = "IR16") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_10 <- presentation_inter_reg_df_10 %>% mutate(sub = "IR10") %>% mutate(epoch = "presentation")
  presentation_inter_reg_df_9 <- presentation_inter_reg_df_9 %>% mutate(sub = "IR9") %>% mutate(epoch = "presentation")
  
  
  # choice dfs #
  choice_inter_reg_df_57 <- choice_inter_reg_df_57 %>% mutate(sub = "IR57") %>% mutate(epoch = "choice")
  choice_inter_reg_df_35 <- choice_inter_reg_df_35 %>% mutate(sub = "IR35") %>% mutate(epoch = "choice")
  choice_inter_reg_df_28 <- choice_inter_reg_df_28 %>% mutate(sub = "IR28") %>% mutate(epoch = "choice")
  choice_inter_reg_df_26 <- choice_inter_reg_df_26 %>% mutate(sub = "IR26") %>% mutate(epoch = "choice")
  choice_inter_reg_df_16 <- choice_inter_reg_df_16 %>% mutate(sub = "IR16") %>% mutate(epoch = "choice")
  choice_inter_reg_df_10 <- choice_inter_reg_df_10 %>% mutate(sub = "IR10") %>% mutate(epoch = "choice")
  choice_inter_reg_df_9 <- choice_inter_reg_df_9 %>% mutate(sub = "IR9") %>% mutate(epoch = "choice")
  
  
  
  
  
  merge_with_elec_info <- function(df, localizations_sub) {
    ## merge df with elec info, stored in localizations_sub
    ## creates new var of form ROI1-ROI2  
    ## also cleans elec names
    
    # clean the names
    df <- df %>%
      mutate(electrode = gsub("POL ", "", electrode)) %>%
      mutate(electrode = gsub(" POL", "", electrode)) %>%
      mutate(electrode = gsub("-Ref", "", electrode)) %>%
      mutate(electrode = gsub("-Ref-", "-", electrode)) 
    
    # separate two electrodes
    df <- df %>%
      mutate(first_elec= gsub("-.*", "", electrode)) %>%
      mutate(second_elec = gsub(".*-", "", electrode))
    # join on first electrode
    df <- left_join(df,
                    localizations_sub,
                    by = c("first_elec" = "Electrode"))
    # rename to `first_region`
    df <- df %>%
      rename(first_region = `Loc Meeting`)
    # join on second electrode and rename
    df <-  left_join(df,
                     localizations_sub,
                     by = c("second_elec" = "Electrode")) %>%
      rename(second_region = `Loc Meeting`) 
    # paste into useful format
    df <- df %>% mutate(ROI = paste0(first_region, "-", second_region))
    
    return(df)
    
  }
  
  
  ## get region names ##
  localizations <- read_csv(path(here(), "munge", "combined_electrode_info_with_bob.csv"))
  localizations_28 <- localizations %>% filter(subject == "IR28") %>% select(Electrode, `Loc Meeting`)
  localizations_35 <- localizations %>% filter(subject == "IR35") %>% select(Electrode, `Loc Meeting`)
  localizations_57 <- localizations %>% filter(subject == "IR57") %>% select(Electrode, `Loc Meeting`)
  localizations_9 <- localizations %>% filter(subject == "IR9") %>% select(Electrode, `Loc Meeting`)
  localizations_10 <- localizations %>% filter(subject == "IR10") %>% select(Electrode, `Loc Meeting`)
  localizations_10 <- localizations_10 %>%
    mutate(Electrode = gsub("POL ", "", Electrode)) %>%
    mutate(Electrode = gsub(" POL", "", Electrode)) %>%
    mutate(Electrode = gsub("-Ref", "", Electrode)) %>%
    mutate(Electrode = gsub("-Ref-", "-", Electrode)) 
  localizations_26 <- localizations %>% filter(subject == "IR26") %>% select(Electrode, `Loc Meeting`)
  localizations_16 <- localizations %>% filter(subject == "IR16") %>% select(Electrode, `Loc Meeting`)
  localizations_cp34 <- localizations %>% filter(subject == "CP34") %>% select(Electrode, `Loc Meeting`)
  
  
  ## merge dfs ##
  presentation_inter_reg_df_9 <- merge_with_elec_info(presentation_inter_reg_df_9, localizations_9)
  presentation_inter_reg_df_10 <- merge_with_elec_info(presentation_inter_reg_df_10, localizations_10)
  presentation_inter_reg_df_16 <- merge_with_elec_info(presentation_inter_reg_df_16, localizations_16)
  presentation_inter_reg_df_26 <- merge_with_elec_info(presentation_inter_reg_df_26, localizations_26)
  presentation_inter_reg_df_28 <- merge_with_elec_info(presentation_inter_reg_df_28, localizations_28)
  presentation_inter_reg_df_35 <- merge_with_elec_info(presentation_inter_reg_df_35, localizations_35)
  presentation_inter_reg_df_57 <- merge_with_elec_info(presentation_inter_reg_df_57, localizations_57)
  
  choice_inter_reg_df_9 <- merge_with_elec_info(choice_inter_reg_df_9, localizations_9)
  choice_inter_reg_df_10 <- merge_with_elec_info(choice_inter_reg_df_10, localizations_10)
  choice_inter_reg_df_16 <- merge_with_elec_info(choice_inter_reg_df_16, localizations_16)
  choice_inter_reg_df_26 <- merge_with_elec_info(choice_inter_reg_df_26, localizations_26)
  choice_inter_reg_df_28 <- merge_with_elec_info(choice_inter_reg_df_28, localizations_28)
  choice_inter_reg_df_35 <- merge_with_elec_info(choice_inter_reg_df_35, localizations_35)
  choice_inter_reg_df_57 <- merge_with_elec_info(choice_inter_reg_df_57, localizations_57)
  
  
  
  
  all_pres_inter_reg_results <- rbind(presentation_inter_reg_df_9, 
                                      presentation_inter_reg_df_10,
                                      presentation_inter_reg_df_16,
                                      presentation_inter_reg_df_26,
                                      presentation_inter_reg_df_28,
                                      presentation_inter_reg_df_35,
                                      presentation_inter_reg_df_57)
  
  all_choice_inter_reg_results <- rbind(choice_inter_reg_df_9, 
                                        choice_inter_reg_df_10,
                                        choice_inter_reg_df_16,
                                        choice_inter_reg_df_26,
                                        choice_inter_reg_df_28,
                                        choice_inter_reg_df_35,
                                        choice_inter_reg_df_57)
  
  # only main interaction
  all_pres_inter_reg_results <- all_pres_inter_reg_results %>% filter(predictor  == pred)
  all_choice_inter_reg_results <- all_choice_inter_reg_results %>% filter(predictor == pred)
  
  if(epoch == "presentation"){
    return(all_pres_inter_reg_results)
  } else {
    return(all_choice_inter_reg_results)
  }
  
}