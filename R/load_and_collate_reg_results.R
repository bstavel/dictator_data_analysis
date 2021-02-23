load_and_collate_reg_results <- function(tag, epoch, folder) {

    ## first read in just self_payoff regressions ##
    ### any subject general vars ###
    spec_vars <- c("*", "*")
    regions_to_combine <- c("All")
    
    ## IR57 ##
    choice_var_reg_df_57 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR57",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR57", folder))
    
    
    pres_var_reg_df_57 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR57",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR57", folder))
    
    
    ## IR35 ##
    choice_var_reg_df_35 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR35",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR35", folder))
    
    
    pres_var_reg_df_35 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR35",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR35", folder))
    
    
    ## IR28 ##
    choice_var_reg_df_28 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR28",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR28", folder))
    
    
    pres_var_reg_df_28 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR28",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR28", folder))
    
    ## IR26 ##
    choice_var_reg_df_26 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR26",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR26", folder))
    
    
    pres_var_reg_df_26 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR26",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR26", folder))
    
    
    ## IR16 ##
    choice_var_reg_df_16 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR16",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR16", folder))
    
    
    pres_var_reg_df_16 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR16",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR16", folder))
    
    
    ## IR10 ##
    choice_var_reg_df_10 <- compile_mult_reg_results(regions = regions_to_combine,
                                                          type = tag,
                                                          sub = "IR10",
                                                          spec_vars = spec_vars,
                                                          path = path(here(), "results", "IR10", folder))
    
    
    pres_var_reg_df_10 <- compile_mult_reg_results(regions = regions_to_combine,
                                                        type = tag,
                                                        sub = "IR10",
                                                        spec_vars = spec_vars,
                                                        path = path(here(), "results", "IR10", folder))
    
    
    ## IR9 ##
    choice_var_reg_df_9 <- compile_mult_reg_results(regions = regions_to_combine,
                                                         type = tag,
                                                         sub = "IR9",
                                                         spec_vars = spec_vars,
                                                         path = path(here(), "results", "IR9", folder))
    
    
    pres_var_reg_df_9 <- compile_mult_reg_results(regions = regions_to_combine,
                                                       type = tag,
                                                       sub = "IR9",
                                                       spec_vars = spec_vars,
                                                       path = path(here(), "results", "IR9", folder))
    
    
    # presentation dfs #
    pres_var_reg_df_57 <- pres_var_reg_df_57 %>% mutate(sub = "IR57") %>% mutate(epoch = "presentation")
    pres_var_reg_df_35 <- pres_var_reg_df_35 %>% mutate(sub = "IR35") %>% mutate(epoch = "presentation")
    pres_var_reg_df_28 <- pres_var_reg_df_28 %>% mutate(sub = "IR28") %>% mutate(epoch = "presentation")
    pres_var_reg_df_26 <- pres_var_reg_df_26 %>% mutate(sub = "IR26") %>% mutate(epoch = "presentation")
    pres_var_reg_df_16 <- pres_var_reg_df_16 %>% mutate(sub = "IR16") %>% mutate(epoch = "presentation")
    pres_var_reg_df_10 <- pres_var_reg_df_10 %>% mutate(sub = "IR10") %>% mutate(epoch = "presentation")
    pres_var_reg_df_9 <- pres_var_reg_df_9 %>% mutate(sub = "IR9") %>% mutate(epoch = "presentation")
    
    
    # choice dfs #
    choice_var_reg_df_57 <- choice_var_reg_df_57 %>% mutate(sub = "IR57") %>% mutate(epoch = "choice")
    choice_var_reg_df_35 <- choice_var_reg_df_35 %>% mutate(sub = "IR35") %>% mutate(epoch = "choice")
    choice_var_reg_df_28 <- choice_var_reg_df_28 %>% mutate(sub = "IR28") %>% mutate(epoch = "choice")
    choice_var_reg_df_26 <- choice_var_reg_df_26 %>% mutate(sub = "IR26") %>% mutate(epoch = "choice")
    choice_var_reg_df_16 <- choice_var_reg_df_16 %>% mutate(sub = "IR16") %>% mutate(epoch = "choice")
    choice_var_reg_df_10 <- choice_var_reg_df_10 %>% mutate(sub = "IR10") %>% mutate(epoch = "choice")
    choice_var_reg_df_9 <- choice_var_reg_df_9 %>% mutate(sub = "IR9") %>% mutate(epoch = "choice")
    
    # bind presentation #
    all_pres_var_reg_results <- rbind(pres_var_reg_df_9, 
                                        pres_var_reg_df_10,
                                        pres_var_reg_df_16,
                                        pres_var_reg_df_26,
                                        pres_var_reg_df_28,
                                        pres_var_reg_df_35,
                                        pres_var_reg_df_57)
    
    # clean electrode var #
    all_pres_var_reg_results <- all_pres_var_reg_results %>%
      mutate(electrode = gsub("POL ", "", electrode)) %>%
      mutate(electrode = gsub(" POL", "", electrode)) %>%
      mutate(electrode = gsub("-Ref", "", electrode)) %>%
      mutate(electrode = gsub("-Ref-", "-", electrode)) 
    
    # bind choice #
    all_choice_var_reg_results <- rbind(choice_var_reg_df_9, 
                                             choice_var_reg_df_10,
                                             choice_var_reg_df_16,
                                             choice_var_reg_df_26,
                                             choice_var_reg_df_28,
                                             choice_var_reg_df_35,
                                             choice_var_reg_df_57)
    # clean electrode var #
    all_choice_var_reg_results <- all_choice_var_reg_results %>%
      mutate(electrode = gsub("POL ", "", electrode)) %>%
      mutate(electrode = gsub(" POL", "", electrode)) %>%
      mutate(electrode = gsub("-Ref", "", electrode)) %>%
      mutate(electrode = gsub("-Ref-", "-", electrode)) 
    
    
    if(epoch == "presentation"){
      return(all_pres_var_reg_results)
    } else {
      return(all_choice_var_reg_results)
    }
    
    
    


}