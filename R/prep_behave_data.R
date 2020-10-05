prep_behave_data = function(path_to_behave_munge){
  # read in data #
  df <- read.csv(path_to_behave_munge)
  # fix some of pre calc vars #
  new_df <- df %>%
    mutate(side_chosen_numeric = side.chosen) %>%
    mutate(side_chosen = if_else(side.chosen == 76, "Left", if_else(side.chosen == 82, "Right", "Other"))) %>%
    select(-X, -side_chosen_numeric, -side.chosen, -ineq)
  # change names #
  colnames(new_df) <- gsub("\\.", "_", colnames(new_df))
  new_df <- new_df %>%
    # create derived vars based on github readme vars#
    mutate(self_var_payoff = self_payoff + self_foregone - 10) %>%
    mutate(other_var_payoff = other_payoff + other_foregone - 10) %>%
    mutate(ineq_var_abs = abs(self_var_payoff - other_var_payoff)) %>%
    mutate(ineq_choice_abs = abs(self_payoff - other_payoff)) %>%
    mutate(ineq_foregone_abs = abs(self_foregone - other_foregone)) %>%
    mutate(ineq_disadvent = (as.numeric(other_var_payoff > self_var_payoff)) * (other_var_payoff - self_var_payoff)) %>%
    mutate(ineq_advent = (as.numeric(other_var_payoff < self_var_payoff)) * (other_var_payoff - self_var_payoff)) %>%
    mutate(ineq_disadvent_choice = (as.numeric(other_payoff > self_payoff)) * (other_payoff - self_payoff)) %>%
    mutate(ineq_advent_choice = (as.numeric(other_payoff < self_payoff)) * (other_payoff - self_payoff)) %>%
    mutate(ineq_var = self_var_payoff - other_var_payoff) %>%
    mutate(ineq_choice = self_payoff - other_payoff) %>%
    mutate(chose_equality = if_else(self_payoff == 10 & other_payoff == 10, "10-10", if_else(self_payoff > other_payoff, "Advent", if_else(other_payoff > self_payoff, "Disadvent", "Equal")))) %>%
    mutate(self_diff = self_payoff - self_foregone) %>%
    mutate(other_diff = other_payoff - other_foregone) %>%
    # trial design variables
    mutate(L_self = if_else(side_chosen == "Left", self_payoff, self_foregone)) %>%
    mutate(R_self = if_else(side_chosen == "Right", self_payoff, self_foregone)) %>%
    mutate(L_other = if_else(side_chosen == "Left", other_payoff, other_foregone)) %>%
    mutate(R_other = if_else(side_chosen == "Right", other_payoff, other_foregone))
  # write csv to munge #
  write.csv(new_df, fs::path(here(), "munge", "clean_behavioral_data.csv"))
    
}