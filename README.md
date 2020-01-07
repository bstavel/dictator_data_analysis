# Dictator Behavioral Data Analysis

These scripts process the behavioral data for a modified Dictator Task in ECoG patients to investiagte neural implementations of decision making during presentations of advantageous and inadvenategous inequity.

This document is designed to help us document tough choices that we made throughout the analysis stream.

## Participant Exclusion

How much missingness is okay for us to include subjects? 10%, 15%, 20%?

Of the 60 sanity check trials, when would we exclude a subject for selecting too many of the worse for everyone trials? 

## Counterbalancing

We identified an issue with counterbalancing (or lack thereof) in the intial set of subjects. See the `behavior_brs.Rmd` file for details. To cope with this we made two decisions:

* collect some additional datasets with a properly counterbalanced trial file (created in `behavior_brs.Rmd`)
* Run all our regression models with a vector specifying `side of choice` and an interaction between `side of choice` and `side value`.

## Regions

Will we run eletrodes that are partly in the wm in their notes? What percentage cut-off should we use? Some subjects look like they were run on an old system and don't have the full brain atlas information, can we rerun them?

Here is the first pass numbers for ofc and insula across subjects:

*OFC: ~ 42 electrodes to test*

| subject 	| loc_meeting 	| AAL 	| AFNI 	| BrainWeb 	| Brainnetome 	| JuBrain 	| VTPM 	|
|---------	|-------------	|-----	|------	|----------	|-------------	|---------	|------	|
| IR9     	| 0           	| 10  	| 13   	| 0        	| 9           	| 3       	| 0    	|
| IR57    	| 14          	| 14  	| 17   	| 0        	| 16          	| 5       	| 0    	|
| IR39    	| 6           	| 4   	| 5    	| 0        	| 3           	| 1       	| 0    	|
| IR35    	| 7           	| 7   	| 10   	| 0        	| 8           	| 4       	| 0    	|
| IR28    	| 0           	| 0   	| 0    	| 0        	| 1           	| 0       	| 0    	|
| IR26    	| 0           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR19    	| 11          	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR16    	| 5           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR10    	| 3           	| 5   	| 1    	| 0        	| 3           	| 0       	| 0    	|
| ST40    	| 4           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|


*insula: ~64 electrodes to test*

| subject 	| loc_meeting 	| AAL 	| AFNI 	| BrainWeb 	| Brainnetome 	| JuBrain 	| VTPM 	|
|---------	|-------------	|-----	|------	|----------	|-------------	|---------	|------	|
| IR9     	| 0           	| 10  	| 13   	| 0        	| 9           	| 3       	| 0    	|
| IR57    	| 14          	| 14  	| 17   	| 0        	| 16          	| 5       	| 0    	|
| IR39    	| 6           	| 4   	| 5    	| 0        	| 3           	| 1       	| 0    	|
| IR35    	| 7           	| 7   	| 10   	| 0        	| 8           	| 4       	| 0    	|
| IR28    	| 0           	| 0   	| 0    	| 0        	| 1           	| 0       	| 0    	|
| IR26    	| 0           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR19    	| 11          	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR16    	| 5           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|
| IR10    	| 3           	| 5   	| 1    	| 0        	| 3           	| 0       	| 0    	|
| ST40    	| 4           	| 0   	| 0    	| 0        	| 0           	| 0       	| 0    	|


## Referencing

I thought we did pairwise referencing, not CAR? Ignacio's scripts seem to be using CAR.

## High Gamma

### Morelt Wavelets or the filter-Hilbert method?

I&Z used the filter + the hilbert transform for their high frequency analysis. This is the the absolute value of the hilbert transform which is called the analytical amplitude. This is the value that is averaged in the moving time window.

### Filtering

Z&I look like they bandpassed at 70-200Hz/250Hz.

## Time binning (for regressions and high gamma)

Z & I were looking at the following bin sizes: `200, 160, 120, 80, 40` They all had overlaps of 25% regardless of size. It also looks like they were epoching based on trial start, presentation time, and choice time.

## Regressions-- Modeling Choices

| Possible Predictors | Priority? |
| ------------------- | :-------: |
| self_choice_payoff | |
| other_choice_payoff | |
| self_foregone | |
| other_foregone | |
| self_var_payoff | |
| other_var_payoff | |
| ineq_var_abs | |
| ineq_choice_abs | |
| ineq_foregone_abs | |
| ineq_var_disadvant | |
| ineq_var_advant | |
| self_diff_ch_foregone | |
| other_diff_ch_foregone | |
| side_chosen | |
| RT | |

## Regressions-- Corrections, Significance Thresholds, CV, etc
