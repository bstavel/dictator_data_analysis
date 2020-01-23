# Dictator Data Analysis

These scripts process the behavioral data and the raw High Gamma data for a modified Dictator Task in ECoG patients to investiagte neural implementations of decision making during presentations of advantageous and inadvenategous inequity. These scripts include the main analyses of interest, as well as parallelized version that were run to the SCF cluster in Winter/Spring of 2020.

This document is designed to help us document tough choices that we made throughout the analysis stream.

_*An earlier version of this readme was only for the behavioral analyses, but it has since expanded_

## Participant Exclusion

How much missingness is okay for us to include subjects? 10%, 15%, 20%?

Of the 60 sanity check trials, when would we exclude a subject for selecting too many of the worse for everyone trials? 

## Counterbalancing

We identified an issue with counterbalancing (or lack thereof) in the intial set of subjects. See the `behavior_brs.Rmd` file for details. To cope with this we made two decisions:

* collect some additional datasets with a properly counterbalanced trial file (created in `behavior_brs.Rmd`)
* Run all our regression models with a vector specifying `side of choice` and an interaction between `side of choice` and `side value`.

## Regions

Will we run eletrodes that are partly in the wm in their notes?

*?*

What percentage cut-off should we use?

- While the parcellations are useful, at the end of the day we use expert judgement. This means using the loc meeting variable, where Bob & Co have determined the true location of the electrode.

Some subjects look like they were run on an old system and don't have the full brain atlas information, can we rerun them?

- These have since been run.

Should I add the cingulate?

- Yes, added the cingulate for perlim analysis with IR35.

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


The superior temporal sulcus is also very likely an important region given that includes the tempoparietal junction. Should create a table for this.


## Referencing

I thought we did pairwise referencing, not CAR? Ignacio's scripts seem to be using CAR. 

- This is likely because he was using grid electrodes.

We are using pairwise referencing, so each electrode that we have as a dependent variable, is actually made up of two electrodes.

## High Gamma

### Morlet Wavelets or the filter-Hilbert method?

I&Z used the filter + the hilbert transform for their high frequency analysis. This is the the absolute value of the hilbert transform which is called the analytical amplitude. This is the value that is averaged in the moving time window.

We have decided to use the hilbert method as well, but we also separate the filtering in a two-step process (see below). We run the hilbert transform in windows of ~12 seconds to avoid edge artifacts.

### Filtering

Z&I look like they bandpassed at 70-200Hz/250Hz. We chose to do filtering in at least two parts within the high band, run the hilbert transformation, and then recombine to have better accuracy for the higher requencies.

## Time binning (for regressions and high gamma)

Z & I were looking at the following bin sizes: `200, 160, 120, 80, 40` They all had overlaps of 25% regardless of size. It also looks like they were epoching based on trial start, presentation time, and choice time.

*We decided on two sets of analyses: 1) presentation locked continueing for 1 second and 2) choice locked with .5 before to 1 second after*

## Regressions-- Modeling Choices

### Predictors

| Possible Predictors | Definition | Priority? |
| ------------------- | :-------: | :-------: |
| self_choice_payoff | self payoff | *Yes* |
| other_choice_payoff | other payoff | *Yes* |
| self_foregone | self foregone | *Yes* |
| other_foregone | other foregone  | *Yes* |
| self_var_payoff | self payoff + self foregone - 10 | |
| other_var_payoff | other payoff + other foregone - 10  | |
| ineq_var_abs | abs(self_var_payoff - other_var_payoff) | |
| ineq_choice_abs | abs(self payoff - other payoff) | |
| ineq_foregone_abs | abs(self foregone - other foregone) | |
| ineq_var_disadvant |  (other_var_payoff > self_var_payoff) * (other_var_payoff - self_var_payoff) |  *Yes* |
| ineq_var_advant | (other_var_payoff < self_var_payoff) * (other_var_payoff - self_var_payoff) |  *Yes* |
| self_diff_ch_foregone | self payoff - self foregone (?) | |
| other_diff_ch_foregone | other payoff - other foregone (?) | |
| side_chosen | side chosen | |
| RT | reaction time | |

We are wary of the overlap and correlations between these variables. All the variables marked as a priority are being sued in the prelim analyses to establish our pipeline. The two pipelines are explained below. For `Pipeline 1` we are using only `ineq_var_advant` and `ineq_vardisadvant`. For `pipeline 2` we are using all variables marked as a priority.

_(?): couldn't find definition in the code so I took my best guess_

### Dependent Variables

We z score the high gamma (I think check on this). Then for the binning we as a first pass in subject IR35 are using bins with size 100ms and overalp of 50ms. In presentation locked analyses we take 200ms before presentation and 3 seconds afterwards (though will change them). The 200ms before presentation onset is used to baseline the data. This process results in ~60 bins, which corresponds to 60 dependent variables of interest.

## Regressions-- Corrections, Significance Thresholds, CV, etc

### Individual regressions

The OFC-DM paper said there had to be a 5 bin stretch to count as an active electrode-- do we want to keep this scheme? Wasn't in the code.

Use permuted null distribution?

### Model selection

`Pipeline 1:` 
[1] Run linear models for each bin for each electrode for each inequity measure
[2] Find the longest stretch of significant bins for each inequite measure for each electrode
[3] Sum the F-statistic for the longest stretch
[4] Shuffle the predictor 10,000 times and repreat steps 1-3.
[5] Calcualte the permutation p val for the sum of f statistics.
[6] Mark all electrodes with the permuted p val under .05 as active electrodes (should we correct here?)
[7] For all active electrodes and for all signigicant bins run an anova between two regressions where the base model is the consitituent parts of either advant. or disadvant. inequality and the second region is the base model with either apporpriate inequality mesure included as a predictor.
[8] FDR correct these anovas.

We might want to alter how we do the anovas to leverage the time course component-- perhaps following a similar procedure as the summed f stat and permutation testing.

`Pipeline 2` 
_This follow exactly the anova procedure layed out in the OFC-DM paper_
[1] Run linear models for each bin for each electrode for each inequity measure and each self_payoof, self_foregone, other_payoff, other_foregone.
[2] Find the longest stretch of significant bins for each inequite measure for each electrode
[3] Sum the F-statistic for the longest stretch
[4] Shuffle the predictor 10,000 times and repreat steps 1-3.
[5] Run ...


