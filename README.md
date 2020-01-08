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

Will we run eletrodes that are partly in the wm in their notes? What percentage cut-off should we use? Some subjects look like they were run on an old system and don't have the full brain atlas information, can we rerun them? Should I add the cingulate?

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

### Morlet Wavelets or the filter-Hilbert method?

I&Z used the filter + the hilbert transform for their high frequency analysis. This is the the absolute value of the hilbert transform which is called the analytical amplitude. This is the value that is averaged in the moving time window.

### Filtering

Z&I look like they bandpassed at 70-200Hz/250Hz.

## Time binning (for regressions and high gamma)

Z & I were looking at the following bin sizes: `200, 160, 120, 80, 40` They all had overlaps of 25% regardless of size. It also looks like they were epoching based on trial start, presentation time, and choice time.

## Regressions-- Modeling Choices

| Possible Predictors | Definition | Priority? |
| ------------------- | :-------: | :-------: |
| self_choice_payoff | self payoff | |
| other_choice_payoff | other payoff | |
| self_foregone | self foregone | |
| other_foregone | other foregone  | |
| self_var_payoff | self payoff + self foregone - 10 | |
| other_var_payoff | other payoff + other foregone - 10  | |
| ineq_var_abs | abs(self_var_payoff - other_var_payoff) | |
| ineq_choice_abs | abs(self payoff - other payoff) | |
| ineq_foregone_abs | abs(self foregone - other foregone) | |
| ineq_var_disadvant |  (other_var_payoff > self_var_payoff) * (other_var_payoff - self_var_payoff) | |
| ineq_var_advant | (other_var_payoff < self_var_payoff) * (other_var_payoff - self_var_payoff) | |
| self_diff_ch_foregone | self payoff - self foregone (?) | |
| other_diff_ch_foregone | other payoff - other foregone (?) | |
| side_chosen | side chosen | |
| RT | reaction time | |

_(?): couldn't find definition in the code so I took my best guess_

Do we want to z score or baseline HG? 

## Regressions-- Corrections, Significance Thresholds, CV, etc

### Individual regressions

Some cut off for number of epochs an electrode has to be active for it to count as an encoding electrode. In OFC-DM paper it was 5 epochs.

Use permuted null distribution?

### Model selection

Looks like I&Z were using stepwise methods for their regressions, but only as a confirmation that findings were not due to collinearity. The scripts are easy to adapt for our purposes, but if we choose to deviate from his methods here are some other options:

LASSO, Ridge, etc. The metric of interest would then be the beta weights, where "high" beta weights would indicate a relationship w/ HG. How we quantify "high" would require some null model/bootstrapping something since coef testing is inappropriate for these methods. It also would require some hold test to identify a high lambda and we would have to decide on a error metric (or just go with l2). But these methods are very robust to noise, it is very interpretable, and we would have a single model for each of our elctrodes. With ~200 datapoints we do probably have enough for train-test splits. There is some caution about using LASSO to identify "true models" rather than models with low prediction error. But this is so exploratory and we don't have a chance so for the true model, so I think it probably okay?

Alternatively, we could have small models that we test in each model. We could group the variables by what we are invetisgating-- for inequality we might start with a model with the difference for self and the difference for other. Then, if we get a hit, we can parse down into the other ways of quanitfying inequality listed above. If we have real preferences/thoughts about what it is we are looking for that could be a good solution. But, since it is in the p-val framework, I do not know how we can apply corrections in particularly rigorous way. We would also probably want to lay out the hierarchy of tests in the beginning, but we do not have all that much to go on.






