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

## Time binning (for Regressions)

## High Gamma Filtering

## Regressions-- Modeling Choices

## Regressions-- Corrections, Significance Thresholds, CV, etc
