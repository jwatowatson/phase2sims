# Optimal design for phase 2 studies of SARS-CoV-2 antiviral drugs


This does a bunch of simulations to help design phase 2 trials of SARS-CoV-2 antivirals.

The model and data are taken from https://github.com/gradlab/CtTrajectories
The model has been slightly re-parameterised so that clearance slopes are inferred rather than clearance times. This helps make the posterior predictive behave better (otherwise clearance slope is dependent on peak viral load).

The RMarkdown script *Comparing_time_vs_rate.Rmd* does a series of power calculations for different effect sizes and sample sizes, comparing the use of time-to-clearance and rate-of-clearance. This produces Figures 1-4 in the paper.

The RMarkdown script *Adaptive_randomisation.Rmd* simulates an adaptive randomisation trials (the randomisation ratios of intervention arms are made proportional to the probability of superiority compared to control) under two scenarios:

* Scenario 1: there are 4 intervention arms, 1 control arm. One of the 4 intervention arms increases viral clearance by 10\% relative to control.
* Scenario 1: there are 4 intervention arms, 1 control arm. None of the 4 intervention arms increase viral clearance relative to control.

Trials can stop early for success (>99\% probability that the viral clearance is increased by at least 1\% for at least 1 drug); or for futility (<10\% probability that the viral clearance is increased by 1\% for all drugs).


