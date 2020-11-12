# Optimal design for phase 2 studies of SARS-CoV-2 antiviral drugs


This does a bunch of simulations to help design phase 2 trials of SARS-CoV-2 antivirals.

The model and data are taken from https://github.com/gradlab/CtTrajectories
The model has been slightly re-parameterised so that clearance slopes are inferred rather than clearance times. This helps make the posterior predictive behave better (otherwise clearance slope is dependent on peak viral load).

The R script run_power_simulations.R does a whole of power calculations comparing the use of time to clearance and rate of clearance.

The script adaptive_randomisation_simulation.R simulates an adaptive randomisation trial whereby the randomisation ratios of intervention arms are made proportional to the probability of superiority compared to control.


