# Markov_simulations
The goal of this collection is to simulate the specific class of Markov chains and calculate relevant information theory-based measures. The code is oriented around `markov_process` S4 class object. 

## Current Functionalities
1. Markov chain simulation
2. Stationary density calculation
3. Markov chain simulation of selected nodes only (so called Markov Filtration)
4. Calculate transfer entropy between target nodes and all other. 

## How to run
All required functions and procedures are present at `simulations.R` file. 

### Comparison to Ay Polani information flow using an example from their article.
Script for it may be found in `mutual_info_comparison.R` followed by some other visualisations. 

## Bulid version and dependancies
- R version 4.2.2 (2022-10-31 ucrt)
- data.table_1.14.6 
- dplyr_1.0.10  




