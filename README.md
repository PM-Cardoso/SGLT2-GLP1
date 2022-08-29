# SGLT2-GLP1
Collection of functions and scripts to investigate how clinical features can be used for prediction.


---

Files:
- 0.0: Functions
    - .1: Functions for plotting and some calculations.
    - .2: Functions for _bartMachine_ tree analysis. Slightly modified _bartMan_ R package (not kept up to date) (package requirements were modified).
    
- 1.0: Detailed explanation of steps taken in the selection of patients for our cohorts.

- 2.0: Collection of plots demonstrating specific details/quirks of the dataset.

- 3.0: R packages to model causal treatment effect.
    - .1: Fitting of a causal model using _grf_ R package. This includes an evaluation of model fit.
    - .2: Fitting of a causal model using _bcf_ R package. This includes an evaluation of model fit.
    
- 4.0: _bartMachine_ models for treatment heterogeneity. 
    - .1: Fitting a BART model with variable selection for propensity score and outcome model. This includes an evaluation of model fit.
    - .2: Fitting a BART model with routine variables in propensity score model and biomarkers in outcome model. This includes an evaluation of model fit.
    - .3: Fitting a BART model with variable selection for propensity score and variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
    - .4: Fitting a BART model with variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
    
- 5.0: _bartMachine_ models using no methodical procedure.
    - .1: Fitting a collection of naive Bart models for HbA1c outcome using routine clinical variables / all variables / propensity scores, alternating between them.
    
- 6.0: Comparing models.
    - .1: Collection of plots comparing naive BART models in 5.0.
    - .2: Collection of plots comparing _bcf_ and _bartMachine_ with the same variables. Head-to-head comparisons of treatment effect for 3.2. vs 5.1. model 1 (Complete/Routine)
    - .3: Collection of plots comparing 4.1-4.4 models.

- 7.0: Simulations.
    - .1: 