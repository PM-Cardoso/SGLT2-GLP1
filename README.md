# SGLT2-GLP1
Collection of functions and scripts to investigate how clinical features can be used for prediction.


---

Files:
- 0.0:
    - .1: Functions for plotting and some calculations.
- 1.0: Detailed explanation of steps taken in the selection of patients for our cohorts.
- 2.0: Collection of plots demonstrating specific details/quirks of the dataset.
- 3.0: 
    - .1: Fitting of a causal model using _grf_ R package. This includes an evaluation of model fit.
    - .2: Fitting of a causal model using _bcf_ R package. This includes an evaluation of model fit.
- 4.0:
    - .1: Fitting a BART model with variable selection for propensity score and outcome model. This includes an evaluation of model fit.
    - .2: Fitting a BART model with routine variables in propensity score model and biomarkers in outcome model. This includes an evaluation of model fit.
    - .3: Fitting a BART model with variable selection for propensity score and variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
    - .4: Fitting a BART model with variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
