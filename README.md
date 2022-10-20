# SGLT2-GLP1
Collection of functions and scripts to investigate how clinical features can be used for prediction.

## Final structure for analysis: (Model 4.7)

1. Fit a BART propensity score model (function _bartMachine::bartMachine_) with all variable available.
2. Perform BART variable selection (function _bartMachine::var_selection_by_permute_) to choose variables for the PS model.
3. Re-fit a BART propensity score model (function _bartMachine::bartMachine_) with the selected variables.
4. Fit a BART response model (function _bartMachine::bartMachine_):
    - Include missing data feature (use_missing_data = FALSE)
    - Include all available variables + propensity score values
5. Perform BART variable selection (function _bartMachine::var_selection_by_permute_cv_) to choose variables for the response model.
6. Re-fit a BART response model (function _bartMachine::bartMachine_) with the selected variables:
    - Include missing data feature (use_missing_data = FALSE)
7. Check model fit for the outcome:
    - Plot predicted vs observed to check for any structure in the values.
    - Plot standardised residuals to check for any structure in the residuals.
8. Check model fit for treatment effects: (model fitted in observational data)
    - Plot predicted CATE vs ATE.


---

Files:
- 0.0: Functions
    - .1: Functions for plotting and some calculations.
    - .2: Functions for _bartMachine_ tree analysis. Slightly modified _bartMan_ R package (not kept up to date) (package requirements were modified).
    
- 1.0: Detailed explanation of steps taken in the selection of patients for our cohorts.

- 2.0: Descriptive analysis of data
    - .1: Collection of plots demonstrating specific details/quirks of the dataset.
    - .2: Table description of Development and Validation data
    - .3: Treatment Effects/Variable description table/plot (Model 4.5)

- 3.0: R packages to model causal treatment effect.
    - .1: Fitting of a causal model using _grf_ R package. This includes an evaluation of model fit.
    - .2: Fitting of a causal model using _bcf_ R package. This includes an evaluation of model fit.
    
- 4.0: _bartMachine_ models for treatment heterogeneity. 
    - .1: Fitting a BART model with variable selection for propensity score and outcome model. This includes an evaluation of model fit.
    - .2: Fitting a BART model with routine variables in propensity score model and biomarkers in outcome model. This includes an evaluation of model fit.
    - .3: Fitting a BART model with variable selection for propensity score and variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
    - .4: Fitting a BART model with variable selection using BART + _grf_ for the outcome model. This includes an evaluation of model fit.
    - .5: Fitting a BART model with variable selection for using BART + _grf_ for the outcome model. This includes an evaluation of model fit. (Change from 4.4 - instead of 'score', we use 'score.excl.mi')
    - .6: Fitting a BART propensity score model, variable selection, matching individuals, BART model with all variables.
    - .7: **Fitting a BART propensity score model, variable selection, refit propensity score model, BART HbA1c model + propensity score as covariate, variable selection, refit BART HbA1c model.** (Model chosen) 
    
- 5.0: _bartMachine_ models using no methodical procedure.
    - .1: Fitting a collection of naive Bart models for HbA1c outcome using routine clinical variables / all variables / propensity scores, alternating between them.
    
- 6.0: Comparing models.
    - .1: Collection of plots comparing naive BART models in 5.0.
    - .2: Collection of plots comparing _bcf_ and _bartMachine_ with the same variables. Head-to-head comparisons of treatment effect for 3.2. vs 5.1. model 1 (Complete/Routine)
    - .3: Collection of plots comparing 4.1-4.4 models.
    - .4: Differential treatment effect.

- 7.0: Sensitivity analysis.
    - .1: Exclusion of GLP1 patients before 2013.
    - .2: Variable importance model 4.4.
    - .3: Modelling predicted treatment effect against model variables.
    
- 8.0: Presentations/Slides
    - .1: MRC: 29th September London
    - .2: SGLT2 vs GLP1 paper for publish
    
- 9.0: Shiny App
    - .1: Model 4.4 - probability of achieving target HbA1c.
    
- 10.0: Modelling accompanying data
    - .1: Weight reduction
    - .2: Discontinuation
    
---

- 11.0: Aurum download modelling
    - .1: Functions used specifically for this portion.
    - .2: Detailed explanation of the selection of cohorts.
    - .3: Descriptive analysis of datasets.
    - .4: Model heterogeneity.


    

