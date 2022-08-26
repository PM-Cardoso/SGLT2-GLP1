####################
## Description:
##  - In this file we make a collected of plots comparing assessment
##      measurements of several models: 3.2. bcf vs 5.1. Complete/Routine model
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")


library(tidyverse)





#### BCF vs 5.1. model1 and model2.

## Plot x-axis Predicted tau from mean/median BCF vs prediction interval of treatment effect for the models.

## Plot x-axis Predited tau from 5.1. model vs bcf tau - predicted tau (iteration by iteration)