# Annotation Models: Statistical Models Developed For Estimating Annotator and Item Parameters in Human Annotation Data

`{AnnotationModelsR}` implements a number of well-known statistical models developed to estimate annotator and item parameters in human annotation data. 

## Currently featured models

- `em()`: the Dawid-Skene EM annotator-effect model for categorical annotation data proposed in Dawid, A. P. and A. M. Skene (1979) "Maximum likelihood estimation of observer error-rates using the EM algorithm." *Applied statistics*, 28(1), pp. 20--28. DOI: [10.2307/2346806](https://doi.org/10.2307/2346806)

## To be done

- implement tests for `em()`
- implement Methods for `em.fit` plot functions for categorical data (extend plots for binary label classes)