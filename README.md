# K-Fold Cross Validation Helper Functions

Date: 2017-06-08

Most Recent Update: 2017-06-08

Author: Connor Lenio

Dependencies: dplyr, caret, MASS, glmnet, parallel, foreach

Description: Functions to conduct k-fold cross validation, especially for with-sample error estimation of root mean squared error

Currently includes:

* klist - Creates a "klist" or a list of cross validation folds for a given data set - this list can then be used to quickly validate a model by estimate its average RMSE over the k-split data contained in the list. 

* validate_lm - Validates a linear model for a given klist and returns the RMSE values

* validate_glm_nb - Validates a Negative Binomial GLM using the splits contained in the klist

* validate_lasso - Validates a Lasso linear regression using the splits contained in the klist

* validate_lasso_nb - Validates a Lasso negative binomial regression using the splits contained in the klist
