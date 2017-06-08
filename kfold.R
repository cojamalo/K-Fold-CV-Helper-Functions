
# Title: K-Fold Cross Validation Helper Functions
# Date: 2017-06-08
# Most Recent Update: 2017-06-08
# Author: Connor Lenio
# Dependencies: dplyr, caret, MASS, glmnet
# Description: Functions to conduct k-fold cross validation, especially for with-sample error estimation of root mean squared error


# Creates a "klist" or a list of cross validation folds for a given data set - this list can then be used to quickly validate a model by estimate its average RMSE over the k-split data contained in the list. 
klist <- function(data, k = 10, repeats = 10) {
    for (i in 1:repeats) {
        set.seed(122+i)
        fold_size =  nrow(data) * 1/k
        # Generate data for folds
        final_list <- NULL
        test <- sample_n(data, size = fold_size)
        train <- suppressMessages(anti_join(data, test))
        final_list[[1]] <- list(run = i, train = train, test = test, rmse = NA)
        used <- NULL
        for (j in 2:k) {
            for (n in 1:(j-1)) {
                used <- rbind(used, final_list[[n]]$test)   
            }
            leftover <- suppressMessages(anti_join(data, used))
            test <- sample_n(leftover, size = fold_size)
            train <- suppressMessages(anti_join(data, test))
            new_row <- list(run = j, train = train, test = test, rmse = NA)
            final_list[[j]] <- new_row
        }
    }
    return(final_list)
}

# Validates a linear model for a given klist and returns the average RMSE value
validate_lm <- function(kfold_list, model_formula, method= "lm",k = 10, parallel = TRUE, ...) {
    pb <- txtProgressBar(min = 0, max = 10, style = 3)
    rmse_vec <- c()
    if (parallel) {
        output <- foreach(m = 1:k, .packages = c("caret"),  .combine='rbind')  %dopar%
        {
            model <- train(model_formula, method = method, data = kfold_list[[m]]$train, ...) # caret
            pred <- predict(model, kfold_list[[m]]$test)  
            actual <- kfold_list[[m]]$test[,all.vars(model_formula)[1]] 
            resid <- pred - actual
            rmse <- sqrt(mean(resid^2))
            kfold_list[[m]]$rmse <- rmse
            output <- data.frame(rmse = rmse)
            return(output)
        }
        return(output)
    }
    else {
        for (m in 1:k) {
            setTxtProgressBar(pb, m)
            model <- train(model_formula, method = method, data = kfold_list[[m]]$train, ...)
            pred <- predict(model, kfold_list[[m]]$test)   
            resid <- pred - kfold_list[[m]]$test[,all.vars(model_formula)[1]]
            rmse <- sqrt(mean(resid^2))
            kfold_list[[m]]$rmse <- rmse
            rmse_vec <- c(rmse_vec, rmse)
        }   
        return(rmse_vec)
    }
    
}

## Validates a Negative Binomial GLM using the splits contained in the klist
validate_glm_nb <- function(kfold_list, model_formula, method= "lm",k = 10, parallel = TRUE, ...) {
    rmse_vec <- c()
    if (parallel) {
        output <- foreach(m = 1:k, .packages = c("MASS"),  .combine='rbind')  %dopar%
        {
            model <- glm.nb(model_formula,data = kfold_list[[m]]$train) # MASS
            pred <- predict(model, kfold_list[[m]]$test, type = "response")  
            actual <- kfold_list[[m]]$test[,all.vars(model_formula)[1]] 
            resid <- pred - actual
            rmse <- sqrt(mean(resid^2))
            kfold_list[[m]]$rmse <- rmse
            output <- data.frame(rmse = rmse)
            return(output)
        }
        return(output)
    }
    else {
        for (m in 1:k) {
            setTxtProgressBar(pb, m)
            model <- glm.nb(model_formula,data = kfold_list[[m]]$train) # MASS
            pred <- predict(model, kfold_list[[m]]$test, type = "response")   
            resid <- pred - kfold_list[[m]]$test[,all.vars(model_formula)[1]]
            rmse <- sqrt(mean(resid^2))
            kfold_list[[m]]$rmse <- rmse
            rmse_vec <- c(rmse_vec, rmse)
        }   
        return(rmse_vec)
    }
    
}

## Validates a Lasso linear regression using the splits contained in the klist
validate_lasso <- function(kfold_list, response_str, k = 10, s = "lambda.1se", ...) {
    pb <- txtProgressBar(min = 0, max = 10, style = 3)
    rmse_vec <- c()
    for (m in 1:k) {
        setTxtProgressBar(pb, m) # utils
        registerDoMC() # doMC
        xval <- as.matrix(kfold_list[[m]]$train[,!colnames(kfold_list[[m]]$train) %in% response_str])
        yval <- as.matrix(kfold_list[[m]]$train[,response_str])
        newx <- as.matrix(kfold_list[[m]]$test[,!colnames(kfold_list[[m]]$test) %in% response_str])
        lasso_cv <- cv.glmnet(x = xval, y = yval, type.measure = "mse", nfolds = 10, parallel = TRUE, ...) #glmnet, dplyr
        pred <- predict(lasso_cv, newx = newx, s = s)
        actual <- kfold_list[[m]]$test[,response_str]
        resid <-  pred - actual
        rmse <- sqrt(mean(resid^2))
        kfold_list[[m]]$rmse <- rmse
        rmse_vec <- c(rmse_vec, rmse)
    }
    return(rmse_vec)
    
}

## Validates a Lasso negative binomial regression using the splits contained in the klist
validate_lasso_nb <- function(kfold_list, response_str, k = 10, s = "lambda.optim", theta = 1, alpha = 1, ...) {
    pb <- txtProgressBar(min = 0, max = 10, style = 3)
    rmse_vec <- c()
    for (m in 1:k) {
        setTxtProgressBar(pb, m) # utils
        registerDoMC() # doMC
        newx <- as.matrix(kfold_list[[m]]$test[,!colnames(kfold_list[[m]]$test) %in% response_str])
        lasso_cv <- cv.glmreg(count ~ ., data = train, family = "negbin", theta = theta, alpha = alpha) 
        if (s == "lambda.optim") { 
            which <- lasso_cv$lambda.optim 
            }
        else {
            which <- lasso_cv$lambda.which 
            }
        pred <- predict(lasso_cv, newx = newx, which = which)
        actual <- kfold_list[[m]]$test[,response_str]
        resid <-  pred - actual
        rmse <- sqrt(mean(resid^2))
        kfold_list[[m]]$rmse <- rmse
        rmse_vec <- c(rmse_vec, rmse)
    }
    return(rmse_vec)
    
}