## Original author: Vincent Aranzana-Climent, Associate professor, University of Poitiers
## Modified by: Michael BURTIN, PhD student, University of Poitiers

## Generic function to calculate R-square values
calculate_rsq <- function(model)
{
  if(anyNA(model))
  {
    return(NA)
  }
  else {
    summary(model)
    coefficients(model)
    
    pred <- predict(model)
    n <- length(pred)
    res <- resid(model)
    w <- weights(model)
    if (is.null(w))
      w <- rep(1, n)
    rss <- sum(w * res ^ 2)
    resp <- pred + res
    center <- weighted.mean(resp, w)
    r.df <- summary(model)$df[2]
    int.df <- 1
    tss <- sum(w * (resp - center) ^ 2)
    r.sq <- 1 - rss / tss
    adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
    out <- list(pseudo.R.squared = r.sq,
                adj.R.squared = adj.r.sq)
  }
}

## Model definition
Emax_model <-  function(IndexValue, I0, Imax, IC50, H)
{
  I0 - (Imax * IndexValue ^ H) / (IC50 ^ H + IndexValue ^ H)
}

## Main function
index_PKPD_fit_curve <- function(data) {
  
  # Require check if the package is loadead. If not, the package is load.
  require(dplyr)        # for data manipulation
  require(tidyr)        # for nested tibbles
  require(purrr)        # for map
  require(minpack.lm)   # for nlsLM
  require(car)          # DeltaMethod
  
  ## Format the data
  nestedData <- data%>%mutate(deltaLog10CFU = MaxCFU - Log10CFU) |>  group_by(PKPD_Index)%>%nest()
  
  ## Model fitting with NLS
  nestedData$nls <-
    map(nestedData$data,
        ~ tryCatch(
          nlsLM(
          deltaLog10CFU ~ Emax_model(value, I0, Imax, IC50, H),
          data = .x,
          start = list(
            I0 = 4,
            Imax = 7,
            IC50 = 50,
            H = 2
          ),
          control = list(maxiter=100)
        ), error = function(cond){NA}))
  
  ## Calculate Rsq & RsqAdj
  nestedData$Rsq <-
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            }
            else{
              calculate_rsq(.x)[[1]]
            })
  
  ## Estimate each parameters
  nestedData$I0 <-
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            }
            else{
              coefficients(.x)[[1]]
            })
  
  nestedData$Imax <-
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            }
            else{
              coefficients(.x)[[2]]
            })
  
  nestedData$IC50 <-
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            }
            else{
              coefficients(.x)[[3]]
            })
  
  nestedData$H <- 
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            }
            else{
              coefficients(.x)[[4]]
            })
  
  nestedData$Sd_Residuals <- 
    map_dbl(nestedData$nls,
            ~ if (anyNA(.x)) {
              return(NA)
            } else {
              sd(residuals(.x))
            })
  
  return(nestedData)
}