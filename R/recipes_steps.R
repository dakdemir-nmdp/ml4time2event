#'
#'
t2edata_split<-function(data, ...){
  train_test_split <- rsample::initial_split(data, ...)
  list(
    Train = rsample::training(train_test_split),
    Test = rsample::testing(train_test_split),
    Split = train_test_split
  )
}




#'
#'
t2emodel_data_recipe_init<-function(timevar, eventvar, expvar,idvars,  traindata){

  frml<-as.formula(paste(paste(timevar, eventvar, sep="+"),paste(expvar, collapse="+"), sep="~") )
  # Ensure subsetting returns a data.frame, not a vector
  recipe_vars <- c(timevar, eventvar, expvar, idvars)
  # Initialize recipe with the full training data
  out<-recipes::recipe(frml,
                       data = traindata) # Use full traindata for initialization
  out$vars <- recipe_vars # Restore custom attribute assignment
  return(out)
}





#'
#'
minimal_data_recipe<-function(model_recipe, pmiss=.3, pother=.05,dummy=TRUE, onehot=FALSE){
  vars<-model_recipe$vars
  if (dummy){
  out<-model_recipe %>%

    recipes::step_filter_missing(threshold = pmiss) %>%
  # mean impute numeric variables
    recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
  # convert the additional ingredients variable to dummy variables
    recipes::step_other(threshold = pother)%>%
    recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = onehot) %>%
  # remove predictor variables that are almost the same for every entry
    recipes::step_nzv(recipes::all_predictors()) %>% # Move nzv before range
    recipes::step_zv(recipes::all_predictors()) %>% # Add step_zv
  # rescale all numeric variables except for vanilla, salt and baking powder to lie between 0 and 1
    recipes::step_range(recipes::all_numeric_predictors(), min = 0, max = 1)

  } else {
   out <- model_recipe %>%
     recipes::step_filter_missing(threshold = pmiss) %>%
      # mean impute numeric variables
     recipes::step_impute_mean(recipes::all_numeric_predictors()) %>%
      # convert the additional ingredients variable to dummy variables
     recipes::step_other(threshold = pother)%>%
     recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
      # remove predictor variables that are almost the same for every entry
     recipes::step_nzv(recipes::all_predictors()) %>% # Move nzv before range
     recipes::step_zv(recipes::all_predictors()) %>% # Add step_zv
      # rescale all numeric variables except for vanilla, salt and baking powder to lie between 0 and 1
     recipes::step_range(recipes::all_numeric_predictors(), min = 0, max = 1)
  }

  out$vars<-vars
  return(out)
}




#'
#'
data_recipe_corstep<-function(model_recipe, treshold=.9){
  vars<-model_recipe$vars

    out<-model_recipe %>%
      recipes::step_corr(recipes::all_numeric_predictors(), threshold = treshold)

  out$vars<-vars
  return(out)
}



#'
#'
rfimp_data_recipe<-function(model_recipe, dummy=TRUE, onehot=FALSE, pmiss=.3, pother=.05, trees=10, options = list(nbagg=10,keepX = TRUE)){
  vars<-model_recipe$vars

  if (dummy){
  out<-model_recipe %>%
    recipes::step_filter_missing(threshold = pmiss) %>%
    recipes::step_other(threshold = pother)%>%
    # mean impute numeric variables
    recipes::step_impute_bag(recipes::all_predictors(), trees=trees, options=options)%>%
    recipes::step_dummy(recipes::all_nominal_predictors(),one_hot = onehot) %>%
    # remove predictor variables that are almost the same for every entry
    recipes::step_nzv(recipes::all_predictors()) %>% # Move nzv before range
    recipes::step_zv(recipes::all_predictors()) %>% # Add step_zv
    # rescale all numeric variables except for vanilla, salt and baking powder to lie between 0 and 1
    recipes::step_range(recipes::all_numeric_predictors(), min = 0, max = 1)
  } else {
    out<-model_recipe %>%
      recipes::step_filter_missing(threshold = pmiss) %>%
      recipes::step_other(threshold = pother)%>%
      # mean impute numeric variables
      recipes::step_impute_bag(recipes::all_predictors(), trees=trees, options=options)%>%
      # remove predictor variables that are almost the same for every entry
      recipes::step_nzv(recipes::all_predictors()) %>% # Move nzv before range
      recipes::step_zv(recipes::all_predictors()) %>% # Add step_zv
      # rescale all numeric variables except for vanilla, salt and baking powder to lie between 0 and 1
      recipes::step_range(recipes::all_numeric_predictors(), min = 0, max = 1)
  }
  out$vars<-vars
  return(out)
}



#'
#'
prep_data_recipe<-function(model_recipe,training){
  # Simplify: Directly return the output of prep, remove vars handling
  return(recipes::prep(model_recipe, training = training))
}
###############


#'
#'
bake_data_recipe<-function(model_recipe,data){
  # Pass the full new data to bake()
  out<-recipes::bake(model_recipe, new_data = data) # Use new_data argument
  return(out)
}
