## The following script reads the prepared mangabey data, conducts the model comparisons for all four interaction types (aggression, supplants, grooming per second, 1m proximity)
## Please keep in mind that these analyses can take a very long time to run (about 1h for all models each of the interaction types with 8 cores)
## Raw data are not provided due to data sharing policies of the field site. Individuals are anonymised. For questions or raw data, please contact mielke.alexand@gmail.com

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(brms)
library(bayesplot)
library(rstanarm)
library(rstan)


# Parameters --------------------------------------------------------------
cores <- 8
chains <- 3
iterations <- 3000
warmup <- 1500


# Load Data ---------------------------------------------------------------

mangabey.data <- readRDS('mangabey_data.RDS')

# Model Specifications ----------------------------------------------------
## please see paper for questions regarding the models
## all rank variables are z-standardised

## Indices
### 'elo' = Elo ranking sensu Foerster et al 2016, Model 3
### 'david' = Normalised David Scores, using steepness package

## Standardisation
### 'raw' = raw Elo or David scores, not equidistant
### 'ordinal' = ordinal ranks (highest-ranking individual is 1, lowest-ranking is N[subjects])
### 'stan' = proportional ranks (highest-ranking individual is 1, lowest-ranking is 0)

## Specification
### 'main' = main effects only (rank receiver + rank sender)
### 'difference' = rank difference variable and main effect sender (rank difference + rank sender)
### 'absolute.diff' = absolute rank difference variable and main effect sender (absolute rank difference + rank sender)
### 'interaction' = interaction term of sender and receiver ranks (rank receiver * rank sender)
### 'square.interaction' = squared interaction term of sender and receiver ranks ((rank receiver^2) * rank sender)

form.list <- list(
  elo.main.raw = "behaviour ~ 1 + z.elo.receiver.raw + z.elo.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.main.ordinal = "behaviour ~ 1 + z.elo.receiver.ordinal + z.elo.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.main.stan = "behaviour ~ 1 + z.elo.receiver.stan + z.elo.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.difference.raw = "behaviour ~ 1 + z.elo.difference.raw + z.elo.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.difference.stan = "behaviour ~ 1 + z.elo.difference.stan + z.elo.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.difference.ordinal = "behaviour ~ 1 + z.elo.difference.ordinal + z.elo.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.absolute_diff.raw = "behaviour ~ 1 + z.abs.elo.difference.raw * z.elo.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.absolute_diff.stan = "behaviour ~ 1 + z.abs.elo.difference.stan * z.elo.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.absolute_diff.ordinal = "behaviour ~ 1 + z.abs.elo.difference.ordinal * z.elo.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.interaction.raw = "behaviour ~ 1 + z.elo.sender.raw * z.elo.receiver.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.interaction.stan = "behaviour ~ 1 + z.elo.sender.stan * z.elo.receiver.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.interaction.ordinal = "behaviour ~ 1 + z.elo.sender.ordinal * z.elo.receiver.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.square_interaction.raw = "behaviour ~ 1 + z.elo.sender.raw * (z.elo.receiver.raw + I(z.elo.receiver.raw^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.square_interaction.stan = "behaviour ~ 1 + z.elo.sender.stan * (z.elo.receiver.stan + I(z.elo.receiver.stan^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  elo.square_interaction.ordinal = "behaviour ~ 1 + z.elo.sender.ordinal * (z.elo.receiver.ordinal + I(z.elo.receiver.ordinal^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  higher_lower = "behaviour ~ 1 + z.elo.sender.stan + higher_lower + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.main.raw = "behaviour ~ 1 + z.david.receiver.raw + z.david.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.main.stan = "behaviour ~ 1 + z.david.receiver.stan + z.david.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.main.ordinal = "behaviour ~ 1 + z.david.receiver.ordinal + z.david.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.difference.raw = "behaviour ~ 1 + z.david.difference.raw + z.david.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.difference.stan = "behaviour ~ 1 + z.david.difference.stan + z.david.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.difference.ordinal = "behaviour ~ 1 + z.david.difference.ordinal + z.david.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.absolute_diff.raw = "behaviour ~ 1 + z.abs.david.difference.raw * z.david.sender.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.absolute_diff.stan = "behaviour ~ 1 + z.abs.david.difference.stan * z.david.sender.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.absolute_diff.ordinal = "behaviour ~ 1 + z.abs.david.difference.ordinal * z.david.sender.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.interaction.raw = "behaviour ~ 1 + z.david.sender.raw * z.david.receiver.raw + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.interaction.stan = "behaviour ~ 1 + z.david.sender.stan * z.david.receiver.stan + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.interaction.ordinal = "behaviour ~ 1 + z.david.sender.ordinal * z.david.receiver.ordinal + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.square_interaction.raw = "behaviour ~ 1 + z.david.sender.raw * (z.david.receiver.raw + I(z.david.receiver.raw^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.square_interaction.stan = "behaviour ~ 1 + z.david.sender.stan * (z.david.receiver.stan + I(z.david.receiver.stan^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)",
  david.square_interaction.ordinal = "behaviour ~ 1 + z.david.sender.ordinal * (z.david.receiver.ordinal + I(z.david.receiver.ordinal^2)) + (1|year) + (1|sender) + (1|receiver) + offset(log.observation.time)"
)


# Priors ------------------------------------------------------------------
## only weakly informative priors: fixed effects are set to a normal distribution around 0 with a standard deviation of 2.5
m1priors <- c(prior(normal(0, 1), class = "b"))


# Aggression Models -------------------------------------------------------

## select aggression as behaviour of interest
mangabey.data$behaviour <- mangabey.data$aggression.sent

# run all specified model specifications using brm()
results.agg <- lapply(1:length(form.list), function(x) {
  res <- brm(
    data = mangabey.data,
    prior = m1priors,
    formula = form.list[[x]],
    family = zero_inflated_poisson(),
    warmup = warmup,
    iter = iterations,
    cores = cores,
    chains = chains,
    seed = 123, 
    control = list(adapt_delta = 0.99)
  )

  # save information and models
  return(list(
    c(
      loglik = mean(as.numeric(logLik(res))),
      loo_R2 = loo_R2(res)[1],
      R2 = bayes_R2(res)[1],
      loo = loo(res)$estimates[3]
    ),
    model = res
  ))
})
# extract summary information
results.aggression <-
  bind_rows(lapply(results.agg, function(x) {
    x[[1]]
  }))
# define model name
results.aggression <- results.aggression %>% 
  mutate(model = names(form.list),
# define behaviour
  behaviour = "aggression",
# define delta LOO
  deltaLOO = loo - min(loo)) %>%
  #arrange by LOO
  arrange(results.aggression$loo)

formula.meta <- str_split(results.aggression$model, pattern = '[.]') %>% 
  data.frame() %>% 
  t()
colnames(formula.meta) <- c('rank_variable', 'model_specification', 'standardisation')

results.aggression <- bind_cols(results.aggression, formula.meta)




# Grooming Models -------------------------------------------------------

## select grooming as behaviour of interest
mangabey.data$behaviour <- round(mangabey.data$grooming.seconds/60)

# run all specified model specifications using brm()
results.groom <- lapply(1:length(form.list), function(x) {
  res <- brm(
    data = mangabey.data,
    prior = m1priors,
    formula = form.list[[x]],
    family = zero_inflated_negbinomial(),
    warmup = warmup,
    iter = iterations,
    cores = cores,
    chains = chains,
    seed = 123, 
    control = list(adapt_delta = 0.99)
  )
  
  # save information and models
  return(list(
    c(
      loglik = mean(as.numeric(logLik(res))),
      loo_R2 = loo_R2(res)[1],
      R2 = bayes_R2(res)[1],
      loo = loo(res)$estimates[3]
    ),
    model = res
  ))
})
# extract summary information
results.grooming <-
  bind_rows(lapply(results.groom, function(x) {
    x[[1]]
  }))
# define model name
results.grooming <- results.grooming %>% 
  mutate(model = names(form.list),
         # define behaviour
         behaviour = "grooming",
         # define delta LOO
         deltaLOO = loo - min(loo)) %>%
  #arrange by LOO
  arrange(results.grooming$loo)

formula.meta <- str_split(results.grooming$model, pattern = '[.]') %>% 
  data.frame() %>% 
  t()
colnames(formula.meta) <- c('rank_variable', 'model_specification', 'standardisation')

results.grooming <- bind_cols(results.grooming, formula.meta)



# Proximity Models -------------------------------------------------------

## select proximity as behaviour of interest
mangabey.data$behaviour <- mangabey.data$proximity.sent

# run all specified model specifications using brm()
results.prox <- lapply(1:length(form.list), function(x) {
  res <- brm(
    data = mangabey.data,
    prior = m1priors,
    formula = form.list[[x]],
    family = poisson(),
    warmup = warmup,
    iter = iterations,
    cores = cores,
    chains = chains,
    seed = 123, 
    control = list(adapt_delta = 0.99)
  )
  
  # save information and models
  return(list(
    c(
      loglik = mean(as.numeric(logLik(res))),
      loo_R2 = loo_R2(res)[1],
      R2 = bayes_R2(res)[1],
      loo = loo(res)$estimates[3]
    ),
    model = res
  ))
})
# extract summary information
results.proximity <-
  bind_rows(lapply(results.prox, function(x) {
    x[[1]]
  }))
# define model name
results.proximity <- results.proximity %>% 
  mutate(model = names(form.list),
         # define behaviour
         behaviour = "proximity",
         # define delta LOO
         deltaLOO = loo - min(loo)) %>%
  #arrange by LOO
  arrange(results.proximity$loo)

formula.meta <- str_split(results.proximity$model, pattern = '[.]') %>% 
  data.frame() %>% 
  t()
colnames(formula.meta) <- c('rank_variable', 'model_specification', 'standardisation')

results.proximity <- bind_cols(results.proximity, formula.meta)



# Supplant Models -------------------------------------------------------

## select supplant as behaviour of interest
mangabey.data$behaviour <- mangabey.data$supplant.sent

# run all specified model specifications using brm()
results.supp <- lapply(1:length(form.list), function(x) {
  res <- brm(
    data = mangabey.data,
    prior = m1priors,
    formula = form.list[[x]],
    family = zero_inflated_poisson(),
    warmup = warmup,
    iter = iterations,
    cores = cores,
    chains = chains,
    seed = 123, 
    control = list(adapt_delta = 0.99)
  )
  
  # save information and models
  return(list(
    c(
      loglik = mean(as.numeric(logLik(res))),
      loo_R2 = loo_R2(res)[1],
      R2 = bayes_R2(res)[1],
      loo = loo(res)$estimates[3]
    ),
    model = res
  ))
})
# extract summary information
results.supplant <-
  bind_rows(lapply(results.supp, function(x) {
    x[[1]]
  }))
# define model name
results.supplant <- results.supplant %>% 
  mutate(model = names(form.list),
         # define behaviour
         behaviour = "supplant",
         # define delta LOO
         deltaLOO = loo - min(loo)) %>%
  #arrange by LOO
  arrange(results.supplant$loo)

formula.meta <- str_split(results.supplant$model, pattern = '[.]') %>% 
  data.frame() %>% 
  t()
colnames(formula.meta) <- c('rank_variable', 'model_specification', 'standardisation')

results.supplant <- bind_cols(results.supplant, formula.meta)

