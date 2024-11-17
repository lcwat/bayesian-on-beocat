## EXERCISE 11
## BAYESIAN ANALYSIS
## learning how to use brm versions of analyses I have already done this semester
## 
## Luke Watson
## 11/12/24
## 
## code written on machine then analyses ran and saved as .rds on 4 beocat cores
## tidying tables and plotting then done on my machine


# load libraries ----------------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(marginaleffects)
library(broom.mixed)
library(splines)

# load data ---------------------------------------------------------------

# using the usarrests dataset
usarrests <- USArrests

# clean names
usarrests <- usarrests |> 
  rename(
    murder = Murder, assault = Assault, rape = Rape, urban_pop = UrbanPop
  )

# add mean centered predictors
usarrests <- usarrests |> 
  mutate(
    c.assault = assault - mean(assault),
    c.rape = rape - mean(rape),
    c.urban_pop = urban_pop - mean(urban_pop)
  )

# SKIP, FILE ALREADY CLEANED ----------------------------------------------

# clean names and factor
causality <- causality |> 
  rename(
    subject = Subject, 
    detection_prob = `Detection Probability`,
    sex = Sex,
    correct = Correct,
    median_latency = `Median(Latency)`
  )

# recode sex variable to introduce na so . is not treated as sep category
# also median latency should be numeric not character
causality <- causality |> 
  mutate(
    sex = case_when(
      sex == "m" ~ "male",
      sex == "f" ~ "female",
      sex == "\\." ~ NA
    ),
    median_latency = as.numeric(median_latency) # turn into numeric, . = na
  )

# make subject and sex into factors
causality <- causality |> 
  mutate(
    subject = factor(subject),
    sex = factor(sex)
  )

# looks good now, nas were introduced and are more appropriate than .
str(causality)
summary(causality)

# overwrite file to load in correct version
write_csv(causality, "data/causality-exp.csv")


# -------------------------------------------------------------------------

# load in causality data
causality <- read_csv("data/causality-exp.csv")

# factor
causality <- causality |> 
  mutate(
    subject = factor(subject),
    sex = factor(sex)
  )

# effect code sex and center detection_prob
causality <- causality |> 
  mutate(
    c.detection_prob = detection_prob - mean(detection_prob)
  )

contrasts(causality$sex) = contr.sum(2)

str(causality)

# load in citation data
citations <- read.csv("data/citations.csv")

# tidy names
citations <- citations |> 
  rename(
    scientist = Scientist,
    index = Index,
    citations = Citations
  )

# done

# compare freq and bayesian lm --------------------------------------------

# run traditional lm with assault rape and urban pop full factorial pred murder
lm_murder <- lm(murder ~ assault*rape*urban_pop, data = usarrests)

### BEOCAT
# run bayesian version of this analysis and take note of the results
brm_murder <- brm(
  murder ~ assault*rape*urban_pop, data = usarrests, cores = 4
)

# save the results
write_rds(brm_murder, "rds/brm-a.rds")
###


# load in the results
brm_murder <- read_rds("rds/brm-a.rds")

# plot
plot(brm_murder)

# small dataset and collinear predictors, check in on how the chains are converging

# rerun model with proper centering to address multicollinearity
lm_murder_c <- lm(murder ~ c.assault * c.rape * c.urban_pop, data = usarrests)

### BEOCAT
brm_murder_c <- brm(
  murder ~ c.assault * c.rape * c.urban_pop, data = usarrests,
  cores = 4
)

# save and read in results
write_rds(brm_murder_c, "rds/brm-d.rds")
###


brm_murder_c <- read_rds("rds/brm-d.rds")

# compare results for centered predictors between trad and brm
table_lm <- tidy(lm_murder_c)

table_brm <- tidy(brm_murder_c)

### BEOCAT
# run a main effects only bayesian model to compare to full factorial model 
# uninformative prior, avoid using BF to compare
brm_main_only <- brm(
  murder ~ c.assault + c.rape + c.urban_pop, data = usarrests,
  cores = 4
)

# save it 
write_rds(brm_main_only, "brm-e.rds")

###

# read it in
brm_main_only <- read_rds("brm-e.rds")

# instead can use waic()
waic(
  brm_main_only, # main effects only
  brm_murder_c # full factorial
)

# plot full density plots of each main effect parameter in each model and compare

# need some code here, mcmc_dens() within the hw for each pars = c(b_c.assault, ...)
# or posterior_draws() like in the heiss blog post

# use conditional_effects() to include plots of three main effects 

# again, need some code, some offered by Mike in hw, also could look at heiss blog


# mixed effects bayesian regression model ---------------------------------

### BEOCAT
# run brm with mixed effects
brm_causality <- brm(
  correct ~ c.detection_prob * sex + (c.detection_prob | subject),
  family = bernoulli(), # bernoulli error dist.
  data = causality, cores = 4
)

# save and read it in
write_rds(brm_causality, "rds/brm-2a.rds")
brm_causality <- read_rds("rds/brm-2a.rds")

# check out results and convergence for parameters
summary(brm_causality)

# tidy report
tidy_causality <- tidy(brm_causality)

# plot out each par
plot(brm_causality, variable = "^b_", regex = T)

###

# check out conditional effects, check heiss blog post about interpretation of this
# I think it is the partial deriv of effect when random effects are held constant
causality_95 <- conditional_effects(brm_causality)

# show 68% credible interval, check on how to do this
causality_68 <- conditional_effects(brm_causality)

### BEOCAT
# run a second model in with an informed prior for sex difference
get_prior(
  correct ~ c.detection_prob * sex + (c.detection_prob | subject),
  family = bernoulli(), # bernoulli error dist.
  data = causality
)

# looking to set prior of the sex1 parameter, using just sex will not work bc
# priors are created for each level of sex

brm_informed_causality <- brm(
  correct ~ c.detection_prob * sex + (c.detection_prob | subject),
  family = bernoulli(), # bernoulli error dist. 
  prior = c(set_prior("normal(.5, .1)", class = "b", coef = "sex1")),
  sample_prior = T, save_all_pars = T,
  data = causality, cores = 4
)
summary(brm_informed_causality)
plot(brm_informed_causality, variable = "^b_", regex = T)
# save and read it in
write_rds(brm_informed_causality, "rds/brm-2c.rds")
###

brm_informed_causality <- read_rds("rds/brm-2c.rds")

# run a hypothesis test to see if posterior is different than 0
hypoth_causality <- hypothesis(brm_informed_causality, "sex = 0")

# superimpose the prior and posterior 
plot(hypoth_causality)

# see if there is a ggplot version of this, may be posterior draws


# mixed effects nls ------------------------------------------------------

### BEOCAT
# create a brm version of a mixed effects nonlinear model 
brm_citations <- brm(
  bf(citations ~ A * exp(index*B), A+B ~ 1 + (1|scientist), nl = T),
  data = citations, family = poisson(),
  prior = c(
    set_prior("normal(400,400)", nlpar = "A"),
    set_prior("normal(-1, 1)", nlpar = "B")
  ),
  iter = 5000, warmup = 4000, # more samples to address converg. issues
  cores = 4
)

summary(brm_citations)

plot(brm_citations, variable = "^b_", regex = T)

# save it
write_rds(brm_citations, "rds/brm-3a.rds")
###

# read in model
brm_citations <- read_rds("rds/brm-3a.rds")

# some code from mike to plot individual scientists
conditions <- data.frame(scientist = unique(citations$scientist))

rownames(conditions) <- unique(citations$scientist)

ce <- conditional_effects(
  brm_citations, conditions = conditions, re_formula = NULL
)

ce[[1]] |> 
  ggplot(aes(x = index, y = estimate__, col = scientist)) +
  
  geom_line() +
  
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = .3) +
  
  geom_point(data = citations, aes(x = index, y = citations)) +
  
  facet_wrap(~scientist)

### BEOCAT
# try a spline fit 
spline_citations <- brm(
  citations ~ ns(index) + (ns(index) | scientist), data = citations,
  family = poisson()
)

# verify convergence
plot(spline_citations, variable = "^b_", regex = T)

# save it
write_rds(spline_citations, "rds/brm-3d.rds")
###

# read it in
spline_citations <- read_rds("rds/brm-3d.rds")

# plot the spline to visualize fit
# does not provide code to do this, could use predict or sim. function for brm?


