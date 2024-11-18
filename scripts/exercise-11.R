## EXERCISE 11
## BAYESIAN ANALYSIS
## learning how to use brm versions of analyses I have already done this semester
## 
## Luke Watson
## 11/12/24 - 11/18/24
## 
## code written on machine then analyses ran and saved as .rds on 4 beocat cores
## tidying tables and plotting then done on my machine


# load libraries ----------------------------------------------------------

library(tidyverse)
library(showtext) # font import and graphics output
library(brms)
library(bayesplot)
library(tidybayes)
library(marginaleffects)
library(broom.mixed)
library(splines)


# theme set ---------------------------------------------------------------

# palette set, acadia national forest
clrs <- NatParksPalettes::natparks.pals("Acadia")

# take a look at the pretty colors!
df <- tibble(
  color = seq(1:9), 
  name = c(as.character(clrs))
)

ggplot(df, aes(x = color, fill = name)) +
  geom_bar(position = "stack") + 
  scale_fill_manual(values = c(clrs))

# add font from google
font_add_google("Palanquin Dark")

# ensure font is displayed in plot viewer
showtext_auto()

# ensure showtext renders to ggsave, will make text look huge in viewer!
showtext_opts(dpi = 300)

new_theme <- function() {
  theme_bw() + 
    theme(
      text = element_text(family = "Palanquin Dark"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      plot.title = element_text(face = "bold"), 
      strip.text = element_text(face = "bold"), 
      strip.background = element_rect(fill = "grey80", color = NA), 
      #legend.position = "none",
      panel.border = element_blank(),
      axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
    ) 
}

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

summary(brm_murder)

# terrible convergence, rhats are huge!!

# plot
plot(brm_murder)

# no fuzzy caterpillars either! chains are not showing good convergence at all!!

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

summary(brm_murder_c)
# much better rhat and ESS

plot(brm_murder_c)

# chains are converging much better than before

# compare results for centered predictors between trad and brm
table_lm <- tidy(lm_murder_c)
table_lm <- table_lm |> 
  mutate(
    model = "linear", 
    .before = "term"
  )

table_lm <- table_lm |> 
  mutate(
    across(where(is.numeric), ~round(., digits = 3))
  )

table_brm <- tidy(brm_murder_c)
table_brm <- table_brm |> 
  mutate(
    model = "bayesian lm", 
    .before = "term",
    across(where(is.numeric), ~round(., digits = 3))
  ) |> 
  select(c(4:9))

# save as csv
write.table(table_lm, "tables/lm-murder-full.csv", sep = ", ")
write.table(table_brm, "tables/brm-murder-full.csv", sep = ", ")


### BEOCAT
# run a main effects only bayesian model to compare to full factorial model 
# uninformative prior, avoid using BF to compare
brm_main_only <- brm(
  murder ~ c.assault + c.rape + c.urban_pop, data = usarrests,
  cores = 4
)

# save it 
write_rds(brm_main_only, "rds/brm-e.rds")

###

# read it in
brm_main_only <- read_rds("rds/brm-e.rds")

summary(brm_main_only)
plot(brm_main_only)

# instead can use waic()
waic(
  brm_main_only, # main effects only
  brm_murder_c # full factorial
)

# appears to be in favor of the main effects only model, I think I will go with
# this as my preferred model due to parsimony and slightly improved fit

# plot full density plots of each main effect parameter in each model and compare

# need some code here, mcmc_dens() within the hw for each pars = c(b_c.assault, ...)
mcmc_dens(
  as.array(brm_main_only), pars = c("b_c.assault", "b_c.rape", "b_c.urban_pop")
)

color_scheme_set("viridisA")

# plot central posterior estimates for main eff parameters, 
# shows density distrib. with shaded region for 68% cred interval
mcmc_areas(
  as.array(brm_main_only), pars = c("b_c.assault", "b_c.rape", "b_c.urban_pop"), 
  
) +
  labs(
    x = "b",
    y = "Predictor"
  ) +
  scale_y_discrete(labels = c("Assault", "Rape", "Urban Pop.")) +
  new_theme()

# save it
ggsave(
  "plots/usa-main-eff-area-plot.png", device = "png",
  width = 6.5, height = 6, units = "in"
)

# use conditional_effects() to include plots of three main effects 

# again, need some code, some offered by Mike in hw, also could look at heiss blog
cond_eff_us_arrests <- conditional_effects(brm_main_only)

# plot assaults me
cond_eff_us_arrests[[1]] |> 
  ggplot(aes(x = c.assault + mean(usarrests$assault), y = estimate__)) +
  
  # plot fitted line
  geom_line(color = clrs[8], linewidth = 1) + 
  
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), 
    fill = clrs[5], alpha = .5
  ) +
  
  # plot orig data
  geom_point(
    data = usarrests, aes(x = assault, y = murder),
    color = clrs[6], alpha = .7
  ) +
  
  labs(
    x = "Assault", y = "Murder"
  ) + 
  
  new_theme()

# save it
ggsave(
  "plots/assault-me-cond-eff-plot.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# plot rape
cond_eff_us_arrests[[2]] |> 
  ggplot(aes(x = c.rape + mean(usarrests$rape), y = estimate__)) +
  
  # plot fitted line
  geom_line(color = clrs[8], linewidth = 1) + 
  
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), 
    fill = clrs[5], alpha = .5
  ) +
  
  # plot orig data
  geom_point(
    data = usarrests, aes(x = rape, y = murder),
    color = clrs[6], alpha = .7
  ) +
  
  labs(
    x = "Rape", y = "Murder"
  ) + 
  
  new_theme()

# save it
ggsave(
  "plots/rape-me-cond-eff-plot.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# plot urban pop
cond_eff_us_arrests[[3]] |> 
  ggplot(aes(x = c.urban_pop + mean(usarrests$urban_pop), y = estimate__)) +
  
  # plot fitted line
  geom_line(color = clrs[8], linewidth = 1) + 
  
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), 
    fill = clrs[5], alpha = .5
  ) +
  
  # plot orig data
  geom_point(
    data = usarrests, aes(x = urban_pop, y = murder),
    color = clrs[6], alpha = .7
  ) +
  
  labs(
    x = "Urban Pop.", y = "Murder"
  ) + 
  
  new_theme()

# save it
ggsave(
  "plots/urbanpop-me-cond-eff-plot.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# or could try slopes 
conf_slope_assault <- slopes(
  brm_main_only
)

# seems to give a LOT smaller errors, which I'm not experienced enough to know 
# where this difference is coming from, likely should just save this function
# for plotting more complex relationships within model like when I need to 
# hold random effects constant or display an interaction 

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
tidy_causality <- tidy(brm_causality) |> 
  filter(effect == "fixed")

# plot out each par, both rhats and fuzzy caterpillars show good convergence 
plot(brm_causality, variable = "^b_", regex = T)

# save it, idk what else to do for ggplot of this type of plot 
ggsave(
  "plots/causality-flat-param-plot.png", device = "png",
  width = 8, height = 8, units = "in"
)

###

# check out conditional effects
causality_95 <- conditional_effects(brm_causality)

# plot out prob of detection main effect 
causality_95[[1]] |> 
  ggplot(
    aes(x = c.detection_prob + mean(causality$detection_prob), y = estimate__)
  ) +
  
  geom_line(col = clrs[1], linewidth = 1) + 
  
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), 
    alpha = .4, fill = clrs[3]
  ) +
  
  labs(
    x = "P(Detection)", 
    y = "P(Correct)"
  ) +
  
  new_theme()

# save it 
ggsave(
  "plots/causality-flat-prob-me.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# try epreds
pd_me_dist <- brm_causality |> 
  epred_draws(
    newdata = expand_grid(
      sex = "female", 
      c.detection_prob = seq(0 - mean_dp, 1 - mean_dp, .05)
    ),
    re_formula = NA
  )

# plot out with ribbon
pd_me_dist |> 
  ggplot(aes(x = c.detection_prob + mean_dp, y = .epred)) +
  
  # plot line with ribbons
  stat_lineribbon() + 
  
  scale_fill_manual(values = c(clrs[4], clrs[3], clrs[2])) + 
  
  labs(
    x = "P(Detection)", y = "P(Correct)", fill = "Credible Interval"
  ) +
  
  new_theme() +
  
  theme(legend.position = "bottom")

# looks good! save it over the other one
ggsave(
  "plots/causality-flat-prob-me-ribbon.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# check out the me of sex with conditional effects vs predictions
# would just be a bar chart w error bars
causality_95[[2]] |> 
  ggplot(aes(x = sex, y = estimate__, fill = sex)) +
  
  # bar
  geom_bar(
    stat = "identity", 
    width = .5
  ) +
  
  geom_errorbar(
    aes(ymin = lower__, ymax = upper__), 
    col = "black", linewidth = 1, width = .25
  ) +
  
  labs(
    x = "Sex", 
    y = "P(Correct)"
  ) +
  
  scale_fill_manual(values = c(clrs[3], clrs[5])) +
  
  new_theme()

# save it 
ggsave(
  "plots/causality-sex-bar-plot.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# try predictions to get cool halfeye comparison that is a lot more bayesian
sex_pred <- predictions(
  brm_causality, 
  newdata = datagrid(sex = c("female", "male")), 
  by = "sex", 
  # conf_level = .68, change conf level, doesn't change halfeye interpret.
  re_formula = NA # conditional effect where random offsets set to 0
) |> 
  posterior_draws()

sex_pred |> 
  ggplot(
    aes(x = draw, fill = sex)
  ) +
  
  # halfeye plot combines interval with slab density plot
  stat_halfeye() + 
  
  scale_fill_manual(values = c(clrs[3], clrs[5])) + 
  
  labs(
    x = "P(Correct)", y = "Density", fill = "Sex"
  ) + 
  
  new_theme()

# looks great!! save it
ggsave(
  "plots/causality-flat-sex-halfeyes.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# now need to see what we can do with the interaction
# can try this code for conditional effect for each level of sex across detect. prob
mean_dp <- mean(causality$detection_prob)

sex_detect_ixn_dist <- brm_causality |> 
  epred_draws(
    newdata = expand_grid(
      c.detection_prob = seq(0 - mean_dp, 1 - mean_dp, by = .05), 
      sex = c("female", "male")
    ), 
    re_formula = NA # conditional effect
  )

sex_detect_ixn_dist |> 
  ggplot(aes(x = c.detection_prob + mean_dp, y = .epred)) +
  
  # plot line with ribbons
  stat_lineribbon() + 
  
  scale_fill_manual(values = c(clrs[4], clrs[3], clrs[2])) + 
  
  labs(
    x = "P(Detection)", y = "P(Correct)", fill = "Credible Interval"
  ) +
  
  facet_wrap(~sex) +
  
  new_theme() +
  
  theme(legend.position = "bottom")

# looks good as well! save it
ggsave(
  "plots/causality-flat-ixn-ribbon-plot.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# show 68% credible interval, check on how to do this
causality_68 <- conditional_effects(brm_causality, prob = .68)

# plot out prob of detection main effect 
causality_68[[1]] |> 
  ggplot(
    aes(x = c.detection_prob + mean(causality$detection_prob), y = estimate__)
  ) +
  
  geom_line(col = clrs[1], linewidth = 1) + 
  
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), 
    alpha = .4, fill = clrs[3]
  ) +
  
  ylim(.7, .9) +
  
  labs(
    x = "P(Detection)", 
    y = "P(Correct)"
  ) +
  
  new_theme()

# save it
ggsave(
  "plots/causality-flat-prob-me-68.png", device = "png", 
  width = 7, height = 5, units = "in"
)

# bar
causality_68[[2]] |> 
  ggplot(aes(x = sex, y = estimate__, fill = sex)) +
  
  # bar
  geom_bar(
    stat = "identity", 
    width = .5
  ) +
  
  geom_errorbar(
    aes(ymin = lower__, ymax = upper__), 
    col = "black", linewidth = 1, width = .25
  ) +
  
  labs(
    x = "Sex", 
    y = "P(Correct)"
  ) +
  
  scale_fill_manual(values = c(clrs[3], clrs[5])) +
  
  new_theme()

# save it
ggsave(
  "plots/causality-sex-bar-plot-68.png", device = "png", 
  width = 7, height = 5, units = "in"
)

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
hypoth_causality <- hypothesis(brm_informed_causality, "sex1 = 0")

# superimpose the prior and posterior 
plot(hypoth_causality)

# looks good, but not a ggplot object! hard to format
# can combine samples/draws from prior and posterior to recreate same plot
prior_caus <- hypoth_causality[[3]] |> 
  mutate(
    dist = "prior"
  )
post_caus <- hypoth_causality[[2]] |> 
  mutate(
    dist = "posterior"
  )

both_dist <- rbind(prior_caus, post_caus)

both_dist |> 
  ggplot(aes(x = H1, fill = dist)) +
  
  # halfeyes
  stat_halfeye(alpha = .6) +
  
  scale_fill_manual(values = c(clrs[3], clrs[6])) +
  
  labs(
    x = "Difference in P(Correct) between F/M", 
    y = "Density", 
    fill = ""
  ) + 
  
  new_theme()
  
# looks good! save it
ggsave(
  "plots/hypoth-test-dist-halfeyes.png", device = "png", 
  width = 7, height = 5, units = "in"
)


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

# both rhats and fuzzy caterpillars show good convergence for both parameters

# save it
write_rds(brm_citations, "rds/brm-3a.rds")
###

# read in model
brm_citations <- read_rds("rds/brm-3a.rds")

# report individ scientist coef
table_AB <- coef(brm_citations)

# some code from mike to plot individual scientists
conditions <- data.frame(scientist = unique(citations$scientist))

rownames(conditions) <- unique(citations$scientist)

ce <- conditional_effects(
  brm_citations, conditions = conditions, re_formula = NULL
)

ce[[1]] |> 
  ggplot(aes(x = index, y = estimate__, col = scientist)) +
  
  # fitted model line from conditional effects
  geom_line() +
  
  scale_color_manual(values = c(clrs[3:8])) +
  
  # error ribbons
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), alpha = .3
  ) +
  
  scale_fill_manual(values = c(clrs[3:8])) +
  
  # orig data points
  geom_point(
    data = citations, aes(x = index, y = citations)
  ) +
  
  labs(
    x = "Index", y = "Citations"
  ) +
  
  facet_wrap(~scientist) +
  
  new_theme() + 
  
  theme(legend.position = "none")

# save it
ggsave(
  "plots/citation-individ-curves.png", device = "png", 
  width = 8, height = 6, units = "in"
)

# try with epred draws to get ribbon of overall curve 
sci_ce_dist <- brm_citations |> 
  epred_draws(
    newdata = expand_grid(
      index = seq(1, 40, 1)
    ),
    re_formula = ~ (1 | scientist)
  )

# plot
sci_ce_dist |> 
  ggplot(aes(x = index, y = .epred)) +
  
  # plot line with ribbons
  stat_lineribbon() + 
  
  scale_fill_manual(values = c(clrs[6], clrs[7], clrs[8])) + 
  
  labs(
    x = "Index", y = "Citations", fill = "Credible Interval"
  ) +
  
  new_theme() +
  
  theme(legend.position = "bottom")

# save it!
ggsave(
  "plots/citations-overall-nls-mixed.png", device = "png",
  width = 7, height = 5, units = "in"
)

# looks like the same ribbon and line for each! not sure how to troubleshoot atm

### BEOCAT
# try a spline fit 
spline_citations <- brm(
  citations ~ ns(index) + (ns(index) | scientist), data = citations,
  family = poisson()
)

# verify convergence
plot(spline_citations, variable = "^b_", regex = T)

# seems to have converged well with fuzzy caterpillar

# save it
write_rds(spline_citations, "rds/brm-3d.rds")
###

# read it in
spline_citations <- read_rds("rds/brm-3d.rds")

# plot the spline to visualize fit
# does not provide code to do this, could use predict or sim. function for brm?

conditions <- data.frame(scientist = unique(citations$scientist))

rownames(conditions) <- unique(citations$scientist)

ce <- conditional_effects(
  spline_citations, conditions = conditions, re_formula = NULL
)

ce[[1]] |> 
  ggplot(aes(x = index, y = estimate__, col = scientist)) +
  
  # fitted model line from conditional effects
  geom_line() +
  
  scale_color_manual(values = c(clrs[3:8])) +
  
  # error ribbons
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__), alpha = .3
  ) +
  
  scale_fill_manual(values = c(clrs[3:8])) +
  
  # orig data points
  geom_point(
    data = citations, aes(x = index, y = citations)
  ) +
  
  labs(
    x = "Index", y = "Citations"
  ) +
  
  facet_wrap(~scientist) +
  
  new_theme() + 
  
  theme(legend.position = "none")

# save it, looks really similar to the mixed eff nls ran
ggsave(
  "plots/citations-spline-individ-curves.png", device = "png",
  width = 8, height = 6, units = "in"
)

# done!