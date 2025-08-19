# ATNr workflow?
# Note the below use random number generation, so need to be careful with seeds.
# It doesn't look like ATNr implements errors or testing procedures
# (e.g. checking if there is a binarized matrix provided when one is asked for).
# Some initial problems regarding misreading some things, but should be fixed.

# ATNr also requires the data set "schneider" be findable, but it only loads
# on attaching the package.
library(ATNr)

# 0. Parameters: ##############################################################
species <- 200
connectance <- 0.1 # not sure about value. Chosen from:
# wiki.santafe.edu/images/f/f7/Dunne_-_Encyclopedia_of_Complex_Systems%2c_Food_Webs.pdf

# 1. Create a Matrix: #########################################################
## 1a: Niche Matrix: ##########################################################

foodwebNiche <- ATNr::create_niche_model(S = species, C = connectance)
# Is sufficient in and of itself for a foodweb matrix.

## 1b: Allometric Matrix: #####################################################

# The alternative with allometric scaling instead is to create an "L matrix".
# To create the equivalent:
speciesBasal <- sum(colSums(foodwebNiche) == 0)
# I think this is supposed to be body size as opposed to biomass.
# That is, big eats small.
masses <- 1e-2 * 10 ^ (ATNr::TroLev(foodwebNiche) - 1)

foodwebL <- ATNr::create_Lmatrix(masses, speciesBasal)
# Note Ropt : optimal body mass ratio.
#      gamma: width of things willing to eat.
#      th   : threshold for removal of links.
foodwebLBinarized <- foodwebL
foodwebLBinarized[foodwebLBinarized > 0] <- 1

# 2. Create the model with parameters: ########################################
## 2a: Schneider et al. 2016 Nutrient Model: ##################################
# I didn't see guidance on how to choose the number of nutrients.
# ATNr uses 3 in their vignette.
nutrients <- 3

# What about "unscaled"?
# "Biological rates controlling the growth rate are [not] normalised by the
#  growth rate of the smallest basal species."

# Note that the models are now C++ objects with active bindings.

modelNicheUnscaledNutrients <- ATNr::create_model_Unscaled_nuts(
  species, speciesBasal, nutrients, masses, foodwebNiche
)
modelLUnscaledNutrients <- ATNr::create_model_Unscaled_nuts(
  species, speciesBasal, nutrients, masses, foodwebLBinarized
)

# Note that the NicheUnscaledNutrients might be an abuse.
# The vignette notes we need to initialise the models, and the
# unscaled nutrients version specifically requires the L matrix with a
# specific interpretation of what the matix means.
# It isn't immediately clear if the niche version can accommodate that
# interpretation as well, or if it should only be used for the next two models.
# (Note, foodwebNiche is binary {0, 1}, foodwebL is probabilistic [0, 1].)
# No error is encountered, however.

modelNicheUnscaledNutrients <- ATNr::initialise_default_Unscaled_nuts(
  modelNicheUnscaledNutrients, foodwebNiche
) # also includes a temperature parameter? unclear how temperature works.
# code inspection suggests a reproductive rate where smaller than 20 decreases
# reproduction (and maybe interactions?) while bigger increases said rates.

modelLUnscaledNutrients <- ATNr::initialise_default_Unscaled_nuts(
  modelLUnscaledNutrients, foodwebL
)

## 2b: Delmas et al. 2016 Scaled Model: #######################################
# This does *NOT* have a temperature parameter and *DOES* scale so that the
# smallest basal has unit growth rate.
modelNicheScaled <- ATNr::create_model_Scaled(
  species, speciesBasal, masses, foodwebNiche
)
modelLBinarizedScaled <- ATNr::create_model_Scaled(
  species, speciesBasal, masses, foodwebLBinarized
)
# Unclear if the full L matrix should work.
# Results seem to indicate that this is the only combination that truly fails.
modelLScaled <- ATNr::create_model_Scaled(
  species, speciesBasal, masses, foodwebL
)


modelNicheScaled <- ATNr::initialise_default_Scaled(
  modelNicheScaled
)
modelLBinarizedScaled <- ATNr::initialise_default_Scaled(
  modelLBinarizedScaled
)
modelLScaled <- ATNr::initialise_default_Scaled(
  modelLScaled
)

## 2c: Binzer et al. 2016 Unscaled Model: #####################################
# This *DOES* have a temperature parameter and does *NOT* scale so that the
# smallest basal has unit growth rate.
modelNicheUnscaled <- ATNr::create_model_Unscaled(
  species, speciesBasal, masses, foodwebNiche
)
modelLBinarizedUnscaled <- ATNr::create_model_Unscaled(
  species, speciesBasal, masses, foodwebLBinarized
)
# Unclear if the full L matrix should work.
modelLUnscaled <- ATNr::create_model_Unscaled(
  species, speciesBasal, masses, foodwebL
)

# Form of the temperature parameter looks to be similar, but not identical, to
# the scaling used in 2a.

modelNicheUnscaled <- ATNr::initialise_default_Unscaled(
  modelNicheUnscaled
)
modelLBinarizedUnscaled <- ATNr::initialise_default_Unscaled(
  modelLBinarizedUnscaled
)
modelLUnscaled <- ATNr::initialise_default_Unscaled(
  modelLUnscaled
)

## 3: Run manually using deSolve: #############################################
# Why not use the provided wrappers?
# We want to be able to put these into a spatial context.
# That will require modifying how the system runs.

initialConditionsSpecies <- runif(species, 1, 10)
initialConditionsNutrients <- runif(nutrients, 1, 100) # 10 * > species biomass.
times <- c(seq(0, 1.0, .1), seq(1.1, 10.0, .3), seq(10.1, 200.0, 1.0))
rtol <- 1e-5; atol <- 1e-5

# Note that for custom code, they "REQUIRE" initialisations be done.
modelNicheUnscaledNutrients$initialisations()
modelLUnscaledNutrients$initialisations()
modelNicheScaled$initialisations()
modelLBinarizedScaled$initialisations()
modelLScaled$initialisations()
modelNicheUnscaled$initialisations()
modelLBinarizedUnscaled$initialisations()
modelLUnscaled$initialisations()

evalODE <- function(t, y, parms) {
  # Basic:
  return(list(parms$ODE(y, t)))
}
eventFunc <- list(
  func = function(t, y, parms) {y[y < 1e-6] <- 0;y},
  time = times
)

solNUNpre <- deSolve::lsoda( # Trying to get closer to steady state for nutr.
  y = c(initialConditionsNutrients, rep(0, species)), # Nutrients FIRST.
  times = 1:100,
  func = evalODE,
  modelNicheUnscaledNutrients,
  rtol = rtol, atol = atol
)
solNUN <- deSolve::lsoda( # Maybe no consumers???
  y = c(solNUNpre[nrow(solNUNpre), 1+1:nutrients],# Nutrients FIRST.
    initialConditionsSpecies),
  times = times,
  func = evalODE,
  events = eventFunc,
  modelNicheUnscaledNutrients,
  rtol = rtol, atol = atol
)
solLUN <- deSolve::lsoda(       # Single Basal? -> Double Basal?
  c(initialConditionsNutrients, # Nutrients FIRST.
    initialConditionsSpecies),  # Despite large number of inputs, almost all
  times,                        # die out before the end of the time, leaving
  evalODE,                      # essentially only basals.
  events = eventFunc,
  modelLUnscaledNutrients,
  rtol = rtol, atol = atol
)
solNS <- deSolve::lsoda(       # Less exciting than solBS, but seems to give
  c(initialConditionsSpecies), # more consistent results? Consumers still seem
  times,                       # to be far more common than the basal species.
  evalODE,                     # (Curves are generally flatter, slower.)
  events = eventFunc,
  modelNicheScaled,
  rtol = rtol, atol = atol
)
solBS <- deSolve::lsoda(       # Sometimes Single Basal? Seems the best
  c(initialConditionsSpecies), # performer. Non-obvious what to do with the
  times,                       # scaling though. As well, basal species seem
  evalODE,                     # to have substantially less biomass than the
  events = eventFunc,          # consumer species, but labelling is unclear.
  modelLBinarizedScaled,       # Also need to be careful with the elimination
  rtol = rtol, atol = atol     # threshold.
)
try(solLS <- deSolve::lsoda( # FAILED x2
  c(initialConditionsSpecies),
  times,
  evalODE,
  events = eventFunc,
  modelLScaled,
  rtol = rtol, atol = atol
))
solNU <- deSolve::lsoda( # These 3 seem uninteresting? (Flat)
  c(initialConditionsSpecies),
  times,
  evalODE,
  events = eventFunc,
  modelNicheUnscaled,
  rtol = rtol, atol = atol
)
solBU <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  events = eventFunc,
  modelLBinarizedUnscaled,
  rtol = rtol, atol = atol
)
solLU <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  events = eventFunc,
  modelLUnscaled,
  rtol = rtol, atol = atol
)

par(mfrow = c(3, 3))
ATNr::plot_odeweb(solNUN, species)
title(main = "NUN")
plot(x = 1, y = 1)
ATNr::plot_odeweb(solLUN, species)
title(main = "LUN")
ATNr::plot_odeweb(solNS, species)
title(main = "NS")
ATNr::plot_odeweb(solBS, species)
title(main = "BS")
if(exists("solLS")) {ATNr::plot_odeweb(solLS, species)} else {plot(1,1)}
title(main = "LS")
ATNr::plot_odeweb(solNU, species)
title(main = "NU")
ATNr::plot_odeweb(solBU, species)
title(main = "BU")
ATNr::plot_odeweb(solLU, species)
title(main = "LU")

