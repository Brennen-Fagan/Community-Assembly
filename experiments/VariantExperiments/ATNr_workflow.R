# ATNr workflow?
# Note the below use random number generation, so need to be careful with seeds.
# It doesn't look like ATNr implements errors or testing procedures
# (e.g. checking if there is a binarized matrix provided when one is asked for).
# Note also that I've occasionally had bad runs where nothing goes well.
# (Plots that are just blocks, empty, or crashed early.)

# ATNr also requires the data set "schneider" be findable, but it only loads
# on attaching the package.
library(ATNr)

# 0. Parameters: ##############################################################
species <- 200
connectance <- 0.1 # err, not sure of the reliability, but
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
initialConditionsNutrients <- runif(nutrients, 1, 100)
times <- c(seq(0, 100, 1), seq(101, 1000, 3), seq(1010, 5000, 10))

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
  # return(list(parms$ODE(y, t)))

  # with loss of < 1 individual per unit area:
  population <- parms$ODE(y, t)
  # If biomass / (body mass) < 1, remove
  popspecies <- population[1:species]
  popspecies[popspecies/masses < 1] <- 0
  population[1:species] <- popspecies
  return(list(population))
  # Somewhat surprisingly the *unscaled* version doesn't have any losses
  # under this rule.
}

solNUN <- deSolve::lsoda(
  c(initialConditionsSpecies, initialConditionsNutrients),
  times,
  evalODE,
  modelNicheUnscaledNutrients
)
solLUN <- deSolve::lsoda(
  c(initialConditionsSpecies, initialConditionsNutrients),
  times,
  evalODE,
  modelLUnscaledNutrients
)
solNS <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelNicheScaled
)
solBS <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelLBinarizedScaled
)
solLS <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelLScaled
)
solNU <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelNicheUnscaled
)
solBU <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelLBinarizedUnscaled
)
solLU <- deSolve::lsoda(
  c(initialConditionsSpecies),
  times,
  evalODE,
  modelLUnscaled
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
ATNr::plot_odeweb(solLS, species)
title(main = "LS")
ATNr::plot_odeweb(solNU, species)
title(main = "NU")
ATNr::plot_odeweb(solBU, species)
title(main = "BU")
ATNr::plot_odeweb(solLU, species)
title(main = "LU")
