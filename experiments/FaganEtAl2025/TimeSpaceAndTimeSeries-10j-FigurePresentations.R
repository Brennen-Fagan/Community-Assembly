# Presentation variations:

variant <- c("Networks")[1]

if (variant == "Networks") {
  ### Networks oriented: ######################################################
  #### Setup: #################################################################
  
  
  source("TimeSpaceAndTimeSeries-10h-PlotPreparations.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsRichness.R")
  source("TimeSpaceAndTimeSeries-10i-PreparationsAbund.R")
  source(file.path("R", "flattenDiversity.R")) # Req'd by below
  source(file.path("R", "generateNetworks.R")) # To create inset graphs.
  
  figureNetworks <- list(
    graph = list(
      seed = "2", # "11", "17", "2"!,
      time = 25000
    )
  )
  
  #### Cluster of Single Adaptation Type Figures: #############################
  ##### Figure 2: Richness, Networks through Time #############################
  
  ##### Figure 3: Intervention, Richness, Networks through Time ###############
  
  ##### Figure 4: Richness, Abundance, Turnover, Complexity (RATC) ############
  
  #### 2 and Multiple Adaptation Type Figures: ################################
  ##### Figure 5: Richness w/Time, Abundance, Turnover, Complexity ############
  
  ##### Figure 6: Short Term Int. RATC ########################################
  
  #### Summary Images: ########################################################
  ##### Figure 7a: Parameters Cause RATC ######################################
  
  ##### Figure 7b: Network reorganisation over short time scales ##############
  
}