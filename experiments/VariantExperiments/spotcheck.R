load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/Data_2024-01-29/MNA-ExampleExtProp-Result-Env10-Ring-5-1-1-ExtProp1-Intervention.RData")
resultIntervention <- result
load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/Data_2024-01-29/MNA-ExampleExtProp-Result-Env10-Ring-5-1-1-ExtProp1.RData")
resultNull <- result
rm(result)

all.equal(resultIntervention$NumEnvironments, resultNull$NumEnvironments) %>% print()
all.equal(resultIntervention$ReactionTime, resultNull$ReactionTime)%>% print()
all.equal(resultIntervention$HistorySeed, resultNull$HistorySeed)%>% print()
all.equal(resultIntervention$Parameters, resultNull$Parameters)%>% print()
all.equal(resultIntervention$Ellipsis, resultNull$Ellipsis)%>% print()

all.equal(
  resultNull$Abundance[
    1:(which.max(resultNull$Abundance[, 1] > 
                   resultNull$Ellipsis$Intervention$Time) - 1
    ),], 
  resultIntervention$Abundance[
    1:(which.max(resultIntervention$Abundance[, 1] > 
                   resultNull$Ellipsis$Intervention$Time) - 1),])%>% print()

all.equal(
  resultNull$Events[
    1:(which.max(resultNull$Events$Times > 
                   resultNull$Ellipsis$Intervention$Time) - 1),], 
  resultIntervention$Events[
    1:(which.max(resultIntervention$Events$Times > 
                   resultIntervention$Ellipsis$Intervention$Time) - 1),])%>% print()