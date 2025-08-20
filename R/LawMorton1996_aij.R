LawMorton1996_aij <- function(
  Species_i, Species_j, k = c(0.01, 10, 0.5, 0.2, 100, 0.1), # Table 2 values.
  CompetitionBasal = 0, Connectance = 1, DiagParam = 0 # Suggested by Jon.
) {
  # returns 'p's for the effect of j on i, aij
  
  # i == j
  if (Species_i$ID == Species_j$ID) {
    if (Species_i$Type == "Basal") {
      return(- Species_i$ReproductionRate * Species_i$Size / k[5] + DiagParam)
    } else if(Species_i$Type == "Consumer") {
      return(0 + DiagParam)
    } else {
      return(NA)
    }
  }
  
  # i != j
  # This is not quite fair; this implements amensalism and commensalism
  # rather than proper connectance reduction in the undirected graph.
  if (runif(1) > Connectance) {
    return(0)
  }
  
  if (Species_i$Type == "Basal") {
    if (Species_j$Type == "Basal") {
      
      if (CompetitionBasal > 0) {
        # Without a good idea of what they otherwise would look like, I am
        # just going to borrow the predation effects.
        # This is not without precedent; Coyte et al. 2015 used the same dist.
        # but with different sign pairs for different interaction types.
        # We will use the CompetitionBasal parameter to scale this effect though.
        if (Species_i$Size < Species_j$Size) {
          # Predator on Prey
          # Note negativity.
          return(
            - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2) * CompetitionBasal
          )
          
        } else if (Species_i$Size > Species_j$Size) {
          # Prey on Predator
          # Note double negativity in the original, we remove one for competition.
          aji <-
            - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
          return(
            aji * k[4] * Species_j$Size / Species_i$Size * CompetitionBasal
          )
        } else {
          return(0)
        }
        
      } else {
        # "Basal species are assumed to be independent of one another."
        return(0)
      }
      
    } else if (Species_j$Type == "Consumer") {
      
      # average effect of consumer j on victim i.
      # aij = -k1 exp( -[ log10 (k2 si / sj) / k3]^2 ) if si < sj,
      # 0 otherwise
      if (Species_i$Size < Species_j$Size) {
        return(
          - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2)
        )
      } else {
        return(0)
      }
      
    } else {
      return(NA)
    }
  } else if (Species_i$Type == "Consumer") {
    if (Species_j$Type == "Basal") {
      
      # Then average effect of victim j on consumer i
      # aij = - aji k4 sj / si if sj < si, 0 otherwise.
      if (Species_j$Size < Species_i$Size) {
        aji <-
          - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
        return(
          - aji * k[4] * Species_j$Size / Species_i$Size
        )
      } else {
        return(0)
      }
      
    } else if (Species_j$Type == "Consumer") {
      # Then it comes down to their sizes.
      
      if (Species_i$Size < Species_j$Size) {
        # Predator on Prey
        return(
          - k[1] * exp(-(log10(k[2] * Species_i$Size / Species_j$Size) / k[3]) ^ 2)
        )
        
      } else if (Species_i$Size > Species_j$Size) {
        # Prey on Predator
        aji <-
          - k[1] * exp(-(log10(k[2] * Species_j$Size / Species_i$Size) / k[3]) ^ 2)
        return(
          - aji * k[4] * Species_j$Size / Species_i$Size
        )
      } else {
        return(0)
      }
      
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}