# if (calculationsBII) {
# Adjusting
# https://adrianadepalma.github.io/BII_tutorial/bii_example.html#Compositional_Similarity
# to match our formatting.
# Note we get both the old and new versions.

getJacAbSym <- function(s1, s2, s1Control, s2Control, data){

  # "get the list of species that are present in site 1 (i.e., their abundance was greater than 0)"
  s1species <- data %>%

    # "filter out[sic] the SSBS that matches s1" (pristine)
    # (Note extra filter to make sure pristine/control because timeseries.)
    dplyr::filter(Patch == s1,
                  Control == s1Control ) %>%

    # "filter out[sic] the species where the Measurement (abundance) is greater than 0"
    dplyr::filter(Abundance > 0) %>%

    # "get the unique species from this dataset"
    dplyr::distinct(ID) %>%

    # "pull out the column into a vector"
    dplyr::pull(ID)

  # "for site 2, get the total abundance of species that are also present in site 1"

  s2abundance_s1species <- data %>%

    # "filter out[sic] the SSBS that matches s2"
    # (If timeseries, s1 might be s2, in which case grab not control data.)
    dplyr::filter(Patch == s2, Control == s2Control ) %>%

    # "filter out[sic] the species that are also present in site 1"
    dplyr::filter(ID %in% s1species) %>%

    # "pull out the Measurement into a vector"
    dplyr::pull(Abundance) %>%

    # "calculate the sum"
    sum()

  # "calculate the total abundance of all species in site 2"
  s2_sum <- data %>%

    # "filter out[sic] the SSBS that matches s2"
    dplyr::filter(Patch == s2) %>%

    # "pull out the measurement column (the abundance)"
    dplyr::pull(Abundance) %>%

    # "calculate the sum"
    sum()


  # "Now calculate the compositional similarity
  # this is the number of individuals of species also found in s1, divided by the total abundance in s2
  # so that equates to the proportion of individuals in s2 that are of species also found in s1"

  sor <- s2abundance_s1species / s2_sum


  # "if there are no taxa in common, then sor = 0
  # if abundances of all taxa are zero, then similarity becomes NaN."
  return(sor)

}

# NOTE: SIMILARITIES
get_bray <- function(s1, s2, s1Control, s2Control, data){

  require(betapart)

  sp_data <- data %>%

    # filter patches to the pair we're interested in
    dplyr::filter(paste(Patch, Control) %in%
                    c(paste(s1, s1Control),
                      paste(s2, s2Control))
    ) %>%

    dplyr::mutate(Patch = paste(Patch, Control)) %>%

    # pull out the site name, species name and abundance information
    dplyr::select(Patch, ID, Abundance) %>%

    # pivot so that each column is a species and each row is a Patch
    tidyr::pivot_wider(names_from = ID, values_from = Abundance,
                       values_fill = 0 # For some reason missing in orig. func?
    ) %>%

    # set the rownames to the Patch and then remove that column
    tibble::column_to_rownames("Patch")

  # if one of the sites doesn't have any individuals in it
  # i.e. the row sum is 0
  if(sum(rowSums(sp_data) == 0, na.rm = TRUE) == 1){
    # then the similarity between sites should be 0
    bray <- 0
    # if both sites have no individuals
  }else if(sum(rowSums(sp_data) == 0, na.rm = TRUE) == 2){
    # then class the similarity as NA
    bray <- NA
    # otherwise if both sites have individuals, calculate the balanced bray-curtis
    # as similarity (1-bray)
  }else{
    bray <- 1 -
      betapart::bray.part(sp_data) %>%
      purrr::pluck("bray.bal") %>%
      purrr::pluck(1)
  }
  return(bray)
}

# NOTE: DISSIMILARITIES
get_bray_all <- function(s1, s2, s1Control, s2Control, data){
  require(betapart)
  sp_data <- data %>%
    # filter patches to the pair we're interested in
    dplyr::filter(paste(Patch, Control) %in%
                    c(paste(s1, s1Control),
                      paste(s2, s2Control))
    ) %>%
    dplyr::mutate(Patch = paste(Patch, Control)) %>%
    # pull out the site name, species name and abundance information
    dplyr::select(Patch, ID, Abundance) %>%
    # pivot so that each column is a species and each row is a Patch
    tidyr::pivot_wider(names_from = ID, values_from = Abundance,
                       values_fill = 0 # For some reason missing in orig. func?
    ) %>%
    # set the rownames to the Patch and then remove that column
    tibble::column_to_rownames("Patch")

  bray <- betapart::bray.part(sp_data) %>% unlist() %>% t() %>% as.data.frame()
  return(bray)
}

inv_logit <- function(f, a) {
  a <- (1 - 2*a)
  (a * (1 + exp(f) ) + (exp(f) - 1)) / (2 * a * (1 + exp(f) ))
}
# }
