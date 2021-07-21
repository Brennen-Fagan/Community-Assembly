Microbiome_NumericalAssembly <- function(
  Species, CarryingCapacity,
  FracExploit = 0.5,
  FracCompete = 0.25,
  FracMutual = 0.25
) {
  stopifnot(Species > 0)
  stopifnot(CarryingCapacity > 0)
  stopifnot(FracExploit > 0)
  stopifnot(FracCompete > 0)
  stopifnot(FracMutual > 0)

  # Normalise parameters. This allows ratio specification as well.
  if ((FracSum <- FracExploit + FracCompete + FracMutual) > 1) {
    FracExploit <- FracExploit / FracSum
    FracCompete <- FracCompete / FracSum
    FracMutual <- FracMutual / FracSum
  }

  # Create a pool.
  # Create interaction matrix.
  # While community non-empty, community invadable, and community not steady...
  # Draw species from pool and add to community.
  # Run community.
  # Note the carrying capacity addresses the mutual benefaction problem.
  # Evaluate community size, invadability, and steady-state.

}

Microbiome_PermanenceAssembly <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # While community non-empty, community invadable, and community not steady...
  # Draw species from pool and add to community.
  # Evaluate community permanence.
  # If not permanent, compute what it collapses to.
  # Evaluate community size, invadability, and steady-state.
}

Microbiome_SolveForStability <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # Add all of pool to community.
  # Solve growth rates of community.
  # Check boundedness of community; if unbounded, we discard.
  # Check repulsion of boundary; if not repelling, we discard.
}

Microbiome_CollapseCommunity <- function(

) {
  # Create a pool.
  # Create interaction matrix.
  # Add all of pool to community.
  # Run community.
  # (Do we need to check for extinction? Arguable.)
  # Result is whatever remains after running.
}
