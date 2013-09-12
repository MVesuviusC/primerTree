\name{search_primer}
\alias{search_primer}
\title{Automatic primer searching'}
\usage{
  search_primer(forward, reverse, name = NULL,
    num_aligns = 500, num_permutations = 25,
    simplify = TRUE, ..., .parallel = FALSE)
}
\arguments{
  \item{forward}{forward primer to search by 5'-3'}

  \item{reverse}{reverse primer to search by 3'-5'}

  \item{name}{name to give to the primer pair}

  \item{num_aligns}{name to give to the primer pair}

  \item{simplify}{use simple names for primer hit results
  or complex}

  \item{the}{number of primer permutations to search, if
  the degenerate bases cause more than this number of
  permutations to exist, this number will be sampled from
  all possible permutations.}

  \item{...}{additional arguments for primer_search, run
  primer_search with no arguments to see all available
  options}

  \item{.parallel}{use the parallel foreach, see foreach
  and plyr}
}
\description{
  Automatic primer searching'
}

