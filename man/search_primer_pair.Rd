\name{search_primer_pair}
\alias{search_primer_pair}
\title{Automatic primer searching'}
\usage{
  search_primer_pair(forward, reverse, name = NULL,
    num_aligns = 500, num_permutations = 25,
    simplify = TRUE, ..., .parallel = FALSE,
    .progress = "none")
}
\arguments{
  \item{name}{name to give to the primer pair}

  \item{simplify}{use simple names for primer hit results
  or complex}

  \item{.progress}{name of the progress bar to use, see
  'plyr::create_progress_bar'}

  \item{forward}{forward primer to search by 5'-3' on plus
  strand}

  \item{reverse}{reverse primer to search by 5'-3' on minus
  strand}

  \item{num_aligns}{number of alignment results to keep}

  \item{num_permutations}{the number of primer permutations
  to search, if the degenerate bases cause more than this
  number of permutations to exist, this number will be
  sampled from all possible permutations.}

  \item{.parallel}{if 'TRUE', perform in parallel, using
  parallel backend provided by foreach}
}
\description{
  Automatic primer searching'
}

