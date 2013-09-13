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

  \item{num_aligns}{number of alignment results to keep}

  \item{simplify}{use simple names for primer hit results
  or complex}

  \item{.progress}{name of the progress bar to use, see
  'plyr::create_progress_bar'}

  \item{forward}{forward primer to search by 5'-3' on plus
  strand}

  \item{reverse}{reverse primer to search by 5'-3' on minus
  strand}

  \item{.parallel}{if 'TRUE', perform in parallel, using
  parallel backend provided by foreach}
}
\description{
  Automatic primer searching'
}

