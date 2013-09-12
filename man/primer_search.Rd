\name{primer_search}
\alias{primer_search}
\title{Query a pair of primers using ncbi's Primer-BLAST, if primers contain iupac}
\usage{
  primer_search(forward, reverse, ..., .parallel = FALSE,
    .progress = "none")
}
\arguments{
  \item{forward}{forward primer to search.}

  \item{reverse}{reverse primer to search.}

  \item{.parallel}{if 'TRUE', perform in parallel, using
  parallel backend provided by foreach}

  \item{.progress}{name of the progress bar to use, see
  'plyr::create_progress_bar'}
}
\value{
  data.frame of primer hits
}
\description{
  ambiguity codes, enumerate all possible combinations and
  combine the results.
}

