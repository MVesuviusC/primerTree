\name{primers_search}
\alias{primers_search}
\title{Query multiple pairs of primers using ncbi's Primer-BLAST, if primers contain iupac}
\usage{
  primers_search(forward = NULL, reverse = NULL, ...,
    .parallel = FALSE)
}
\arguments{
  \item{forward}{forward primer to search.}

  \item{reverse}{reverse primer to search.}
}
\value{
  list of data.frames of primer hits
}
\description{
  ambiguity codes, enumerate all possible combinations and
  combine the results.
}

