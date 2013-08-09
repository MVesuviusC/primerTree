\name{primer_search}
\alias{primer_search}
\title{Query a pair of primers using ncbi's Primer-BLAST, if primers contain iupac}
\usage{
  primer_search(forward = NULL, reverse, ...)
}
\arguments{
  \item{forward}{forward primer to search.}

  \item{reverse}{reverse primer to search.}
}
\value{
  data.frame of primer hits
}
\description{
  ambiguity codes, enumerate all possible combinations and
  combine the results.
}

