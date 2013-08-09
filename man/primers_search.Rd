\name{primers_search}
\alias{primers_search}
\title{Query multiple pairs of primers using ncbi's Primer-BLAST, if primers contain iupac}
\usage{
  primers_search(forward = NULL, reverse = NULL,
    name = paste(forward, reverse, sep = "_"), ...)
}
\arguments{
  \item{name}{name to use for the primer pair, if not given
  is simply the forward and reverse primers seperated by an
  underscore.}

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

