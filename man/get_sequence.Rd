\name{get_sequence}
\alias{get_sequence}
\title{Retrieves a fasta sequence from NCBI nucleotide database.}
\usage{
  get_sequence(gi, start = NULL, end = NULL)
}
\arguments{
  \item{gi}{nucleotide gi to retrieve.}

  \item{start}{start base to retrieve, numbered beginning
  at 1.  If NULL the beginning of the sequence.}

  \item{end}{end base to retrieve, numbered beginning at 1.
  if NULL the end of the sequence.}
}
\value{
  an ape::DNAbin object.
}
\description{
  Retrieves a fasta sequence from NCBI nucleotide database.
}
\seealso{
  \code{\link{ape::DNAbin}}
}

