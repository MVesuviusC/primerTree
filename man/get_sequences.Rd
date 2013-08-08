\name{get_sequences}
\alias{get_sequences}
\title{Retrieves fasta sequences from NCBI nucleotide database.}
\usage{
  get_sequences(gi, start, end)
}
\arguments{
  \item{start}{start bases to retrieve, numbered beginning
  at 1.  If NULL the beginning of the sequence.}

  \item{end}{end base to retrieve, numbered beginning at 1.
  if NULL the end of the sequence.}
}
\value{
  an ape::DNAbin object.
}
\description{
  Retrieves fasta sequences from NCBI nucleotide database.
}
\seealso{
  \code{\link{ape::DNAbin}}
}

