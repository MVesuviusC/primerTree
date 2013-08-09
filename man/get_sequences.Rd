\name{get_sequences}
\alias{get_sequences}
\title{Retrieves fasta sequences from NCBI nucleotide database.}
\usage{
  get_sequences(gi, start = NULL, stop = NULL,
    simplify = TRUE)
}
\arguments{
  \item{start}{start bases to retrieve, numbered beginning
  at 1.  If NULL the beginning of the sequence.}

  \item{stop}{stop bases to retrieve, numbered beginning at
  1. if NULL the stop of the sequence.}

  \item{simplify}{simplify the FASTA headers to include
  only the genbank accession.}
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

