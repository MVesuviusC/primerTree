\name{search_primer_pair}
\alias{search_primer_pair}
\title{Automatic primer searching Search a given primer pair, retrieving the alignment
results, their product sequences, the taxonomic information for the sequences,
a multiple alignment of the products}
\usage{
  search_primer_pair(forward, reverse, name = NULL,
    num_aligns = 500, num_permutations = 25,
    simplify = TRUE, clustal_options = list(), ...,
    .parallel = FALSE, .progress = "none")
}
\arguments{
  \item{name}{name to give to the primer pair}

  \item{simplify}{use simple names for primer hit results
  or complex}

  \item{.progress}{name of the progress bar to use, see
  \code{\link{create_progress_bar}}}

  \item{clustal_options}{a list of options to pass to
  clustal omega, run \code{link{clustalo}} for a list of
  options}

  \item{forward}{forward primer to search by 5'-3' on plus
  strand}

  \item{reverse}{reverse primer to search by 5'-3' on minus
  strand}

  \item{num_aligns}{number of alignment results to keep}

  \item{num_permutations}{the number of primer permutations
  to search, if the degenerate bases cause more than this
  number of permutations to exist, this number will be
  sampled from all possible permutations.}

  \item{...}{additional arguments passed to Primer-Blast}

  \item{.parallel}{if 'TRUE', perform in parallel, using
  parallel backend provided by foreach}
}
\value{
  A list with the following elements, \item{name
  name}{description name of the primer pair} \item{name
  BLAST_result}{description html blast results from
  Primer-BLAST as 'a \code{\link{response}}} object.
  \item{name taxonomy}{description taxonomy for the primer
  products from NCBI} \item{name sequence}{description
  sequence of the primer products} \item{name
  alignment}{description multiple alignment of the primer
  products} \item{name tree}{description phylogenetic tree
  of the reconstructed from the 'multiple alignment}
}
\description{
  Automatic primer searching Search a given primer pair,
  retrieving the alignment results, their product
  sequences, the taxonomic information for the sequences, a
  multiple alignment of the products
}
\examples{
\dontrun{
#simple search
mammals_16S = search_primer_pair(name='Mammals 16S',
 'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT')
#returning 1000 alignments, allow up to 3 mismatches in primer
mammals_16S = search_primer_pair(name='Mammals 16S',
 'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT',
 num_aligns=1000, total_primer_specificity_mismatch=3)
}
}
\seealso{
  \code{\link{primer_search}}, \code{\link{clustalo}}
}

