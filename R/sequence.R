#' Retrieves a fasta sequence from NCBI nucleotide database.
#'
#' @param gi nucleotide gi to retrieve.
#' @param start start base to retrieve, numbered beginning at 1.  If NULL the
#'        beginning of the sequence.

#' @param end end base to retrieve, numbered beginning at 1. if NULL the end of
#'        the sequence.
#' @return an ape::DNAbin object.
#' @seealso \code{\link{ape::DNAbin}}
#' @export get_sequence

get_sequence = function(gi, start=NULL, stop=NULL){

  fetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

  query=list(db='nuccore', rettype='fasta', retmode='text', id=gi)

  if(!is.null(start))
    query$seq_start = start

  if(!is.null(stop))
    query$seq_stop = stop

  response = POST(fetch_url, body=query)

  #stop if response failed
  stop_for_status(response)

  content = content(response, as='raw')

  #from ape package read.FASTA
  res <- .Call("rawStreamToDNAbin", content, PACKAGE = "ape")
  names(res) <- sub("^ +", "", names(res))
  class(res) <- "DNAbin"
  res
}

#' Retrieves fasta sequences from NCBI nucleotide database.
#'
#' @param start start bases to retrieve, numbered beginning at 1.  If NULL the
#'        beginning of the sequence.

#' @param stop stop bases to retrieve, numbered beginning at 1. if NULL the stop of
#'        the sequence.
#' @param simplify simplify the FASTA headers to include only the genbank
#'        accession.
#' @param .parallel perform loops in parallel
#' @return an ape::DNAbin object.
#' @seealso \code{\link{ape::DNAbin}}
#' @export get_sequences

get_sequences = function(gi, start=NULL, stop=NULL, simplify=TRUE, ..., .parallel=FALSE){
  size = length(gi)
  get_sequence_itr = function(i){
    sequence = get_sequence(gi[i], start[recycle(i,length(start))], stop[recycle(i, length(stop))])
  }
  sequences = alply(seq_along(gi), .margins=1, .parallel=.parallel, get_sequence_itr)
  names = if(simplify) gi else laply(sequences, names)
  sequences = llply(sequences, `[[`, 1)
  names(sequences) = names
  class(sequences) = 'DNAbin'
  sequences
}
recycle = function(x, length){
  ((x-1) %% length) + 1
}
#' Construct a neighbor joining tree from a dna alignment
#'
#' @param dna fasta dna object the tree is to be constructed from
#' @param ... furthur arguments to dist.dna
#' @seealso \code{\link{ape::dist.dna}}, \code{\link{ape::nj}}
#' @export tree_from_alignment
tree_from_alignment = function(dna, ...){
  nj(dist(dna, ...))
}
