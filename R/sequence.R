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

get_sequence = function(gi, start=NULL, end=NULL){
  library(httr)
  library(ape)

  fetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

  query=list(db='nuccore', rettype='fasta', retmode='text', id=gi)

  if(!is.null(start))
    query$seq_start = start

  if(!is.null(end))
    query$seq_stop = end

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

#' @param end end base to retrieve, numbered beginning at 1. if NULL the end of
#'        the sequence.
#' @return an ape::DNAbin object.
#' @seealso \code{\link{ape::DNAbin}}
#' @export get_sequences

#TODO fix this list to be one DNAbin object
get_sequences = function(gi, start, end){
  size = length(gis)
  sequences = vector('list', size)
  for(i in seq_along(gis)){
    sequences[[i]] = fetch_sequence(gis[i], starts[recycle(i,length(starts))], ends[recycle(i, length(ends))])
  }
  sequences
}
recycle = function(x, length){
  ((x-1) %% length) + 1
}
#TODO document this function
tree_from_fasta = function(dna){
  library(ape)
  nj(dist(dna))
}
