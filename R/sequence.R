#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Retrieves a fasta sequence from NCBI nucleotide database.
#'
#' @param gi nucleotide gi to retrieve.
#' @param start start base to retrieve, numbered beginning at 1.  If NULL the
#'        beginning of the sequence.

#' @param stop last base to retrieve, numbered beginning at 1. if NULL the end of
#'        the sequence.
#' @return an DNAbin object.
#' @seealso \code{\link{DNAbin}}
#' @export

get_sequence = function(gi, start=NULL, stop=NULL){

  fetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

  query=list(db='nuccore', rettype='fasta', retmode='text', id=gi)

  if(!is.null(start))
    query$seq_start = start

  if(!is.null(stop))
    query$seq_stop = stop

  response = POST_retry(fetch_url, body=query)

  #stop if response failed
  stop_for_status(response)

  content = content(response, as='raw')

  #from ape package read.FASTA
  res <- .Call("rawStreamToDNAbin", content)
  names(res) <- sub("^ +", "", names(res))
  class(res) <- "DNAbin"
  res
}

#' Retrieves fasta sequences from NCBI nucleotide database.
#'
#' @param gi the gi number of the sequence to retrieve
#' @param start start bases to retrieve, numbered beginning at 1.  If NULL the
#'        beginning of the sequence.

#' @param stop stop bases to retrieve, numbered beginning at 1. if NULL the stop of
#'        the sequence.
#' @param simplify simplify the FASTA headers to include only the genbank
#'        accession.
#' @param .parallel if 'TRUE', perform in parallel, using parallel backend
#'        provided by foreach
#' @param .progress name of the progress bar to use, see 'create_progress_bar'
#' @return an DNAbin object.
#' @seealso \code{\link{DNAbin}}
#' @export

get_sequences = function(gi, start=NULL, stop=NULL, simplify=TRUE, .parallel=FALSE, .progress='none'){
  #expand arguments by recycling
  args = expand_arguments(gi=gi, start=start, stop=stop)
  #assign expanded arguments to actual arguments
  lapply(seq_along(args), function(i) names(args)[i] <<- args[i])

  size = length(gi)
  get_sequence_itr = function(i){
    sequence = get_sequence(gi[i], start[i], stop[i])
  }
  sequences = alply(seq_along(gi), .margins=1, .parallel=.parallel, .progress=.progress, failwith(NA, f=get_sequence_itr))
  names = if(simplify) gi else laply(sequences, names)
  sequences = llply(sequences, `[[`, 1)
  names(sequences) = names
  class(sequences) = 'DNAbin'
  sequences
}
# from http://stackoverflow.com/questions/9335099/implementation-of-standard-recycling-rules
expand_arguments <- function(...){
  dotList <- list(...)
  max.length <- max(sapply(dotList, length))
  lapply(dotList, rep, length=max.length)
}
#' Construct a neighbor joining tree from a dna alignment
#'
#' @param dna fasta dna object the tree is to be constructed from
#' @param pairwise.deletion a logical indicating if the distance matrix should 
#' be constructed using pairwise deletion
#' @param ... furthur arguments to dist.dna
#' @seealso \code{\link{dist.dna}}, \code{\link{nj}}
#' @export
tree_from_alignment = function(dna, pairwise.deletion=TRUE, ...){
  nj(dist.dna(dna, model="N", pairwise.deletion=pairwise.deletion, ...))
}
#' Multiple sequence alignment with clustal omega
#'
#' Calls clustal omega to align a set of sequences of class DNAbin.  Run
#' without any arguments to see all the options you can pass to the command
#' line clustal omega.
#' @param x an object of class 'DNAbin'
#' @param exec a character string with the name or path to the program
#' @param quiet whether to supress output to stderr or stdout
#' @param original.ordering use the original ordering of the sequences
#' @param ... additional arguments passed to the command line clustalo
#' @export
clustalo = function (x, exec = 'clustalo', quiet = TRUE, original.ordering = TRUE, ...)
{
    help_text = system(paste(exec, '--help'), intern=TRUE)
    all_options = get_command_options(help_text)

    inf <- tempfile(fileext='.fas')
    outf <- tempfile(fileext='.aln')

    options = c(infile=inf, outfile=outf, list(...))
    match_args = pmatch(names(options), names(all_options), duplicates.ok=TRUE)
    bad_args = is.na(match_args)

    if (missing(x)){
        message(paste(help_text, collapse="\n"))
        stop('No input')
    }
    if(any(bad_args)){
      stop(paste(names(options)[bad_args], collapse=','), ' not valid option\n')
    }

    write.dna(x, inf, "fasta")
    args = paste(paste(all_options[match_args], options, collapse=' '))
    system2(exec, args=args, stdout = ifelse(quiet, FALSE, ''), stderr = ifelse(quiet, FALSE, ''))
    res <- read.dna(outf, "fasta")
    if (original.ordering)
        res <- res[labels(x), ]
    res
}

#parses the usage and enumerates the commands
get_command_options = function(usage){
  m = gregexpr('-+\\w+', usage)
  arguments = unlist(regmatches(usage, m))
  names(arguments) = gsub('-+', '', arguments)
  arguments
}
