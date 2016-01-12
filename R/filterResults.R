#PrimerTree
#Copyright (C) 2013 Jim Hester


#' Filter out sequences retrieved by search_primer_pair() that are either too
#' short or too long. The alignment and tree will be recalculated after removing
#' unwanted reads.
#' @param x a primerTree object
#' @param min_length the minimum sequence length to keep
#' @param max_length the maximum sequence length to keep
#' @param ... additional arguments passed to methods.
#' @return a primerTree object
#' @examples
#' \dontrun{
#' # filter out sequences longer or shorter than desired:
#' mammals_16S_filtered <- filter_seqs(mammals_16S, min_length=131, max_length=156)
#' }
#' @export
filter_seqs = function(x, ...) UseMethod("filter_seqs")

#' @describeIn filter_seqs Method for primerTree objects
#' @export
filter_seqs.primerTree = function(x, min_length = 0, max_length = Inf, ...) {
  #if(is.null(primerTreeObj$sequence)) {
  #  print "No sequences in primerTree object"
  #}
  
  # calculate how many removed and print out to screen
  lengths <- unlist(lapply(primertree_obj$sequence, length))
  above_max <- lengths >= max_length
  below_min <- lengths <= min_length
  message(sum(below_min), " sequences below ", min_length) # same for min
  message(sum(above_max), " sequences above ", max_length) # same for min
  
  #filter the sequences
  primertree_obj$sequence <- primertree_obj$sequence[lengths > min_length & lengths < max_length]
      
  # realign filtered sequences
  primertree_obj$alignment <- clustalo(primertree_obj$sequence)
  
  # remake tree
  primertree_obj$tree <- tree_from_alignment(primertree_obj$alignment)
  primertree_obj
}

#' Get a summary of sequence lengths from a primerTree object
#' @param x a primerTree object. 
#' @param summarize a logical indicating if a summary should be displayed
#' @return a table of sequence length frequencies
#' @examples 
#' 
#' # Show the counts for each length
#' seq_lengths(mammals_16S)
#' 
#' # Plot the distribution of lengths
#' seqLengths <- seq_lengths(mammals_16S)
#' barplot(seqLengths, 
#'  main = "Frequency of sequence lengths for 16S mammal primers", 
#'  xlab="Amplicon length (in bp)", 
#'  ylab=("Frequency"))
#' @export
seq_lengths <- function(x,summarize = TRUE) UseMethod("seq_lengths")

#' Method for primerTree objects
#' @inheritParams seq_lengths
#' @export
seq_lengths.primerTree <- function(x, summarize = TRUE) {
  lengths <- unlist(lapply(x$sequence, length))
  message("sequence length distribution: ")
  table(as.factor(lengths))
}
