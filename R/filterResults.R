#PrimerTree
#Copyright (C) 2013 Jim Hester


#' Filter out sequences retrieved by search_primer_pair() that are either too
#' short or too long. The alignment and tree will be recalculated after removing
#' unwanted reads.
#' @param x a primerTree object
#' @return a primerTree object
#' @examples
#' 
#' # filter out sequences longer or shorter than desired:
#' mammals_16S_filtered <- filter_seqs(mammals_16S, min_length=131, max_length=156)
#' 
#' @export

filter_seqs = function(x, ...) UseMethod("filter_seqs")

#' Filter out sequences retrieved by search_primer_pair() that are either too
#' short or too long. The alignment and tree will be recalculated after removing
#' unwanted reads.
#' @param primertree_obj a primerTree object
#' @param min_length the minimum sequence length to keep
#' @param max_length the maximum sequence length to keep
#' @return a primerTree object
#' @export
 

filter_seqs.primerTree = function(primertree_obj, min_length = 0, max_length = Inf ) {
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

#' Get sequence lengths from a primerTree object
#' Also returns a summary of the sequence lengths 
#' @param primertree_obj a primerTree object
#' @param summarize a logical indicating if a summary should be displayed
#' @export
 
table <- function(...) UseMethod("table")
table.default <- base::table

#' Get sequence lengths from a primerTree object
#' Also returns a summary of the sequence lengths 
#' @param primertree_obj a primerTree object
#' @param summarize a logical indicating if a summary should be displayed
#' @return sequence lengths
#' @examples 
#' 
#' # Show the counts for each length
#' table(mammals_16S)
#' 
#' # Plot the distribution of lengths
#' seqLengths <- table(mammals_16S)
#' barplot(seqLengths)
#' 
#' @export

table.primerTree <- function(primertree_obj, summarize = TRUE) {
  lengths <- unlist(lapply(primertree_obj$sequence, length))
  message("sequence length distribution: ")
  table(as.factor(lengths))
}
