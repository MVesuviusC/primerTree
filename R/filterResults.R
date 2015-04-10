#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Filter out sequences retrieved by search_primer_pair() that are either too
#' short or too long. The alignment and tree will be recalculated after removing
#' unwanted reads.
#' @param primerTreeObj a primerTree object
#' @param maxLength the maximum sequence length to keep
#' @param minLength the minimum sequence length to keep
#' @return a primerTree object
#' @export
filter_PrimerTree_results = function(primerTreeObj, maxLength=99999, minLength=1) {
  #if(is.null(primerTreeObj$sequence)) {
  #  print "No sequences in primerTree object"
  #}
  
  # calculate how many removed and print out to screen
  numKickedOutBig <- sum(as.logical(lapply(primerTreeObj$sequence,function(x) length(x)<maxLength))=="FALSE")
  message(paste(numKickedOutBig, "sequences greater than", maxLength, "bases in length" ) ) 
  numKickedOutSmall <- sum(as.logical(lapply(primerTreeObj$sequence,function(x) length(x)>minLength))=="FALSE")
  message(paste(numKickedOutSmall, "sequences less than", minLength, "bases in length" ) )
    
  # get rid of long sequences that may mess up the alignment
  primerTreeObj$sequence <- primerTreeObj$sequence[as.logical(lapply(primerTreeObj$sequence,function(x) length(x)<maxLength))] 
  
  # get rid of short sequences that may mess up the alignment
  primerTreeObj$sequence <- primerTreeObj$sequence[as.logical(lapply(primerTreeObj$sequence,function(x) length(x)>minLength))] 
  # realign filtered sequences
  primerTreeObj$alignment <- clustalo(primerTreeObj$sequence)
  # remake tree
  primerTreeObj$tree <- tree_from_alignment(primerTreeObj$alignment)
  return(primerTreeObj)
}

#' Get sequence lengths from a primerTree object
#' Also returns a summary of the sequence lengths 
#' @param primerTreeObj a primerTree object
#' @param summarize a logical indicating if a summary should be displayed
#' @return numeric values of sequence lengths
#' @export
get_PrimerTree_Seq_Lengths <- function(primerTreeObj, summarize = TRUE) {
  lengths<-as.numeric(lapply(primerTreeObj$sequence,function(x) length(x)))
  message("sequence length distribution: ")
  print(table(as.factor(lengths)))
}
