#PrimerTree
#Copyright (C) 2013 Jim Hester

#' \pkg{primerTree} Visually Assessing the Specificity and Informativeness of Primer Pairs
#'
#' \code{primerTree} has two main commands:
#' \code{\link{search_primer_pair}} which takes a primer pair and returns an
#' primerTree object of the search results
#' \code{\link{plot.primerTree}} a S3 method for plotting the primerTree object
#' obtained using \code{\link{search_primer_pair}}
#' @name primerTree
#' @docType package
NULL

#' PrimerTree results for the mammalian 16S primers
#' @name mammals_16S
#' @docType data
NULL

#' PrimerTree results for the bryophyte trnL primers
#' @name bryophytes_trnL
#' @docType data
NULL

#' @S3method print primerTree
print.primerTree = function(x, ...){
  cat("Name: ", x$name, "\n",
      "  Arguments: ", paste(names(x$arguments), x$arguments, sep=":", collapse=' '), "\n")
  cat("\nHTTP Response\n")
      print(x$response[[1]])
  cat("\nPrimer Products\n")
      print(x$sequence)
  cat("\nAligned Products\n")
      print(x$alignment)
  cat("\nPhylogenetic Tree\n")
      print(x$tree)
}
#' plot function for a primerTree object, calls plot_tree_ranks
#' @param x primerTree object to plot
#' @param ranks The ranks to include, defaults to all common ranks, if NULL
  #' print all ranks.  If 'none' just print the layout.
#' @param ... additional arguments passed to plot_tree_ranks
#' @method plot primerTree
#' @export plot.primerTree
#' @seealso \code{\link{plot_tree_ranks}}, \code{\link{plot_tree}}
#' @examples
#' #plot with all common ranks
#' plot(mammals_16S)
#'
#' #plot only the class
#' plot(mammals_16S, 'class')
#'
#' #plot the layout only
#' plot(mammals_16S, 'none')
plot.primerTree = function(x, ranks=NULL, ...){
  if(is.null(ranks)){
    plot_tree_ranks(x$tree, x$taxonomy, x$name, ...)
  }
  else if(length(ranks) > 1){
    plot_tree_ranks(x$tree, x$taxonomy, x$name, ranks=ranks, ...)
  }
  else {
    if(ranks == 'none') {
      plot_tree(x$tree, ...)
    }
    else{
      plot_tree(x$tree, taxonomy=x$taxonomy, rank=ranks, ...)
    }
  }
}
#' Automatic primer searching Search a given primer pair, retrieving the alignment
#' results, their product sequences, the taxonomic information for the sequences,
#' a multiple alignment of the products
#' @param name name to give to the primer pair
#' @param simplify use simple names for primer hit results or complex
#' @param .progress name of the progress bar to use, see
#' \code{\link{create_progress_bar}}
#' @param clustal_options a list of options to pass to clustal omega, run
#'    \code{link{clustalo}} for a list of options
#' @inheritParams primer_search
#' @return A list with the following elements,
#' \item{name name}{description name of the primer pair}
#' \item{name BLAST_result}{description html blast results from Primer-BLAST as
#'  'a \code{\link{response}}} object.
#' \item{name taxonomy}{description taxonomy for the primer products from NCBI}
#' \item{name sequence}{description sequence of the primer products}
#' \item{name alignment}{description multiple alignment of the primer products}
#' \item{name tree}{description phylogenetic tree of the reconstructed from the
#' 'multiple alignment}
#' @seealso \code{\link{primer_search}}, \code{\link{clustalo}}
#' @export search_primer_pair
#' @examples
#' \dontrun{
#' #simple search
#' mammals_16S = search_primer_pair(name='Mammals 16S',
#'  'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT')
#' #returning 1000 alignments, allow up to 3 mismatches in primer
#' mammals_16S = search_primer_pair(name='Mammals 16S',
#'  'CGGTTGGGGTGACCTCGGA', 'GCTGTTATCCCTAGGGTAACT',
#'  num_aligns=1000, total_primer_specificity_mismatch=3)
#' }
search_primer_pair = function(forward, reverse, name=NULL, num_aligns=500,
    num_permutations=25, simplify=TRUE, clustal_options=list(), ...,
    .parallel=FALSE, .progress='none'){

  #HACK, primerTree is an environment rather than a list so we can treat it as
  #a pointer, I could make it a reference class, but that seems to be overkill
  #as I am converting to a list at the end of the function anyway...


  primer_search = new.env(parent=globalenv())
  #list all primers used to search
  primer_search = env2list(
    try_default({
      primer_search$name = if(!is.null(name)) name
                    else name=paste(forward, reverse, sep='_')

      primer_search$arguments =
        c(forward=forward, reverse=reverse, name=name,
          num_aligns=num_aligns, num_permutations = num_permutations,
          simplify=simplify, clustal_options=clustal_options, list(...))

      primer_search$response = primer_search(forward, reverse,
                                             num_permutations=num_permutations,
                                             .progress=.progress,
                                             .parallel=.parallel,
                                             num_aligns=num_aligns,
                                             ...)
      start_time = now()
      primer_search$BLAST_result =
        filter_duplicates(ldply(primer_search$response, parse_primer_hits, .parallel=.parallel))

      message(nrow(primer_search$BLAST_result), ' BLAST alignments parsed in ', seconds_elapsed_text(start_time))

      start_time = now()
      primer_search$taxonomy = get_taxonomy(primer_search$BLAST_result$gi)
      message('taxonomy retrieved in ', seconds_elapsed_text(start_time))

      start_time = now()
      primer_search$sequence = get_sequences(primer_search$BLAST_result$gi,
                                             primer_search$BLAST_result$product_start,
                                             primer_search$BLAST_result$product_stop,
                                             simplify=simplify,
                                             .parallel=.parallel)

      lengths = laply(primer_search$sequence, length)
      message(length(primer_search$sequence), ' sequences retrieved from NCBI',
              ' in ', seconds_elapsed_text(start_time), ', product length',
              ' min:', min(lengths),
              ' mean:', round(mean(lengths),2),
              ' max:', max(lengths))

      start_time = now()
      primer_search$alignment = do.call(clustalo, c(list(primer_search$sequence, threads=getDoParWorkers()), clustal_options))
      message(nrow(primer_search$alignment), ' sequences aligned in ',
              seconds_elapsed_text(start_time),
              ' length:', ncol(primer_search$alignment))

      start_time = now()
      primer_search$tree = tree_from_alignment(primer_search$alignment)
      message('constructed neighbor joining tree in ', seconds_elapsed_text(start_time))

      primer_search
    }, default=primer_search)
  )
  class(primer_search) = 'primerTree'
  primer_search
}

#given a start time print out the number of seconds which have elapsed
seconds_elapsed_text = function(start_time){
  paste((start_time %--% now()) %/% seconds(1), 'seconds')
}

#fast way to convert a environment to a list
env2list = function(env){
  names = ls(env)
  mget(names, env)
}
