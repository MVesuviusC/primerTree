#' accessors for various objects
#' @export as.DNAbin.primerTree as.character.primerTree length.primerTree content.primerTree
names.primerTree = function(x){
  x$name
}
as.DNAbin.primerTree = function(x){
  x$dna
}
as.character.primerTree = function(x){
  as.character(x$dna)
}
length.primerTree = function(x){
  length(x$dna)
}
content.primerTree = function(x, ...){
  content(x$response, ...)
}
#' plot function for a primerTree object, calls plot_tree_ranks
#' @param primerTree object to plot
#' @param type The type of tree to plot, default unrooted.
#' @param ranks The ranks to include, defaults to all common ranks, if null print all ranks.
#' @param size The size of the colored points
#' @param guide_size The size of the length guide.  If NULL auto detects a
#'        reasonable size.
#' @param legend_cutoff The number of different taxa names after which the
#'        names are no longer printed.
#' @param ... additional arguments passed to plot_tree_ranks
#' @export plot.primerTree
plot.primerTree = function(x, ...){
  plot_tree_ranks(x$tree, x$taxonomy, x$name, ...)
}
#' Automatic primer searching'
#' @param name name to give to the primer pair
#' @param simplify use simple names for primer hit results or complex
#' @param .progress name of the progress bar to use, see 'plyr::create_progress_bar'
#' @inheritParams primer_search
#' @export search_primer_pair
search_primer_pair = function(forward, reverse, name=NULL, num_aligns=500, num_permutations=25, simplify=TRUE, ...,
                         .parallel=FALSE, .progress='none'){

  #HACK, primerTree is an environment rather than a list so we can treat it as a pointer

  primer_search = new.env(parent=globalenv())
  #list all primers used to search
  primer_search = env2list(
    try_default({
      primer_search$name = if(!is.null(name)) name
                    else name=paste(forward, reverse, sep='_')

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
              ' in ', seconds_elapsed_text(start_time), ', ',
              ' min:', min(lengths),
              ' mean:', round(mean(lengths),2),
              ' max:', max(lengths))

      start_time = now()
      primer_search$alignment = clustalo(primer_search$sequence, threads=getDoParWorkers())
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

seconds_elapsed_text = function(start_time){
  paste((start_time %--% now()) %/% seconds(1), 'seconds')
}

env2list = function(env){
  names = ls(env)
  mget(names, env)
}
