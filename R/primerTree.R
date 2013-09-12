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
#' @param forward forward primer to search by 5'-3'
#' @param reverse reverse primer to search by 3'-5'
#' @param name name to give to the primer pair
#' @param num_aligns name to give to the primer pair
#' @param simplify use simple names for primer hit results or complex
#' @param the number of primer permutations to search, if the degenerate bases
#'        cause more than this number of permutations to exist, this number will be
#'        sampled from all possible permutations.
#' @param ... additional arguments for primer_search, run primer_search with no
#'        arguments to see all available options
#' @param .parallel use the parallel foreach, see foreach and plyr
#' @export search_primer
search_primer = function(forward, reverse, name=NULL, num_aligns=500, num_permutations=25, simplify=TRUE, ...,
                         .parallel=FALSE){

  #HACK, primerTree is an environment rather than a list so we can treat it as a pointer

  primer_search = new.env(parent=globalenv())
  class(primer_search) = 'primerTree'
  #list all primers used to search
  env2list(
    try_default({
      primer_search$name = if(!is.null(name)) name
                    else name=paste(forward, reverse, sep='_')

      primer_search$primers = enumerate_primers(forward, reverse)
      if(nrow(primer_search$primers) > num_permutations){
        message('sampling ', num_permutations, ' primers from ', nrow(primer_search$primers), ' possible combinations')
        primer_search$primers = primer_search$primers[ sample.int(nrow(primer_search$primers), num_permutations, replace=F), ]
      }
      message('blasting ', nrow(primer_search$primers), ' primer combinations')
      primer_search$response = primer_search(forward, reverse, ...,
                                      .parallel=.parallel,
                                      num_targets_with_primers=num_aligns %/%
                                      nrow(primer_search$primers))

      primer_search$blast_result =
        filter_duplicates(ldply(primer_search$response, parse_primer_hits, .parallel=.parallel))

      message(nrow(primer_search$blast_result), ' blast alignments parsed')

      primer_search$sequence = get_sequences(primer_search$blast_result$gi,
                                             primer_search$blast_result$product_start,
                                             primer_search$blast_result$product_stop,
                                             simplify=simplify,
                                             .parallel=.parallel)

      lengths = laply(primer_search$sequence, length)
      message(length(primer_search$sequence), ' sequences retrieved from NCBI,',
              ' min:', min(lengths),
              ' mean:', mean(lengths),
              ' max:', max(lengths))

      primer_search$alignment = clustalo(primer_search$sequence, threads=getDoParWorkers())
      message(nrow(primer_search$alignment), ' sequences aligned',
              ' length:', ncol(primer_search$alignment))

      primer_search$tree = tree_from_alignment(primer_search$alignment)
      message('constructed neighbor joining tree')

      primer_search$taxonomy = get_taxonomy(primer_search$blast_result$gi)
      primer_search
    }, default=primer_search)
  )
}

env2list = function(env){
  names = ls(env)
  mget(names, env)
}
