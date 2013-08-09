as.DNAbin.primer = function(x){
  x$dna
}
as.character.primer = function(x){
  as.character(x$dna)
}
length.primer = function(x){
  length(x$dna)
}
content.primer = function(x, ...){
  content(x$response, ...)
}
plot.primer = function(x, ...){
  plot_tree_ranks(x$tree, x$taxonomy, x$name, ...)
}
#' @export search_primer
search_primer = function(forward, reverse, name=NULL, ..., .parallel=FALSE, simplify=TRUE){
  primer = list()
  primer$response = primer_search(forward, reverse, ..., .parallel=.parallel)
  primer$blast_result = deduplicate(ldply(primer$response, parse_primer_hits, .parallel=.parallel))
  primer$blast_result
  message(nrow(primer$blast_result), ' blast alignments parsed')
  primer$sequence = get_sequences(primer$blast_result$gi, primer$blast_result$start, primer$blast_result$stop, simplify=simplify, .parallel=.parallel)
  lengths = laply(primer$sequence, length)
  message(length(primer$sequence), ' sequences retrieved from NCBI,',
          ' min:', min(lengths),
          ' mean:', mean(lengths),
          ' max:', max(lengths))

  primer$alignment = clustal(primer$sequence)
  message(nrow(primer$alignment), ' sequences aligned',
          ' length:', ncol(primer$alignment))

  primer$tree = tree_from_alignment(primer$alignment)
  message('constructed neighbor joining tree')

  primer$taxonomy = get_taxonomy(primer$blast_result$gi)
  primer$name = if(!is.null(name)) name else name=paste(forward, reverse, sep='_')
  class(primer) = 'primer'
  primer
}
