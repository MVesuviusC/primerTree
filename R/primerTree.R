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
#' @import ggplot2 XML ape httr plyr directlabels gridExtra
#'   stringr foreach
#' @importFrom lubridate %--% seconds now
#' @importFrom grid grid.locator
#' @importFrom scales expand_range
#' @importFrom grDevices dev.cur dev.off dev.set pdf
#' @importFrom stats na.omit quantile
#' @importFrom utils capture.output
#' @useDynLib primerTree
NULL

#' PrimerTree results for the mammalian 16S primers
#' @name mammals_16S
#' @docType data
NULL

#' PrimerTree results for the bryophyte trnL primers
#' @name bryophytes_trnL
#' @docType data
NULL

#' @export
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
#'   print all ranks.  If 'none' just print the layout.
#' @param main an optional title to display, if NULL displays the name as the title
#' @param ... additional arguments passed to plot_tree_ranks
#' @export
#' @seealso \code{\link{plot_tree_ranks}}, \code{\link{plot_tree}}
#' @examples
#' library(gridExtra)
#' library(directlabels)
#' #plot with all common ranks
#' plot(mammals_16S)
#'
#' #plot only the class
#' plot(mammals_16S, 'class')
#'
#' #plot the layout only
#' plot(mammals_16S, 'none')
plot.primerTree = function(x, ranks=NULL, main=NULL, ...){
  if(is.null(ranks)){
    if(is.null(main))
      main = x$name
    plot_tree_ranks(x$tree, x$taxonomy, main=main, ...)
  }
  else if(length(ranks) > 1){
    if(is.null(main))
      main = x$name
    plot_tree_ranks(x$tree, x$taxonomy, ranks=ranks, main=main, ...)
  }
  else {
    if(ranks == 'none') {
      plot_tree(x$tree, main=main, ...)
    }
    else{
      plot_tree(x$tree, taxonomy=x$taxonomy, rank=ranks, main=main, ...)
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
#' @param clustal_options a list of options to pass to clustal omega, see
#'    \code{link{clustalo}} for a list of options
#' @param distance_options a list of options to pass to dist.dna, see
#'    \code{link{dist.dna}} for a list of options
#' @inheritParams primer_search
#' @return A list with the following elements,
#' \item{name}{name of the primer pair}
#' \item{BLAST_result}{html blast results from Primer-BLAST as
#'  'a \code{\link{response}}} object.
#' \item{taxonomy}{taxonomy for the primer products from NCBI}
#' \item{sequence}{sequence of the primer products}
#' \item{alignment}{multiple alignment of the primer products}
#' \item{tree}{phylogenetic tree of the reconstructed from the
#' 'multiple alignment}
#' @seealso \code{\link{primer_search}}, \code{\link{clustalo}}
#' @export
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
    num_permutations=25, simplify=TRUE, clustal_options=list(), distance_options=list(model="N", pairwise.deletion=T),
    ..., .parallel=FALSE, .progress='none'){

  #HACK, primerTree is an environment rather than a list so we can treat it as
  #a pointer, I could make it a reference class, but that seems to be overkill
  #as I am converting to a list at the end of the function anyway...

  if(missing(forward) || missing(reverse)){
    BLAST_primer()
    return()
  }

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
      primer_search$distances = do.call(dist.dna, c(list(primer_search$alignment, distance_options)))
      message('pairwise DNA distances calculated in ',
              seconds_elapsed_text(start_time))

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
#' identify the point closest to the mouse click
#' only works on single ranks
#' @param x the plot to identify
#' @param ... additional arguments passed to annotate
#' @export
identify.primerTree_plot = function(x, ...) {
  point <- gglocator(x$layers[[4]])
  distances <- distance(point, x$layers[[4]]$data[,c('x','y')])
  closest <- which(distances == min(distances))[1]
  point$label <- x$layers[[4]]$data[closest,deparse(x$layers[[4]]$mapping$colour)]
  x + annotate("text", label=point$label, x=point$x, y=point$y, ...)
}
gglocator = function(object) {
  loc <-  as.numeric(grid.locator("npc"))

  xrng <- with(object, range(data[,deparse(mapping$x)]))
  yrng <- with(object, range(data[,deparse(mapping$y)]))

  point <- data.frame(xrng[1] + loc[1]*diff(xrng), yrng[1] + loc[2]*diff(yrng))
  names(point) <- with(object, c(deparse(mapping$x), deparse(mapping$y)))
  point
}

#returns the distance from a point in point to the points in points
distance <- function(point,points){
  sqrt((point$x-points$x)^2 + (point$y-points$y)^2)
}

#' Summarize a primerTree result, printing quantiles of sequence length and
#' pairwise differences.

#' @param object the primerTree object to summarise
#' @param ... Ignored options
#' @param probs quantile probabilities to compute, defaults to 0, 5, 50, 95,
#' and 100 probabilities.
#' @param ranks ranks to show unique counts for, defaults to the common ranks
#' @return invisibly returns a list containing the printed results
#' @export
summary.primerTree <- function(object, ..., probs=c(0, .05, .5, .95, 1), ranks = common_ranks) {

  res = list()
  res[['lengths']] = t(data.frame('Sequence lengths'=labeled_quantile(laply(object$sequence, length), sprintf('%.0f%%', probs*100), probs=probs), check.names=F))
  print(res[['lengths']])

  res[['distances']] = t(data.frame('Pairwise differences'=labeled_quantile(object$distance, sprintf('%.0f%%', probs*100), probs=probs), check.names=F))
  cat('\n')
  print(res[['distances']])
  
  res[['rankDistances']] = calc_rank_dist_ave(object, common_ranks)
  #cat('\n', 'Median pairwise differences for taxonomic levels')
  print(res[['rankDistances']])

  res[['ranks']] = laply(object$taxonomy[common_ranks], function(x) length(unique(x)))
  cat('\n', 'Unique taxa out of ', length(object$sequence), ' sequences\n', sep='')
  names(res[['ranks']]) = ranks
  print(res[['ranks']])

  invisible(res)
}

labeled_quantile = function(x, labels, ...){
  res = quantile(x, ...)
  names(res) = labels
  res
}

#' @title calc_rank_dist_ave
#' Summarize 
#' pairwise differences.

#' @param x a primerTree object 
#' @param ... Ignored options
#' @param ranks ranks to show unique counts for, defaults to the common ranks
#' @return returns a data frame of results
#' @details 
#' The purpose of this function is to calculate the average number
#' of nucleotide differences between species within a given taxonomic 
#' level. 
#' 
#' There are several key assumptions and calculations made in this 
#' function.
#' 
#' First, the function randomly selects one sequence from each species
#' in the primerTree results. This is to keep any one species (e.g.
#' human, cow, etc.) with many hits from skewing the results.
#' 
#' Second, for each taxonomic level tested, the function divides the
#' sequences by each taxon at that level and calculates the mean 
#' number of nucleotide differences within that taxa, then returns the
#' mean of those values. 
#' 
#' Third, when calculating the average distance, any taxa for which 
#' there is only one species is omitted, as the number of nucleotide
#' differences will always be 0.
#' 
#' @examples
#' \dontrun{
#' calc_rank_dist_ave(mammals_16S)
#' 
#' calc_rank_dist_ave(bryophytes_trnL)
#' 
#' # Note that the differences between the results from these two primers
#' # the mean nucleotide differences is much higher for the mammal primers
#' # than the byrophyte primers. This suggests that the mammal primers have
#' # better resolution to distinguish individual species.
#' }
#' @export

# using tree data format info from http://www.phytools.org/eqg/Exercise_3.2/ 
calc_rank_dist_ave <- function(x, ranks = common_ranks, ...) {
  usedRanks <- grep("species", ranks, invert = T, value = T)
  rankDistMean <- data.frame(matrix(nrow=1,ncol=0))
  for(rank in usedRanks) {
    
    # Raw taxonomy data
    taxa <- as.data.frame(x$taxonomy)
    
    # Randomize the order of the taxa data frame for the next step
    taxa <- taxa[ sample(1:nrow(taxa), size = nrow(taxa), replace = F) ,]
    rownames(taxa) <- taxa$gi
    
    # Pick random example per species and add back in the taxa info
    uniqueFactors <- as.data.frame(unique( taxa$species ))
    colnames(uniqueFactors) <- "species"
    uniqueFactors <- join(uniqueFactors, taxa, type = "left",match="first", by = "species")
    uniqueFactors <- uniqueFactors[, colnames(uniqueFactors) %in% c("gi","species",rank)]
    
    # Drop any row in (rank) where there is only one species represented
    # Any instance of this leads to a distance within that rank of 0, skewing the results downward
    counts <- as.data.frame(table(uniqueFactors[[rank]]))
    colnames(counts) <- c(rank,"count")
    uniqueFactors <- join(uniqueFactors, counts, by = rank)
    uniqueFactors <- subset(uniqueFactors, count > 1)
    
    # Get sequences for randomly selected species
    seqs <- x$sequence
    seqs <- seqs[names(seqs) %in% uniqueFactors$gi ]
    
    # Align and calculate pairwise distances and convert dists to dataframe
    align <- clustalo(seqs)
    dists <- as.data.frame(as.matrix(dist.dna(align, model="N", pairwise.deletion=T)))
    dists$gi <- row.names(dists)
    
    # Melt the dists dataframe so I can drop any distance that isn't within the (rank)
    melted <- melt(dists,id="gi",variable.name = "gi2", value.name = "dist")
    
    # Replace the rank1 gi with the rank1 taxa
    #melted2 <- join(melted, uniqueFactors[,colnames(uniqueFactors) %in% c("gi", rank)],by = "gi")
    melted <- join(melted, uniqueFactors,by = "gi")
    melted$rank1 <- melted[[rank]]
    
    # Drop all columns except the three needed so the next join doesn't get messed up
    melted <- melted[,colnames(melted) %in% c("gi2", "dist", "rank1", "species")]
    colnames(melted)[1] <- "gi"
    
    # Replace the rank2 gi with the rank2 taxa
    melted <- join(melted, uniqueFactors, by = "gi")
    melted$rank2 <- melted[[rank]]
    
    # Drop all columns except the three needed
    melted <- melted[,colnames(melted) %in% c("rank2", "dist", "rank1","species")]
    
    # We only want distances within a taxa, so drop all comparisons between taxa 
    # We also want to drop any comparisons of a species to itself, which will have dist == 0
    melted <- subset(melted, rank1 == rank2 & species != species.1)
    
    # Calculate the median distance for each taxa compared
    #   We calculate each separately to avoid any one taxa with lots of hits (like human seqs) 
    #   from skewing the median
    melted$group <- paste(melted$rank1,melted$rank2)
    perGroup <- ddply(melted,.(group), summarize, median = median(dist))
    
    # Plug the medians into the storage dataframe
    rankDistMean[[rank]] <- mean(melted$dist)
  }
  message("\nAverage number of nucleotide differences between sequences within a given taxonomic group")
  message("See function description for further details")
  return(rankDistMean)
}

