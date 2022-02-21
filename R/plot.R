#PrimerTree
#Copyright (C) 2013 Jim Hester

#' plots a tree along with a series of taxonomic ranks
#' @param ranks The ranks to include, defaults to all common ranks, if null
#' print all ranks.
#' @inheritParams plot_tree
#' @seealso \code{\link{plot_tree}} to plot only a single rank or the just the
#' tree layout.
#' @export
#' @examples
#' library(gridExtra)
#' library(directlabels)
#' #plot all the common ranks
#' plot_tree_ranks(mammals_16S$tree, mammals_16S$taxonomy)
#' #plot specific ranks, with a larger dot size
#' plot_tree_ranks(mammals_16S$tree, mammals_16S$taxonomy,
#'   ranks=c('kingdom', 'class', 'family'), size=3)

plot_tree_ranks = function(tree, taxonomy, main=NULL, type='unrooted',
                      ranks=common_ranks,  size=2,
                                    guide_size=NULL, legend_cutoff=25, ...){
  if(is.null(ranks))
    ranks = setdiff(names(taxonomy), c('accession', 'taxId'))

  plots = list()
  plots$structure = plot_tree(tree, main=main, guide_size=guide_size, type=type, ...)

  for(rank in intersect(ranks, names(taxonomy))){
    if(length(na.omit(taxonomy[rank])) > 0){
      plots[[rank]]= plot_tree(tree, guide_size=guide_size, type=type, rank=rank, taxonomy=taxonomy, size=size, legend_cutoff=legend_cutoff, ...)
    }
  }
  p = do.call(arrangeGrob, plots)
  class(p) <- c("primerTree_plot_multi", class(p))

  p
}
common_ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#' plots a tree, optionally with colored and labeled points by taxonomic rank
#'
#' @param tree to be plotted, use layout_tree to layout tree.
#' @param taxonomy A data.frame with an accession field corresponding to the
#'        tree tip labels.
#' @param main An optional title for the plot
#' @param type The type of tree to plot, default unrooted.
#' @param rank The rank to include, if null only the tree is plotted
#' @param size The size of the colored points
#' @param guide_size The size of the length guide.  If NULL auto detects a
#'        reasonable size.
#' @param legend_cutoff The number of different taxa names after which the
#'        names are no longer printed.
#' @param ... additional arguments passed to \code{\link{layout_tree_ape}}
#' @return plot to be printed.
#' @export

plot_tree = function(tree, type='unrooted', main=NULL, guide_size=NULL,
                     rank=NULL, taxonomy=NULL,  size=2, legend_cutoff=25, ...){

  x = layout_tree_ape(tree, type=type, ...)

  range_x = range(x$edge$x, x$tip$x)
  range_y = expand_range(range(x$edqe$y, x$tip$y), mul=.1)

  if(is.null(guide_size)){
    guide_size = 10**(round_any(log10(range_x[2]-range_x[1]), 1)-1)
  }
  p = ggplot() +
   geom_segment(data=x$edge, aes_string(x='x', y='y', xend='xend', yend='yend')) +
    theme_noaxis() +
    annotate('segment', x=range_x[1], xend=range_x[1]+guide_size,
             y=range_y[1], yend=range_y[1],
             arrow=arrow(ends="both",angle=90, length=unit(.2,"cm"))) +
    annotate('text', x=range_x[1]+(guide_size/2), y=range_y[1], label=guide_size, vjust=-.5)

  if(!is.null(rank)){
    if(is.null(taxonomy))
      stop('Must provide a taxonomy if plotting a rank')

    if(is.null(main))
      main = rank

      x$tip = merge(x$tip, taxonomy, by.x='label', by.y='accession', all.x=T)

      rows = na.omit(x$tip[, c('x','y',rank)])
      p = p+geom_point(data=rows, aes_string(x='x', y='y', color=rank),
                     size=size, na.rm=T) + theme(legend.position='none')

        if(length(unique(rows[[rank]])) < legend_cutoff){
          smart.grid2 = list('get.means', 'calc.boxes', 'empty.grid')
          p = p + geom_dl(data=rows, method=smart.grid2,
                    aes_string(x='x', y='y', color=rank, label=rank))
        }
  }
  p = p + ggtitle(main)
  class(p) = c(class(p), 'primerTree_plot')
  p
}

#' layout a tree using ape, return an object to be plotted by
#' \code{\link{plot_tree}}
#' @param tree The \code{\link{phylo}} tree to be plotted
#' @param ... additional arguments to \code{\link{plot.phylo}}
#' @return \item{edge}{list of x, y and xend, yend coordinates
#' as well as ids for the edges}
#' \item{tips}{list of x, y, label and id for the tips}
#' \item{nodes}{list of x, y and id for the nodes}
layout_tree_ape = function(tree, ...){
  #hack to write no output
  cur_dev = dev.cur() #store previous dev
  temp_file = tempfile()
  pdf(file=temp_file)
  plot.phylo(tree, plot=F, ...)
  dev.off()
  unlink(temp_file)
  dev.set(cur_dev) #restore previous dev

  last = .PlotPhyloEnv$last_plot.phylo
  new = list()
  new$edge$x = last$xx[last$edge[,1]]
  new$edge$xend = last$xx[last$edge[,2]]

  new$edge$y = last$yy[last$edge[,1]]
  new$edge$yend = last$yy[last$edge[,2]]
  new$edge$id = tree$edge

  new$edge = data.frame(new$edge, stringsAsFactors=F)

  new$tip$x = last$xx[1:last$Ntip]
  new$tip$y = last$yy[1:last$Ntip]
  new$tip$label = tree$tip.label
  new$tip$id = 1:last$Ntip

  new$tip = data.frame(new$tip, stringsAsFactors=F)

  new$node$x = last$xx[(last$Ntip + 1):length(last$xx)]
  new$node$y = last$yy[(last$Ntip + 1):length(last$yy)]
  new$node$label = tree$node.label
  new$node = data.frame(new$node, stringsAsFactors=F)

  new
}
theme_noaxis = function(){
  theme(panel.border=element_blank(), panel.grid=element_blank(),
        axis.line=element_blank(), axis.text=element_blank(),
        axis.title=element_blank(), axis.ticks=element_blank(),
        plot.margin=unit(c(0, 0, -1, -1), 'lines'))
}

#' @export
print.primerTree_plot_multi <- function(x, ...) {
   grid.arrange(x, ...)
}
