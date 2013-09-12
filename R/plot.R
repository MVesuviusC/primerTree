#' plots a tree with colored and labeled points
#'
#' @param tree to be plotted, use layout_tree to layout tree.
#' @param taxonomy A data.frame with an accession field corresponding to the
#'        tree tip labels.
#' @param main An optional title for the plot
#' @param type The type of tree to plot, default unrooted.
#' @param ranks The ranks to include, defaults to all common ranks, if null print all ranks.
#' @param size The size of the colored points
#' @param guide_size The size of the length guide.  If NULL auto detects a
#'        reasonable size.
#' @param legend_cutoff The number of different taxa names after which the
#'        names are no longer printed.
#' @param ... additional arguments passed to layout_tree_ape
#' @return plot to be printed.
#' @export plot_tree_ranks

plot_tree_ranks = function(tree, taxonomy, main=NULL, type='unrooted',
                      ranks=common_ranks,  size=2,
                                    guide_size=NULL, legend_cutoff=25, ...){
  if(is.null(ranks))
    ranks = setdiff(names(taxonomy), c('accession', 'gi', 'taxId'))

  x = layout_tree_ape(tree, type=type, ...)

  x$tip = merge(x$tip, taxonomy, by.x='label', by.y='gi', all.x=T)

  range_x = range(x$edge$x, x$tip$x)
  min_y = expand_range(range(x$edge$y, x$tip$y), mul=.1)[1]

  if(is.null(guide_size)){
    guide_size = 10**(round_any(log10(range_x[2]-range_x[1]), 1)-1)
  }
  p = ggplot() +
   geom_segment(data=x$edge, aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_noaxis() +
    annotate('segment', x=range_x[1], y=min_y, xend=range_x[1]+guide_size, yend=min_y, vjust=1, arrow=arrow(ends="both",angle=90,length=unit(.2,"cm"))) +
    annotate('text', x=range_x[1], y=min_y, label=guide_size, vjust=-1)

  plots = list()
  plots$structure = p + ggtitle(main)

  smart.grid2 = list('get.means', 'calc.boxes', 'empty.grid')
  for(rank in intersect(ranks, names(taxonomy))){
    rows = na.omit(x$tip[, c(rank, 'x', 'y')])
    if(nrow(rows) > 0){
      plots[[rank]]=
      p+geom_point(data=rows,
                   aes_string(x='x', y='y', color=rank), size=size, na.rm=T) +
        ggtitle(rank) + scale_x_continuous(expand=c(.1,0)) +
        scale_y_continuous(expand=c(.1, 0)) + theme(legend.position='none')

      if(length(unique(rows[[rank]])) < legend_cutoff){
        plots[[rank]] = plots[[rank]] +
        geom_dl(data=rows, method=smart.grid2,
                aes_string(x='x', y='y', color=rank, label=rank))
      }
    }
  }
  p = do.call(arrangeGrob, plots)
  p
}

common_ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")

layout_tree_ape = function(tree, ...){
  pdf(file='/dev/null') #hack to write no output
  plot.phylo(tree, plot=F, ...)
  dev.off()
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
        axis.title=element_blank(), axis.ticks=element_blank())
}

round_any = function(x, accuracy, f=round){
  f(x/accuracy) * accuracy
}
