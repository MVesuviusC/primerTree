\name{plot_tree_ranks}
\alias{plot_tree_ranks}
\title{plots a tree with colored and labeled points}
\usage{
  plot_tree_ranks(tree, taxonomy, main = NULL,
    type = "unrooted", ranks = common_ranks, size = 2,
    guide_size = NULL, legend_cutoff = 25, ...)
}
\arguments{
  \item{tree}{to be plotted, use layout_tree to layout
  tree.}

  \item{taxonomy}{A data.frame with an accession field
  corresponding to the tree tip labels.}

  \item{main}{An optional title for the plot}

  \item{type}{The type of tree to plot, default unrooted.}

  \item{ranks}{The ranks to include, defaults to all common
  ranks, if null print all ranks.}

  \item{size}{The size of the colored points}

  \item{size}{The size of the colored points}

  \item{guide_size}{The size of the length guide.  If NULL
  auto detects a reasonable size.}

  \item{legend_cutoff}{The number of different taxa names
  after which the names are no longer printed.}
}
\value{
  plot to be printed.
}
\description{
  plots a tree with colored and labeled points
}

