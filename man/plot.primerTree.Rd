\name{plot.primerTree}
\alias{plot.primerTree}
\title{plot function for a primerTree object, calls plot_tree_ranks}
\usage{
  plot.primerTree(x, ranks = NULL, ...)
}
\arguments{
  \item{primerTree}{object to plot}

  \item{type}{The type of tree to plot, default unrooted.}

  \item{ranks}{The ranks to include, defaults to all common
  ranks, if null print all ranks.}

  \item{size}{The size of the colored points}

  \item{guide_size}{The size of the length guide.  If NULL
  auto detects a reasonable size.}

  \item{legend_cutoff}{The number of different taxa names
  after which the names are no longer printed.}

  \item{...}{additional arguments passed to
  plot_tree_ranks}
}
\description{
  plot function for a primerTree object, calls
  plot_tree_ranks
}

