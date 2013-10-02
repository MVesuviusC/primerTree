\name{layout_tree_ape}
\alias{layout_tree_ape}
\title{layout a tree using ape, return an object to be plotted by
\code{\link{plot_tree}}}
\usage{
  layout_tree_ape(tree, ...)
}
\arguments{
  \item{tree}{The \code{\link{phylo}} tree to be plotted}

  \item{...}{additional arguments to
  \code{\link{plot.phylo}}}
}
\value{
  \item{name edge}{description list of x, y and xend, yend
  coordinates as well as ids for the edges} \item{name
  tips}{description list of x, y, label and id for the
  tips} \item{name nodes}{description list of x, y and id
  for the nodes}
}
\description{
  layout a tree using ape, return an object to be plotted
  by \code{\link{plot_tree}}
}

