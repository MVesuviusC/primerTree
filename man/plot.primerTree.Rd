\name{plot.primerTree}
\alias{plot.primerTree}
\title{plot function for a primerTree object, calls plot_tree_ranks}
\usage{
  \method{plot}{primerTree} (x, ranks = NULL, ...)
}
\arguments{
  \item{x}{primerTree object to plot}

  \item{ranks}{The ranks to include, defaults to all common
  ranks, if NULL print all ranks.  If 'none' just print the
  layout.}

  \item{...}{additional arguments passed to
  plot_tree_ranks}
}
\description{
  plot function for a primerTree object, calls
  plot_tree_ranks
}
\examples{
library(gridExtra)
library(directlabels)
#plot with all common ranks
plot(mammals_16S)

#plot only the class
plot(mammals_16S, 'class')

#plot the layout only
plot(mammals_16S, 'none')
}
\seealso{
  \code{\link{plot_tree_ranks}}, \code{\link{plot_tree}}
}

