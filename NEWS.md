# primerTree 1.1.0

* Updated code to comply with current CRAN standards
* Fixed missing package anchors in roxygen \link{} sections for external packages
* All code formatted to tidyverse style guide standards using the air package
* Added package prefix to all non-base R functions
* Updated internal function to use httr2 instead of httr due to NCBI eutils API bug, which adds httr2 as a dependency

# primerTree 1.0.7

* Couple of internal changes to fix CRAN issue
* Also updated documentation for filter_seqs

# primerTree 1.0.6

* Changed how matches are handled internally to use accession numbers instead of
  gi numbers.

# primerTree 1.0.5

* `search_primer_pair()` gains a `api_key` parameter, or you can set `NCBI_API_KEY` as an environment variable (@MVesuviusC, #36)
* Add RCurl to Imports, it is an implicit dependency due to a S4 Class being
  included in the data for the package and CRAN is now checking that packages
  load without Suggested packages.

# primerTree 1.0.4

* Add RCurl to Suggests, it is an implicit dependency due to a S4 Class being
  included in the data for the package.

* Added a `NEWS.md` file to track changes to the package.
