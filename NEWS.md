# primerTree 1.0.5

* `search_primer_pair()` gains a `api_key` parameter, or you can set `NCBI_API_KEY` as an environment variable (@MVesuviusC, #36)
* Add RCurl to Imports, it is an implicit dependency due to a S4 Class being
  included in the data for the package and CRAN is now checking that packages
  load without Suggested packages.

# primerTree 1.0.4

* Add RCurl to Suggests, it is an implicit dependency due to a S4 Class being
  included in the data for the package.

* Added a `NEWS.md` file to track changes to the package.
