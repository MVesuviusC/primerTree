#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Retrieve the taxonomy information from NCBI for a set of nucleotide gis.
#'
#' @param accessions a character vector of the accessions to retrieve
#' @return data.frame of the 'accessions, taxIds, and taxonomy
#' @export

get_taxonomy <- function(accessions) {
    accessions <- unique(as.character(accessions))
    taxids <- accession2taxid(accessions)

    taxonomy <- fetch_taxonomy(unique(taxids))
    merge(
        data.frame(
            accession = names(taxids),
            taxId = taxids,
            stringsAsFactors = FALSE
        ),
        taxonomy
    )
}
#' Maps a nucleotide database accession to a taxonomy database taxId
#'
#' @param accessions accessions character vector to lookup.
#' @return named vector of taxIds.
#' @export

accession2taxid <- function(accessions) {
    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    taxids <-
        httr2::request(url) |>
        httr2::req_method("POST") |>
        httr2::req_body_form(
            db = "nuccore",
            id = paste0(accessions, collapse = ",")
        ) |>
        httr2::req_retry(max_tries = 5) |>
        httr2::req_perform() |>
        httr2::resp_body_string() |>
        XML::xmlParse() |>
        XML::xpathSApply(
            "//DocSum",
            parse_docsum
        )

    taxids
}
parse_docsum <- function(docsum) {
    gid <- XML::xpathSApply(
        docsum,
        ".//Item[@Name='AccessionVersion']", XML::xmlValue
    )
    taxid <- XML::xpathSApply(
        docsum,
        ".//Item[@Name='TaxId']", XML::xmlValue
    )
    if (length(taxid) != 1) {
        taxid <- NA
    }
    names(taxid) <- gid
    taxid
}

fetch_taxonomy <- function(taxid) {
    fetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    query <- list(
        db = "taxonomy",
        rettype = "null",
        retmode = "xml",
        id = paste(taxid, collapse = ",")
    )

    response <- POST_retry(fetch_url, body = query)

    #stop if response failed
    stop_for_status(response)

    parse_taxonomy_xml(XML::xmlParse(content(response, as = "text")))
}
parse_taxonomy_xml <- function(xml) {
    plyr::rbind.fill(XML::xpathApply(xml, "//TaxaSet/Taxon", parse_taxon))
}
parse_taxon <- function(taxon) {
    tax_id <- XML::xpathSApply(taxon, "./TaxId", xmlValue)
    ranks <- XML::xpathSApply(taxon, ".//Rank", xmlValue)
    names <- XML::xpathSApply(taxon, ".//ScientificName", xmlValue)
    names(names) <- ranks
    names <- names[ranks != "no rank"]
    names["taxId"] <- tax_id
    as.data.frame(t(names), stringsAsFactors = FALSE)
}
