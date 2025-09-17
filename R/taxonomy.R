#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Retrieve the taxonomy information from NCBI for a set of nucleotide gis.
#'
#' @param accessions a character vector of the accessions to retrieve
#' @return data.frame of the 'accessions, taxIds, and taxonomy
#' @export

get_taxonomy <- function(accessions, api_key = Sys.getenv("NCBI_API_KEY")) {
    accessions <- unique(as.character(accessions))
    taxids <- accession2taxid(accessions, api_key)

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

accession2taxid <- function(accessions, api_key = Sys.getenv("NCBI_API_KEY")) {
    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

    names(accessions) <- rep("id", times = length(accessions))
    # query <- c(
    #     list(db = "taxonomy", dbfrom = "nuccore", idtype = "acc"),
    #     accessions[41:50]
    # )

    request_base <-
        httr2::request(url) |>
        httr2::req_method("POST")

    acc_chunks <- split(accessions, ceiling(seq_along(accessions) / 100))
    # this keeps the names from populating taxids
    names(acc_chunks) <- NULL

    taxids <-
        lapply(
            acc_chunks,
            function(acc) {
                request_base |>
                    httr2::req_body_form(
                        db = "taxonomy",
                        dbfrom = "nuccore",
                        idtype = "acc",
                        id = acc,
                        .multi = "explode"
                    ) |>
                    httr2::req_retry(max_tries = 5) |>
                    httr2::req_perform() |>
                    httr2::resp_body_string() |>
                    XML::xmlParse() |>
                    XML::xpathSApply(
                        "//LinkSet",
                        parse_LinkSet
                    )
            }
        ) |>
        unlist()

    taxids
}
parse_LinkSet <- function(LinkSet) {
    gid <- XML::xpathSApply(LinkSet, ".//IdList/Id", xmlValue)
    taxid <- XML::xpathSApply(LinkSet, ".//LinkSetDb/Link/Id", xmlValue)
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
