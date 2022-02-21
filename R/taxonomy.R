#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Retrieve the taxonomy information from NCBI for a set of nucleotide gis.
#'
#' @param accessions a character vector of the accessions to retrieve
#' @return data.frame of the 'accessions, taxIds, and taxonomy
#' @export

get_taxonomy = function(accessions){

  accessions = unique(as.character(accessions))
  taxids = accession2taxid(accessions)

  taxonomy=fetch_taxonomy(unique(taxids))
  merge(
    data.frame(accession=names(taxids), taxId=taxids, stringsAsFactors=FALSE),
    taxonomy
  )
}
#' Maps a nucleotide database accession to a taxonomy database taxId
#'
#' @param accessions accessions character vector to lookup.
#' @return named vector of taxIds.
#' @export

accession2taxid = function(accessions){
  url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'

  names(accessions) = rep('id', times=length(accessions))
  query=c(list(db='taxonomy', dbfrom='nuccore', idtype='acc'), accessions)

  response=POST_retry(url, body=query)

  #stop if response failed
  stop_for_status(response)

  parsed = xmlParse(content(response, as="text"))

  xpathSApply(parsed, '//LinkSet', parse_LinkSet)
}
parse_LinkSet = function(LinkSet){
  gid = xpathSApply(LinkSet, './/IdList/Id', xmlValue)
  taxid = xpathSApply(LinkSet, './/LinkSetDb/Link/Id', xmlValue)
  if(length(taxid) != 1)
    taxid = NA
  names(taxid) = gid
  taxid
}

fetch_taxonomy = function(taxid) {

  fetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

  query=list(db='taxonomy', rettype='null', retmode='xml', id=paste(taxid, collapse=','))

  response = POST_retry(fetch_url, body=query)

  #stop if response failed
  stop_for_status(response)

  parse_taxonomy_xml(xmlParse(content(response, as="text")))
}
parse_taxonomy_xml = function(xml){
  rbind.fill(xpathApply(xml, '//TaxaSet/Taxon', parse_taxon))
}
parse_taxon = function(taxon){
  tax_id = xpathSApply(taxon, './TaxId', xmlValue)
  ranks = xpathSApply(taxon, './/Rank', xmlValue)
  names = xpathSApply(taxon, './/ScientificName', xmlValue)
  names(names) = ranks
  names = names[ranks != 'no rank']
  names['taxId'] = tax_id
  as.data.frame(t(names), stringsAsFactors=FALSE)
}
