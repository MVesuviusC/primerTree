#PrimerTree
#Copyright (C) 2013 Jim Hester

#' Retrieve the taxonomy information from NCBI for a set of nucleotide gis.
#'
#' @param gis a character vector of the gis to retrieve
#' @return data.frame of the 'gis, taxIds, and taxonomy
#' @export

get_taxonomy = function(gis){

  gis = unique(as.character(gis))
  taxids = gi2taxid(gis)

  taxonomy=fetch_taxonomy(unique(taxids))
  merge(
    data.frame(gi=names(taxids), taxId=taxids, stringsAsFactors=FALSE),
    taxonomy
  )
}
#' Maps a nucleotide database gi to a taxonomy database taxId
#'
#' @param gi gi character vector to lookup.
#' @return named vector of taxIds.
#' @export

gi2taxid = function(gi){
  url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'

  names(gi) = rep('id', times=length(gi))
  query=c(list(db='taxonomy', dbfrom='nuccore'), gi)

  response=POST_retry(url, body=query)

  #stop if response failed
  stop_for_status(response)

  parsed = content(response, type='text/xml')

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
#' Maps a genbank accession to a nuclotide database gi.
#'
#' @param accession accession character vector to lookup.
#' @return named vector of gis.
#' @export

accession2gi = function(accession){
  url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
  query=list(db='nuccore', rettype='seqid', id=paste(collapse=',', accession))

  response=POST_retry(url, body=query)

  #stop if response failed
  stop_for_status(response)

  mapping = gsub('.*?accession "([^"]+)"[^v]+version (\\d+).*?gi (\\d+)\n+',
                 '"\\1.\\2"="\\3",', content(response, 'text'))
  #remove trailing , from the substitution
  mapping = substr(mapping, 0, nchar(mapping)-1)

  #convert to named vector
  eval(parse(text=paste('c(', mapping, ')', sep='')))
}
fetch_taxonomy = function(taxid) {

  fetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

  query=list(db='taxonomy', rettype='null', retmode='xml', id=paste(taxid, collapse=','))

  response = POST_retry(fetch_url, body=query)

  #stop if response failed
  stop_for_status(response)

  parse_taxonomy_xml(content(response))
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
