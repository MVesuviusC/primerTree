#' Query a pair of primers using ncbi's Primer-BLAST, if primers contain iupac
#'
#' ambiguity codes, enumerate all possible combinations and combine the
#' results.
#' @param forward forward primer to search.
#' @param reverse reverse primer to search.
#' @param name name to use for the primer pair, if not given is simply the
#'        forward and reverse primers seperated by an underscore.
#' @return data.frame of primer hits
#' @export primer_search
primer_search = function(forward, reverse,
                          name=paste(forward, reverse, sep='_'), ...){
  primers_search(forward=forward, reverse=reverse, name=name, ...)[[1]]
}

#' Query multiple pairs of primers using ncbi's Primer-BLAST, if primers contain iupac
#'
#' ambiguity codes, enumerate all possible combinations and combine the
#' results.
#' @inheritParams primer_search
#' @return list of data.frames with primer hits
#' @export primers_search

primers_search = function(forward, reverse,
                          name=paste(forward, reverse, sep='_'), ...){
  #enumerate all combinations to handle ambiguity codes
  forward_primers = enumerate_ambiguity(forward)
  primers = data.frame(forward=forward, 
                       reverse=rep(enumerate_ambiguity(reverse), 
                                   each=length(forward_primers)))
  alply(primers, .margins=1, expand=F, blast_primer, ...)
}
iupac = list( "M" = list("A", "C"),
              "R" = list("A", "G"),
              "W" = list("A", "T"),
              "S" = list("C", "G"),
              "Y" = list("C", "T"),
              "K" = list("G", "T"),
              "V" = list("A", "C", "G"),
              "H" = list("A", "C", "T"),
              "D" = list("A", "G", "T"),
              "B" = list("C", "G", "T"),
              "N" = list("A", "C", "G", "T"),
              "I" = list("A", "T", "C"))

enumerate_ambiguity = function(sequence){
  library(stringr)
  search_regex = paste(names(iupac), collapse='|')
  locs = str_locate_all(sequence, search_regex)
  sequences = list()
  count = 1
  for (i in seq_len(nrow(locs[[1]]))){
    loc = locs[[1]][i,]
    ambiguity = str_sub(sequence, loc[1], loc[2])
    for(type in iupac[[ambiguity]]){
      new_seq = sequence
      str_sub(new_seq, loc[1], loc[2]) <- type
      sequences[[count]] = enumerate_ambiguity(new_seq)
      count = count + 1
    }
    return(unlist(sequences))
  }
  return(sequence)
}


blast_primer = function(forward, reverse, ..., organism='',
  primer_specificity_database='nt', exclude_env='on'){
  library(httr)
  library(plyr)

  url = 'http://www.ncbi.nlm.nih.gov/tools/primer-blast/'
  form = GET(url)

  stop_for_status(form)

  content = parsable_html(form)

  all_options = get_options(content)

  options = list(..., primer_left_input=forward, primer_right_input=reverse,
                 organism=organism,
                 primer_specificity_database=primer_specificity_database,
                 search_specific_primer='on')

  names(options) = toupper(names(options))

  match_args = pmatch(names(options), all_options$name)
  bad_args = is.na(match_args)

  if(any(bad_args)){
    message(subset(all_options, is.na(type) | type != 'hidden',  select=c(name, type, defval)))
    stop(paste(names(options)[bad_args], collapse=','), ' not valid option\n')
  }

  options = get_defaults(options, all_options)

  response = POST(paste(url, 'primertool.cgi', sep=''), body=options)

  stop_for_status(response)

  values = get_refresh_from_meta(response)

  while(length(values) > 0){
    message('blast alignment processing, refreshing in ', values[1], ' seconds...')
    Sys.sleep(values[1])
    response = GET(values[2])
    stop_for_status(response)
    values = get_refresh_from_meta(response)
  }
  message('blast alignment complete')
  return(response)
}

parse_primer_hits = function(response){
  content = parsable_html(response)

  rbind.fill(xpathApply(content, '//pre', parse_pre))
}
parse_a = function(a){
  gid = gsub('.*id=(\\d+)', '\\1', xmlAttrs(a)['href'])
  data.frame(gid=as.character(gid), accession=as.character(xmlValue(a)))
}

parse_pre = function(pre){
  library(stringr)
  pre_text = xmlValue(pre)

  a = getNodeSet(pre, './preceding-sibling::a[1]')
  if(length(a) <= 0)
    stop('Parsing failed for ', pre_text)

  ids = parse_a(a[[1]])

  product_length_regex = 'product length = (\\d+)'
  template_regex = 'Template[^\\d]+(\\d+)[^.ACGT]+([.ACGT]+)[^\\d]+(\\d+)'
  full_regex = paste('[\\S\\W]*', product_length_regex, '[\\S\\W]*?',
                     template_regex, '[\\S\\W]*', template_regex, '[\\S\\W]*', sep='')
  values = str_split(gsub(full_regex, paste('\\', 1:8, sep='', collapse='|'),
                           pre_text, perl=T), '[|]')[[1]]
  data.frame(ids,
             product_length=as.numeric(values[1]),
             mismatch_forward=str_count(values[3], '[ACGT]'),
             mismatch_reverse=str_count(values[6], '[ACGT]'),
             forward_start = as.numeric(values[2]),
             forward_stop = as.numeric(values[4]),
             reverse_start = as.numeric(values[5]),
             reverse_stop = as.numeric(values[7]),
             start=min(values[c(2,4,5,7)]),
             stop=max(values[c(2,4,5,7)])
             )
}
get_refresh_from_meta = function(response){
  library(stringr)
  content = parsable_html(response)
  meta = content['//meta[@http-equiv="Refresh"]']
  if(length(meta) > 0){
    values = str_split(xmlAttrs(meta[[1]])['content'], '; URL=')[[1]]
    return(values)
  }
  return()
}

get_defaults = function(set_options, options){
  #only look at values with set defaults
  options = options[ !is.na(options$defval), ]
  unchanged_options = setdiff(options$name, names(set_options))
  default_values = as.character(options[ match(unchanged_options, options$name), 'defval' ])
  names(default_values) = unchanged_options
  c(set_options, default_values)
}

get_options = function(content){
  options = rbind.fill(xpathApply(content, '//form//input | //form//select', parse_attributes))
  options$type = as.character(options$type)

  #add dropdown type if they are NA
  options$type[is.na(options$type)] <- 'dropdown'

  #make default values for checkboxes on or off rather than checked or unchecked
  options$defval = as.character(options$defval)
  check_map = c('checked' = 'on', 'unchecked' = '')
  checkboxes = which(options$type == 'checkbox')
  options$defval[ checkboxes ] = check_map[ options$defval[ checkboxes ] ]

  subset(options, type != 'hidden',  select=c(name, type, defval))
}

parse_attributes = function(x){
  as.data.frame(t(xmlAttrs(x)))
}
parsable_html = function(response){
  library(XML)
  content = htmlParse(content(response, as='text'))
}
