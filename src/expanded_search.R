options(stringsAsFactors=F)
library(tidyverse)
library(janitor)
library(stringr)
if(interactive()) {
  setwd('~/d/sci/src/nhp_models')
}


already_checked = read_tsv('dataset/studies.tsv')
hits = read_csv('data/advanced-search-2022-06-02-c-english.csv') %>% clean_names()
reviews = read_csv('data/advanced-search-2022-06-02-c-reviews.csv') %>% clean_names()

#sum(hits$pmid %in% already_checked$pmid)
#View(hits[!(hits$pmid %in% already_checked$pmid) & !(hits$pmid %in% reviews$pmid),])
#View(already_checked[!(already_checked$pmid %in% hits$pmid),])


hits %>%
  filter(!pmid %in% c(reviews$pmid, already_checked$pmid)) %>%
  mutate(website=paste0('https://pubmed.ncbi.nlm.nih.gov/',pmid),
         pmid,
         study=title,
         first_author,
         year=publication_year,
         date_searched='2022-06-02',
         search_terms_used='see Methods',
         review_status='awaiting',
         comments='',
         detailed_comments='',
         nih_non_nih='',
         year5bin=5*floor(year/5)) %>%
  select(website, pmid, study, first_author, year, date_searched, 
         search_terms_used, review_status, comments, detailed_comments, 
         nih_non_nih, year5bin) -> for_spreadsheet

write_tsv(for_spreadsheet, 'data/search_c_for_spreadsheet.tsv')

to_curate = read_tsv('data/search_c_for_spreadsheet.tsv')

html_raw = read_lines('data/comoy-2022.html')
pubmed_uids = str_extract(html_raw, 'list_uids=[0-9]+')
pmids = gsub('list_uids=','',pubmed_uids[!is.na(pubmed_uids)])
sum(pmids %in% already_checked$pmid | pmids %in% to_curate$pmid)
sum(!(pmids %in% already_checked$pmid | pmids %in% to_curate$pmid))

possible_new_pmids = pmids[!(pmids %in% already_checked$pmid) & !(pmids %in% to_curate$pmid) & !(pmids %in% reviews$pmid)]
search_string = paste0('https://pubmed.ncbi.nlm.nih.gov/?term=', paste(possible_new_pmids, collapse='+or+'))

comoy_2022_new_refs = read_csv('data/comoy-2022-new-refs.csv') %>% clean_names()
comoy_2022_new_refs %>%
  mutate(website=paste0('https://pubmed.ncbi.nlm.nih.gov/',pmid),
         pmid,
         study=title,
         first_author,
         year=publication_year,
         date_searched='2022-06-10',
         search_terms_used='Comoy 2022',
         review_status='awaiting',
         comments='',
         detailed_comments='',
         nih_non_nih='',
         year5bin=5*floor(year/5)) %>%
  select(website, pmid, study, first_author, year, date_searched, 
         search_terms_used, review_status, comments, detailed_comments, 
         nih_non_nih, year5bin) %>%
  arrange(desc(year)) -> comoy_for_spreadsheet

write_tsv(comoy_for_spreadsheet, 'data/comoy_for_spreadsheet.tsv')
