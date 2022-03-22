options(stringsAsFactors = F)
if(interactive()) {
  setwd('~/d/sci/src/nhp_models')
}
library(openxlsx)


#### FUNCTIONS

fill_blank = function(df, col_to_fill, cols_to_require_same = character(0)) {
  last_filled = 1
  for (i in 2:nrow(df)) {
    if (is.na(df[i,col_to_fill]) | df[i,col_to_fill] == '') {
      
      if (length(cols_to_require_same) > 0) {
        if (df[i,cols_to_require_same] == df[last_filled,cols_to_require_same]) {
          df[i,col_to_fill] = df[last_filled,col_to_fill]
        } else {
          next
        }
      } else {
        df[i,col_to_fill] = df[last_filled,col_to_fill]
      }
    } else {
      last_filled = i
    }
    
  }
  return (df)
}


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}


#### CONSTANTS


roa_params = data.frame(roa=c('ic','iv','oral','multiple','other'),
                        x=1:5,
                        color=c('#FFA500','#008B8B','#91219E','#734A12','#ACACAC'))

strain_params = data.frame(strain=c('other','CJD','Kuru','BSE/vCJD','CWD'),
                           x=1:5,
                           color=c('#ACACAC','#8DD3C7','#DAA520','#BEBADA','#FB8072'))

outcome_params = data.frame(outcome=c('endpoint','intercurrent','censored'),
                            x=1:3, 
                            color = c('#386cb0','#beaed4','#C0C0C0'))

endpoint_params = data.frame(endpoint=c('histology','symptom onset','terminal prion disease'),
                             disp=c('histology','symptom onset','terminal disease'),
                             x=1:3,
                             color=c('#FED98E','#FE9929','#CC4C02'))

species_meta = data.frame(grp=c('ape','old world monkey','new world monkey','prosimian'),
                          grp_x=1:4)


path = 'data/full_spreadsheet_2022-03-17-1348.xlsx'

studies = read.xlsx(path,sheet=2)
colnames(studies) = gsub('[^0-9a-z_]','_',tolower(colnames(studies)))
studies = studies[!is.na(studies$website),]

studies$year5bin = 5*floor(studies$year/5)

species = read.xlsx(path,sheet='Notes',startRow = 15,rows=15:43)
colnames(species)[1:6] = c('combined_name','class1','brain_size_notes','body_size','website','comments')
species$grp = tolower(species$class1)
species$grp[species$class1=='Strepsirrhine primate'] = 'prosimian'
species$combined_name[grepl('Capuchin',species$combined_name)] = 'Capuchin monkey (Cebus sp.)'

species$common_name = gsub(' \\(.*','',species$combined_name)
species$scientific_name = gsub('\\)','',gsub('.*\\(','',species$combined_name))

species$grp_x = species_meta$grp_x[match(species$grp, species_meta$grp)]
species = species[with(species, order(grp_x, common_name)),]


#gg = read.table('data/gibbs_gajdusek.tsv',sep='\t',header=T,quote='',comment.char='',nrows = 188)
gg = read.xlsx(path,sheet='NIH')
colnames(gg) = gsub('[^0-9a-z_]','_',tolower(colnames(gg)))

#ao = read.table('data/final_all_other.tsv',sep='\t',header=T,quote='',comment.char='')
ao = read.xlsx(path,sheet='Non-NIH')
colnames(ao) = gsub('[^0-9a-z_]','_',tolower(colnames(ao)))

b9 = read.xlsx(path,sheet='Brown1994')
colnames(b9) = gsub('[^0-9a-z_]','_',tolower(colnames(b9)))

# colnames(gg)
# colnames(ao)
# colnames(b9)

gg$primate_id = ''
b9$primate_id = ''

final_colnames = intersect(intersect(colnames(gg), colnames(ao)), colnames(b9))

data = rbind(gg[,final_colnames], ao[,final_colnames], b9[,final_colnames])

while (any(grepl('__',colnames(data)))) {
  colnames(data) = gsub('__','_',colnames(data))
}

colnames(data)[colnames(data)=='website'] = 'url'
colnames(data)[colnames(data)=='inoculation_route'] = 'roa'
colnames(data)[colnames(data)=='inoculation_route_details'] = 'roa_details'
colnames(data)[colnames(data)=='strain_used'] = 'strain'
colnames(data)[colnames(data)=='species_of_primate'] = 'species'
colnames(data)[colnames(data)=='n_prion_infected_animals'] = 'n'
data$n = suppressWarnings(as.integer(data$n))
colnames(data)[colnames(data)=='n_prion_animals_that_reached_endpoint'] = 'n_endpoint'
data$n_endpoint = suppressWarnings(as.integer(data$n_endpoint))
colnames(data)[colnames(data)=='n_prion_animals_lost_to_intercurrent_illness'] = 'n_intercurrent'
data$n_intercurrent = suppressWarnings(as.integer(data$n_intercurrent))

data$pmid = as.character(data$pmid)

data$cohort_id = paste0(data$pmid, '-', data$cohort)

data$year5bin = studies$year5bin[match(data$pmid, studies$pmid)]

data = fill_blank(data, 'url', 'pmid')
data = fill_blank(data, 'study', 'pmid')
data = fill_blank(data, 'first_author', 'pmid')
data = fill_blank(data, 'year', 'pmid')

data = data[,-which(colnames(data) %in% c('date_searched','search_terms_used'))]

data = fill_blank(data, 'strain', 'cohort_id')
data = fill_blank(data, 'species', 'pmid')
data = fill_blank(data, 'roa', 'cohort_id')
data = fill_blank(data, 'roa_details', 'cohort_id')
data = fill_blank(data, 'n', 'cohort_id')
data$n_endpoint[is.na(data$n_endpoint)] = 0
data$n_intercurrent[is.na(data$n_intercurrent)] = 0
data$n_censored = data$n - data$n_endpoint - data$n_intercurrent

data$mean_mpi = suppressWarnings(as.numeric(data$mean_mpi))
data$sd_mpi = suppressWarnings(as.numeric(data$sd_mpi))

data$endpointx = endpoint_params$x[match(data$endpoint, endpoint_params$endpoint)]

# species within two genuses need to be grouped together b/c authors used them interchangeably or did not specify - Cebus sp. and Saguinus sp.
data$species[grepl('Capuchin',data$species)] = 'Capuchin monkey (Cebus sp.)'
data$species[data$species=='Saddle-backed tamarin (Saguinus fuscicollis), White-fronted tamarin (Saguinus nigricollis), or White-lipped tamarin (Saguinus labiatus)'] = 'Tamarins (Saguinus sp.)'

data$roa_group = 'other'
data$roa_group[!(data$roa %in% c('other','unspecified'))] = data$roa[!(data$roa %in% c('other','unspecified'))]


data$strain_group = 'other'
data$strain_group[grepl('BSE|vCJD',data$strain)] = 'BSE/vCJD'
data$strain_group[grepl('CWD',data$strain)] = 'CWD'
data$strain_group[grepl('Kuru',data$strain)] = 'Kuru'
data$strain_group[grepl('CJD',data$strain) & !grepl('vCJD',data$strain)] = 'CJD'

# QC checks:
# which(!(data$species_of_primate %in% species$name))
# table(data$strain_group)
# table(data$roa_group)
# table(data[,c('year5bin','roa_group')])
# table(data[,c('year5bin','strain_group')])

data = data[with(data, order(cohort_id, -endpointx)),]
data$row = 1:nrow(data)

data$dup_status = 'unique'
data$dup_status[duplicated(data$cohort_id)] = 'duplicate' # different endpoint, same animals, same paper
data$dup_status[data$exclude != 'include'] = 'duplicate' # same or overlapping animals, different paper
data$dup_status[data$pmid=='8179297'] = 'indeterminate'


write.table(studies,'dataset/studies.tsv',sep='\t',row.names=F,col.names=T,quote=F,na='')
write.table(species,'dataset/species.tsv',sep='\t',row.names=F,col.names=T,quote=F,na='')
write.table(data,'dataset/data.tsv',sep='\t',row.names=F,col.names=T,quote=F,na='')


write.xlsx(list(table_s1_studies=studies,
                table_s2_data=data,
                table_s3_species=species),
           file='display_items/supplementary-file.xlsx',
           overwrite = T,
           firstRow=T,
           headerStyle=createStyle(textDecoration='bold'))



