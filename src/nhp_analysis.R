options(stringsAsFactors = F)
if(interactive()) {
  setwd('~/d/sci/src/nhp_models')
}

start_time = Sys.time()
cat(file=stderr(), 'Loading packages, constants, and datasets...')

suppressMessages(library(stats))
suppressMessages(library(sqldf))
suppressMessages(library(survival))
suppressMessages(library(DiagrammeR))
suppressMessages(library(DiagrammeRsvg))
suppressMessages(library(magrittr))
suppressMessages(library(rsvg))

#### FUNCTIONS

# apply numeric transparency value to a hex color
alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # default degree of transparency for shading confidence intervals in plot


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

year5bins = data.frame(year5bin=seq(1965,2015,5))

#### READ IN ANALYTICAL DATASET

studies = read.table('dataset/studies.tsv', sep='\t', header=T, quote='', comment.char='')
species = read.table('dataset/species.tsv', sep='\t', header=T, quote='', comment.char='')
data = read.table('dataset/data.tsv', sep='\t', header=T, quote='', comment.char='')

data_dedup = data[data$dup_status=='unique',]






### FIGURE 1: FLOWCHART

cat(file=stderr(), 'done.\nCreating Figure 1...')

studies$flowchart = ''
studies$flowchart[studies$review_status != 'excluded'] = 'Included'
studies$flowchart[studies$comments %in% c('all animals included in other publications','review article','no nhp endpoint data')] = '1 No original NHP endpoint data*'
studies$flowchart[studies$comments %in% c('full text not found')] = '4 Full text not found'
studies$flowchart[studies$comments %in% c('not enough details')] = '3 Not enough detail'
studies$flowchart[studies$comments %in% c('no nhp data')] = '2 No NHP data'
studies$flowchart[studies$comments %in% c('drug interventions')] = '5 Drug intervention'

exclusion_table = sqldf("select flowchart, count(*) n from studies where review_status = 'excluded' group by 1;")

exclusion_text = paste0(paste0(gsub('^[0-9] ','',exclusion_table$flowchart), ' n=', exclusion_table$n), collapse='\\n')

n_coh_epts = nrow(data)
n_coh_epts_indet = sum(data$dup_status=='indeterminate')
n_coh_epts_dup = sum(data$dup_status=='duplicate')
n_coh_epts_uniq = sum(data$dup_status=='unique')

n_cohorts_uniq1 = length(unique(data$cohort_id[data$dup_status=='unique']))
n_cohorts_uniq2 = nrow(data_dedup)
stopifnot(n_cohorts_uniq2 == n_cohorts_uniq1)

n_animals_dedup = sum(data_dedup$n)

resx=300
grViz(paste0("digraph flowchart {
             # node definitions with substituted label text
             node [fontname = Helvetica, shape = rectangle]        
             tab1 [label = '@@1']
             tab2 [label = '@@2']
             tab3 [label = '@@3']
             tab4 [label = '@@4']
             tab5 [label = '@@5']
             tab6 [label = '@@6']
             tab7 [label = '@@7']
             tab8 [label = '@@8']
             
             # edge definitions with the node IDs
             tab1 -> tab2;
             tab1 -> tab3;
             tab3 -> tab4;
             tab4 -> tab6;
             tab4 -> tab5;
             tab4 -> tab7;
             tab5 -> tab8
             }
             
             [1]: 'Articles examined: n=",nrow(studies),"'
             [2]: 'Articles excluded: n=",sum(studies$review_status == 'excluded'),"\\n\\n",exclusion_text,"'
             [3]: 'Articles included: n=",sum(studies$review_status != 'excluded'),"'
             [4]: 'Total rows: n=",n_coh_epts,"'
             [5]: 'Unique cohorts: n=",n_coh_epts_uniq,"'
             [6]: 'Duplicated cohorts: n=",n_coh_epts_dup,"'
             [7]: 'Cohorts of indeterminate\\nduplicate status: n=",n_coh_epts_indet,"'
             [8]: 'Unique animals: n=",n_animals_dedup,"'
             ")) %>% export_svg %>% charToRaw %>% rsvg_png("display_items/figure-1.png",width=3.5*resx)

### END FIGURE 1













### FIGURE 2: CHARACTERIZATION OF THE DATASET

cat(file=stderr(), 'done.\nCreating Figure 2...')

resx=300
png('display_items/figure-2.png',width=5.3*resx,height=5*resx,res=resx)
#layout_matrix = matrix(c(1,1,2,2,3,3,4,5,6,6,6,6),byrow=T,nrow=6)
layout_matrix = matrix(c(1,1,1,4,4,
                         2,2,2,5,5,
                         3,3,3,6,6,
                         7,7,7,7,7,
                         7,7,7,7,7),byrow=T,nrow=5)
layout(layout_matrix)
panel = 1


year_study = sqldf("
                   select   y.year5bin, sum(case when review_status != 'excluded' then 1 else 0 end) n_studies
                   from     year5bins y left outer join studies s
                   on       y.year5bin = s.year5bin
                   group by 1
                   order by 1
                   ;")
data_dedup$year5bin = studies$year5bin[match(data_dedup$pmid, studies$pmid)]
year_animal = sqldf("
                    select   y.year5bin, sum(case when exclude = 'include' then d.n else 0 end) n_animals
                    from     year5bins y left outer join data_dedup d
                    on       y.year5bin = d.year5bin
                    group by 1
                    order by 1
                    ;")

par(mar=c(2,3,3,6))
study_col = '#0071BC'
animal_col = '#FF8E43'
barwidth=4
xlims = c(1965, 2020)
ylims = c(0,11)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(1960, 2020, 10), lwd=0, line=-0.5)
axis(side=1, at=seq(1960, 2020, 5), labels=NA, lwd=1, lwd.ticks=0.5, tck=-0.1)
axis(side=1, at=seq(1960, 2020, 1), labels=NA, lwd=0, lwd.ticks=0.5, tck=-0.05)
axis(side=2, at=0:2*5, tck=-0.1, labels=NA)
axis(side=2, at=0:10, tck=-0.05, labels=NA)
axis(side=2, at=0:2*5, las=2, line=-0.5, lwd=0)
mtext(side=2, line=1.75, text='N studies', col='black')
year_study$mid = year_study$year5bin + 2.5
year_animal$mid = year_animal$year5bin + 2.5
rect(xleft=year_study$mid-barwidth/2, xright=year_study$mid+barwidth/2, ybottom=rep(0, nrow(year_study)), ytop=year_study$n_studies, col=study_col, border=NA)
mtext(side=3, line=1, text='studies included', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.1, text=LETTERS[panel])
panel = panel + 1


year_roa = sqldf("
select   d.year5bin, r.x, r.color, r.roa, sum(d.n) n_animals
from     data_dedup d, roa_params r
where    d.roa_group = r.roa
and      d.exclude = 'include'
group by 1, 2, 3, 4
;")

strain_data = sqldf("
select   d.year5bin, d.strain_group, sum(d.n) n_animals
from     data_dedup d
where    d.exclude = 'include'
group by 1, 2
;")
year_strain_cp = sqldf("select y.year5bin, p.* from year5bins y, strain_params p;")

year_strain_complete = sqldf("
select   c.year5bin, c.x, c.color, c.strain, n_animals
from     year_strain_cp c left outer join strain_data d
on       c.year5bin = d.year5bin and c.strain = d.strain_group
order by 1, 2
;")
year_strain_complete$n_animals[is.na(year_strain_complete$n_animals)] = 0
year_strain_complete$cumn = 0
for (i in 1:nrow(year5bins)) {
  rows = year_strain_complete$year5bin == year5bins$year5bin[i]
  year_strain_complete$cumn[rows] = cumsum(year_strain_complete$n_animals[rows])
}
year_strain_complete$mid = year_strain_complete$year5bin+2.5

barwidth=4
xlims = c(1965, 2020)
ylims = c(0,120)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(1960, 2020, 10), lwd=0, line=-0.5)
axis(side=1, at=seq(1960, 2020, 5), labels=NA, lwd=1, lwd.ticks=0.5, tck=-0.1)
axis(side=1, at=seq(1960, 2020, 1), labels=NA, lwd=0, lwd.ticks=0.5, tck=-0.05)
axis(side=2, at=0:3*50, tck=-0.1, labels=NA)
axis(side=2, at=0:3*50, las=2, line=-0.5, lwd=0)
axis(side=2, at=0:15*10, tck=-0.05, labels=NA)
mtext(side=2, line=1.75, text='N animals', col='black')
rect(xleft=year_strain_complete$mid-barwidth/2, xright=year_strain_complete$mid+barwidth/2, ybottom=year_strain_complete$cumn-year_strain_complete$n_animals, ytop=year_strain_complete$cumn, col=year_strain_complete$color, border=NA)
par(xpd=T)
legend(x=max(xlims)+2,y=max(ylims)*1.75,cex=1.0,legend=strain_params$strain,pch=15,col=strain_params$color,text.col=strain_params$color,bty='n')
par(xpd=F)
mtext(side=3, line=1, text='prion strain', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.1, text=LETTERS[panel])
panel = panel + 1




roa_data = sqldf("
                 select   d.year5bin, d.roa_group, sum(d.n) n_animals
                 from     data_dedup d
                 where    d.exclude = 'include'
                 group by 1, 2
                 ;")
year_roa_cp = sqldf("select y.year5bin, p.* from year5bins y, roa_params p;")

year_roa_complete = sqldf("
                          select   c.year5bin, c.x, c.color, c.roa, n_animals
                          from     year_roa_cp c left outer join roa_data d
                          on       c.year5bin = d.year5bin and c.roa = d.roa_group
                          order by 1, 2
                          ;")
year_roa_complete$n_animals[is.na(year_roa_complete$n_animals)] = 0
year_roa_complete$cumn = 0
for (i in 1:nrow(year5bins)) {
  rows = year_roa_complete$year5bin == year5bins$year5bin[i]
  year_roa_complete$cumn[rows] = cumsum(year_roa_complete$n_animals[rows])
}
year_roa_complete$mid = year_roa_complete$year5bin+2.5

barwidth=4
xlims = c(1965, 2020)
ylims = c(0,120)
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(1960, 2020, 10), lwd=0, line=-0.5)
axis(side=1, at=seq(1960, 2020, 5), labels=NA, lwd=1, lwd.ticks=0.5, tck=-0.1)
axis(side=1, at=seq(1960, 2020, 1), labels=NA, lwd=0, lwd.ticks=0.5, tck=-0.05)
axis(side=2, at=0:3*50, tck=-0.1, labels=NA)
axis(side=2, at=0:3*50, lwd=0, las=2, line=-0.5)
axis(side=2, at=0:15*10, tck=-0.05, labels=NA)
mtext(side=2, line=1.75, text='N animals', col='black')
rect(xleft=year_roa_complete$mid-barwidth/2, xright=year_roa_complete$mid+barwidth/2, ybottom=year_roa_complete$cumn-year_roa_complete$n_animals, ytop=year_roa_complete$cumn, col=year_roa_complete$color, border=NA)
par(xpd=T)
legend(x=max(xlims)+2,y=max(ylims)*1.75,cex=1.0,legend=roa_params$roa,pch=15,col=roa_params$color,text.col=roa_params$color,bty='n')
par(xpd=F)
mtext(side=3, line=1, text='route of administration', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.1, text=LETTERS[panel])
panel = panel + 1

par(mar=c(2,4,3,3))
endpoint_n = sqldf("
select   p.x, p.color, p.disp, sum(dd.n) n_animals
from     data_dedup dd, endpoint_params p
where    dd.endpoint = p.endpoint
group by 1, 2, 3
order by 1 ,2, 3
;")
barwidth=0.4
xlims = c(0,400)
ylims = range(endpoint_n$x) + c(-0.5, 0.5)
endpoint_n$y = max(endpoint_n$x)-endpoint_n$x+1
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4*100, tck=-0.1, labels=NA)
axis(side=1, at=0:40*10, tck=-0.05, labels=NA)
axis(side=1, at=0:4*100, lwd=0, line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, at=endpoint_n$y, text=endpoint_n$disp, line=0.25, las=2, cex=0.6)
rect(xleft=rep(0,nrow(endpoint_n)), xright=endpoint_n$n_animals, ybottom=endpoint_n$y-barwidth, ytop=endpoint_n$y+barwidth, col=endpoint_n$color,border=NA)
mtext(side=3, line=1, text='endpoint assessed', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.3, text=LETTERS[panel])
panel = panel + 1

passage_params = data.frame(passage=c('primary','2nd','3rd+'),
                            color=c('#FFAA00','#9ACD32','#4393C3'),
                            x=1:3)
data_dedup$passaged = grepl('\u2192',data_dedup$strain) # one right arrow means 2nd passage
data_dedup$passage3 = grepl('\u2192.*\u2192',data_dedup$animal_comments) # 2+ right arrows means 3rd or higher
data_dedup$passage3 = grepl('3rd',data_dedup$animal_comments) # "3rd" also means 3rd or higher
data_dedup$passage = 'primary'
data_dedup$passage[data_dedup$passaged] = '2nd'
data_dedup$passage[data_dedup$passage3] = '3rd+'

passage_n = sqldf("
select   p.x, p.color, p.passage, sum(dd.n) n_animals
from     data_dedup dd, passage_params p
where    dd.passage = p.passage
group by 1, 2, 3
order by 1, 2, 3
;")
par(mar=c(2,4,3,3))
xlims = c(0,100*ceiling(max(passage_n$n_animals)/100))
ylims = range(passage_n$x) + c(-0.5, 0.5)
passage_n$y = max(passage_n$x)-passage_n$x+1
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4*100, tck=-0.15, labels=NA)
axis(side=1, at=0:40*10, tck=-0.075, labels=NA)
axis(side=1, at=0:4*100, lwd=0, line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, at=passage_n$y, text=passage_n$passage, line=0.25, las=2, cex=0.6)
rect(xleft=rep(0,nrow(passage_n)), xright=passage_n$n_animals, ybottom=passage_n$y-barwidth, ytop=passage_n$y+barwidth, col=passage_n$color, border=NA)
mtext(side=3, line=1, text='passage', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.3, text=LETTERS[panel])
panel = panel + 1


outcome_params$n = 0
outcome_params$n[outcome_params$outcome=='endpoint'] = sum(data_dedup$n_endpoint)
outcome_params$n[outcome_params$outcome=='intercurrent'] = sum(data_dedup$n_intercurrent)
outcome_params$n[outcome_params$outcome=='censored'] = sum(data_dedup$n_censored)
par(mar=c(2,4,3,3))
xlims = c(0,400)
ylims = range(outcome_params$x) + c(-0.5, 0.5)
outcome_params$y = max(outcome_params$x)-outcome_params$x+1
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4*100, tck=-0.15, labels=NA)
axis(side=1, at=0:40*10, tck=-0.075, labels=NA)
axis(side=1, at=0:4*100, lwd=0, line=-0.5)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, at=outcome_params$y, text=outcome_params$outcome, line=0.25, las=2, cex=0.6)
rect(xleft=rep(0,nrow(outcome_params)), xright=outcome_params$n, ybottom=outcome_params$y-barwidth, ytop=outcome_params$y+barwidth, col=outcome_params$color, border=NA)
mtext(side=3, line=1, text='outcome', cex=0.8)
mtext(side=3, line=0.25, cex=1.5, adj=-0.3, text=LETTERS[panel])
panel = panel + 1

# double check that no species will drop out on the join due to spelling mismatches:
# table(data_dedup$species[which(!(data_dedup$species %in% species$combined_name))])

species_count = sqldf("
select   s.grp_x, s.grp, s.common_name, s.scientific_name, dd.species, sum(dd.n) n_animals, sum(dd.n_endpoint) n_endpoint
from     data_dedup dd, species s
where    dd.species = s.combined_name
group by 1, 2, 3, 4, 5
order by 1, 2, 3
;")

# color gradients for reference:
# ['#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506']
# ['#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026']

species_count$color[species_count$grp_x == 1] = c('#238b45','#005a32')
species_count$color[species_count$grp_x == 2] = c('#fe9929','#ec7014','#cc4c02','#993404','#662506','#e31a1c','#bd0026','#800026') #rev(c('#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'))
species_count$color[species_count$grp_x == 3] = c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
species_count$color[species_count$grp_x == 4] = c('#6a51a3')

species_count$y = nrow(species_count):1

species_grps = sqldf("
select   grp, max(y) maxy, min(y) miny, avg(y) midy
from     species_count
group by 1
order by 4 desc
;")

par(mar=c(3,21,2,3))
ylims = range(species_count$y) + c(-0.6, 0.6)
xlims = c(0,50*ceiling(max(species_count$n_animals)/50))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, at=species_count$y, line=8, las=2, text=species_count$common_name, cex=0.5, font=1)
mtext(side=2, at=species_count$y, line=0.25, las=2, text=species_count$scientific_name, cex=0.5, font=3)
#mtext(side=2, at=species_count$y, line=0.25, las=2, text=species_count$common_name, cex=0.5, font=1)
axis(side=1, at=0:5*50, tck=-0.05, labels=NA)
axis(side=1, at=0:5*50, lwd=0, line=-0.5)
axis(side=1, at=0:25*10, tck=-0.025, labels=NA)
rect(xleft=rep(0,nrow(species_count)), xright=species_count$n_animals, ybottom=species_count$y-barwidth, ytop=species_count$y+barwidth, col=alpha(species_count$color,ci_alpha), border=NA)
rect(xleft=rep(0,nrow(species_count)), xright=species_count$n_endpoint, ybottom=species_count$y-barwidth, ytop=species_count$y+barwidth, col=species_count$color, border=NA)

bracket_dist = 15
overshoot = 0.3
for (i in 1:nrow(species_grps)) {
  axis(side=2, line=bracket_dist + 0.5, at=species_grps[i,c('miny','maxy')] + c(-overshoot,overshoot), labels=NA, tck=0.025)
  axis(side=2, line=bracket_dist - 0.25, at=species_grps$midy[i], labels=gsub(' ','\n',species_grps$grp[i]), las=2, lwd=0)
}

mtext(side=3, line=0.25, cex=1.5, adj=-0.3, text=LETTERS[panel])

unnecessary_message = dev.off()

### END FIGURE 2







### FIGURE 3: INCUBATON TIME IN COHORTS MEETING POWER CALCULATION CRITERIA

cat(file=stderr(), 'done.\nCreating Figure 3...')

resx=300
png('display_items/figure-3.png',width=6.5*resx,height=4.5*resx,res=resx)

power_cohorts = sqldf("
                      select   *
                      from     data_dedup
                      where    n >= 3
                      and      n_endpoint >= 3
                      and      n = n_endpoint + n_intercurrent
                      and      endpoint in ('terminal prion disease','symptom onset')
                      and      mean_mpi is not null and sd_mpi is not null
                      ;")

power_cohorts$species_y = species_count$y[match(power_cohorts$species, species_count$species)]
power_cohorts$color = species_count$color[match(power_cohorts$species, species_count$species)]

power_cohorts$max_mpi = as.numeric(gsub('.+-','',power_cohorts$range_mpi))
power_cohorts$min_mpi = as.numeric(gsub('-.*','',power_cohorts$range_mpi))

power_cohorts$l95 = power_cohorts$mean_mpi - 1.96 * power_cohorts$sd_mpi/sqrt(power_cohorts$n_endpoint)
power_cohorts$u95 = power_cohorts$mean_mpi + 1.96 * power_cohorts$sd_mpi/sqrt(power_cohorts$n_endpoint)

power_cohorts = power_cohorts[with(power_cohorts, order(-species_y, mean_mpi)),]
power_cohorts$y = nrow(power_cohorts):1

studies$shortname = paste0(gsub(',.*','',studies$first_author),' ',studies$year)
power_cohorts$shortname = studies$shortname[match(power_cohorts$pmid, studies$pmid)]

power_cohorts$common_name = species_count$common_name[match(power_cohorts$species, species_count$species)]

power_cohorts$roa_disp = power_cohorts$roa_details
power_cohorts$roa_disp[power_cohorts$roa_details == ''] = power_cohorts$roa[power_cohorts$roa_details == '']
power_cohorts$details_disp = paste0(power_cohorts$strain, ' ', power_cohorts$roa_disp)
power_cohorts$details_disp = gsub('ic','i.c.',gsub('iv','i.v.',gsub('ip','i.p.',power_cohorts$details_disp)))
power_cohorts$details_disp[power_cohorts$details_disp=='BSE → Cynomolgus macaque i.v.'] = 'BSE→Cyno i.v.'
power_cohorts$details_disp[power_cohorts$details_disp=='vCJD i.c./intratonsilar'] = 'vCJD i.c./tonsil'

power_cohort_species = sqldf("select common_name, color, max(y) maxy from power_cohorts group by 1, 2 order by 1, 2;")

ylims = range(power_cohorts$y) + c(-0.6, 0.6)
xlims = c(0, max(power_cohorts$max_mpi, na.rm=T)*1.05)

par(mar=c(3,12,1,6))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
mtext(side=2, at=power_cohorts$y, text=power_cohorts$shortname, las=2, line=0.25)
axis(side=1, at=seq(0,108,12), labels=seq(0,108,12)/12)
axis(side=1, at=seq(0,108,1), labels=NA, lwd=0, lwd.ticks=0.5, tck=-0.01)
mtext(side=1, line=2, text='years to endpoint (mean, range)')
points(x=power_cohorts$mean_mpi, y=power_cohorts$y, pch=19, col=power_cohorts$color)
arrows(x0=power_cohorts$min_mpi, x1=power_cohorts$max_mpi, y0=power_cohorts$y, col=power_cohorts$color, lwd=3, angle=90, length=0.05, code=3)
abline(h=power_cohort_species$maxy + 0.5, lwd=0.25)
text(x=4*12, y=power_cohort_species$maxy, pos=4, col=power_cohort_species$color, labels=power_cohort_species$common_name, font=2, cex=0.7)
mtext(side=4, line=0.25, at=power_cohorts$y, text=power_cohorts$details_disp, las=2)
unnecessary_message = dev.off()

### END FIGURE 3














### TABLE 1

cat(file=stderr(), 'done.\nCreating Table 1...')

data_dedup$roa_ic = ''
data_dedup$roa_ic[grepl('ic',data_dedup$roa_group)] = 'i.c.'
data_dedup$roa_ic[grepl('ic',data_dedup$roa_details)] = 'i.c.'

power_cohorts$roa_ic = ''
power_cohorts$roa_ic[grepl('ic',power_cohorts$roa_group)] = 'i.c.'
power_cohorts$roa_ic[grepl('ic',power_cohorts$roa_details)] = 'i.c.'
power_cohorts$strain_roa = paste0(power_cohorts$strain, ' ', power_cohorts$roa_ic)

power_cohorts$mean_sd_bc = paste0(power_cohorts$mean_mpi,'±',formatC(power_cohorts$sd_mpi,format='f',digits=1))

tbl1_cols = c('common_name','strain','roa_group','strain_roa','shortname','mean_mpi','sd_mpi','mean_sd_bc','cohort_id')
table1_bestcase = power_cohorts[c(1,3,14,16,19),tbl1_cols]
data_dedup$common_name = species$common_name[match(data_dedup$species,species$combined_name)]

table1 = table1_bestcase
table1$attack_rate = ''
table1$mean_sd = ''
table1$mean = as.numeric(NA)
table1$sd = as.numeric(NA)
table1$expected_duration_oth = as.integer(NA)
table1$power_oth = as.numeric(NA)

row = data$cohort_id=='8179297-1' & data$endpoint=='symptom onset'
# table1$attack_rate[1] = paste0(data$n_endpoint[row],'/',data$n[row])
table1$attack_rate[1] = '24/29' # account for the 4 outliers they excluded from their mean±sd calculation - explained in text
table1$mean[1] = data$mean_mpi[row]
table1$sd[1] = data$sd_mpi[row]
table1$mean_sd[1] = paste0(data$mean_mpi[row],'±',data$sd_mpi[row],'*')

# this study was only 2 animals - we should require at least 3, like in Figure 3..
# row = data$cohort_id == '11259641-9' & data$endpoint=='terminal prion disease'
# table1$attack_rate[2] = paste0(data$n_endpoint[row],'/',data$n[row])
# table1$mean[2] = data$mean_mpi[row]
# table1$sd[2] = data$sd_mpi[row]
# table1$mean_sd[2] = paste0(data$mean_mpi[row],'±',data$sd_mpi[row])

row = data$cohort_id=='8179297-3' & data$endpoint=='symptom onset'
table1$attack_rate[3] = paste0(data$n_endpoint[row],'/',data$n[row])
table1$mean[3] = data$mean_mpi[row]
table1$sd[3] = data$sd_mpi[row]
table1$mean_sd[3] = paste0(data$mean_mpi[row],'±',data$sd_mpi[row])

row = data$cohort_id=='8179297-2' & data$endpoint=='symptom onset'
table1$attack_rate[4] = paste0(data$n_endpoint[row],'/',data$n[row])
table1$mean[4] = data$mean_mpi[row]
table1$sd[4] = data$sd_mpi[row]
table1$mean_sd[4] = paste0(data$mean_mpi[row],'±',data$sd_mpi[row])

table1$p_attack[1] = 24/29 
table1$p_attack[2] = 2/2
table1$p_attack[3] = 30/31
table1$p_attack[4] = 196/211
table1$p_attack[5] = NA

n_iter = 1000
set.seed(1)
for (i in 1:nrow(table1)) {
  n_integer = 6
  maxes = numeric(0)
  pvals = numeric(0)
  for (j in 1:n_iter) {
    placebo_survival = rnorm(n=n_integer, m=table1$mean_mpi[i]*1.0, s=table1$sd_mpi[i]*1.0)
    treated_survival = rnorm(n=n_integer, m=table1$mean_mpi[i]*1.5, s=table1$sd_mpi[i]*1.5)
    survival_data = data.frame(time=c(placebo_survival,treated_survival),event=rep(1,2*n_integer),grp=rep(c('p','t'),each=n_integer))
    sdiff = survdiff(Surv(time, event) ~ grp, data=survival_data)
    pval = 1-pchisq(q=sdiff$chisq,df=1)
    pvals = c(pvals, pval)
    max_survival = max(c(placebo_survival, treated_survival))
    maxes = c(maxes, max_survival)
  }
  table1$expected_duration_bc[i] = round(mean(maxes),0)
  table1$power_bc[i] = mean(pvals < 0.05)
}

# now do the power calculation under the "other" scenario where available = rows 1, 3, 4
for (i in c(1,3,4)) {
  n_integer = 6
  maxes = numeric(0)
  pvals = numeric(0)
  for (j in 1:n_iter) {
    placebo_survival = rnorm(n=n_integer, m=table1$mean[i]*1.0, s=table1$sd[i]*1.0)
    treated_survival = rnorm(n=n_integer, m=table1$mean[i]*1.5, s=table1$sd[i]*1.5)
    survival_data = data.frame(time=c(placebo_survival,treated_survival),event=sample(c(0,1),prob=c(1-table1$p_attack[i],table1$p_attack[i]),size=2*n_integer,replace=T),grp=rep(c('p','t'),each=n_integer))
    sdiff = survdiff(Surv(time, event) ~ grp, data=survival_data)
    pval = 1-pchisq(q=sdiff$chisq,df=1)
    pvals = c(pvals, pval)
    max_survival = max(c(placebo_survival, treated_survival))
    maxes = c(maxes, max_survival)
  }
  table1$expected_duration_oth[i] = round(mean(maxes),0)
  table1$power_oth[i] = mean(pvals < 0.05)
}


disp_cols = c('common_name','strain_roa','shortname','mean_sd_bc','expected_duration_bc','power_bc','attack_rate','mean_sd','expected_duration_oth','power_oth')
table1_disp = table1[,disp_cols]

write.table(table1_disp, 'display_items/table-1.tsv', sep='\t', col.names=T, row.names=F, quote=F, na='')
# clipcopy(table1_disp) # useful in interactive mode






### STATS FOR TEXT

cat(file=stderr(), 'done.\nWriting out statistics for text...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

write(paste0('Proportion of animals IC inoculated: ',sum(data_dedup$n[data_dedup$roa=='ic']),'/',sum(data_dedup$n[!is.na(data_dedup$roa)])),text_stats_path,append=T)

write(paste0('Animals censored versus reaching endpoint: ',sum(outcome_params$n[outcome_params$outcome %in% c('censored','intercurrent')]),' vs. ',outcome_params$n[outcome_params$outcome=='endpoint']),text_stats_path,append=T)

write(paste0('Number of cohorts with N < 4 animals: ',sum(data_dedup$n < 4)),text_stats_path,append=T)
write(paste0('Number of cohorts with N = 1 animal: ',sum(data_dedup$n == 1)),text_stats_path,append=T)
write(paste0('Total number of cohorts: ',sum(!is.na(data_dedup$n))),text_stats_path,append=T)

write(paste0('\n'),text_stats_path,append=T)


### STATS FOR TEXT DONE


elapsed_time = Sys.time() - start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',formatC(as.numeric(elapsed_time), format='f', digits=1),' ',units(elapsed_time),'.\n'))

