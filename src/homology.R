options(stringsAsFactors=F)

if (interactive()) {
  setwd('~/d/sci/src/nhp_models')
}

result = data.frame(species=c('chimpanzee','cynomolgus','spider monkey','squirrel monkey','gray mouse lemur'),
                    filename=c('human_chimpanzee.aln','human_cynomolgus_rcomp.aln','human_spidermonkey_rcomp.aln','human_squirrelmonkey.aln','human_graymouselemur.aln'),
                    identity_numerator=integer(5),identity_denominator=integer(5),identity_proportion=numeric(5),
                    tile_numerator=integer(5),tile_denominator=integer(5),tile_proportion=numeric(5))

for (r in 1:nrow(result)) {
  
  path = paste0('homology/alignments/',result$filename[r])

  aln = read.fwf(path, skip=38, header=F, widths=c(13,8,50,8), col.names=c('content','start','alignment','end'))
  
  aln = aln[!is.na(aln$content),]
  
  leftseq = aln$content[1]
  rightseq = aln$content[3]
  
  lefttbl = aln[aln$content==leftseq,]
  righttbl = aln[aln$content==rightseq,]
  alntbl = aln[aln$content=="             ",]
  
  alnfull = paste0(alntbl$alignment,collapse='')
  
  # start index is same for left/right
  start_index = min(which(unlist(strsplit(alnfull, split='')) != ' '))
  
  # end index is more complicated - have to figure out what base in left seq it corresponds to
  end_index_in_pair = max(which(unlist(strsplit(alnfull, split='')) != ' '))
  end_row = ceiling(end_index_in_pair / 50)
  end_index_in_left = lefttbl$start[end_row] + end_index_in_pair %% 50
  
  aligned_length_in_left = end_index_in_left - start_index + 1
  
  # trim missing bases at beginning/end - may just be due to different annotation of where transcript starts
  alnsplit = unlist(strsplit(trimws(alnfull, which='both'), split=''))
  
  # identity
  result$identity_numerator[r] = sum(alnsplit=='|')
  result$identity_denominator[r] =aligned_length_in_left
  result$identity_proportion[r] = sum(alnsplit=='|')/aligned_length_in_left
  
  # 20 bp tiling identity
  tile_size = 20
  matched_tiles = 0
  for (i in 1:(length(alnsplit)-tile_size)) {
    if (all(alnsplit[i:(i+tile_size-1)] == '|')) {
      matched_tiles = matched_tiles + 1
    }
  }
  
  potential_tiles = aligned_length_in_left - tile_size
  
  result$tile_numerator[r] = matched_tiles
  result$tile_denominator[r] = potential_tiles
  result$tile_proportion[r] = matched_tiles/potential_tiles
}

write.table(result, 'output/homology_calculations.tsv', sep='\t', col.names=T, row.names=F, quote=F, na='')
