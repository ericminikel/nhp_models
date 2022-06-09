options(stringsAsFactors=F)

if (interactive()) {
  setwd('~/d/sci/src/nhp_models')
}

result = data.frame(species=c('chimpanzee','cynomolgus','spider monkey','squirrel monkey','gray mouse lemur'),
                    dna_filename=c('human_chimpanzee.aln','human_cynomolgus_rcomp.aln','human_spidermonkey_rcomp.aln','human_squirrelmonkey.aln','human_graymouselemur.aln'),
                    prot_filename=c('human_chimpanzee.aln','human_cynomolgus.aln','human_spidermonkey.aln','human_squirrelmonkey.aln','human_graymouselemur.aln'),
                    identity_numerator=integer(5),identity_denominator=integer(5),identity_proportion=numeric(5),
                    tile_numerator=integer(5),tile_denominator=integer(5),tile_proportion=numeric(5),
                    protein_numerator=integer(5),protein_denominator=integer(5),protein_proportion=numeric(5))

for (r in 1:nrow(result)) {
  
  ## dna alignments
  
  dna_path = paste0('homology/dna/alignments/',result$dna_filename[r])

  aln = read.fwf(dna_path, skip=38, header=F, widths=c(13,8,50,8), col.names=c('content','start','alignment','end'))
  
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
  
  ## protein alignments
  
  protein_path = paste0('homology/protein/alignments/',result$prot_filename[r])
  aln = read.fwf(protein_path, skip=38, header=F, widths=c(13,8,50,8), col.names=c('content','start','alignment','end'))
  
  aln = aln[!is.na(aln$content),]
  
  leftseq = aln$content[1]
  rightseq = aln$content[3]
  
  lefttbl = aln[aln$content==leftseq,]
  righttbl = aln[aln$content==rightseq,]
  alntbl = aln[aln$content=="             ",]
  
  alnleft = paste0(lefttbl$alignment, collapse='')
  alnright = paste0(righttbl$alignment, collapse='')
  alnfull = paste0(alntbl$alignment,collapse='')
  alnsymbols = strsplit(alnfull,split='')[[1]]
  alnleftchars = strsplit(alnleft,split='')[[1]]
  alnindices = numeric(nchar(alnfull))
  human_codon_number = 1
  for (i in 1:length(alnindices)) {
    if (substr(alnleft,i,i) == '-') {
      alnindices[i] = human_codon_number
    } else {
      alnindices[i] = human_codon_number
      human_codon_number = human_codon_number + 1
    }
  }
  min_aa = 23
  max_aa = 230
  # amino acids counted for human
  denominator = max_aa - min_aa + 1
  # exact matches in pairwise alignment
  matches = sum(alnsymbols[alnindices >= min_aa & alnindices <= max_aa]=='|')
  # insertions in nhp relative to human - count each as just one
  insertions = str_count(alnleft,'-+')
  # count a whole 9-residue deletion of first OPRI as just 1 mismatch
  if (grepl(paste0(rep(' ',9),collapse=''),alnfull)) {
    matches = matches + 8
  }
  
  result$protein_numerator[r] = matches - insertions
  result$protein_denominator[r] = denominator
  result$protein_proportion[r] = matches/denominator
  
}

write.table(result, 'output/homology_calculations.tsv', sep='\t', col.names=T, row.names=F, quote=F, na='')
