library(tidyverse)

setwd('~/src/skeeters') 
contig_df = read_tsv('~/src/skeeters/data/contig_quality_df.tsv')

contig_df %>% mutate(cov = sub('.*\\_','',qseqid)) %>% filter((cov > 2 & qlength > 250 & kingdom != "Arthropoda") | 
                                                                (cov > 3 & qlength > 1000)) %>% select(sample, qseqid) %>%
  write_csv('~/src/skeeters/data/reblast_contigs.csv')

contig_df %>% mutate(cov = sub('.*\\_','',qseqid)) %>% filter((cov > 2 & qlength > 250 & kingdom != "Arthropoda")) %>%
  write_csv('~/src/skeeters/data/reblast_contigs_no_arthopods.csv')