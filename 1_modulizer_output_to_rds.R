library(dplyr)
library(tidyr)

binaryEvents <- list(cassette=read.table("data/modulize/cassette.tsv", header=T, sep="\t"),
              alternative_intron=read.table("data/modulize/alternative_intron.tsv", header=T, sep="\t"),
              alt5prime=read.table("data/modulize/alt5prime.tsv", header=T, sep="\t"),
              alt3prime=read.table("data/modulize/alt3prime.tsv", header=T, sep="\t"),
              alt3and5prime=read.table("data/modulize/alt3and5prime.tsv", header=T, sep="\t"),
              alternate_first_exon=read.table("data/modulize/alternate_first_exon.tsv", header=T, sep="\t"),
              alternate_last_exon=read.table("data/modulize/alternate_last_exon.tsv", header=T, sep="\t"),
              p_alt5prime=read.table("data/modulize/p_alt5prime.tsv", header=T, sep="\t"),
              p_alt3prime=read.table("data/modulize/p_alt3prime.tsv", header=T, sep="\t"),
              p_alternate_first_exon=read.table("data/modulize/p_alternate_first_exon.tsv", header=T, sep="\t"),
              p_alternate_last_exon=read.table("data/modulize/p_alternate_last_exon.tsv", header=T, sep="\t"),
              mutually_exclusive=read.table("data/modulize/mutually_exclusive.tsv", header=T, sep="\t"),
              tandem_cassette=read.table("data/modulize/tandem_cassette.tsv", header=T, sep="\t"),
              multi_exon_spanning=read.table("data/modulize/multi_exon_spanning.tsv", header=T, sep="\t"))

saveRDS(binaryEvents, "rds_files/binaryEventsByModulizer.rds")


