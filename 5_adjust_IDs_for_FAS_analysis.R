library(tidyverse)
library(Biostrings)

fasDF <- read.table("tsv_files/FAS_pairing.tsv",
                    header=T, sep="\t")
fasDF$Distal_header <- fasDF$Distal_header %>%
    strsplit(., " = ", fixed=T) %>%
    sapply(., "[[", 2) %>%
    strsplit(., " | ") %>%
    sapply(., "[[", 1)

fasDF$Distal_header <- gsub("Distal", "D", fasDF$Distal_header, fixed=TRUE)

fasDF$Proximal_header <- fasDF$Proximal_header %>%
    strsplit(., " = ", fixed=T) %>%
    sapply(., "[[", 2) %>%
    strsplit(., " | ") %>%
    sapply(., "[[", 1)
fasDF$Proximal_header <- gsub("Proximal", "P", fasDF$Proximal_header, fixed=TRUE)

write.table(fasDF,
            "tsv_files/FAS_pairing_adjusted_IDs.tsv",
            col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)


regulatedFASTA <- Biostrings::readAAStringSet(
    "fasta_files/regulated_isoforms.fasta")

names(regulatedFASTA) <- names(regulatedFASTA) %>%
    strsplit(., " = ", fixed=T) %>%
    sapply(., "[[", 2) %>%
    strsplit(., " | ") %>%
    sapply(., "[[", 1)

names(regulatedFASTA) <- gsub("Distal", "D", names(regulatedFASTA), fixed=TRUE)
names(regulatedFASTA) <- gsub("Proximal", "P", names(regulatedFASTA), fixed=TRUE)

writeXStringSet(regulatedFASTA, "fasta_files/regulated_isoforms_adjusted_IDs.fasta")

dups <- names(regulatedFASTA) %>% duplicated 
regulatedFASTA <- regulatedFASTA[!dups]

writeXStringSet(regulatedFASTA, "fasta_files/regulated_isoforms_adjusted_IDs_no_redundant_IDs.fasta")