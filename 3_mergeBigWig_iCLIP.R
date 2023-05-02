library(dplyr)
library(GenomicRanges)
library(rtracklayer)

#Load the BigWig Files as Rle

rep1_p <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep1L.v4.duprm.plus.bw", as = "Rle")
rep1_m <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep1L.v4.duprm.minus.bw", as = "Rle")

rep2_p <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep2L.v4.duprm.plus.bw", as = "Rle")
rep2_m <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep2L.v4.duprm.minus.bw", as = "Rle")

rep3_p <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep3L.v4.duprm.plus.bw", as = "Rle")
rep3_m <- import("data/iCLIP/uf_muellermcnicoll_2022_02_HEK293T_USP39_rep3L.v4.duprm.minus.bw", as = "Rle")

rep1_p <- rep1_p[names(rep1_p) %in% standardChromosomes(rep1_p)]
rep1_m <- rep1_m[names(rep1_m) %in% standardChromosomes(rep1_m)]

rep2_p <- rep2_p[names(rep2_p) %in% standardChromosomes(rep2_p)]
rep2_m <- rep2_m[names(rep2_m) %in% standardChromosomes(rep2_m)]

rep3_p <- rep3_p[names(rep3_p) %in% standardChromosomes(rep3_p)]
rep3_m <- rep3_m[names(rep3_m) %in% standardChromosomes(rep3_m)]

#IMPORTANT: I checked the order of chrs in the RleList

combined_p <- rep1_p + rep2_p + rep3_p
combined_m <- rep1_m + rep2_m + rep3_m

export(combined_p, con = "data/iCLIP/HEK293T_USP39_combined_L_plus.bw", format = "BigWig")
export(combined_m, con = "data/iCLIP/HEK293T_USP39_combined_L_minus.bw", format = "BigWig")
