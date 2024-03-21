library(tidyverse)
library(GenomicRanges)

binaryEvents <- readRDS("rds_files/binaryEventsByModulizer.rds")

# Filter for junctions that match the criterion and check if all four
#   junctions of a CE event survived the filtering
nonReg <- binaryEvents$alt5prime %>%
    dplyr::filter(abs(KD.CT_median_dpsi) <= 0.025) %>%
    dplyr::count(event_id, name="Frequency") %>%
    dplyr::filter(Frequency==2)

# Add the module_id to the data.frame
nonRegulatedA5SS <- left_join(nonReg,
                              binaryEvents$alt5prime %>%
                                  dplyr::select(module_id, event_id) %>%
                                  dplyr::distinct(),
                              by="event_id")

# Select only one non-regulated CE per AS module
nonRegulatedA5SS <- nonRegulatedA5SS %>%
    group_by(module_id) %>%
    arrange(event_id) %>%
    dplyr::slice(1) %>% 
    pull(event_id)

# Create data.frame
nonRegulatedA5SS <-  binaryEvents$alt5prime %>% 
    dplyr::filter(event_id %in% nonRegulatedA5SS)


nonRegulatedA5SSranges <- lapply(split(nonRegulatedA5SS, nonRegulatedA5SS$event_id), function(df){
    
    df <- df %>%
        dplyr::rename(chr=seqid, exon=spliced_with)
    
    tmpE1 <- df$reference_exon_coord %>% unique %>% strsplit(., "-", fixed=T) %>% unlist %>% as.integer
    tmpE2 <- df$spliced_with_coord %>% unique %>% strsplit(., "-", fixed=T) %>% unlist %>% as.integer
    
    if(c(tmpE1, tmpE2) %>% is.na %>% any){
        return(GRanges())
    }
    
    
    
    # Incorrectly reported coordinates by MAJIQ
    if(df$event_id[1] == "ENSG00000161547.17_1_A5_1"){
        return(GRanges())
    }
    
    gr <- data.frame(chr=df$chr[1],
                     strand=df$strand[1],
                     start=c(tmpE1[1], tmpE2[1]),
                     end=c(tmpE1[2], tmpE2[2])) %>%
        makeGRangesFromDataFrame() %>%
        sort
    
    if(df$strand[1] == "+"){
        gr$exon <- c("E1P", "E2")
    } else {
        gr$exon <- c("E2", "E1P")
    }
    
    gr <- c(gr, gr[gr$exon == "E1P"])
    
    if(df$strand[1] == "+"){
        end(gr)[3] <- df %>% dplyr::filter(junction_name == "Distal") %>% pull(junction_coord) %>% strsplit(., "-", fixed=T) %>% sapply(., "[[", 1) %>% as.integer
        gr$exon[3] <- "E1D"
    } else {
        start(gr)[3] <- df %>% dplyr::filter(junction_name == "Distal") %>% pull(junction_coord) %>% strsplit(., "-", fixed=T) %>% sapply(., "[[", 2) %>% as.integer
        gr$exon[3] <- "E1D"
    }
    
    gr <- gr %>% sort()
    
    gr$module_id <- df$module_id[1]
    gr$event_id <- df$event_id[1]
    gr$lsv_id <- df$lsv_id[1]
    gr$gene_id <- df$gene_id[1]
    gr$gene_name <- df$gene_name[1]
    gr$event_size <- df$event_size[1]
    
    gr$CT_median_psi <- NA
    gr$CT_median_psi[gr$exon == "E1P"] <- df$CT_median_psi[df$junction_name=="Proximal"]
    gr$CT_median_psi[gr$exon == "E1D"] <- df$CT_median_psi[df$junction_name=="Distal"]
    
    gr$KD_median_psi <- NA
    gr$KD_median_psi[gr$exon == "E1P"] <- df$KD_median_psi[df$junction_name=="Proximal"]
    gr$KD_median_psi[gr$exon == "E1D"] <- df$KD_median_psi[df$junction_name=="Distal"]
    
    gr$KD.CT_median_dpsi <- NA
    gr$KD.CT_median_dpsi[gr$exon == "E1P"] <- df$KD.CT_median_dpsi[df$junction_name=="Proximal"]
    gr$KD.CT_median_dpsi[gr$exon == "E1D"] <- df$KD.CT_median_dpsi[df$junction_name=="Distal"]
    
    gr$KD.CT_probability_changing <- NA
    gr$KD.CT_probability_changing[gr$exon == "E1P"] <- df$KD.CT_probability_changing[df$junction_name=="Proximal"]
    gr$KD.CT_probability_changing[gr$exon == "E1D"] <- df$KD.CT_probability_changing[df$junction_name=="Distal"]
    
    gr$denovo <- NA
    gr$denovo[gr$exon == "E1P"] <- df$denovo[df$junction_name=="Proximal"]
    gr$denovo[gr$exon == "E1D"] <- df$denovo[df$junction_name=="Distal"]
    
    return(gr)
}) %>% as(., "GRangesList")
nonRegulatedA5SSranges <- nonRegulatedA5SSranges[lengths(nonRegulatedA5SSranges) == 3]
saveRDS(nonRegulatedA5SSranges, "rds_files/nonRegA5SSranges.rds")

