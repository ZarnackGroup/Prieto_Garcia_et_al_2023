library(dplyr)
library(GenomicRanges)

con <- DBI::dbConnect(RSQLite::SQLite(), "data/splicegraph.sql")


junctionReadsTable <- dplyr::tbl(con, "junction_reads") %>%
    as.data.frame %>%
    dplyr::filter(experiment_name %in% c("H_ctrl1Aligned.sortedByCoord.out",
                                         "H_ctrl2Aligned.sortedByCoord.out")) %>%
    group_by(junction_gene_id, junction_start, junction_end) %>%
    summarise(sum_reads = sum(reads, na.rm = TRUE),
              n_experiments = dplyr::n()) %>%
    ungroup()


A5SSranges <- readRDS("rds_files/A5SSranges.rds")

A5SSjunctions <- A5SSranges %>%
    lapply(., function(gr){
        if(strand(gr[1]) %>% as.character == "+"){
            start1 <- end(gr[gr$exon=="E1D"])
            end1 <- start(gr[gr$exon=="E2"])
            start2 <- end(gr[gr$exon=="E1P"])
            end2 <- start(gr[gr$exon=="E2"])
        } else {
            start1 <- end(gr[gr$exon=="E2"])
            end1 <- start(gr[gr$exon=="E1D"])   
            start2 <- end(gr[gr$exon=="E2"])
            end2 <- start(gr[gr$exon=="E1P"])
        }

        return(data.frame(
            event_id = gr$event_id[1],
            start = c(start1, start2),
            end = c(end1, end2),
            gene_id = gr$gene_id[1]
        ))
        
        
    }) %>%
    bind_rows()

A5SSjunctions <- left_join(
    A5SSjunctions,
    junctionReadsTable,
    by=c("start" = "junction_start",
         "end" = "junction_end",
         "gene_id" = "junction_gene_id")
)

A5SSreadSupport <- A5SSjunctions %>% 
    group_by(event_id) %>%
    summarize(readSupport = sum(sum_reads, na.rm = TRUE))


nonRegA5SSranges <- readRDS("rds_files/nonRegA5SSranges.rds")

nonRegA5SSjunctions <- nonRegA5SSranges %>%
    lapply(., function(gr){
        if(strand(gr[1]) %>% as.character == "+"){
            start1 <- end(gr[gr$exon=="E1D"])
            end1 <- start(gr[gr$exon=="E2"])
            start2 <- end(gr[gr$exon=="E1P"])
            end2 <- start(gr[gr$exon=="E2"])
        } else {
            start1 <- end(gr[gr$exon=="E2"])
            end1 <- start(gr[gr$exon=="E1D"])   
            start2 <- end(gr[gr$exon=="E2"])
            end2 <- start(gr[gr$exon=="E1P"])
        }
        
        return(data.frame(
            event_id = gr$event_id[1],
            start = c(start1, start2),
            end = c(end1, end2),
            gene_id = gr$gene_id[1]
        ))
        
        
    }) %>%
    bind_rows()

nonRegA5SSjunctions <- left_join(
    nonRegA5SSjunctions,
    junctionReadsTable,
    by=c("start" = "junction_start",
         "end" = "junction_end",
         "gene_id" = "junction_gene_id")
)

nonRegA5SSreadSupport <- nonRegA5SSjunctions %>% 
    group_by(event_id) %>%
    summarize(readSupport = sum(sum_reads, na.rm = TRUE))


saveRDS(A5SSreadSupport, "rds_files/A5SSreadSupport.rds")
saveRDS(nonRegA5SSreadSupport, "rds_files/nonRegA5SSreadSupport.rds")