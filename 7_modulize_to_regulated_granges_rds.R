
# libraries
library(GenomicRanges)
library(dplyr)


# create a list of datasets that should be analyzed
HeLa_regulated <- readRDS("./data/rds_objects/HeLa_regulated_a5ss_events.rds")
RPE_regulated <- readRDS("./data/rds_objects/RPE_regulated_a5ss_events.rds")
GSE213633_regulated <- readRDS("./data/rds_objects/GSE213633_regulated_a5ss_events.rds")

datasets_regulated <- list(
                           GSE213633=GSE213633_regulated,
                           HeLa=HeLa_regulated,
                           RPE=RPE_regulated
                           )

for (i in 1:length(datasets_regulated)){
  dataset <- datasets_regulated[[i]]
  dataset_name <- names(datasets_regulated)[i]
  
  
  column_name <- grep("_median_dpsi$", colnames(dataset), value = TRUE)
  comps <- sub("_median_dpsi$", "", column_name)
  
  
  # Split the string using the dot as the delimiter
  result <- strsplit(comps, "\\.")
  
  # Extract the individual words
  word1 <- result[[1]][1]
  word2 <- result[[1]][2]
  
  control_median_psi <- paste0(word2,"_median_psi")
  knockout_median_psi <- paste0(word1,"_median_psi")
  dpsiCol <- paste0(comps,"_median_dpsi")
  probCol <- paste0(comps,"_probability_changing")
  
  
  RegulatedA5SSranges <- lapply(split(dataset, dataset$event_id), function(df){
    
    
    df <- df %>%
      dplyr::rename(chr=seqid, exon=spliced_with)
    
    tmpE1 <- df$reference_exon_coord %>% unique %>% strsplit(., "-", fixed=T) %>% unlist %>% as.integer
    tmpE2 <- df$spliced_with_coord %>% unique %>% strsplit(., "-", fixed=T) %>% unlist %>% as.integer
    
    if(c(tmpE1, tmpE2) %>% is.na %>% any){
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
    
    mcols(gr)[[control_median_psi]] <- NA
    mcols(gr)[[control_median_psi]][gr$exon == "E1P"] <- df[[control_median_psi]][df$junction_name=="Proximal"]
    mcols(gr)[[control_median_psi]][gr$exon == "E1D"] <- df[[control_median_psi]][df$junction_name=="Distal"]
    
    mcols(gr)[[knockout_median_psi]] <- NA
    mcols(gr)[[knockout_median_psi]][gr$exon == "E1P"] <- df[[knockout_median_psi]][df$junction_name=="Proximal"]
    mcols(gr)[[knockout_median_psi]][gr$exon == "E1D"] <- df[[knockout_median_psi]][df$junction_name=="Distal"]
    
    mcols(gr)[[dpsiCol]] <- NA
    mcols(gr)[[dpsiCol]][gr$exon == "E1P"] <- df[[dpsiCol]][df$junction_name=="Proximal"]
    mcols(gr)[[dpsiCol]][gr$exon == "E1D"] <- df[[dpsiCol]][df$junction_name=="Distal"]
    
    mcols(gr)[[probCol]] <- NA
    mcols(gr)[[probCol]][gr$exon == "E1P"] <- df[[probCol]][df$junction_name=="Proximal"]
    mcols(gr)[[probCol]][gr$exon == "E1D"] <- df[[probCol]][df$junction_name=="Distal"]
    
    gr$denovo <- NA
    gr$denovo[gr$exon == "E1P"] <- df$denovo[df$junction_name=="Proximal"]
    gr$denovo[gr$exon == "E1D"] <- df$denovo[df$junction_name=="Distal"]
    
    return(gr)
  }) %>% as(., "GRangesList")
  
  filename <- paste0("./data/rds_objects/",dataset_name, "_A5SSranges.rds")
  
  RegulatedA5SSranges <- RegulatedA5SSranges[lengths(RegulatedA5SSranges) == 3]
  saveRDS(RegulatedA5SSranges, filename)
  
 
  print("iteration finished")
}