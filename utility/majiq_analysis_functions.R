# libraries required 
library(dplyr)
library(tidyr)

# Author: Mario Keller
# argument: path to a folder containing .tsv files from a majiq modulize output 
# return: a list containing the contents of the .tsv files of the different events as data frames
loadTSVs <- function(tsvDir){
  return(
    list(cassette=read.table(paste0(tsvDir,"/cassette.tsv"), header=T, sep="\t"),
         alternative_intron=read.table(paste0(tsvDir,"/alternative_intron.tsv"), header=T, sep="\t"),
         alt5prime=read.table(paste0(tsvDir,"/alt5prime.tsv"), header=T, sep="\t"),
         alt3prime=read.table(paste0(tsvDir,"/alt3prime.tsv"), header=T, sep="\t"),
         alt3and5prime=read.table(paste0(tsvDir,"/alt3and5prime.tsv"), header=T, sep="\t"),
         alternate_first_exon=read.table(paste0(tsvDir,"/alternate_first_exon.tsv"), header=T, sep="\t"),
         alternate_last_exon=read.table(paste0(tsvDir,"/alternate_last_exon.tsv"), header=T, sep="\t"),
         p_alt5prime=read.table(paste0(tsvDir,"/p_alt5prime.tsv"), header=T, sep="\t"),
         p_alt3prime=read.table(paste0(tsvDir,"/p_alt3prime.tsv"), header=T, sep="\t"),
         p_alternate_first_exon=read.table(paste0(tsvDir,"/p_alternate_first_exon.tsv"), header=T, sep="\t"),
         p_alternate_last_exon=read.table(paste0(tsvDir,"/p_alternate_last_exon.tsv"), header=T, sep="\t"),
         mutually_exclusive=read.table(paste0(tsvDir,"/mutually_exclusive.tsv"), header=T, sep="\t"),
         tandem_cassette=read.table(paste0(tsvDir,"/tandem_cassette.tsv"), header=T, sep="\t"),
         multi_exon_spanning=read.table(paste0(tsvDir,"/multi_exon_spanning.tsv"), header=T, sep="\t")
    )  
  )
}


# Author: Mario Keller, Felix Haidle
# arguments: 
#           binaryEvents: list with splicing event (output from loadTSVs)
#           dpsi_threshold: threshold of dpsi value to be filtered, default: 0.1
#           propability_threshold: threshold of propability that should be filtered, default: 0.9
# return:   dataframe containing the event class, the gene name, module_id,event_id of each regulated event


identifySignificantEvents <- function(binaryEvents, dpsi_threshold=0.1, propability_threshold=0.9){
  
  significantEvents <- lapply(names(binaryEvents), function(eventClass){
    
    eventsClassDf <- binaryEvents[[eventClass]]
    
    # Extract the part before the first underscore in the last column header
    last_col_name <- names(eventsClassDf)[ncol(eventsClassDf)]
    
    
    # building data set specific column names for comparison
    comps <- sub("^(.*?)_.*", "\\1", last_col_name)
    
    lapply(comps, function(comp){
      probCol <- paste0(comp,"_probability_changing")
      dpsiCol <- paste0(comp,"_median_dpsi")
      
      df <- eventsClassDf[,c("gene_name", "module_id", "event_id", "lsv_id", probCol, dpsiCol)]
      
      # filteringthe events by the conditions (found in supplementary methods part of the paper in the setciton "Identification of differentially regulated splicing events")
      df <- df %>%
        dplyr::filter(., !is.na(get(probCol)) & !is.na(get(dpsiCol))) %>%
        dplyr::filter(., !lsv_id == "") %>% 
        group_by(event_id, lsv_id) %>%
        summarise(., gene_name = unique(gene_name), module_id = unique(module_id),
                  significant= sum(get(probCol) >= propability_threshold) == 2,
                  change = sum(abs(get(dpsiCol)) >= dpsi_threshold) == 2,
                  opposite = sum(sign(get(dpsiCol))) == 0,
                  fraction = sum(abs(get(dpsiCol)/max(abs(get(dpsiCol)))) >= 0.5) == 2) %>%
        summarise(gene_name=unique(gene_name), module_id=unique(module_id),
                  significant = any(significant),
                  change=all(change),
                  opposite=all(opposite),
                  fraction=all(fraction))
      
      if(eventClass == "multi_exon_spanning"){
        df <- df %>% dplyr::filter(significant, change, fraction) %>%
          dplyr::select(gene_name, module_id, event_id)
      } else {
        df <- df %>% dplyr::filter(significant, change, opposite, fraction) %>%
          dplyr::select(gene_name, module_id, event_id)
      }
      
      if(nrow(df) == 0){
        df <- data.frame(eventClass=character(), gene_name=character(), module_id=character(), event_id=character(), placeholder=logical())
        colnames(df)[5] <- comp
      } else {
        df$placeholder <- TRUE
        colnames(df)[4] <- comp
        df <- cbind(eventClass = eventClass, df)
      }
      
      return(df)
    }) %>%
      purrr::reduce(full_join, by=c("eventClass", "gene_name", "module_id", "event_id")) %>%
      replace(is.na(.), FALSE)
  }) %>% bind_rows()
  return(significantEvents)
}
  