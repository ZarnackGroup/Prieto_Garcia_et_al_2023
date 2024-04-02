# filtering for regulated events

#reading in functions
source("./utility/majiq_analysis_functions.R")


#libraries
library(dplyr)


# create a list of datasets that should be analyzed
directories <- c(HeLa="./data/modulize/HeLa",
                 RPE= "./data/modulize/RPE",
                 Zebrafish="./data/modulize/Zebrafish",
                 GSE213633="./data/modulize/GSE213633")


for (i in 1:length(directories)) {
  
  # reading in the data 
  dataset_name <- names(directories)[i]
  unfiltered <- loadTSVs(directories[i])


  # find the regulated events
  filtered <- identifySignificantEvents(unfiltered)
  ids <- filter(filtered, eventClass == "alt5prime")$event_id
  
  
  last_col_name <- tail(names(unfiltered[[3]]),n=1)
    
    comps <- sub("^(.*?)_.*", "\\1", last_col_name)
    
    # Split the string using the dot as the delimiter
    result <- strsplit(comps, "\\.")
    
    # Extract the individual words
    word1 <- result[[1]][1]
    word2 <- result[[1]][2]
    dpsiCol <- paste0(comps,"_median_dpsi")
    
  # save wether the junction is up or down regulated
  regulated_full_data <- filter(unfiltered[[3]], event_id %in% ids) %>%
    mutate(
      regulation = case_when(
        !!sym(dpsiCol) > 0 ~ "up",
        !!sym(dpsiCol) < 0 ~ "down",
        TRUE ~ "not"
      ))%>%
    mutate(annotated=ifelse(denovo == "True", "cryptic", "annotated"))
  
  # adding extra columns for junction coordinates to use with granges ater on
  regulated_full_data$start <- as.integer(sub("-.*", "", regulated_full_data$junction_coord))
  regulated_full_data$end <- as.integer(sub(".*-", "", regulated_full_data$junction_coord))
  
  # saving which coordinate is the canonic and which is the cryptic splice site
  regulated_full_data$junction_splice_position_changing <- ifelse(regulated_full_data$strand == '+',
                                                                  regulated_full_data$start, regulated_full_data$end)
  regulated_full_data$junction_splice_position_same <- ifelse(regulated_full_data$strand == '-',
                                                              regulated_full_data$start, regulated_full_data$end)
  
  # saving annotation status for cryptic splice site
  regulated_full_data <- regulated_full_data %>%
    arrange(module_id, row_number()) %>%
    group_by(module_id) %>%
    mutate(alternative_splice_site = ifelse(row_number() %% 2 == 1, lead(junction_splice_position_changing), lag(junction_splice_position_changing)))%>%
    mutate(alternative_splice_site_annotated = ifelse(row_number() %% 2 == 1, lead(annotated), lag(annotated)))
   
   
   
  filename <- paste0("./data/rds_objects/",dataset_name, "_regulated_a5ss_events.rds")
  loaded_df <- saveRDS(regulated_full_data,file = filename)
  }
  
 
  

  
  
  

