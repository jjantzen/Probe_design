#script for calculating percent similarity

library(stringr)
library(ape)
library(dplyr)
library(tidyverse)
library(reshape2)

#Percent Identity = (Matches x 100)/Length of aligned region (with gaps)

#make dataframe
stat_df_full_PA <- data.frame(probe = as.character(), sample = as.character(), id_fraction = as.character(), id_percent = as.character(), sim_fraction = as.character(), sim_percent = as.character(), stringsAsFactors = FALSE)
stat_df_full_JJ <- data.frame(probe = as.character(), sample = as.character(), id_fraction = as.character(), id_percent = as.character(), sim_fraction = as.character(), sim_percent = as.character(), stringsAsFactors = FALSE)

#alignment in wrong format?
format_sim_files <- function(text_file){
  
  stats <- text_file[22:23,1]
  
  name <- as.character(text_file[4,])
  
  name_sep <- strsplit(name, "\\s+")[[1]][3]
  
  name_sep2 <- strsplit(name_sep, "/")
  
  probe <- name_sep2[[1]][1]
  sample <- gsub(".fa", "", name_sep2[[1]][2])
  
  ident <- as.character(stats[1])
  sim <- as.character(stats[2])
  
  sim_sep <- strsplit(sim, "\\s+")
  id_sep <- strsplit(ident, "\\s+")
  
  stat_df <- data.frame(probe = as.character(probe), sample = as.character(sample), id_fraction = as.character(id_sep[[1]][3]), id_percent = as.character(id_sep[[1]][4]), sim_fraction = as.character(sim_sep[[1]][3]), sim_percent = as.character(sim_sep[[1]][4]), stringsAsFactors = FALSE)
  
  stat_df[1,] <- gsub("\\(", "", x = stat_df[1,])
  stat_df[1,] <- gsub("\\)", "", x = stat_df[1,])
  stat_df[1,] <- gsub("%", "", x = stat_df[1,])
  
  stat_df$id_percent <- as.numeric(stat_df$id_percent)
  stat_df$sim_percent <- as.numeric(stat_df$sim_percent)
  
  return(stat_df)
}


#read files

list_of_files_JJ <- list.files("Data/Percent_similarity/Needle_JJ/", full.names = TRUE)
list_of_files_PA <- list.files("Data/Percent_similarity/Needle_PA/", full.names = TRUE)

#run function for all files
for (i in 1:length(list_of_files_JJ)){
  text_file <- read.delim(list_of_files_JJ[i])
  stat_df <- format_sim_files(text_file)
  stat_df_full_JJ <- rbind(stat_df_full_JJ, stat_df)
}

write.csv(stat_df_full_JJ, "Output_files/All_files/JJ_percents_full.csv", row.names = TRUE)


for (i in 1:length(list_of_files_PA)){
  text_file <- read.delim(list_of_files_PA[i])
  stat_df <- format_sim_files(text_file)
  stat_df_full_PA <- rbind(stat_df_full_PA, stat_df)
}

write.csv(stat_df_full_PA, "Output_files/All_files/PA_percents_full.csv", row.names = TRUE)


#get percent similarities by probe and by locus etc 

PA_percents <- read.csv("Output_files/All_files/PA_percents_full.csv", row.names =1, stringsAsFactors = FALSE)
JJ_percents <- read.csv("Output_files/All_files/JJ_percents_full.csv", row.names =1, stringsAsFactors = FALSE)

head(JJ_percents)

sort(unique(JJ_percents$probe))

#get right probe source columns
get_columns <- function(file){
  for (i in 1:nrow(file)){
    file$split[i] <- str_split(file$probe[i], "_")
    file$source[i] <- file$split[[i]][2]
    file$locus[i] <- file$split[[i]][1]
  }
  
  file2 <- subset(file, select = -split)
  return (file2)
}


PA_percent_columns <- get_columns(PA_percents)

JJ_percent_columns <- get_columns(JJ_percents)


#add grouping columns
grouping_by_method <- function(file){
  
  file$method <- NA
  
  file <- file[order(file$probe),]
  
  file$probe <- as.character(file$probe)
  
  for (i in 1:nrow(file)){
    if (startsWith(file$probe[i], "AJ") | startsWith(file$probe[i], "NM")){
      file$method[i] <- "FN"
    } else if (startsWith(file$probe[i], "AT" ) | startsWith(file$probe[i], "TC")) {
      file$method[i] <- "MM"
    } else if (startsWith(file$probe[i], "KT")) {
      file$method[i] <- "SCN"
    } else {
      file$method[i] <- "Ang353"
    }
  }
  return (file)
}

grouping_by_source <- function(file){
  
  file$genome <- NA
  
  file <- file[order(file$probe),]
  
  for (i in 1:nrow(file)){
    if (startsWith(file$source[i], "Tib") | startsWith(file$source[i], "Brachy")){
      file$genome[i] <- "Tibouchina"
    } else if (startsWith(file$source[i], "Meme" ) | startsWith(file$source[i], "Torr") | startsWith(file$source[i], "Aff")) {
      file$genome[i] <- "Memecylon"
    } else if (startsWith(file$source[i], "Trans") | startsWith(file$source[i], "Tetra")) {
      file$genome[i] <- "Transcriptome"
    } else if (startsWith(file$source[i], "Micon")) {
      file$genome[i] <- "Miconia"
    } else {
      file$genome[i] <- "Other_Rosid"
    }
  }
  return (file)
}


grouped_PA_percents <- grouping_by_source(PA_percent_columns)
grouped_JJ_percents <- grouping_by_source(JJ_percent_columns)

grouped_PA_percents2 <- grouping_by_method(grouped_PA_percents)
grouped_JJ_percents2 <- grouping_by_method(grouped_JJ_percents)


grouped_JJ_percents2[which(grouped_JJ_percents2$probe == "KT377070.1_Tibouchina" | grouped_JJ_percents2$probe == "KT377070_Tibouchina"),]  %>% 
  group_by(sample) %>% 
  select(probe, id_percent) %>% 
  tail() 

#add distinguishing column
grouped_PA_percents2$clade <- "Memecylon"
grouped_JJ_percents2$clade <- "Tibouchina"

combo_percents <- rbind(grouped_PA_percents2, grouped_JJ_percents2)

avg_by_probes <- combo_percents %>% 
  group_by(genome, clade) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_probes_combo <- combo_percents %>% 
  group_by(probe, clade) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_locus_combo <- combo_percents %>% 
  group_by(locus, clade) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_probes_genome_combo <- combo_percents %>% 
  group_by(probe, clade, genome) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_probes_method_combo <- combo_percents %>% 
  group_by(probe, clade, method) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_locus_genome_combo <- combo_percents %>% 
  group_by(locus, clade, genome) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

avg_by_locus_method_combo <- combo_percents %>% 
  group_by(locus, clade, method) %>% 
  summarize(avg_ID_pct = mean(id_percent, na.rm = TRUE), min_ID_pct = min(id_percent, na.rm = TRUE), max_ID_pct = max(id_percent, na.rm = TRUE), avg_sim_pct = mean(sim_percent, na.rm = TRUE), min_sim_pct = min(sim_percent, na.rm = TRUE), max_sim_pct = max(sim_percent, na.rm = TRUE))

#and then reformat table to get 
#first gather then spread
gather_probe_avg <- avg_by_probes_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_locus_avg <- avg_by_locus_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_probe_genome_avg <- avg_by_probes_genome_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_locus_genome_avg <- avg_by_locus_genome_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_probe_method_avg <- avg_by_probes_method_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_locus_method_avg <- avg_by_locus_method_combo %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)

gather_genome_avg <- avg_by_probes %>% 
  gather(Measurement, Value, avg_ID_pct:max_sim_pct) %>% 
  unite(Variable, clade, Measurement, sep = ".") %>% 
  spread(Variable, Value)


probe_avg_df <- as.data.frame(gather_probe_avg)
locus_avg_df <- as.data.frame(gather_locus_avg)
locus_genome_avg_df <- as.data.frame(gather_locus_genome_avg)
locus_method_avg_df <- as.data.frame(gather_locus_method_avg)
probe_genome_avg_df <- as.data.frame(gather_probe_genome_avg)
probe_method_avg_df <- as.data.frame(gather_probe_method_avg)

genome_avg_df <- as.data.frame(gather_genome_avg)


probe_avg_df[which(probe_avg_df$probe == "KT377070.1_Tibouchina" | probe_avg_df$probe == "KT377070_Tibouchina"),]
  


head(probe_genome_avg_df)

#write files

write.csv(combo_percents, "Output_files/Tables/For_formatting/Percent_sim_id_combo_expanded.csv", row.names = TRUE)
write.csv(probe_avg_df, "Output_files/paper_scripts_output/Revision/Percent_sim_id_combo_by_probes1.csv", row.names = TRUE)
write.csv(locus_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_locus.csv", row.names = TRUE)
write.csv(probe_method_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_probe_method.csv", row.names = TRUE)
write.csv(probe_genome_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_probe_genome.csv", row.names = TRUE)
write.csv(locus_method_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_locus_method.csv", row.names = TRUE)
write.csv(locus_genome_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_locus_genome.csv", row.names = TRUE)

write.csv(genome_avg_df, "Output_files/Tables/For_formatting/Percent_sim_id_combo_by_genome.csv", row.names = TRUE)


