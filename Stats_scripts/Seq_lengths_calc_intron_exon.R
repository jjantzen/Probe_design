#Calc intron, exon and supercontig lengths

library(tidyverse)

#import lengths
file_list_JJ <- list.files("Data/Seq_lengths_calcs/JJ/", full.names = TRUE)
file_list_PA <- list.files("Data/Seq_lengths_calcs/PA/", full.names = TRUE)


sorting_seq_counts <- function(file_list){
  file_df <- data.frame(sample = as.character(), length = as.character(), probe = as.character(), stringsAsFactors = FALSE)
  
  for (i in 1:length(file_list)){
    file <- read.delim(file_list[i], header = FALSE)
    colnames(file) = c("sample", "length")
    namestring <- strsplit(file_list[i], "\\/")[[1]][4]
    probe <- gsub("_counts.txt", "", namestring)
    file$probe <- probe
    file_df <- rbind(file_df, file)
  } 
  
  file_df$sample <- as.character(file_df$sample)
  
  file_df$split <- strsplit(file_df$sample, split = "-")  
  
  for (i in 1:nrow(file_df)){
    file_df$sample[i] <- file_df$split[[i]][1]
  }
  
  file_df$split <- strsplit(file_df$probe, split = "_")
  
  for (i in 1:nrow(file_df)){
    file_df$seq_source[i] <- file_df$split[[i]][length(file_df$split[[i]])]
    file_df$locus[i] <- file_df$split[[i]][1]
  }
  

  file_df_complete <- subset(file_df, select=-c(split))
  
  file_df_complete$probe <- gsub("_introns", "", file_df_complete$probe)
  file_df_complete$probe <- gsub("_supercontig", "", file_df_complete$probe)
  
  return(file_df_complete)
  
  
}

file_df_complete_JJ <- sorting_seq_counts(file_list_JJ)
file_df_complete_PA <- sorting_seq_counts(file_list_PA)

file_df_complete_JJ$seq_source[which(file_df_complete_JJ$seq_source != "supercontig" & file_df_complete_JJ$seq_source != "introns")] <- "exon"
file_df_complete_PA$seq_source[which(file_df_complete_PA$seq_source != "supercontig" & file_df_complete_PA$seq_source != "introns")] <- "exon"


write.csv(file_df_complete_JJ, "Output_files/seq_lengths_by_type_JJ.csv", row.names = FALSE)
write.csv(file_df_complete_PA, "Output_files/seq_lengths_by_type_PA.csv", row.names = FALSE)

#combine into one dataframe
#read first

file_df_complete_JJ <- read.csv("Output_files/seq_lengths_by_type_JJ.csv", stringsAsFactors = FALSE)
file_df_complete_PA <- read.csv("Output_files/seq_lengths_by_type_PA.csv", stringsAsFactors = FALSE)

file_df_complete_JJ$clade <- "Tibouchina"
file_df_complete_PA$clade <- "Memecylon"

# #redundant now
# file_df_complete_JJ$probe <- gsub("_introns", "", file_df_complete_JJ$probe)
# file_df_complete_JJ$probe <- gsub("_supercontig", "", file_df_complete_JJ$probe)
# 
# file_df_complete_PA$probe <- gsub("_introns", "", file_df_complete_PA$probe)
# file_df_complete_PA$probe <- gsub("_supercontig", "", file_df_complete_PA$probe)

combo_lengths <- rbind(file_df_complete_PA, file_df_complete_JJ)



#average by probe
probe_averages <- combo_lengths %>% 
  group_by(probe, seq_source, clade) %>% 
  summarize(average_length = mean(length), min_length = min(length), max_length = max(length))

write.csv(probe_averages, "Output_files/Tables/For_formatting/seq_length_avgs.csv", row.names = FALSE)

#fix this part
probe_averages_gathered <- as.data.frame(probe_averages) %>% 
  gather(Measurement, Lengths, average_length:max_length) %>% 
  unite(Variable, seq_source, Measurement, sep = ".") %>% 
  unite(sortable, Variable, clade, sep = ".") %>% 
  spread(sortable, Lengths)


write.csv(probe_averages_gathered, "Output_files/Tables/For_formatting/seq_lengths_avgs_grouped.csv", row.names = FALSE)

get_columns <- function(file){
  for (i in 1:nrow(file)){
    file$split[i] <- str_split(file$probe[i], "_")
    file$source[i] <- file$split[[i]][2]
    file$locus[i] <- file$split[[i]][1]
  }
  
  file2 <- subset(file, select = -split)
  return (file2)
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


str(probe_averages_gathered)

combo_lengths2 <- get_columns(combo_lengths)

str(combo_lengths$probe)

probe_genome_avgs <- grouping_by_source(combo_lengths2)
colnames(probe_genome_avgs)


genome_avgs <- probe_genome_avgs %>% 
  group_by(clade, seq_source, genome) %>% 
  summarize(avg_length = mean(length, na.rm = TRUE), min_length = min(length, na.rm = TRUE), max_length = max(length, na.rm = TRUE))



genome_averages_gathered <- as.data.frame(genome_avgs) %>% 
  gather(Measurement, Lengths, avg_length:max_length) %>% 
  unite(Variable, seq_source, Measurement, sep = ".") %>% 
  unite(sortable, Variable, clade, sep = ".") %>% 
  spread(sortable, Lengths)

write.csv(genome_averages_gathered, "Output_files/Tables/For_formatting/seq_lengths_avgs_grouped.csv", row.names = FALSE)
