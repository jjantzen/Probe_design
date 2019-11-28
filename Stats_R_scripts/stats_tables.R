#load libraries
library(ggplot2)
library(reshape2)
library(data.table)
library(ggthemes)
library(ggalt)
library(scales)
library(viridis)
library(RColorBrewer)
library(ggThemeAssist)
library(stringr)
library(tidyr)
library(forcats)
library(gsubfn) 
library(dplyr)

library(stringr)
library(tidyverse)


#First, calculate sequence length stats:

#read file - output from sequencing success stats file

#for all JJ run at once
seq_lengths_JJ <- read.csv("../Output_files/paper_scripts_output/JJ_all_lengths.csv", row.names = 1, stringsAsFactors = FALSE)

#for all PA run at once
seq_lengths_PA <- read.csv("../Output_files/paper_scripts_output/PA_all_lengths.csv", row.names = 1, stringsAsFactors = FALSE)


# #remove bad samples before running stats

seq_lengths_JJ_2 <- seq_lengths_JJ[which(seq_lengths_JJ$Var1 != "A22_Inopinata_clean" & seq_lengths_JJ$Var1 != "120_Rosanae_clean"),]

#change species names
PA_names <- read.csv("../Data/Revision_files/PA/name_list.csv", stringsAsFactors = FALSE)
JJ_names <- read.csv("../Data/JJ_sample_names.csv", stringsAsFactors = FALSE)


for (i in 1:nrow(PA_names)){
  seq_lengths_PA$species[which(seq_lengths_PA$sample == PA_names$HybPiper_Directory_Name[i])] <- PA_names$Species[i]
}

for (i in 1:nrow(JJ_names)){
  seq_lengths_JJ_2$species[which(seq_lengths_JJ_2$Var1 == JJ_names$HybPiper_Directory_Name[i])] <- JJ_names$Species[i]
}

# #get source, probe info 

#write seq length function

seq_lengths_calc <- function(seq_lengths){
  
  colnames(seq_lengths) <- c("sample", "probe", "value", "source", "species")#, "genome", "method")
  
  seq_lengths_calc_output <- c()
  
  for (i in 1:nrow(seq_lengths)){
    seq_lengths$split[i] <- str_split(seq_lengths$probe[i], "_")
    seq_lengths$source[i] <- seq_lengths$split[[i]][2]
    seq_lengths$locus[i] <- seq_lengths$split[[i]][1]
  }
  
  #get avg by sample across all loci
  samp_avg <- seq_lengths %>% group_by(sample) %>% summarize(mean = mean(value), min = min(value), max = max(value), mean_no_zero = mean(value[value!=0]))
  
  #get avg by locus across all samples
  probe_avg <- seq_lengths %>% group_by(probe) %>% summarize(mean = mean(value), min = min(value), max = max(value), mean_no_zero = mean(value[value!=0]))
  
  locus_avg <- seq_lengths %>% group_by(locus) %>% summarize(mean = mean(value), min = min(value), max = max(value), mean_no_zero = mean(value[value!=0]))
  
  #get avg by species across all loci
  species_avg <- seq_lengths %>% group_by(species) %>% summarize(mean = mean(value), min = min(value), max = max(value), mean_no_zero = mean(value[value!=0]))
  
  #get avg by locus and species - not quite
  Locus_samp_avg <- seq_lengths %>% group_by(species, probe) %>% summarize(mean = mean(value), min = min(value), max = max(value), mean_no_zero = mean(value[value!=0]))
  
  #get taxon counts - make sure dplyr loaded last (after plyr)
  
  taxon_count_byprobe <- seq_lengths %>%
    group_by(probe) %>% 
    filter(value > 0) %>% 
    summarise(n = as.integer(n_distinct(species)))
  
  
  taxon_count_bylocus <- seq_lengths %>%
    group_by(locus) %>% 
    filter(value > 0) %>% 
    summarise(n = as.integer(n_distinct(species)))
  #summarize(locus_count = n_distinct(species)) %>% 
  # summarize(avg_locus_count_r = round(mean(locus_count)))
  
  
  #fill in zeros for missing loci
  taxon_count_bylocus_complete <- rbind(taxon_count_bylocus, data.frame(locus = unique(seq_lengths$locus)[-which(unique(seq_lengths$locus) %in% taxon_count_bylocus$locus)], n = as.integer("0")))
  
  taxon_count_bylocus_complete <- taxon_count_bylocus_complete[order(taxon_count_bylocus_complete$locus),]
  
  
  #get counts of success by species
  locus_count_byspecies <- seq_lengths %>%
    filter(value > 0) %>% 
    group_by(species) %>% 
    summarise(n = n_distinct(locus))
  
  
  probe_count_byspecies <- seq_lengths %>%
    filter(value > 0) %>% 
    group_by(species) %>% 
    summarise(n = n_distinct(probe))
  
  
  #overall averages across species - average of species specific counts
  taxon_count_avg <- data.frame(avg_taxa_byprobe = mean(taxon_count_byprobe$n), min_taxa_byprobe = min(taxon_count_byprobe$n), max_taxa_byprobe = max(taxon_count_byprobe$n), avg_taxa_bylocus = mean(taxon_count_bylocus_complete$n), min_taxa_bylocus = min(taxon_count_bylocus_complete$n), max_taxa_bylocus = max(taxon_count_bylocus_complete$n))
  
  probe_locus_avg <- data.frame(avg_probe_byspecies = mean(probe_count_byspecies$n), min_probe_byspecies = min(probe_count_byspecies$n), max_probe_byspecies = max(probe_count_byspecies$n), avg_locus_byspecies = mean(locus_count_byspecies$n), min_locus_byspecies = min(locus_count_byspecies$n), max_locus_byspecies = max(locus_count_byspecies$n))
  
  
  seq_lengths_calc_output[["sample_avg"]] <- c(samp_avg)
  seq_lengths_calc_output[["probe_avg"]] <- c(probe_avg)
  seq_lengths_calc_output[["locus_avg"]] <- c(locus_avg)
  seq_lengths_calc_output[["species_avg"]] <- c(species_avg)
  #seq_lengths_calc_output[["locus_sample_avg"]] <- c(Locus_samp_avg)
  seq_lengths_calc_output[["taxon_count_avg"]] <- c(taxon_count_avg)
  seq_lengths_calc_output[["probe_locus_count_avg"]] <- c(probe_locus_avg)
  seq_lengths_calc_output[["taxon_count_byprobe"]] <- c(taxon_count_byprobe)
  seq_lengths_calc_output[["taxon_count_bylocus"]] <- c(taxon_count_bylocus_complete)
  
  
  return(seq_lengths_calc_output)
  
}


seq_lengths_calc_output_table_PA <- seq_lengths_calc(seq_lengths_PA)
seq_lengths_calc_output_table_JJ <- seq_lengths_calc(seq_lengths_JJ_2)


#Get stats for probes based on hybpiper_stats file
#run for clades separately

#for Prabha all
seq_stats_PA <- read.delim("../Data/Revision_files/PA/01_All_loci_ntd_regular/hybpiper_stats_ntd.txt", row.names = 1, stringsAsFactors = FALSE)

#for JJ all
seq_stats_JJ <- read.delim("../Data/All_JJ/hybpiper_stats_complete.txt", row.names = 1, stringsAsFactors = FALSE)

PA_names <- read.csv("../Data/Revision_files/PA/name_list.csv", stringsAsFactors = FALSE)
JJ_names <- read.csv("../Data/JJ_sample_names.csv", stringsAsFactors = FALSE)


get_stats_calc <- function(seq_stats, names){
  
  species_output <- c()
  seq_stats$species <- NA
  
  
  #get species names for aggregating
  for (i in 1:nrow(names)){
    seq_stats$species[which(rownames(seq_stats) == names$HybPiper_Directory_Name[i])] <- names$Species[i]
  }
  
  
  
  #get percent on target (by sample across all loci (input) and then avg by species across all loci)
  PctOnTarget_species_avg <- seq_stats %>% group_by(species) %>% summarize(mean_PctOnTarget = mean(PctOnTarget))
  
  #then average across all species
  PctOnTarget_avg <- seq_stats %>% summarize(mean_PctOnTarget = mean(PctOnTarget), min_PctOnTarget = min(PctOnTarget), max_PctOnTarget = max(PctOnTarget))
  
  #get number of loci (actually probes) with success (by species and overall)
  #wherever it says "genes" is actually probes - from hybpiper
  num_loci_totals_species <- seq_stats %>% group_by(species) %>% summarize(mean_gene_contigs = mean(GenesWithContigs), max_gene_contigs = max(GenesWithContigs), min_gene_contigs = min(GenesWithContigs), mean_gene_seqs = mean(GenesWithSeqs), max_gene_seqs = max(GenesWithSeqs), min_gene_seqs = min(GenesWithSeqs), mean_genes_mapped = mean(GenesMapped), max_genes_mapped = max(GenesMapped), min_genes_mapped = min(GenesMapped), mean_genesat25 = mean(GenesAt25pct), max_genesat25 = max(GenesAt25pct), min_genesat25 = min(GenesAt25pct), mean_genesat50 = mean(GenesAt50pct), max_genesat50 = max(GenesAt50pct), min_genesat50 = min(GenesAt50pct), mean_genesat75 = mean(GenesAt75pct), max_genesat75 = max(GenesAt75pct), min_genesat75 = min(GenesAt75pct), mean_genesat150 = mean(Genesat150pct), max_genesat150 = max(Genesat150pct), min_genesat150 = min(Genesat150pct))
  
  
  num_loci_totals_all <- seq_stats %>% summarize(mean_gene_contigs = mean(GenesWithContigs), max_gene_contigs = max(GenesWithContigs), min_gene_contigs = min(GenesWithContigs), mean_gene_seqs = mean(GenesWithSeqs), max_gene_seqs = max(GenesWithSeqs), min_gene_seqs = min(GenesWithSeqs), mean_genes_mapped = mean(GenesMapped), max_genes_mapped = max(GenesMapped), min_genes_mapped = min(GenesMapped), mean_genesat25 = mean(GenesAt25pct), max_genesat25 = max(GenesAt25pct), min_genesat25 = min(GenesAt25pct), mean_genesat50 = mean(GenesAt50pct), max_genesat50 = max(GenesAt50pct), min_genesat50 = min(GenesAt50pct), mean_genesat75 = mean(GenesAt75pct), max_genesat75 = max(GenesAt75pct), min_genesat75 = min(GenesAt75pct), mean_genesat150 = mean(Genesat150pct), max_genesat150 = max(Genesat150pct), min_genesat150 = min(Genesat150pct))
  
  #average across all samples
  readdepth_avg <- seq_stats %>% summarize(avg_readdepth = mean(ReadsMapped), avg_totalreads = mean(NumReads))
  
  #average for each species
  readdepth_avg_sp <- seq_stats %>% group_by(species) %>% summarize(avg_readdepth = mean(ReadsMapped), avg_totalreads = mean(NumReads))
  
  #average for each species
  paralog_sp_avgs <- seq_stats %>% group_by(species) %>% summarize(paralogs = mean(ParalogWarnings, na.rm = TRUE))
  
  #average overall (across all samples)
  paralog_avgs <- seq_stats  %>% summarize(paralogs = mean(ParalogWarnings, na.rm = TRUE))
  
  
  species_output[["PctOnTarget_species_avg"]] <- c(PctOnTarget_species_avg)
  
  species_output[["PctOnTarget_avg"]] <- c(PctOnTarget_avg)
  
  species_output[["Num_loci_totals_species"]] <- c(num_loci_totals_species)
  
  species_output[["Num_loci_totals_all"]] <- c(num_loci_totals_all)
  
  species_output[["Readdepth_avg"]] <- c(readdepth_avg)
  
  species_output[["Readdepth_species_avg"]] <- c(readdepth_avg_sp)
  
  species_output[["Paralog_species_avg"]] <- c(paralog_sp_avgs)
  
  species_output[["Paralog_avg"]] <- c(paralog_avgs)
  
  
  return(species_output)
  
}

#remove bad samples and then run function for datasets (samples removed before stats run)
seq_stats_JJ_2 <- seq_stats_JJ[which(rownames(seq_stats_JJ) != "A22_Inopinata_clean" & rownames(seq_stats_JJ) != "120_Rosanae_clean"),]

stats_calc_output_PA <- get_stats_calc(seq_stats_PA, PA_names)
stats_calc_output_JJ <- get_stats_calc(seq_stats_JJ_2, JJ_names)


#Get output in right format
#format output into species tables, probe tables, locus tables and sample tables?


reformat_tables <- function(stats_calc_output, seq_lengths_calc_output_table, newname){
  
  #get stats for locus averages
  locus_based <- data.frame(locus = seq_lengths_calc_output_table$locus_avg$locus, locus_avg_length = seq_lengths_calc_output_table$locus_avg$mean, locus_avg_no_zero_length = seq_lengths_calc_output_table$locus_avg$mean_no_zero, locus_min_length = seq_lengths_calc_output_table$locus_avg$min, locus_max_length = seq_lengths_calc_output_table$locus_avg$max, locus_taxon_count = seq_lengths_calc_output_table$taxon_count_bylocus$n)
  
  probe_based <- data.frame(probe = seq_lengths_calc_output_table$probe_avg$probe, probe_avg_length = seq_lengths_calc_output_table$probe_avg$mean, probe_avg_no_zero_length = seq_lengths_calc_output_table$probe_avg$mean_no_zero, probe_min_length = seq_lengths_calc_output_table$probe_avg$min, probe_max_length = seq_lengths_calc_output_table$probe_avg$max )
  
  short_probe_based <- data.frame(probe = seq_lengths_calc_output_table$taxon_count_byprobe$probe, probe_taxon_count = seq_lengths_calc_output_table$taxon_count_byprobe$n)
  
  combo_probe_based <- merge(probe_based, short_probe_based, by = "probe", all = TRUE)
  
  
  species_based <- data.frame(species = seq_lengths_calc_output_table$species_avg$species, species_avg_length = seq_lengths_calc_output_table$species_avg$mean, species_avg_no_zero_length = seq_lengths_calc_output_table$species_avg$mean_no_zero, species_min_length = seq_lengths_calc_output_table$species_avg$min, species_max_length = seq_lengths_calc_output_table$species_avg$max, PctOnTarget_avg = stats_calc_output$PctOnTarget_species_avg$mean_PctOnTarget, read_depth_avg = stats_calc_output$Readdepth_species_avg$avg_readdepth, total_reads_avg = stats_calc_output$Readdepth_species_avg$avg_totalreads, Gene_contigs_mean = stats_calc_output$Num_loci_totals_species$mean_gene_contigs, Gene_contigs_min = stats_calc_output$Num_loci_totals_species$min_gene_contigs, Gene_contigs_max = stats_calc_output$Num_loci_totals_species$max_gene_contigs, Gene_seqs_mean = stats_calc_output$Num_loci_totals_species$mean_gene_seqs, Gene_seqs_min = stats_calc_output$Num_loci_totals_species$min_gene_seqs, Gene_seqs_max = stats_calc_output$Num_loci_totals_species$max_gene_seqs, Gene_mapped_mean = stats_calc_output$Num_loci_totals_species$mean_genes_mapped, Gene_mapped_min = stats_calc_output$Num_loci_totals_species$min_genes_mapped, Gene_mapped_max = stats_calc_output$Num_loci_totals_species$max_genes_mapped, Gene_25_pct_mean = stats_calc_output$Num_loci_totals_species$mean_genesat25, Gene_25_pct_min = stats_calc_output$Num_loci_totals_species$min_genesat25, Gene_25_pct_max = stats_calc_output$Num_loci_totals_species$max_genesat25, Gene_50_pct_mean = stats_calc_output$Num_loci_totals_species$mean_genesat50, Gene_50_pct_min = stats_calc_output$Num_loci_totals_species$min_genesat50, Gene_50_pct_max = stats_calc_output$Num_loci_totals_species$max_genesat50, Gene_75_pct_mean = stats_calc_output$Num_loci_totals_species$mean_genesat75, Gene_75_pct_min = stats_calc_output$Num_loci_totals_species$min_genesat75, Gene_75_pct_max = stats_calc_output$Num_loci_totals_species$max_genesat75, Gene_150_pct_mean = stats_calc_output$Num_loci_totals_species$mean_genesat150, Gene_150_pct_min = stats_calc_output$Num_loci_totals_species$min_genesat150, Gene_150_pct_max = stats_calc_output$Num_loci_totals_species$max_genesat150, paralog_count_avg = stats_calc_output$Paralog_species_avg$paralogs)
  
  
  clade_based <- data.frame(clade = newname, avg_PctOnTarget = stats_calc_output$PctOnTarget_avg$mean_PctOnTarget, min_PctOnTarget = stats_calc_output$PctOnTarget_avg$min_PctOnTarget, max_PctOnTarget = stats_calc_output$PctOnTarget_avg$max_PctOnTarget, mean_read_depth = stats_calc_output$Readdepth_avg$avg_readdepth, mean_total_reads = stats_calc_output$Readdepth_avg$avg_totalreads, avg_locus_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$avg_locus_byspecies, min_locus_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$min_locus_byspecies, max_locus_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$max_locus_byspecies, avg_probe_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$avg_probe_byspecies, min_avg_probe_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$min_probe_byspecies, max_avg_probe_count_perspecies = seq_lengths_calc_output_table$probe_locus_count_avg$max_probe_byspecies, avg_taxon_count_per_locus = seq_lengths_calc_output_table$taxon_count_avg$avg_taxa_bylocus, min_taxon_count_per_locus = seq_lengths_calc_output_table$taxon_count_avg$min_taxa_bylocus, max_taxon_count_per_locus = seq_lengths_calc_output_table$taxon_count_avg$max_taxa_bylocus, avg_taxon_count_per_probe = seq_lengths_calc_output_table$taxon_count_avg$avg_taxa_byprobe, min_taxon_count_per_probe = seq_lengths_calc_output_table$taxon_count_avg$min_taxa_byprobe, max_taxon_count_per_probe = seq_lengths_calc_output_table$taxon_count_avg$max_taxa_byprobe, num_probe_totals = stats_calc_output$Num_loci_totals_all)
  
  
  #write output tables - change output names
  write.csv(locus_based, paste0("../Output_files/Tables/locus_based_",newname,".csv"))
  write.csv(combo_probe_based, paste0("../Output_files/Tables/probe_based_",newname,".csv"))
  write.csv(species_based, paste0("../Output_files/Tables/species_based_",newname,".csv"))
  write.csv(clade_based, paste0("../Output_files/Tables/clade_based_",newname,".csv"))
  
  str(locus_based)
  for_output <- list()
  
  for_output[["locus_based"]] <- locus_based
  for_output[["probe_based"]] <- combo_probe_based
  for_output[["species_based"]] <- species_based
  for_output[["clade_based"]] <- clade_based
  
  return(for_output)
  #write.csv(locus_based, "../Output_files/Tables/locus_based_Tibouchina.csv")
  #write.csv(combo_probe_based, "../Output_files/Tables/probe_based_Tibouchina.csv")
  #write.csv(species_based, "../Output_files/Tables/species_based_Tibouchina.csv")
  #write.csv(clade_based, "../Output_files/Tables/clade_based_Tibouchina.csv")
  
}

PA_out <- reformat_tables(stats_calc_output_PA, seq_lengths_calc_output_table_PA, "PA")
JJ_out <- reformat_tables(stats_calc_output_JJ, seq_lengths_calc_output_table_JJ, "JJ")


#group by method and by genome source

get_columns <- function(file){
  for (i in 1:nrow(file)){
    file$split[i] <- str_split(file$probe[i], "_")
    file$source[i] <- file$split[[i]][2]
    file$locus[i] <- file$split[[i]][1]
  }
  
  file2 <- subset(file, select = -split)
  return (file2)
}


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

#read files
probe_for_grouping <- read.csv("../Output_files/Tables/For_formatting/combo_probe_based.csv", row.names = 1)

locus_for_grouping <- read.csv("../Output_files/Tables/For_formatting/combo_locus_based.csv", row.names = 1)

#run for functions

probe_for_grouping <- get_columns(probe_for_grouping)
#locus_for_grouping <- get_columns(locus_for_grouping)

probe_grouped <- grouping_by_source(probe_for_grouping)
#locus_grouped <- grouping_by_source(locus_for_grouping)

probe_grouped2 <- grouping_by_method(probe_grouped)
#locus_grouped2 <- grouping_by_method(locus_grouped)


#get aggregated numbers by method and by genome

genome_probe_nums <- probe_grouped2 %>% 
  group_by(genome) %>% 
  summarize(Memecylon_probe_avg_no_zero_length = mean(Memecylon_probe_avg_no_zero_length, na.rm = TRUE), Memecylon_probe_min_length = mean(Memecylon_probe_min_length, na.rm = TRUE), Memecylon_probe_max_length = mean(Memecylon_probe_max_length, na.rm = TRUE), Memecylon_probe_taxon_count = mean(Memecylon_probe_taxon_count, na.rm = TRUE), Tibouchina_probe_avg_no_zero_length = mean(Tibouchina_probe_avg_no_zero_length, na.rm = TRUE), Tibouchina_probe_min_length = mean(Tibouchina_probe_min_length, na.rm = TRUE), Tibouchina_probe_max_length = mean(Tibouchina_probe_max_length, na.rm = TRUE), Tibouchina_probe_taxon_count = mean(Tibouchina_probe_taxon_count, na.rm = TRUE))


method_probe_nums <- probe_grouped2 %>% 
  group_by(method) %>% 
  summarize(Memecylon_probe_avg_no_zero_length = mean(Memecylon_probe_avg_no_zero_length, na.rm = TRUE), Memecylon_probe_min_length = mean(Memecylon_probe_min_length, na.rm = TRUE), Memecylon_probe_max_length = mean(Memecylon_probe_max_length, na.rm = TRUE), Memecylon_probe_taxon_count = mean(Memecylon_probe_taxon_count, na.rm = TRUE), Tibouchina_probe_avg_no_zero_length = mean(Tibouchina_probe_avg_no_zero_length, na.rm = TRUE), Tibouchina_probe_min_length = mean(Tibouchina_probe_min_length, na.rm = TRUE), Tibouchina_probe_max_length = mean(Tibouchina_probe_max_length, na.rm = TRUE), Tibouchina_probe_taxon_count = mean(Tibouchina_probe_taxon_count, na.rm = TRUE))


colnames(probe_grouped2)
#write output files

write.csv(genome_probe_nums, "../Output_files/Tables/For_formatting/combo_probe_based_grouped_genome.csv", row.names = TRUE)

write.csv(method_probe_nums, "../Output_files/Tables/For_formatting/combo_probe_based_grouped_method.csv", row.names = TRUE)


#write.csv(locus_grouped, "../Output_files/Tables/combo_locus_based_grouped.csv", row.names = TRUE)


#for locus
for (i in 2:6){
  colnames(PA_out[[1]])[i] <- paste0("Memecylon_", colnames(PA_out[[1]])[i])
  colnames(JJ_out[[1]])[i] <- paste0("Tibouchina_", colnames(JJ_out[[1]])[i])
}

locus_out <- as.data.frame(cbind(PA_out$locus_based, JJ_out$locus_based), stringsAsFactors = FALSE)

locus_out <- locus_out[,c(1,2:6,8:12)]

write.csv(locus_out, "../Output_files/Tables/For_formatting/combo_locus_based.csv")

#for probe
for (i in 2:6){
  colnames(PA_out[[2]])[i] <- paste0("Memecylon_", colnames(PA_out[[2]])[i])
  colnames(JJ_out[[2]])[i] <- paste0("Tibouchina_", colnames(JJ_out[[2]])[i])
}

probe_out <- as.data.frame(cbind(PA_out$probe_based, JJ_out$probe_based), stringsAsFactors = FALSE)
colnames(probe_out)
probe_out <- probe_out[,c(1,2:6,8:12)]

write.csv(probe_out, "../Output_files/Tables/For_formatting/combo_probe_based.csv")

#for species
# for (i in 2:30){
#   colnames(PA_out[[3]])[i] <- paste0("Memecylon_", colnames(PA_out[[3]])[i])
#   colnames(JJ_out[[3]])[i] <- paste0("Tibouchina_", colnames(JJ_out[[3]])[i])
# }

species_out <- as.data.frame(rbind(PA_out$species_based, JJ_out$species_based), stringsAsFactors = FALSE)
colnames(species_out)
#species_out <- species_out[,c(1,2:6,8:12)]

write.csv(species_out, "../Output_files/Tables/For_formatting/combo_species_based.csv")

#for clade

clade_out <- as.data.frame(rbind(PA_out$clade_based, JJ_out$clade_based), stringsAsFactors = FALSE)
colnames(clade_out)
#probe_out <- probe_out[,c(1,2:6,8:12)]

write.csv(clade_out, "../Output_files/Tables/For_formatting/combo_clade_based.csv")



  