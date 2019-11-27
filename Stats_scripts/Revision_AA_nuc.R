#Revision comparison of aminoacid vs nucleotide
library(dplyr)


#First, read the files. These files are produced by the get_seq_lengths.py and hybpiper_stats.py scripts from HybPiper.

#Tibouchina samples (assembly run at once)
seq_lengths_JJ <- read.delim("Data/Revision_files/JJ/hybpiper_seq_lengths_reseq_aa.txt", row.names = 1, stringsAsFactors = FALSE)

#Memecylon samples (assembly run at once)
seq_lengths_PA <- read.delim("Data/Revision_files/PA/02_All_loci_aminoacid/hybpiper_seq_lengths_aa.txt", row.names = 1, stringsAsFactors = FALSE)

#Remove ITS and ETS from PA set
#seq_lengths_PA <- seq_lengths_PA[,-which(colnames(seq_lengths_PA) %in% c("ITSKC523259_Mrivulare", "ETS_KC523102_Mrivulare"))]

#Make data into matrix format
seq_lengths_matrix_PA <- as.matrix(seq_lengths_PA) #seq_lengths_PA_2 if used list_remove
seq_lengths_matrix_JJ <- as.matrix(seq_lengths_JJ) #seq_lengths_JJ_2 if used list_remove

#Write matrices
write.csv(seq_lengths_matrix_JJ, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/seq_lengths_matrix_JJ_AA.csv", row.names = TRUE)
write.csv(seq_lengths_matrix_PA, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/seq_lengths_matrix_PA_AA.csv", row.names = TRUE)

##Run functions for datasets

#Run for Memecylon percents
all_percents_PA <- get_seq_percent(seq_lengths_matrix_PA)
sorted_all_p_PA <- melt_sort_PA(all_percents_PA)

#omit failed sample from sorted dataframe
#sorted_all_p_PA_2 <- sorted_all_p_PA[which(sorted_all_p_PA$Var1 != "501_Msp_clean"),]

#write and plot percents
write.csv(sorted_all_p_PA, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/PA_all_percents_AA.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_all_p_PA, "/Revision/Aminoacid_vs_nucleotide/PA_all_percents_AA")

#Run for Memecylon lengths
sorted_all_l_PA <- melt_sort_PA(seq_lengths_matrix_PA)

#omit failed sample from sorted dataframe
#sorted_all_l_PA_2 <- sorted_all_l_PA[which(sorted_all_l_PA$Var1 != "501_Msp_clean"),]

#write and plot lengths
write.csv(sorted_all_l_PA, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/PA_all_lengths_AA.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_all_l_PA, "/Revision/Aminoacid_vs_nucleotide/PA_all_lengths_AA")

#Run for Tibouchina percents
all_percents_JJ <- get_seq_percent(seq_lengths_matrix_JJ)
sorted_all_p_JJ <- melt_sort_JJ(all_percents_JJ)

#omit failed samples
sorted_all_p_JJ_2 <- sorted_all_p_JJ[which(sorted_all_p_JJ$Var1 != "A22_Inopinata_clean" & sorted_all_p_JJ$Var1 != "120_Rosanae_clean"),]

#write and plot percents
write.csv(sorted_all_p_JJ_2, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/JJ_all_percents_AA.csv")
plot_heatmap(seq_lengths_matrix_JJ, sorted_all_p_JJ_2, "/Revision/Aminoacid_vs_nucleotide/JJ_all_percents_AA")

#Run for Tibouchina lengths
sorted_all_l_JJ <- melt_sort_JJ(seq_lengths_matrix_JJ)

#omit failed samples
sorted_all_l_JJ_2 <- sorted_all_l_JJ[which(sorted_all_l_JJ$Var1 != "A22_Inopinata_clean" & sorted_all_l_JJ$Var1 != "120_Rosanae_clean"),]

#write and plot lenghts
write.csv(sorted_all_l_JJ_2, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/JJ_all_lengths_AA.csv")
plot_heatmap(seq_lengths_matrix_JJ, sorted_all_l_JJ_2, "/Revision/Aminoacid_vs_nucleotide/JJ_all_lengths_AA")


###comparing nucleotide vs amino acid lengths for specific loci

#read nucleotide files too

#Tibouchina samples (assembly run at once)
seq_lengths_JJ_nuc <- read.delim("Data/All_JJ/hybpiper_seq_lengths_complete_mod.txt", row.names = 1, stringsAsFactors = FALSE)
seq_lengths_matrix_JJ_nuc <- as.matrix(seq_lengths_JJ_nuc) #seq_lengths_JJ_2 if used list_remove


sorted_lengths_JJ_nuc <- melt_sort_JJ(seq_lengths_matrix_JJ_nuc)

str(seq_lengths_JJ_nuc)

#write and plot percents
write.csv(sorted_lengths_JJ_nuc, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/sorted_lengths_JJ_nuc.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_lengths_JJ_nuc, "/Revision/Aminoacid_vs_nucleotide/sorted_lengths_JJ_nuc")



#Memecylon samples (assembly run at once)
seq_lengths_PA_nuc <- read.delim("Data/Revision_files/PA/01_All_loci_ntd_regular/hybpiper_seq_lengths_ntd.txt", row.names = 1, stringsAsFactors = FALSE)
seq_lengths_matrix_PA_nuc <- as.matrix(seq_lengths_PA_nuc)

sorted_lengths_PA_nuc <- melt_sort_JJ(seq_lengths_matrix_PA_nuc)

#write and plot percents
write.csv(sorted_lengths_PA_nuc, "Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/sorted_lengths_PA_nuc.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_lengths_PA_nuc, "/Revision/Aminoacid_vs_nucleotide/sorted_lengths_PA_nuc")

########################
#read previously run analyses
sorted_all_l_JJ_2 <- read.csv("Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/JJ_all_lengths_AA.csv", stringsAsFactors = FALSE, header = TRUE)
sorted_all_l_PA <- read.csv("Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/PA_all_lengths_AA.csv", stringsAsFactors = FALSE, header = TRUE)
sorted_lengths_JJ_nuc <- read.csv("Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/sorted_lengths_JJ_nuc.csv", stringsAsFactors = FALSE, header = TRUE)
sorted_lengths_PA_nuc <- read.csv("Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/sorted_lengths_PA_nuc.csv", stringsAsFactors = FALSE, header = TRUE)


#combine PA and JJ

AA_data <- rbind(sorted_all_l_JJ_2, sorted_all_l_PA)
AA_data$type <- "AA"
Nuc_data <- rbind(sorted_lengths_JJ_nuc, sorted_lengths_PA_nuc)
Nuc_data$type <- "Nuc"

#make dataframe of data with only angiosperm 353 loci for comparison
comparison_df <- rbind(AA_data, Nuc_data)

comparison_df$split <- strsplit(comparison_df$Var2, "_")
for (i in 1:nrow(comparison_df)){
  comparison_df$locus[i] <- comparison_df$split[[i]][1]
}

str(comparison_df)
#to remove split column: comparison_df <- comparison_df[,-8]


#get only Ang353 loci for comparison 
locus_list <- read.csv("./Data/Revision_files/loci_ang353.csv", stringsAsFactors = FALSE, header = FALSE)

for (i in 1:nrow(locus_list)){
  locus_list$V1[i] <- paste0("X", locus_list$V1[i])
}

Ang353_comparision <- comparison_df[which(comparison_df$locus %in% locus_list$V1),]


str(Ang353_comparision)
unique(Ang353_comparision$locus)

#reshapeloci into rows with different versions of baits in different columns
#eg locus AA_SWGX  AA_WWQZ  Nuc_SWGX  Nuc_WWQZ  Nuc_rosid Nuc_Memecylon  Nuc_Tibouchina

#if source not trans or GS, then say rosid
unique(Ang353_comparision$source)
for (i in 1:nrow(Ang353_comparision)){
  if (Ang353_comparision$source[i] %in% c("SWGX", "WWQZ")) {
    Ang353_comparision$source[i] <- Ang353_comparision$source[i]
  } else if (Ang353_comparision$source[i] %in% c("Torricelli", "Affzelli", "Memecylon_combo")){
    Ang353_comparision$source[i] <- "Memecylon"
  } else if (Ang353_comparision$source[i] %in% c("Tibouchina", "Brachyotum")){
    Ang353_comparision$source[i] <- "Tibouchina"
  } else {
  Ang353_comparision$source[i] <- "rosid"
  }
}

Ang353_comparison_df <- Ang353_comparision[,-c(1,3,8)] %>% 
  group_by(Var1, species, locus, type, source) %>% 
  mutate(ind = row_number()) %>% 
  pivot_wider(names_from= c(source, type), values_from = value) %>% 
  select(-ind) %>% ungroup() %>% data.frame()

str(Ang353_comparison_df)

#do a pairwise comparison of aa vs nuc (for each source) and for overall?
AA_Nuc_SWGX_WWQZ <- Ang353_comparison_df[which(!is.na(Ang353_comparison_df$WWQZ_Nuc) | !is.na(Ang353_comparison_df$SWGX_Nuc)),c(1:5,9,10)]

ggplot(AA_Nuc_SWGX_WWQZ[which(AA_Nuc_SWGX_WWQZ$species == "Aegopogon"),])+
  geom_point(aes(x = locus, y = SWGX_AA), colour = "red")+
  geom_point(aes(x = locus, y = SWGX_Nuc), colour = "green")+
  facet_wrap(~Var1, ncol = 5)

AA_Nuc_Tib <- Ang353_comparison_df[which(!is.na(Ang353_comparison_df$Tibouchina_Nuc)),c(1:5,7,8,9)]


unique(Ang353_comparison_df$locus)

colnames(AA_Nuc_Tib)

# #make comparison table
# 
# pairwise_SWGX_WWQZ <- AA_Nuc_SWGX_WWQZ %>% 
#   group_by(Var1) %>% 
#   mutate(SWGX_diff = SWGX_AA - SWGX_Nuc, WWQZ_diff = WWQZ_AA - WWQZ_Nuc) %>% 
#   select(Var1, species, locus, SWGX_diff, WWQZ_diff)
# 
 # pairwise_TibN_AA <- AA_Nuc_Tib %>% 
 #   group_by(Var1) %>% 
 #   mutate(SWGX_diff_T = SWGX_AA - Tibouchina_Nuc, WWQZ_diff_T = WWQZ_AA - Tibouchina_Nuc, SWGX_diff_M = SWGX_AA - Memecylon_Nuc, WWQZ_diff_M = WWQZ_AA - Memecylon_Nuc) %>% 
 #   select(Var1, species, locus, SWGX_diff_T, WWQZ_diff_T, SWGX_diff_M, WWQZ_diff_M)

pairwise_all <- Ang353_comparison_df %>% 
  #group_by(Var1) %>% 
  mutate(SWGX_AA_nuc = SWGX_AA - SWGX_Nuc, WWQZ_AA_nuc = WWQZ_AA - WWQZ_Nuc, SWGX_AA_T = SWGX_AA - Tibouchina_Nuc, WWQZ_AA_T = WWQZ_AA - Tibouchina_Nuc, SWGX_AA_M = SWGX_AA - Memecylon_Nuc, WWQZ_AA_M = WWQZ_AA - Memecylon_Nuc) %>% 
  select(Var1, species, locus, SWGX_AA_nuc, WWQZ_AA_nuc, SWGX_AA_T, WWQZ_AA_T, SWGX_AA_M, WWQZ_AA_M)

Ang353_comparison_df$SWGX_Nuc
unique(pairwise_all$SWGX_AA_nuc)
unique(summary_pairwise$pos_SWGX_AA_T)

#out of...
length(unique(pairwise_all$locus))
length(unique(pairwise_all$Var1))
nrow(pairwise_all)

#get total number of loci/samples where positive vs negative
summary_pairwise <- Ang353_comparison_df %>% 
  rowwise() %>% 
  mutate(SWGX_AA_nuc = SWGX_AA - SWGX_Nuc, WWQZ_AA_nuc = WWQZ_AA - WWQZ_Nuc, SWGX_AA_T = SWGX_AA - Tibouchina_Nuc, WWQZ_AA_T = WWQZ_AA - Tibouchina_Nuc, SWGX_AA_M = SWGX_AA - Memecylon_Nuc, WWQZ_AA_M = WWQZ_AA - Memecylon_Nuc) %>% 
  select(Var1, species, locus, SWGX_AA_nuc, WWQZ_AA_nuc, SWGX_AA_T, WWQZ_AA_T, SWGX_AA_M, WWQZ_AA_M) %>% 
  ungroup() %>%
  group_by(Var1) %>% 
  summarize(pos_SWGX_AA_T = sum(SWGX_AA_T > 0, na.rm = TRUE), pos_SWGX_AA_nuc = sum(SWGX_AA_nuc > 0, na.rm = TRUE), pos_SWGX_AA_M = sum(SWGX_AA_M > 0,na.rm = TRUE))
 # summarize(num_pos = sum(SWGX_AA_nuc > 0, na.rm = TRUE), num_neg = sum(SWGX_AA_nuc < 0, na.rm = TRUE)) 

summary_pairwise$pos_SWGX_AA_nuc[which(!is.na(summary_pairwise$pos_SWGX_AA_nuc))]

Ang353_comparison_df[1000:1200,]
Ang353_comparison_df$SWGX_Nuc[which(!is.na(Ang353_comparison_df$SWGX_Nuc) & !is.na(Ang353_comparison_df$SWGX_AA))]

comp_summary <- summary_pairwise %>% 
  pivot_longer(cols = c(SWGX_AA_nuc, SWGX_AA_T, SWGX_AA_M, WWQZ_AA_nuc, WWQZ_AA_T, WWQZ_AA_M))

comp_raw <- Ang353_comparison_df %>% 
  pivot_longer(cols = -c(Var1, species, locus))

str(comp_summary)

unique(comp_raw$name)


n_fun <- function(x){
  return(data.frame(y = 0.95*log(4200), label = length(x)))
}

pdf("./Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/lengths_boxplot_by_category.pdf")
ggplot(comp_raw)+
  geom_boxplot(aes(x = factor(name), y = value, group = factor(name)))+
  stat_summary(aes(x = factor(name), y = value), fun.data = n_fun, geom = "text", hjust = 0.5)+
  xlab("Category and type of template sequence")+
  ylab("Sequence length (bp)")+
  theme_bw()
dev.off()  

pdf("./Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/lengths_boxplot_by_category_log.pdf")
ggplot(comp_raw)+
  geom_boxplot(aes(x = factor(name), y = log(value), group = factor(name)))+
  stat_summary(aes(x = factor(name), y = log(value)), fun.data = n_fun, geom = "text", hjust = 0.5)+
  xlab("Category and type of template sequence")+
  ylab("Sequence length (bp)")+
  theme_bw()
dev.off()  

write.csv(summary_pairwise, "./Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/summary_stats.txt")

  


#report this in manuscript: 4916 are longer when AA vs nuc and 345 longer when nuc vs AA
#out of how many loci included?
summary_overall <- summary_pairwise %>% 
  #group_by(species) %>% 
  summarize(num_pos_SWGX_AA_nuc = sum(pos_SWGX_AA_nuc), num_pos_SWGX_AA_M = sum(pos_SWGX_AA_M), num_pos_SWGX_AA_T = sum(pos_SWGX_AA_T))

write.csv(summary_pairwise, "./Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/summary_stats.txt")
write.csv(summary_overall, "./Output_files/paper_scripts_output/Revision/Aminoacid_vs_nucleotide/summary_stats_overall.txt")

