#Script for analysing assembly stats and creating figures for publication

#Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
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
library(gridExtra)
library(gridBase)
library(grid)
library(egg)

#Write function to get sequence length as percent of reference sequence

get_seq_percent <- function(sample_data){
  #get sample seq lengths
  sample.len <- sample_data[2:nrow(sample_data),]
  
  #get reference seq length - mean length
  reference.len <- as.numeric(sample_data[1,])
  
  #calculate the percentage length recovered relative to the reference. 
  percent.len <- sweep(sample.len,2,as.numeric(reference.len),'/')
  
  #if longer than reference, set to 1 
  percent.len <- ifelse(percent.len>1,1,percent.len)
  
}

#Write function to get create new dataframe of just the probe names

get_probe_names <- function(matrix){
  #get probe names into new dataframe
  probe_names <- data.frame(dimnames(matrix)[2], stringsAsFactors = FALSE)
  colnames(probe_names) <- c("total")
  
  #split probe name string based on period (.) 
  #first part of string is probe number (probe)
  #second part of string is source of probe (source)
  
  for (i in 1:nrow(probe_names)){
    probe_names$total[i] <- sub('_', '.', probe_names$total[i])
    probe_names$split[i] <- str_split(probe_names$total[i], "\\.")
    probe_names$source[i] <- probe_names$split[[i]][2]
    probe_names$probe[i] <- probe_names$split[[i]][1]
  }
  
  #sort by source
  probe_names <- probe_names[order(probe_names$source, probe_names$probe),]
  
  #remove split column
  probe_names_c <- probe_names[,-2]
  
}

#Write function for Tibouchina samples to convert dataframe into melted format, set variables as factors, and sort and reformat for plotting 
#Data is the read length/percent dataframe (as opposed to number of read dataframe)
#If data is raw read length data (didn't run get_seq_percent function) - uses if statement to remove mean lengths

melt_sort_JJ <- function(data){
  
  if(max(data) <=1){
    #get sample percents
    data <- data
  } else {
    #omit mean lengths
    data <- data[2:nrow(data),]
  }
  
  #melt dataset
  data.long <- melt(data)
  
  colnames(data.long) <- c("Var1", "Var2", "value")
  
  #set sample as factor for plotting
  data.long$Var1 <- as.factor(data.long$Var1)
  
  #get list of gene names for sorting
  gene_names <- unique(data.long$Var2) %>% 
    as.character() %>% 
    sort()
  
  #run function get_probe_names 
  probe_names_c <- get_probe_names(data)
  
  #add source to dataset
  data.long$source <- NA
  for (i in 1:nrow(data.long)){
    for (j in 1:nrow(probe_names_c)){
      if (sub('_', '.', data.long$Var2[i]) == probe_names_c$total[j]){
        data.long$source[i] <- probe_names_c$source[j]
      }
    }
  }
  
  #correct spelling differences
  data.long$source[which(data.long$source == "Affzeli")] <- "Affzelli"
  data.long$source[which(data.long$source == "Torricelil")] <- "Torricelli"
  data.long$source[which(data.long$source == "Torricellii")] <- "Torricelli"
  data.long$source[which(data.long$source == "Tibouchina1")] <- "Tibouchina"
  data.long$source[which(data.long$source == "Tibouchina2")] <- "Tibouchina"
  data.long$source[which(data.long$source == "Memecyon")] <- "Memecylon"
  data.long$source[which(data.long$source == "1_Tibouchina")] <- "Tibouchina"
  data.long$source[which(data.long$source == "1_Brachyotum")] <- "Brachyotum"
  data.long$source[which(data.long$source == "1")] <- "Miconia"
  data.long$source[which(data.long$source == "1_Miconia")] <- "Miconia"
  data.long$source[which(data.long$source == "Brachyotum_appended")] <- "Brachyotum"
  data.long$source[which(data.long$source == "Tibouchina_appended")] <- "Tibouchina"
  
  sources <- unique(data.long$source) %>% 
    as.character() %>% 
    sort()
  
  #set factor order to be gene_names order
  data.long$Var2 <- factor(data.long$Var2, levels = gene_names)
  data.long$source <- factor(data.long$source, levels = sources)
  
  #sort data.long by sample and by locus name
  sorted <- data.long[order(as.character(data.long$Var1, data.long$Var2)),]
  
  #set variable value as numeric (percent or length)
  sorted$value <- as.numeric(sorted$value)
  
  #get species name from sample name
  species <- data.frame(species = as.character(sorted$Var1))
  unique(species)
  #split sample name
  species <- species %>% 
    separate(species, c("ID", "Species","Modifier", "rest"))
  
  #get species names
  for (i in 1:nrow(species)){
    if (species$Modifier[i] != "clean"){
      species$Name[i] <- paste(species$Species[i], sep = "_", species$Modifier[i])
    } else {
      species$Name[i] <- species$Species[i]
    }
  }
  
  #assign species to new column in dataset
  sorted$species <- species$Name
  sorted$Var1 <- as.character(sorted$Var1)
  
  #sort by species
  sorted_sp <- sorted[order(as.character(sorted$species, sorted$Var1)),]
}

#Write different function for Memecylon samples to melt, reshape, sort and reformat data into correct format for plotting
#Data is the read length/percent data
#If data is raw read length data (didn't run get_seq_percent function) - uses if statement to remove mean lengths

melt_sort_PA <- function(data){
  
  if(max(data) <=1){
    #get sample percents
    data <- data
  } else {
    #omit mean lengths
    data <- data[2:nrow(data),]
  }
  
  #melt dataset
  data.long <- melt(data)
  
  colnames(data.long) <- c("Var1", "Var2", "value")
  
  #set sample as factor for plotting
  data.long$Var1 <- as.factor(data.long$Var1)
  
  #get list of gene names for sorting
  gene_names <- unique(data.long$Var2) %>% 
    as.character() %>% 
    sort()
  
  #run function get_probe_names 
  probe_names_c <- get_probe_names(data)
  
  #add source to dataset
  data.long$source <- NA
  for (i in 1:nrow(data.long)){
    for (j in 1:nrow(probe_names_c)){
      if (sub('_', '.', data.long$Var2[i]) == probe_names_c$total[j]){
        data.long$source[i] <- probe_names_c$source[j]
      }
    }
  }
  
  #correct spelling differences
  data.long$source[which(data.long$source == "Affzeli")] <- "Affzelli"
  data.long$source[which(data.long$source == "Torricelil")] <- "Torricelli"
  data.long$source[which(data.long$source == "Torricellii")] <- "Torricelli"
  data.long$source[which(data.long$source == "Tibouchina1")] <- "Tibouchina"
  data.long$source[which(data.long$source == "Tibouchina2")] <- "Tibouchina"
  data.long$source[which(data.long$source == "Memecyon")] <- "Memecylon"
  data.long$source[which(data.long$source == "1_Tibouchina")] <- "Tibouchina"
  data.long$source[which(data.long$source == "1_Brachyotum")] <- "Brachyotum"
  data.long$source[which(data.long$source == "1")] <- "Miconia"
  data.long$source[which(data.long$source == "1_Miconia")] <- "Miconia"
  
  data.long$source[which(data.long$source == "Brachyotum_appended")] <- "Brachyotum"
  data.long$source[which(data.long$source == "Tibouchina_appended")] <- "Tibouchina"
  
  sources <- unique(data.long$source) %>% 
    as.character() %>% 
    sort()
  
  #set factor order to be gene_names order
  data.long$Var2 <- factor(data.long$Var2, levels = gene_names)
  data.long$source <- factor(data.long$source, levels = sources)
  
  #sort data.long by sample and by locus name
  sorted <- data.long[order(as.character(data.long$Var1, data.long$Var2)),]
  
  #set variable value as numeric (percent or length)
  sorted$value <- as.numeric(sorted$value)
  
  #get species name from sample name
  species_df <- read.csv("Data/Revision_files/PA/name_list.csv", stringsAsFactors = FALSE)
  sorted$species <- NA
  
  for (i in 1: nrow(species_df)){
    #assign species to new column in dataset
    sorted$species[which(sorted$Var1 == species_df$Sample_name[i])] <- species_df$Species_name[i]
  }
  
    #sort by species
  sorted_sp <- sorted[order(as.character(sorted$species, sorted$Var1)),]
}


#Write plotting function for plotting heatmap of length/percent sorted by length for gene and by species for sample
#Data is the percent or length object
#Sorted_data is the sorted data object (output of melt_sort function) to be plotted
#Name is output name for saving (string)

plot_heatmap <- function(data, sorted_data, name){
  
  #get data frame
  sorted_sp <- sorted_data
  
  #Set colour palette - get more distinct palette
  numColors <- length(unique(sorted_sp$species)) # How many colors you need
  divider <- numColors/7
  numColors_unique <- numColors / ceiling(divider) #get fewer colors - limited by number available
  getColors <- brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package) (or try div or seq)
  myPalette <- getColors(round(numColors_unique)) #gets limited number of colours
  
  myPalette_big <- c(rep(myPalette, ceiling(divider)), myPalette[1:(numColors-(trunc(divider)*length(myPalette)))]) #repeats colours to get larger number of colours
  names(myPalette_big) <- unique(sorted_sp$species) #gets names associated with colours
  
  myPalette_big <- myPalette_big[which(is.na(names(myPalette_big)) == FALSE)] #remove colours with NAs as species from dataframe
  
  #make dataframe to refer to with colour and sample for yaxis
  yaxis_colours <- sorted_sp %>% distinct(Var1, species)
  yaxis_colours$colour <- c()
  for (i in 1:length(myPalette_big)){
    colour <- myPalette_big[i][1]
    species <- names(myPalette_big[i])
    for (j in 1:nrow(yaxis_colours)){
      if (as.character(yaxis_colours$species[j]) == species) {
        yaxis_colours$colour[j] <- colour
      } 
    }
  }
  
  getColors <- brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package) (or try div or seq)
  xaxis_colours <- getColors(5) #gets limited number of colours
  myPalette_xaxis <- c(rep(xaxis_colours, 30))#repeats colours to get larger number of colours
  names(myPalette_xaxis) <- unique(sorted_sp$source[order(sorted_sp$Var2, sorted_sp$source)]) #gets names associated with colours
  myPalette_xaxis <- myPalette_xaxis[which(is.na(names(myPalette_xaxis)) == FALSE)]
  
  xaxis_colours_df <- sorted_sp %>% distinct(Var2, source)
  xaxis_colours_df <- xaxis_colours_df[order(xaxis_colours_df$source, xaxis_colours_df$Var2),]
  
  xaxis_colours_df$colour <- ""
  
  for (i in 1:length(myPalette_xaxis)){
    colour <- myPalette_xaxis[i][1]
    source <- names(myPalette_xaxis[i])
    for (j in 1:nrow(xaxis_colours_df)){
      if (as.character(xaxis_colours_df$source[j]) == source) {
        xaxis_colours_df$colour[j] <- colour
      } 
    }
  }
  
  #make species factors for sorting
  sorted_sp$species <- as.factor(sorted_sp$species)
  
  #make sample names factors for plotting
  sorted_sp$Var1 <- factor(sorted_sp$Var1, levels = unique(sorted_sp$Var1))
  sorted_sp$source <- factor(sorted_sp$source, levels = unique(sorted_sp$source))
  
  #change order of probes
  sorted_sp <- sorted_sp[order(as.character(sorted_sp$source), sorted_sp$Var2),]
  sorted_sp$Var2 <- factor(sorted_sp$Var2, levels = unique(sorted_sp$Var2))
  
  #create categories of sequence length for colouring
  sorted_sp$category <- cut(sorted_sp$value, c(0, 1, 250, 500, 1000, 6000), include.lowest = TRUE)
  
  sorted_sp$label <- ""
  
  label_ones <- sorted_sp[!duplicated(sorted_sp$source),]
  
  for (i in 1:nrow(sorted_sp)){
    if (sorted_sp$Var2[i] %in% label_ones$Var2) {
      sorted_sp$label[i] <- as.character(sorted_sp$source[i])
    }
  }
  
  #save plot of percent length for all taxa and all loci
  pdf(paste0("Output_files/paper_scripts_output/",name,".pdf"), height = 4, width = 6)
  
  if (max(sorted_sp$value) <=1) {
    #plot seq percentage heatmap
    
    #Calculate sizes for gene and sample labels, increase the multiplier to make bigger
    gene_size_multiplier <- 0.0008
    sample_size_multiplier <- 0.005
    
    #get the axes to write scale
    gene.size <- dim(data)[2] * gene_size_multiplier
    sample.size <- dim(data)[1] * sample_size_multiplier
    
    #plot1 <- ggplot(data = sorted_sp, aes(x=Var2, y=reorder(Var1, desc(Var1)), fill = value))+
    plot1 <- ggplot(data = sorted_sp, aes(x=factor(Var2), y=factor(Var1), fill = value))+
      geom_raster()+
      #guides(fill=FALSE)+ #remove this line if you want the heat scale to appear
      scale_fill_gradient(high = "#132B43", low = "#FFFFFF")+
      theme(axis.text.x = element_text(angle = 45))+ #, hjust = 1
      ylab(NULL)+xlab(NULL)+
      scale_x_discrete(breaks= sorted_sp$Var2, labels = sorted_sp$label)+
      theme(axis.text.y=element_text(face="italic",size = sample.size, color=yaxis_colours$colour),axis.text.x = element_text(size=gene.size, color = xaxis_colours_df$colour)) #size=gene.size, 
    print(plot1)
  } else {
    #plot seq length heat map
    #get the scale for axes
    gene_size_multiplier <- 0.0008
    gene.size <- dim(data)[2] * gene_size_multiplier
    sample_size_multiplier <- 0.005
    sample.size <- dim(data)[1] * sample_size_multiplier
    
    plot2 <- ggplot(data = sorted_sp, aes(x=factor(Var2), y=factor(Var1)))+
      geom_raster(aes(fill = as.factor(category)))+
      theme(axis.text.x = element_text(size = gene.size, angle = 90, hjust = 1))+ #angle = 45, hjust = 1
      ylab(NULL)+xlab(NULL)+
      theme(axis.text.y=element_text(face="italic",size = sample.size, color=yaxis_colours$colour), legend.text = element_text(size = 3), legend.key.size = unit(0.5, "lines"), legend.position = "top",  legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
      #scale_x_discrete(breaks= sorted_sp$Var2, labels = sorted_sp$label)+
      scale_fill_viridis(name="", discrete=TRUE, labels = c("No Seq", "Short (1-250)", "Medium (250-500)", "Long (500-1000)", "Extra Long (1000-6000)"))
    print(plot2)
  }
  #scale_x_discrete(labels = sorted_sp$source)+size=gene.size, color = xaxis_colours_df$colour ,axis.text.x =element_blank()
  
  dev.off()
  
}



####
#Run functions

#Read output files from HybPiper

#First, read the files. These files are produced by the get_seq_lengths.py and hybpiper_stats.py scripts from HybPiper.

#Memecylon samples (assembly run at once)
seq_lengths_PA <- read.delim("Data/Revision_files/PA/01_All_loci_ntd_regular/hybpiper_seq_lengths_ntd.txt", row.names = 1, stringsAsFactors = FALSE)

#Tibouchina samples (assembly run at once)
seq_lengths_JJ <- read.delim("Data/All_JJ/hybpiper_seq_lengths_complete_mod.txt", row.names = 1, stringsAsFactors = FALSE)

#Remove ITS and ETS from PA set
seq_lengths_PA <- seq_lengths_PA[,-which(colnames(seq_lengths_PA) %in% c("KC523259_M_rivulare_ITS", "KC523102_M_rivulare_ETS"))]

#Make data into matrix format
seq_lengths_matrix_PA <- as.matrix(seq_lengths_PA) #seq_lengths_PA_2 if used list_remove
seq_lengths_matrix_JJ <- as.matrix(seq_lengths_JJ) #seq_lengths_JJ_2 if used list_remove

#Write matrices
write.csv(seq_lengths_matrix_JJ, "Output_files/paper_scripts_output/seq_lengths_matrix_JJ.csv", row.names = TRUE)
write.csv(seq_lengths_matrix_PA, "Output_files/paper_scripts_output/seq_lengths_matrix_PA.csv", row.names = TRUE)

##Run functions for datasets

#Run for Memecylon percents
all_percents_PA <- get_seq_percent(seq_lengths_matrix_PA)
sorted_all_p_PA <- melt_sort_PA(all_percents_PA)

#write and plot percents
write.csv(sorted_all_p_PA, "Output_files/paper_scripts_output/PA_all_percents.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_all_p_PA, "PA_all_percents")

#Run for Memecylon lengths
sorted_all_l_PA <- melt_sort_PA(seq_lengths_matrix_PA)

#write and plot lengths
write.csv(sorted_all_l_PA, "Output_files/paper_scripts_output/PA_all_lengths.csv")
plot_heatmap(seq_lengths_matrix_PA, sorted_all_l_PA, "PA_all_lengths")

#Run for Tibouchina percents
all_percents_JJ <- get_seq_percent(seq_lengths_matrix_JJ)
sorted_all_p_JJ <- melt_sort_JJ(all_percents_JJ)

#omit failed samples
sorted_all_p_JJ_2 <- sorted_all_p_JJ[which(sorted_all_p_JJ$Var1 != "A22_Inopinata_clean" & sorted_all_p_JJ$Var1 != "120_Rosanae_clean"),]

#write and plot percents
write.csv(sorted_all_p_JJ_2, "Output_files/paper_scripts_output/JJ_all_percents.csv")
plot_heatmap(seq_lengths_matrix_JJ, sorted_all_p_JJ_2, "JJ_all_percents")

#Run for Tibouchina lengths
sorted_all_l_JJ <- melt_sort_JJ(seq_lengths_matrix_JJ)

#omit failed samples
sorted_all_l_JJ_2 <- sorted_all_l_JJ[which(sorted_all_l_JJ$Var1 != "A22_Inopinata_clean" & sorted_all_l_JJ$Var1 != "120_Rosanae_clean"),]

#write and plot lenghts
write.csv(sorted_all_l_JJ_2, "Output_files/paper_scripts_output/JJ_all_lengths.csv")
plot_heatmap(seq_lengths_matrix_JJ, sorted_all_l_JJ_2, "JJ_all_lengths")

######
###Combine Memecylon and Tibouchina for plotting together

#get matrices into right format (transpose so that column is samples)
tJJ <- as.data.frame(t(as.data.frame(seq_lengths_matrix_JJ, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
tPA <- as.data.frame(t(as.data.frame(seq_lengths_matrix_PA, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)

#sort by rownames
tJJ_sort <- tJJ[order(rownames(tJJ)),]
tPA_sort <- tPA[order(rownames(tPA)),]

#combine sorted dataframes by columns (samples)
combo_dataframes <- cbind(tJJ_sort[which(rownames(tJJ_sort) %in% rownames(tPA_sort)),], tPA_sort[which(rownames(tPA_sort) %in% rownames(tJJ_sort)),])

#get into matrix
combo_matrices <- as.matrix(t(combo_dataframes))

##combine percent dataframes (sorted_sets_p_JJ and sorted_sets_p_PA)
combo_JJ_PA_p <- rbind(sorted_all_p_JJ_2[which(sorted_all_p_JJ_2$Var2 %in% sorted_all_p_PA$Var2),], sorted_all_p_PA[which(sorted_all_p_PA$Var2 %in% sorted_all_p_JJ_2$Var2),])

plot_heatmap(combo_matrices, combo_JJ_PA_p, "combo_PA_JJ_percent")

#combine length dataframes (sorted_sets_l_JJ and sorted_sets_l_JJ)
combo_JJ_PA_l <- rbind(sorted_all_l_JJ_2[which(sorted_all_l_JJ_2$Var2 %in% sorted_all_l_PA$Var2),], sorted_all_l_PA[which(sorted_all_l_PA$Var2 %in% sorted_all_l_JJ_2$Var2),])

plot_heatmap(combo_matrices, combo_JJ_PA_l, "combo_PA_JJ_length")

#write files
write.csv(combo_JJ_PA_l, "Output_files/paper_scripts_output/JJ_PA_combo_lengths.csv", row.names = TRUE)
write.csv(combo_JJ_PA_p, "Output_files/paper_scripts_output/JJ_PA_combo_percents.csv", row.names = TRUE)
write.csv(combo_matrices, "Output_files/paper_scripts_output/JJ_PA_combo_matrices.csv", row.names = TRUE)

########

#plot the combo lengths sorted correctly
combo_JJ_PA_l$split <- NA
combo_JJ_PA_l$Var2 <- as.character(combo_JJ_PA_l$Var2)

for (i in 1:nrow(combo_JJ_PA_l)){
  combo_JJ_PA_l$split[i] <- str_split(combo_JJ_PA_l$Var2[i], "_")
  combo_JJ_PA_l$probe[i] <- combo_JJ_PA_l$split[[i]][1]
}

combo_JJ_PA_l$group <- NA

combo_JJ_PA_l_2 <- combo_JJ_PA_l[order(combo_JJ_PA_l$probe),]


for (i in 1:nrow(combo_JJ_PA_l_2)){
  if (startsWith(combo_JJ_PA_l_2$probe[i], "AJ") | startsWith(combo_JJ_PA_l_2$probe[i], "NM")){
    combo_JJ_PA_l_2$group[i] <- "FN"
  } else if (startsWith(combo_JJ_PA_l_2$probe[i], "AT" ) | startsWith(combo_JJ_PA_l_2$probe[i], "TC")) {
    combo_JJ_PA_l_2$group[i] <- "MM"
  } else if (startsWith(combo_JJ_PA_l_2$probe[i], "KT")) {
    combo_JJ_PA_l_2$group[i] <- "SCN"
  } else {
    combo_JJ_PA_l_2$group[i] <- "Ang353"
  }
}

#remove split column
combo_JJ_PA_l_3 <- combo_JJ_PA_l_2[,-6]

#write file
write.csv(combo_JJ_PA_l_3, "Output_files/paper_scripts_output/JJ_PA_combo_lengths_groups.csv", row.names = TRUE)

#write function for grouping and axes

#get sources into right number of options
combo_JJ_PA_l_3$source <- as.character(combo_JJ_PA_l_3$source)

combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Memecyon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Memecylon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Memecylon_combo")] <- "Combo"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Torricelli_transcriptome_combo")] <- "Combo"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Memecyon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Torricelli")] <- "Memecylon"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Affzelli")] <- "Memecylon"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Tibouchina")] <- "Melastomateae"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Brachyotum")] <- "Melastomateae"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source == "Tetrazygia")] <- "Transcriptome"
combo_JJ_PA_l_3$source[which(combo_JJ_PA_l_3$source != "Melastomateae" & combo_JJ_PA_l_3$source != "Memecylon" & combo_JJ_PA_l_3$source != "Combo" & combo_JJ_PA_l_3$source != "Miconia" & combo_JJ_PA_l_3$source != "Transcriptome")] <- "Rosid"

#remove Tibouchina sample from Memecylon runs
#combo_JJ_PA_l_3 <- combo_JJ_PA_l_3[-which(combo_JJ_PA_l_3$Var1 == "Tibouchina_clean"),]


#plot the combo percents sorted correctly
combo_JJ_PA_p$split <- NA
combo_JJ_PA_p$Var2 <- as.character(combo_JJ_PA_p$Var2)

for (i in 1:nrow(combo_JJ_PA_p)){
  combo_JJ_PA_p$split[i] <- str_split(combo_JJ_PA_p$Var2[i], "_")
  combo_JJ_PA_p$probe[i] <- combo_JJ_PA_p$split[[i]][1]
}

combo_JJ_PA_p$group <- NA

combo_JJ_PA_p_2 <- combo_JJ_PA_p[order(combo_JJ_PA_p$probe),]


for (i in 1:nrow(combo_JJ_PA_p_2)){
  if (startsWith(combo_JJ_PA_p_2$probe[i], "AJ") | startsWith(combo_JJ_PA_p_2$probe[i], "NM")){
    combo_JJ_PA_p_2$group[i] <- "FN"
  } else if (startsWith(combo_JJ_PA_p_2$probe[i], "AT" ) | startsWith(combo_JJ_PA_p_2$probe[i], "TC")) {
    combo_JJ_PA_p_2$group[i] <- "MM"
  } else if (startsWith(combo_JJ_PA_p_2$probe[i], "KT")) {
    combo_JJ_PA_p_2$group[i] <- "SCN"
  } else {
    combo_JJ_PA_p_2$group[i] <- "Ang353"
  }
}

#remove split column
combo_JJ_PA_p_3 <- combo_JJ_PA_p_2[,-6]

#write file
write.csv(combo_JJ_PA_p_3, "Output_files/paper_scripts_output/JJ_PA_combo_percents_groups.csv", row.names = TRUE)

#write function for grouping and axes

#get sources into right number of options
combo_JJ_PA_p_3$source <- as.character(combo_JJ_PA_p_3$source)

combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Memecyon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Memecylon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Memecylon_combo")] <- "Combo"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Torricelli_transcriptome_combo")] <- "Combo"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Memecyon_transcriptome_combo")] <- "Combo"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Torricelli")] <- "Memecylon"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Affzelli")] <- "Memecylon"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Tibouchina")] <- "Melastomateae"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Brachyotum")] <- "Melastomateae"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source == "Tetrazygia")] <- "Transcriptome"
combo_JJ_PA_p_3$source[which(combo_JJ_PA_p_3$source != "Melastomateae" & combo_JJ_PA_p_3$source != "Memecylon" & combo_JJ_PA_p_3$source != "Combo" & combo_JJ_PA_p_3$source != "Miconia" & combo_JJ_PA_p_3$source != "Transcriptome")] <- "Rosid"

#remove Tibouchina sample from Memecylon runs
#combo_JJ_PA_p_3 <- combo_JJ_PA_p_3[-which(combo_JJ_PA_p_3$Var1 == "Tibouchina_clean"),]

#####
data <- combo_matrices
sorted_data <- combo_JJ_PA_l_3
name <- "combo_PA_JJ_length_grouped"

plot_heatmap_grouped <- function(data, sorted_data, name){
  
  #get data frame
  sorted_sp <- sorted_data
  
  unique(sorted_sp$source)
  
  #Set colour palette - get more distinct palette
  numColors <- length(unique(sorted_sp$species)) # How many colors you need
  divider <- numColors/7
  numColors_unique <- numColors / ceiling(divider) #get fewer colors - limited by number available
  getColors <- brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package) (or try div or seq)
  myPalette <- getColors(round(numColors_unique)) #gets limited number of colours
  
  myPalette_big <- c(rep(myPalette, ceiling(divider)), myPalette[1:(numColors-(trunc(divider)*length(myPalette)))]) #repeats colours to get larger number of colours
  names(myPalette_big) <- unique(sorted_sp$species) #gets names associated with colours
  
  
  myPalette_big <- myPalette_big[which(is.na(names(myPalette_big)) == FALSE)]
  #make dataframe to refer to with colour and sample
  yaxis_colours <- sorted_sp %>% distinct(Var1, species)
  yaxis_colours$colour <- c()
  for (i in 1:length(myPalette_big)){
    colour <- myPalette_big[i][1]
    species <- names(myPalette_big[i])
    for (j in 1:nrow(yaxis_colours)){
      if (as.character(yaxis_colours$species[j]) == species) {
        yaxis_colours$colour[j] <- colour
      } 
    }
  }
  
  getColors <- brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package) (or try div or seq)
  xaxis_colours <- getColors(5) #gets limited number of colours
  myPalette_xaxis <- c(rep(xaxis_colours, 30))#repeats colours to get larger number of colours
  names(myPalette_xaxis) <- unique(sorted_sp$source[order(sorted_sp$Var2, sorted_sp$source)]) #gets names associated with colours
  myPalette_xaxis <- myPalette_xaxis[which(is.na(names(myPalette_xaxis)) == FALSE)]
  
  xaxis_colours_df <- sorted_sp %>% distinct(Var2, source)
  
  xaxis_colours_df <- xaxis_colours_df[order(xaxis_colours_df$source, xaxis_colours_df$Var2),]
  
  xaxis_colours_df$colour <- ""
  
  for (i in 1:length(myPalette_xaxis)){
    colour <- myPalette_xaxis[i][1]
    source <- names(myPalette_xaxis[i])
    for (j in 1:nrow(xaxis_colours_df)){
      if (as.character(xaxis_colours_df$source[j]) == source) {
        xaxis_colours_df$colour[j] <- colour
      } 
    }
  }
  
  #make species and sample names factors for sorting and plotting
  sorted_sp$species <- as.factor(sorted_sp$species)
  sorted_sp$Var1 <- factor(sorted_sp$Var1, levels = unique(sorted_sp$Var1))
  sorted_sp$source <- factor(sorted_sp$source, levels = c("Melastomateae", "Memecylon", "Combo", "Transcriptome", "Miconia", "Rosid"))
  sorted_sp$group <- factor(sorted_sp$group, levels = c("MM", "SCN", "FN", "Ang353"))
  
  unique(sorted_sp$source)
  
  #change order for x axis
  sorted_sp <- sorted_sp[order(sorted_sp$group, sorted_sp$source, sorted_sp$Var2),]
  sorted_sp$Var2 <- factor(sorted_sp$Var2, levels = unique(sorted_sp$Var2))
 
  
  unique(sorted_sp$group)
  #create categories of sequence length for colouring
  sorted_sp$category <- cut(sorted_sp$value, c(0, 1, 250, 500, 1000, 6000), include.lowest = TRUE)
  
  sorted_sp$label <- ""
  
  label_ones <- sorted_sp[!duplicated(sorted_sp$source),]
  
  for (i in 1:nrow(sorted_sp)){
    if (sorted_sp$Var2[i] %in% label_ones$Var2) {
      sorted_sp$label[i] <- as.character(sorted_sp$source[i])
    }
  }
 
 
  #save plot of percent length for all taxa and all loci
  pdf(paste0("Output_files/paper_scripts_output/",name,".pdf"), height = 4, width = 6)
  
  if (max(sorted_sp$value) <=1) {
    #plot seq percentage heatmap
    
    #Calculate sizes for gene and sample labels, increase the multiplier to make bigger
    gene_size_multiplier <- 0.0008
    sample_size_multiplier <- 0.01
    
    #get the axes to write scale
    gene.size <- dim(data)[2] * gene_size_multiplier
    sample.size <- dim(data)[1] * sample_size_multiplier
    
    #plot1 <- ggplot(data = sorted_sp, aes(x=Var2, y=reorder(Var1, desc(Var1)), fill = value))+
    plot1 <- ggplot(data = sorted_sp, aes(x=factor(Var2), y=factor(Var1), fill = value))+
      geom_raster()+
      theme(legend.position = "right",  axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_rect(color = "black", size = 0.5, fill = NA))+
      #guides(fill=FALSE)+ #remove this line if you want the heat scale to appear
      scale_fill_gradient(high = "#132B43", low = "#FFFFFF", guide = guide_colorbar(title = "Proportion", ticks = FALSE, reverse = TRUE))+
     # theme(axis.text.x = element_text(angle = 45))+ #, hjust = 1
      ylab(NULL)+xlab(NULL)
      #labs(fill = "Proportion")+
      #scale_x_discrete(breaks= sorted_sp$Var2, labels = sorted_sp$label)+
      #theme(axis.text.y=element_text(face="italic",size = sample.size),axis.text.x = element_text(size=gene.size)) #size=gene.size, color = xaxis_colours_df$colour, color=yaxis_colours$colour
    print(plot1)
    return(plot1)
  } else {
    #plot seq length heat map
    #get the scale for axes
    gene_size_multiplier <- 0.0008
    gene.size <- dim(data)[2] * gene_size_multiplier
    sample_size_multiplier <- 0.01
    sample.size <- nrow(data) * sample_size_multiplier
    
    
    plot2 <- ggplot(data = sorted_sp, aes(x=factor(Var2), y=factor(Var1)))+
      geom_raster(aes(fill = as.factor(category)))+
      #theme(axis.text.x = element_text(size = gene.size, angle = 90, hjust = 1))+ #angle = 45, hjust = 1
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_rect(color = "black", size = 0.5, fill = NA))+
      ylab(NULL)+xlab(NULL)+
      theme(legend.key.size = unit(0.5, "lines"), legend.position = "right",  legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+ #, color=yaxis_colours$colour, legend.text = element_text(size = 6), axis.text.y=element_text(size = sample.size), 
      #scale_x_discrete(breaks= sorted_sp$Var2, labels = sorted_sp$label)+
      scale_fill_viridis(name="Length (bp)", discrete=TRUE, labels = c("0", "1-250", "251-500", "501-1000", "1001-6000"))
    print(plot2)
  }
  #scale_x_discrete(labels = sorted_sp$source)+size=gene.size, color = xaxis_colours_df$colour ,axis.text.x =element_blank()
 
  dev.off()
  
  return(plot2)
}



#plot grouped and ordered heatmaps
plot2_saved <- plot_heatmap_grouped(combo_matrices, combo_JJ_PA_l_3, "combo_PA_JJ_length_grouped")
plot1_saved <- plot_heatmap_grouped(combo_matrices, combo_JJ_PA_p_3, "combo_PA_JJ_percents_grouped")
dev.off()
grid.arrange(plot2_saved, plot1_saved)

#plot combo of two plots together
heatmap_grouped_prep <- function(data, sorted_data){
  
  #get data frame
  sorted_sp <- sorted_data
  
  #make species and sample names factors for sorting and plotting
  sorted_sp$species <- as.factor(sorted_sp$species)
  sorted_sp$Var1 <- factor(sorted_sp$Var1, levels = unique(sorted_sp$Var1))
  sorted_sp$source <- factor(sorted_sp$source, levels = c("Melastomateae", "Memecylon", "Combo", "Transcriptome", "Miconia", "Rosid"))
  sorted_sp$group <- factor(sorted_sp$group, levels = c("MM", "SCN", "FN", "Ang353"))
  
  #change order for x axis
  sorted_sp <- sorted_sp[order(sorted_sp$group, sorted_sp$source, sorted_sp$Var2),]
  sorted_sp$Var2 <- factor(sorted_sp$Var2, levels = unique(sorted_sp$Var2))
  
  #create categories of sequence length for colouring
  sorted_sp$category <- cut(sorted_sp$value, c(0, 1, 250, 500, 1000, 6000), include.lowest = TRUE)
  
  sorted_sp$label <- ""
  
  label_ones <- sorted_sp[!duplicated(sorted_sp$source),]
  
  for (i in 1:nrow(sorted_sp)){
    if (sorted_sp$Var2[i] %in% label_ones$Var2) {
      sorted_sp$label[i] <- as.character(sorted_sp$source[i])
    }
  }
  
  return(sorted_sp)
 
}

lengths_prepped <- heatmap_grouped_prep(combo_matrices, combo_JJ_PA_l_3)
percents_prepped <- heatmap_grouped_prep(combo_matrices, combo_JJ_PA_p_3)

#save plot of percent length for all taxa and all loci


plot1 <- ggplot(data = percents_prepped, aes(x=factor(Var2), y=factor(Var1)))+
  geom_raster(aes(fill = value))+
  coord_fixed(ratio=2)+
  theme(legend.title = element_text(size = 10), legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), legend.text = element_text(size = 8), legend.key.size = unit(0.5, "lines"),legend.position = "right",  axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_blank())+#panel.border = element_rect(color = "black", size = 0.5, fill = NA)
  scale_fill_gradient(high = "#132B43", low = "#FFFFFF", guide = guide_colorbar(title = "Proportion", ticks = FALSE, reverse = FALSE))+
  ylab(NULL)+xlab(NULL)+
  #y-axis labels
  annotate("segment", x = -15, xend = -15, y = 0, yend = 142, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -15, y = 143, yend = 236, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 143, yend = 143, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 236, yend = 236, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 0, yend = 0, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 142, yend = 142, colour = "black", size = 0.5)+
  annotate("text", x = -30, y = 65, colour = "black", label = "italic(Tibouchina)", angle = 90, size = 2, parse = TRUE)+
  annotate("text", x = -30, y = 190, colour = "black", label = "italic(Memecylon)", angle = 90, size = 2, parse = TRUE)+
  annotate("segment", x = -60, xend = -60, y = 0, yend = 236, colour = "white")+
  annotate("text", x = -45, y = 120, colour = "black", label = "Samples", angle = 90, size = 6)+
  #x-axis labels
  annotate("segment", x = 0, xend = 215, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 216, xend = 238, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 239, xend = 250, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 689, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 0, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 215, xend = 215, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 216, xend = 216, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 238, xend = 238, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 239, xend = 239, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 250, xend = 250, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 251, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 689, xend = 689, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 104, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 105, xend = 114, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 115, xend = 194, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 195, xend = 215, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 371, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 372, xend = 420, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 421, xend = 689, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 0, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 104, xend = 104, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 105, xend = 105, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 114, xend = 114, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 115, xend = 115, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 194, xend = 194, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 195, xend = 195, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 215, xend = 215, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 251, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 371, xend = 371, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 372, xend = 372, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 420, xend = 420, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 421, xend = 421, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 689, xend = 689, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("text", x = 100, y = -20, colour = "black", label = "MarkerMiner", size = 2)+
  annotate("text", x = 225, y = -20, colour = "black", label = "SCN", size = 2)+
  annotate("text", x = 246, y = -20, colour = "black", label = "Fn", size = 2)+
  annotate("text", x = 450, y = -20, colour = "black", label = "Angiosperm 353", size = 2)+
  annotate("text", x = 50, y = -9, colour = "black", label = "italic(Tibouchina)", size = 2, parse = TRUE)+
  annotate("text", x = 110, y = -9, colour = "black", label = "italic(M)", size = 2, parse = TRUE)+
  annotate("text", x = 155, y = -9, colour = "black", label = "Hybrid M+Tr", size = 2)+
  annotate("text", x = 205, y = -9, colour = "black", label = "Tr", size = 2)+
  annotate("text", x = 300, y = -9, colour = "black", label = "italic(Tibouchina)", size = 2, parse = TRUE)+
  annotate("text", x = 395, y = -9, colour = "black", label = "italic(Memecylon)", size = 2, parse = TRUE)+
  annotate("text", x = 550, y = -9, colour = "black", label = "Rosids", size = 2)+
  annotate("segment", x = 0, xend = 689, y = -40, yend = -40, colour = "white")+
  annotate("text", x = 320, y = -30, colour = "black", label = "Loci", size = 6)+
  #annotate("text", x = 770, y = -15, colour = "black", label = "Method of selection", size = 4)+
  #annotate("text", x = 770, y = -5, colour = "black", label = "Genome source for probe", size = 4)+
  #annotate("segment", x = 820, xend = 820, y = 0, yend = -45, colour = "white")+
  #add fig letter  
  annotate("text", x = -40, y = 225, colour = "black", label = "B", size = 8)

print(plot1)

#plot seq length heat map

plot2 <- ggplot(data = lengths_prepped, aes(x=factor(Var2), y=factor(Var1)))+
  geom_raster(aes(fill = as.factor(category)))+
  coord_fixed(ratio=2)+
  theme(legend.title = element_text(size = 10), legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), legend.text = element_text(size = 8), legend.key.size = unit(0.5, "lines"),legend.position = "right",  axis.text = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())+#panel.border = element_rect(color = "black", size = 0.5, fill = NA)
  ylab(NULL)+xlab(NULL)+
  scale_fill_viridis(name="Length (bp)", discrete=TRUE, labels = c("0", "1-250", "251-500", "501-1000", "1001-6000"))+
  #y-axis labels
  annotate("segment", x = -15, xend = -15, y = 0, yend = 142, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -15, y = 143, yend = 236, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 143, yend = 143, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 236, yend = 236, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 0, yend = 0, colour = "black", size = 0.5)+
  annotate("segment", x = -15, xend = -5, y = 142, yend = 142, colour = "black", size = 0.5)+
  annotate("text", x = -30, y = 65, colour = "black", label = "italic(Tibouchina)", angle = 90, size = 2, parse = TRUE)+
  annotate("text", x = -30, y = 190, colour = "black", label = "italic(Memecylon)", angle = 90, size = 2, parse = TRUE)+
  annotate("segment", x = -60, xend = -60, y = 0, yend = 236, colour = "white")+
  annotate("text", x = -45, y = 120, colour = "black", label = "Samples", angle = 90, size = 6)+
  #x-axis labels
  annotate("segment", x = 0, xend = 215, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 216, xend = 238, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 239, xend = 250, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 689, y = -15, yend = -15, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 0, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 215, xend = 215, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 216, xend = 216, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 238, xend = 238, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 239, xend = 239, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 250, xend = 250, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 251, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 689, xend = 689, y = -15, yend = -12, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 104, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 105, xend = 114, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 115, xend = 194, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 195, xend = 215, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 371, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 372, xend = 420, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 421, xend = 689, y = -5, yend = -5, colour = "black", size = 0.5)+
  annotate("segment", x = 0, xend = 0, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 104, xend = 104, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 105, xend = 105, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 114, xend = 114, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 115, xend = 115, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 194, xend = 194, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 195, xend = 195, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 215, xend = 215, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 251, xend = 251, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 371, xend = 371, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 372, xend = 372, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 420, xend = 420, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 421, xend = 421, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("segment", x = 689, xend = 689, y = -5, yend = -2, colour = "black", size = 0.5)+
  annotate("text", x = 100, y = -20, colour = "black", label = "MarkerMiner", size = 2)+
  annotate("text", x = 225, y = -20, colour = "black", label = "SCN", size = 2)+
  annotate("text", x = 246, y = -20, colour = "black", label = "Fn", size = 2)+
  annotate("text", x = 450, y = -20, colour = "black", label = "Angiosperm 353", size = 2)+
  annotate("text", x = 50, y = -9, colour = "black", label = "italic(Tibouchina)", size = 2, parse = TRUE)+
  annotate("text", x = 110, y = -9, colour = "black", label = "italic(M)", size = 2, parse = TRUE)+
  annotate("text", x = 155, y = -9, colour = "black", label = "Hybrid M+Tr", size = 2)+
  annotate("text", x = 205, y = -9, colour = "black", label = "Tr", size = 2)+
  annotate("text", x = 300, y = -9, colour = "black", label = "italic(Tibouchina)", size = 2, parse = TRUE)+
  annotate("text", x = 395, y = -9, colour = "black", label = "italic(Memecylon)", size = 2, parse = TRUE)+
  annotate("text", x = 550, y = -9, colour = "black", label = "Rosids", size = 2)+
  annotate("segment", x = 0, xend = 689, y = -40, yend = -40, colour = "white")+
  annotate("text", x = 320, y = -30, colour = "black", label = "Loci", size = 6)+
  #annotate("text", x = 770, y = -15, colour = "black", label = "Method of selection", size = 4)+
  #annotate("text", x = 770, y = -5, colour = "black", label = "Genome source for probe", size = 4)+
  #annotate("segment", x = 820, xend = 820, y = 0, yend = -45, colour = "white")+
  #add fig letter  
  annotate("text", x = -40, y = 225, colour = "black", label = "A", size = 8)
  


print(plot2)

dev.off()

pdf("Output_files/paper_scripts_output/revision_combo_plots.pdf", height = 8, width = 7.5)
ggarrange(plot2, plot1, newpage = FALSE)
dev.off()


#get dimensions of plot
ggp <- ggplot_build(plot1)
my.ggp.yrange <- ggp$layout$panel_scales_y[[1]]$range$range  # data range!
my.ggp.xrange <- ggp$layout$panel_scales_x[[1]]$range$range  # data range!




#Get mean (reference) length
ref <- seq_lengths_PA[1,]

write.csv(ref, "./Output_files/paper_scripts_output/Revision/ref_lengths.csv")
