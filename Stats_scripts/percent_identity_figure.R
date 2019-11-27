#Plot percent identity figures for manuscript


#plot percent identity by avg exon length

#read length and pct identity file
percent_ID <- read.csv("Output_files/Tables/For_formatting/combo_length_pct_avg_stats_by_probe.csv", header = TRUE, stringsAsFactors = FALSE)

red_percent_ID <- percent_ID[,c(2,5,11,8,14,16,22,28,29,40,41)]

#rename columns by samples
colnames(red_percent_ID)[which(startsWith(colnames(red_percent_ID), "Tibouchina"))] <- c("probe_avg_no_zero_length.Tibouchina", "probe_taxon_count.Tibouchina", "avg_ID_pct.Tibouchina")
colnames(red_percent_ID)[which(startsWith(colnames(red_percent_ID), "Memecylon"))] <- c("probe_avg_no_zero_length.Memecylon", "probe_taxon_count.Memecylon", "avg_ID_pct.Memecylon")

#write functions for splitting locus columns, and grouping by source of probe sequences
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

#run functions for %id data
red_percent_ID_columns <- get_columns(red_percent_ID)
red_percent_ID_grouped <- grouping_by_source(red_percent_ID_columns)

#reshape dataframe
reshaped_percent_ID <- red_percent_ID_grouped %>% 
  gather("clade", "measurements", -probe, -source, -locus, -genome) %>% 
  extract("clade", c("variable", "clade"), "(.*)\\.(.+)$") %>% 
  spread(variable, "measurements") 

#percent of taxa recovered column
reshaped_percent_ID$prop_taxa <- NA
reshaped_percent_ID$prop_taxa[which(reshaped_percent_ID$clade == "Tibouchina")] <- reshaped_percent_ID$probe_taxon_count[which(reshaped_percent_ID$clade == "Tibouchina")]/40*100
reshaped_percent_ID$prop_taxa[which(reshaped_percent_ID$clade == "Memecylon")] <- reshaped_percent_ID$probe_taxon_count[which(reshaped_percent_ID$clade == "Memecylon")]/53*100

#plot data
#plot exon length vs % identity
genome.labs <- c("Other Rosid\n Sequences", "Memecylon","Tibouchina", "Transcriptome", "Miconia")
names(genome.labs) <- unique(reshaped_percent_ID$genome)#[1:4]

str(reshaped_percent_ID$genome)

reshaped_percent_ID$genome <- factor(reshaped_percent_ID$genome, levels = c("Memecylon", "Tibouchina", "Transcriptome", "Other_Rosid", "Miconia"), labels = c("italic(Memecylon)", "italic(Tibouchina)", "Transcriptomes", "Other Rosids", "italic(Miconia)"))
unique(reshaped_percent_ID$genome)

levels(reshaped_percent_ID$genome)= c("Memecylon"=expression(italic("Memecylon")),
                      "Tibouchina"=expression(italic("Tibouchina")),
                      "Transcriptome"=expression("Transcriptome"), 
                      "Other_Rosid"=expression("Rosids"),
                      "Miconia"=expression(italic("Miconia")))



pdf("Output_files/paper_scripts_output/Revision/percent_ID_plots_submission.pdf", width = 7.5, height = 4)
p1 <- ggplot(reshaped_percent_ID[which(reshaped_percent_ID$exon.average_length <= 2000),])+
  geom_point(aes(x = avg_ID_pct, y = exon.average_length, colour = clade, shape = clade), size = 0.75)+
  facet_wrap( ~ genome, labeller=label_parsed, ncol = 5)+#scales = "free_x", labeller = labeller(genome = genome.labs), 
  ylab(label = "Mean total exon \nlength recovered (bp)")+
  xlab(label = "Mean percent identity (%)")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), legend.text = element_text(face = "italic"), axis.text = element_text(size = 3))+
  scale_color_manual(values = c('#601A4A', '#63ACBE'), name = "Sample clade")+
  scale_shape_manual(values = c(19, 15), name = "Sample clade")

p2 <- ggplot(reshaped_percent_ID[which(reshaped_percent_ID$exon.average_length <= 2000),], aes(x = avg_ID_pct, y = prop_taxa))+
  geom_point(aes(colour = clade, shape = clade), size = 0.75)+
  facet_wrap( ~ genome, labeller=label_parsed, ncol = 5)+#scales = "free_x", 
  ylab(label = "Percent of taxa sequenced \nrelative to total (%)")+
  xlab(label = "Mean percent identity (%)")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), legend.text = element_text(face = "italic"), axis.text = element_text(size = 3))+
  scale_color_manual(values = c('#601A4A', '#63ACBE'), name = "Sample clade")+
  scale_shape_manual(values = c(19, 15), name = "Sample clade")

ggarrange(p1, p2, newpage = FALSE)

dev.off()

dev.off()

pdf("output_figures/percent_ID_vs_probe_no_zero_length.pdf")
ggplot(reshaped_percent_ID[which(reshaped_percent_ID$genome != "Miconia" & reshaped_percent_ID$probe_avg_no_zero_length <= 3000),], aes(x = probe_avg_no_zero_length, y = avg_ID_pct))+
  geom_point(aes(colour = clade))+
  facet_wrap( ~ genome, scales = "free")
dev.off()

pdf("Output_figures/percent_ID_vs_taxon_count.pdf")
ggplot(reshaped_percent_ID, aes(x = probe_taxon_count, y = avg_ID_pct))+
  geom_point(aes(colour = clade))
dev.off()

pdf("Output_figures/percent_ID_vs_supercontig.pdf")
ggplot(reshaped_percent_ID[which(reshaped_percent_ID$genome != "Miconia" & reshaped_percent_ID$supercontig.average_length <= 2000),], aes(x = supercontig.average_length, y = avg_ID_pct))+
  geom_point(aes(colour = clade))+
  facet_wrap( ~ genome, scales = "free")
dev.off()
