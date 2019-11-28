#read depth per locus and number of templates

library(tidyr)
library(tidyselect)
library(ggplot2)

JJ_depths <- read.delim("./Data/Revision_files/Read_depth/summary_depth.txt", stringsAsFactors = FALSE, header = FALSE)

JJ_depths <- JJ_depths[c(nrow(JJ_depths),1:(nrow(JJ_depths)-1)),]

rownames(JJ_depths) <- JJ_depths[,1]

JJ_depths_sorted <- JJ_depths[order(rownames(JJ_depths)),]

colnames(JJ_depths_sorted) <- JJ_depths_sorted[81,]

JJ_depths_sorted <- JJ_depths_sorted[-81,]

colnames(JJ_depths_sorted)[1] <- "Samples"

#reformat dataframe to get groups of loci by template

#remove ITS ETS


JJ_depths_sorted <- JJ_depths_sorted[,-c(ncol(JJ_depths_sorted)-1, ncol(JJ_depths_sorted))]

gathered_depths <- JJ_depths_sorted %>% 
  pivot_longer(-Samples, names_to = "Genes")

gathered_depths$split <- strsplit(gathered_depths$Genes, "_")

gathered_depths$locus <- NA

for (i in 1:nrow(gathered_depths)){
  gathered_depths$locus[i] <- gathered_depths$split[[i]][[1]]
}


gathered_depths_templates <- gathered_depths[,-4]

#get number of templates from sum of locus counts

template_data <- gathered_depths_templates %>% 
  group_by(locus, Samples) %>%
  mutate(num_templates = n())  
  

counting <- as.data.frame(template_data) %>% 
  dplyr::select(locus, num_templates) %>% 
  distinct() %>% 
  group_by(num_templates) %>% 
  summarize(count = n())

class(template_data$value) <- "integer"

#results - 1 = 156, 2 = 164, 3 = 55, 4 = 10
pdf("./Output_files/paper_scripts_output/Revision/read_depth_boxplot_log.pdf")
ggplot(template_data)+
  geom_boxplot(aes(x = num_templates, y = log(value), group = num_templates))+
  xlab("Number of template sequences per locus")+
  ylab("Read depth (log)")+
  theme_bw()
dev.off()


