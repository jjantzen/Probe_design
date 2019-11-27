#read depth per locus and number of templates

library(tidyr)
library(tidyselect)
library(ggplot2)

JJ_depths <- read.delim("./Data/Revision_files/Read_depth/summary_depth.txt", stringsAsFactors = FALSE, header = FALSE)

JJ_depths <- JJ_depths[c(nrow(JJ_depths),1:(nrow(JJ_depths)-1)),]

rownames(JJ_depths) <- JJ_depths[,1]

JJ_depths_sorted <- JJ_depths[order(rownames(JJ_depths)),]

JJ_depths_sorted[,c(690:692)]
ncol(JJ_depths_sorted)

colnames(JJ_depths_sorted) <- JJ_depths_sorted[81,]

JJ_depths_sorted <- JJ_depths_sorted[-81,]

colnames(JJ_depths_sorted)[1] <- "Samples"

colnames(JJ_depths_sorted)
rownames(JJ_depths_sorted)

#reformat dataframe to get groups of loci by template

str(JJ_depths_sorted)

JJ_depths_sorted[c(1:10), c(1:10)]

#remove ITS ETS

ncol(JJ_depths_sorted)

JJ_depths_sorted <- JJ_depths_sorted[,-c(ncol(JJ_depths_sorted)-1, ncol(JJ_depths_sorted))]

gathered_depths <- JJ_depths_sorted %>% 
  pivot_longer(-Samples, names_to = "Genes")



gathered_depths$split <- strsplit(gathered_depths$Genes, "_")

gathered_depths$locus <- NA

for (i in 1:nrow(gathered_depths)){
  gathered_depths$locus[i] <- gathered_depths$split[[i]][[1]]
}


unique(gathered_depths$locus)

gathered_depths_templates <- gathered_depths[,-4]

gathered_depths_templates[c(1:10),]

#get number of templates from sum of locus counts

template_data <- gathered_depths_templates %>% 
  group_by(locus, Samples) %>%
  mutate(num_templates = n())  
  

counting <- as.data.frame(template_data) %>% 
  dplyr::select(locus, num_templates) %>% 
  distinct() %>% 
  group_by(num_templates) %>% 
  summarize(count = n())

str(template_data)

class(template_data$value) <- "integer"

#results - 1 = 156, 2 = 164, 3 = 55, 4 = 10
pdf("./Output_files/paper_scripts_output/Revision/read_depth_boxplot_log.pdf")
ggplot(template_data)+
  geom_boxplot(aes(x = num_templates, y = log(value), group = num_templates))+
  xlab("Number of template sequences per locus")+
  ylab("Read depth (log)")+
  theme_bw()
dev.off()

unique(template_data$value)

unique(template_data$Genes[which(template_data$num_templates == 4)])


#try plotting by probe source instead of num of templates
ggplot(template_data)+
  geom_boxplot(aes(x = num_templates, y = log(value), group = num_templates))+
  xlab("Number of template sequences per locus")+
  ylab("Read depth (log)")
