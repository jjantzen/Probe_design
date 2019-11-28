###Calculate mean read depth###

# ((TOTAL ON TARGET READ COUNT) * (READ LENGTH))/(CUMULATIVE LENGTH OF ALL TARGET GENES)

# Cumulative target length is the end-to-end length of the targets. 
# Multiply the mean sequence length by the number of sequences.
# CUMULATIVE LENGTH OF ALL TARGET GENES — for instance if you have 100 genes, 1000 bp on average long, the cumulative length is 100*1000 bp.


#read in file with read numbers

reads_data_JJ <- read.delim("Data/All_JJ/hybpiper_stats_complete.txt")
reads_data_PA <- read.delim("Data/Prabha/All/hybpiper_stats_set5.txt")

length_data_JJ <- read.delim("Data/All_JJ/hybpiper_seq_lengths_complete_mod.txt")
length_data_PA <- read.delim("Data/Prabha/All/hybpiper_seq_lengths_set5_mod.txt")

both_together_reads <- rbind(reads_data_JJ, reads_data_PA)
both_together_lengths <- rbind(length_data_JJ, length_data_PA[which(colnames(length_data_PA) %in% colnames(length_data_JJ))])


#get number of reads (total)
read_num <- data.frame(sample = both_together_reads$Name, num_reads = both_together_reads$ReadsMapped)

read_num_JJ <- data.frame(sample = reads_data_JJ$Name, num_reads = reads_data_JJ$ReadsMapped)
read_num_PA <- data.frame(sample = reads_data_PA$Name, num_reads = reads_data_PA$ReadsMapped)


read_length <- 150

#get cumulative target length

#length for target (from reference length)
mean_lengths_target <- data.frame(both_together_lengths[which(both_together_lengths$Species == "MeanLength"),])

#sum of all target lengths (only once)
c_length <- sum(mean_lengths_target[1,-1])


#calculate read depth

read_num$depth <- read_num$num_reads * read_length / c_length
read_num_JJ$depth <- read_num_JJ$num_reads * read_length / c_length
read_num_PA$depth <- read_num_PA$num_reads * read_length / c_length

average_read_depth <- mean(read_num$depth)
average_read_depth_JJ <- mean(read_num_JJ$depth)
average_read_depth_PA <- mean(read_num_PA$depth)


