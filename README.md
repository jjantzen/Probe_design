# Probe_design
Probe design for target enrichment of Memecylon and Tibouchina (Melastomataceae)

MarkerMiner was used to identify loci for target enrichment (https://bitbucket.org/srikarchamala/markerminer/src/master/).

Prior to assembly, reads were cleaned and trimmed using the secapr_cleaning.txt script (see also https://github.com/AntonelliLab/seqcap_processor).

Assembly was conducted using HybPiper (https://github.com/mossmatters/HybPiper) and assembly stats were calculated using associated scripts.

Custom scripts:

To calculate percent identity of template sequences and recovered sequences, we conducted needle pairwise alignments using the needle_align.txt script. 

To calculate lengths of introns and supercontigs, the calc_length.py script was used by running the seq_calc_sub.txt script. 

Other scripts for manipulating fasta and text files include cat_files.txt and separating_script.txt.

The R scripts included in Stats_R_scripts were used to aggregate and present the statistics on enrichment and assembly success. 
