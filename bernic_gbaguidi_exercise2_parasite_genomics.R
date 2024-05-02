###########################
# MaModAfrica
# Genomics and Bioinformatics course
# Parasite Genomics
#
# Practical exercise 2: DNA sequence alignment and classification of recrudescences and new infections
###########################

# install the necessary package for alignment (done only once)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("QuasR")

# load the necessary packages
library(QuasR) # for alignment
library(tidyr) # for manipulating the data table
library(dplyr) # for manipulating the data table
library(ggplot2) # for plotting
library(stringr) # for string manipulation

# Set the directory to the folder for Day 2 that you downloaded on your computer
setwd("C:/Users/Administrateur/Downloads/Day2- Lecture Monica/Day2")
# File setup and specification:
# Specify the reference genome file which is in fasta format (.fa, .fna).
# This reference will be used for alignment.
genomeFile = "reference_HB37.fna"

# Load the text file containing definitions of sample names and fastq file paths.
# This text file will be used for localizing the sequencing reads.
sampleFile = "samples_info.txt"

# Run the sequence alignment. This operation will take a few seconds.
proj = qAlign(sampleFile, genomeFile)

# Checking how many reads were aligned from each sample
alignmentStats(proj) 

# Generating a report with several diagnostics of the results
qQCReport(proj, "qc_report.pdf")

##############
# PART 2: Classification of infections and calculation of treatment efficacy
##############

# Loading the haplotypes identified for each marker
all_haplotypes = read.csv("haplotypes_in_samples.csv")

# Printing the first lines to see how the table looks like
head(all_haplotypes, 5)

# Observe the naming convention for the samples. Identify the unique patient IDs
sample_pairs = unique(all_haplotypes$SampleName)

# Extract only the patient IDs (remove the "_D0", "_DX")
# There are several ways to do this:
# You can remove them manually by locating the position of the "_"
samples = unique(substr(sample_pairs, 1, str_locate(sample_pairs, "_") - 1))
print(samples)

# Or you can directly replace them with an empty character
samples = str_remove_all(sample_pairs, "_D0")
samples = str_remove_all(samples, "_DX")
samples = unique(samples)
print(samples)


#### Identify for each marker whether we observe a recrudescence or new infection

# Initialise the results table
result_tab = NULL

# Extract the list with the three marker names
marker_list = unique(all_haplotypes$MarkerID)

# Loop over all markers
for (marker in marker_list) {
  # Select only the rows corresponding to the selected marker
  hap_m = all_haplotypes %>% filter(MarkerID==marker)
  
  # Loop through all the samples
  for (s in samples) {
    # Reconstruct the Day 0 and Day X names 
    s_D0 = paste0(s, "_D0")
    s_DX = paste0(s, "_DX")
    # Extract the haplotypes at day 0 and Day X
    hap_D0 = hap_m %>% filter(SampleName == s_D0)
    hap_DX = hap_m %>% filter(SampleName == s_DX)
    
    # Print the intersection
    common_haplotypes = intersect(hap_D0$Haplotype, hap_DX$Haplotype)
    print(common_haplotypes)
    
    # Check if the sets of haplotypes intersect (then recrudescence)
    if(length(common_haplotypes ) > 0) {
      # Recrudescence!
      infection_type = "R"
    } else {
      infection_type = "NI"
    }
    # Append the new result for the considered sample pair and marker
    line = data.frame(SampleName = s, Marker = marker, InfectionType = infection_type)
    result_tab = rbind(result_tab, line)
  }
}

print(result_tab)

# Count how many times you observe recrudescence across the three markers for each sample pair
n_r = result_tab %>% group_by(SampleName) %>% summarise(numberRecrudescence=sum(InfectionType=='R'))
print(n_r)

n_r_1 = result_tab %>% group_by(SampleName) %>% summarise(numberRecrudescence=sum(InfectionType=='R'),Result_3_out_3=ifelse(sum(InfectionType=='R')==3,"R","NI"),Result_3_out_2=ifelse(sum(InfectionType=='R')>=2,"R","NI"))
#n_r_3x2 = result_tab %>% group_by(SampleName) %>% summarise(numberRecrudescence=sum(InfectionType=='R'),ifelse(sum(InfectionType=='R')>=2,"R","NI"))




##### Treatment efficacy

Trteff <- n_r_1 %>% 
  as.data.frame() %>% 
  filter(Result_3_out_2=='NI'| Result_3_out_3=="NI") %>% 
  summarise(TE3X3=sum(Result_3_out_3=='NI')/20,TE3X2=sum(Result_3_out_2=='NI')/20)
