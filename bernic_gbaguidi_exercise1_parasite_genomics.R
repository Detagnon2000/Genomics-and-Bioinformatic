###########################
# MaModAfrica
# Genomics and Bioinformatics course
# Parasite Genomics
#
# Practical exercise: panel design
#
# Visualisation of genetic diversity for the gene cpmp.
# In this exercise you will work with a processed file which contains information
# about genomic variations across the cpmp gene. You will need to follow the 
# instructions in the comments below and fill out the empty fields (marked with "_____")
# 
# Send your script today to monica.golumbeanu@swisstph.ch
# The subject of the email should be "Practical 1 parasite genomics"
###########################

# Loading necessary packages
library(tidyr)
library(dplyr) # for manipulating the data table
library(ggplot2) # for plotting

# Loading the mutation information from MalariaGEN
malaria_gen_mut_file = "cpmp_malariaGEN_mut.txt"
cpmp_mut = read.csv(malaria_gen_mut_file, header = TRUE, sep = "\t")

##?## How many rows does the data table with mutation information have? Fill out your response:
n_rows_cpmp_mut = nrow(cpmp_mut)
print(n_rows_cpmp_mut)

##?## What is the command to print the first 15 lines of the data table:
head(x = cpmp_mut, n = 15)

# The data table contains several alternative sequence variants per row. 
# To ease table processing, we need to make sure we have only one alternative sequence per row
# The function separate_rows can be used for this purpose.
# You can access the help page of this function with the following command:
?separate_longer_delim

##?## Based on the above instructions from the help page, you will need to complete the function call:
cpmp_mut_long = separate_longer_delim(cpmp_mut, c("Alternative", "Frequency"), delim = ",")

##?## Print the first 10 lines of the long dataframe to check the result.
# What do you observe regarding the format of the frequency?
head(cpmp_mut_long,n=10)

# Just making sure that the mutation frequency is numeric
cpmp_mut_long$Frequency = as.double(cpmp_mut_long$Frequency)

# As you may have observed above, the data table contains all types of sequence variation (SNPs, indels, insertions, ...). 
# You can check this by listing the unique entries:
print(sort(unique(cpmp_mut_long$Alternative)))

# We need to remove the non-relevant sequence variations in the Alternative column and keep only the SNPs. 
# Same applies to the Reference column.
# Think of what rule applies to the columns Reference and Alternative.
# You can then use the function filter to select the relevant rows of the data table
?filter

##?## Fill out the necessary information for running the filter command
cpmp_SNPs = cpmp_mut_long %>% filter(nchar(Alternative)==1) %>% filter(nchar(Reference)==1)

# Printing the first 20 lines of the data table to check that the filtering worked
head(cpmp_SNPs, 20)

# The data table contains now at each mutated position on the cpmp gene, the mutation frequency for each observed mutation in the malariaGEN samples.
# This corresponds to the information of how often this mutation is observed in the population of samples form malariaGEN.

##?##Plot the mutation frequencies (column Frequency) at each position (column Position). 
# What do you observe?
ggplot(cpmp_SNPs, aes(x = Position, y =Frequency) )+ 
  geom_point() + theme_minimal()

# Your targeted marker region cannot be longer than 300 nucleotides. 
# You will need to further narrow it down.

##?## Fill out the positions in the command below to define narrower regions. 
# You can try several times until you reach a region less than 300 nucleotides. Note down the coordinates of the final targeted narrow region.
cpmp_SNPs_w = cpmp_SNPs %>% filter(Position >=180000 & Position < 180300)

# Plotting the narrower region
ggplot(cpmp_SNPs_w, aes(x = Position, y = Frequency)) + 
  geom_point() + theme_minimal() + scale_x_continuous(n.breaks = 10)



### Bonus question

bonus <- cpmp_SNPs %>% group_by(Position) %>% summarise(He=1-sum(Frequency^2)-(1-sum(Frequency))^2)

bonus_inter <- merge(cpmp_SNPs,bonus, by='Position')
ggplot(bonus_inter, aes(x = Position, y = Frequency)) + 
  geom_point() +
  geom_line(aes(x=Position,y=He),col='red')+
  theme_minimal() + scale_x_continuous(n.breaks = 10)
