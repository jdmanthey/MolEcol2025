## Estimating Genetic Diversity from Genomic Data - 18 September 2025

Today we will be working with some genomic data in order to estimate population genetic diversity within populations.
First, do the following items:
1. Download the lab4_items.zip file above.
2. Unzip the folder and place it on your desktop (or other working location).
3. Open up RStudio and set your working directory to the newly-unzipped directory.

One optional step (which I often prefer) is to remove scientific notation of large numbers. It is up to you if you want to do 
this, but I often find it helpful. This makes all the numbers always print out fully:

    options(scipen=999)

&nbsp;

### Study System

For today's exercises, there are two potential datasets to use, chromosomes 1 and 19 (the first million bp (Mbp)), of an organism
called the Brown Creeper (_Certhia americana_). To give you a better understanding of where these samples are coming from 
in a taxonomic viewpoint, below is a map showing the distribution of the Brown Creeper (_Certhia americana_). This species has two 
main lineages, denoted here by the blue and red colors:

![distribution](https://github.com/jdmanthey/MolEcol2019/blob/master/05_genetic_differentiation/distribution.png)

For some of the questions on the exercises, you will need to communicate and compare results with one or more other people that have 
analyzed the opposite dataset. Check to make sure that your potential partner(s) chooses the other dataset relative to you for comparison.

Once you have this figured out with your partner(s), read in your respective dataset with the function 'read_vcf.' 

&nbsp;

## 1. Load functions and read in VCF

    # load some functions to read a vcf and calculate diversity
    source("lab4_functions.R")
    
    # choose one or the other (don't run both):
    input_file <- "chr1_1Mbp.recode.vcf"
    
    input_file <- "chr19_1Mbp.recode.vcf"

    # read in vcf
    vcf <- read_vcf(input_file)  
    
    # simply the vcf
    vcf <- simplify_vcf(vcf)

Normally in R, when you just type the name of the object (here we named it simply "input_file"), R displays the whole object. With some large
files, you won't want to do that. Instead, you could look at the structure of an object to see what components it has. Check this out
with the following command:

    str(vcf)

## 2. Define populations 

Next, we will define the populations sampled in the dataset using a text file already in your directory. This file is tab delimited with
the first column giving the number of the sample and the second column giving the population. Read that file in and take a look at the file.

    # read in popmap
    popmap <- read.table("popmap.txt", header=T)

    # identify unique populations
    populations <- unique(popmap[,2])

There are four populations sampled in this dataset:
1. Santa Rita Mountains, Arizona
2. Morelos, Mexico
3. Pinal Mountains, Arizona
4. Utah, USA

&nbsp;

## 3. Estimate nucleotide diversity (pi) of entire region. 

Now that we have our populations defined, we can calculate nucleotide diversity within each population. Remember that the nucleotide 
diversity is the mean number of differences between individuals sampled. 

    # estimate nucleotide diversity for all populations
    pi <- pi_all_pops(vcf, popmap)
    
    # correct for window_size (approximated)
    pi$pi <- pi$pi * pi$n_sites / 1000000

    pi

&nbsp;

## 4. Estimate pi in sliding windows of 25 kbp and 100 kbp.

Now that we've looked at the mean patterns across this whole length of sequence we are analyzing, we will look to see if there
are different patterns in sliding windows across the region. Sliding windows in the case we are setting up are looking at
non-overlapping segments of the chromosome. Here, we will set up two different-sized sliding windows: 25 kbp and 100 kbp. 
These are arbitrary numbers, but are representative of possible window sizes we may investigate across an entire genome. 

    # define window size
    window_size <- 25000

    # estimate pi in sliding windows
    window_output_25000 <- pi_windows(vcf, popmap, window_size)

    # correct for window_size (approximated)
    window_output_25000$pi <- window_output_25000$pi * window_output_25000$n_sites / window_size

### Modify the above code section to make a new object called 'window_output_100000' that has pi estimates in 100 kbp windows.

&nbsp;

## 5. Plot the sliding window estimates.

    # Set up plotting dimensions with 2 rows and 1 column
    par(mfrow=c(2,1))
    # choose plotting colors (feel free to change)
    plot_colors <- c("darkblue", "darkorange2", "deeppink3", "black")
    # 1st color = Santa Rita Mountains, Arizona
    # 2nd color = Morelos, Mexico
    # 3rd color = Pinal Mountains, Arizona
    # 4th color = Utah, USA

    # Set up the plotting dimensions with axis and chart labels for 25kbp windows
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="Nucleotide Diversity", main="25 kbp windows")
    # plot for each population
    for(a in 1:length(populations)) {
      a_rep <- window_output_25000[window_output_25000$population == populations[a],]
      lines((a_rep$end - 12500), a_rep$pi, col=plot_colors[a])
    }

    # Set up the plotting dimensions with axis and chart labels for 100kbp windows
    plot(c(1,1), col="white", ylim=c(0,0.01), xlim=c(0,1000000), xlab="Position(bp)", ylab="Nucleotide Diversity", main="100 kbp windows")
    # plot for each population
    for(a in 1:length(populations)) {
      a_rep <- window_output_100000[window_output_100000$population == populations[a],]
      lines((a_rep$end - 50000), a_rep$pi, col=plot_colors[a])
    }

## 6. Estimate nucleotide diversity in a mixed population.

Make a copy of the popmap.txt file and rename it to: popmap2.txt. Edit the file so that there are 3 additional lines that define a new
population named "mixed" and includes one individual each from three populations. After you have done this and saved the file, redo steps 
2 and 3 above with this new popmap. 


