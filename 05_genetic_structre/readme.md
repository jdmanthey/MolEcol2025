## Estimating Genetic Differentiation from Genomic Data 2 October 2025

Today we will be working with the same genomic data as last week. Refer to the information about populations, etc. from last
week regarding the source of the data. First, do the following items:
1. Download the lab5_items.zip file above.
2. Unzip the folder and place it on your desktop (or other working location).
3. Open up RStudio and set your working directory to the newly-unzipped directory.

Remove scientific notation on numbers if desired:

    options(scipen=999)

&nbsp;

### Study System (reminder; same as last activity)

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

    # load some functions to read a vcf and calculate genetic structure 
    source("lab5_functions.R")
    
    # choose one or the other (don't run both):
    input_file <- "chr1_1Mbp.recode.vcf"
    
    input_file <- "chr19_1Mbp.recode.vcf"

    # read in vcf
    vcf <- read_vcf(input_file)  
    
    # simply the vcf
    vcf <- simplify_vcf(vcf)

Remind yourself of the structure and dimensions of the vcf

    str(vcf)
    dim(vcf)

## 2. Define populations (again, same as last time)

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

## 3. Estimate differentiation (F<sub>ST</sub>) and divergence (D<sub>XY</sub>) of entire region. 

Now that we have our populations defined, we can calculate measures of genetic structure between pairs of populations.

    # estimate fst for all populations
    fst <- differentiation_all_pops(vcf, popmap, "fst")
    fst

    # estimate dxy for all populations
    dxy <- differentiation_all_pops(vcf, popmap, "dxy")
    dxy


What does the relationship between F<sub>ST</sub> and D<sub>XY</sub> look like?

    # plot
    plot(fst$fst, dxy$dxy, xlab="FST", ylab="DXY")

&nbsp;

## 4. Estimate F<sub>ST</sub> and D<sub>XY</sub> in sliding windows of 25 kbp and 100 kbp.

Now that we've looked at the mean patterns across this whole length of sequence we are analyzing, we will look to see if there
are different patterns in sliding windows across the region. Sliding windows in the case we are setting up are looking at
non-overlapping segments of the chromosome. Here, we will set up two different-sized sliding windows: 25 kbp and 100 kbp. 
These are arbitrary numbers, but are representative of possible window sizes we may investigate across an entire genome. 

    # define window size
    window_size <- 25000

    # estimate fst in sliding windows
    fst_window_25kbp <- differentiation_windows(vcf, popmap, window_size, "fst")

    # estimate dxy in sliding windows
    dxy_window_25kbp <- differentiation_windows(vcf, popmap, window_size, "dxy")

### Modify the above code section to make new objects called 'fst_window_100kbp' and 'dxy_window_100kbp' that have F<sub>ST</sub> and D<sub>XY</sub> estimates in 100 kbp windows.

&nbsp;

## 5. Plot the sliding window estimates.

### Set up dimensions, comparisons, and colors of plot

    # define population pairwise comparisons
    all_combinations <- combn(populations, 2)

    # Set up plotting dimensions with 2 rows and 1 column
    par(mfrow=c(2,1))
    # choose plotting colors (feel free to change)
    plot_colors <- c("orange1", "darkgreen", "purple", "black", "gray", "red")
    # 1st color = Santa Rita Mountains, Arizona x Morelos, Mexico
    # 2nd color = Santa Rita Mountains, Arizona x Pinal Mountains, Arizona
    # 3rd color = Santa Rita Mountains, Arizona x Utah, USA
    # 4th color = Morelos, Mexico x Pinal Mountains, Arizona
    # 5th color = Morelos, Mexico x Utah, USA
    # 6th color = Pinal Mountains, Arizona x Utah, USA

### Plot F<sub>ST</sub>

    plot(c(1,1), col="white", ylim=c(-0.1,1), xlim=c(0,1000000), xlab="Position(bp)", ylab="FST", main="25 kbp windows")
    # plot for each population pairwise comparison
    for(a in 1:ncol(all_combinations)) {
      a_rep <- fst_window_25kbp[fst_window_25kbp$population1 == all_combinations[1,a] & fst_window_25kbp$population2 == all_combinations[2,a],]
      lines((a_rep$end - 12500), a_rep$fst, col=plot_colors[a])
    }

    plot(c(1,1), col="white", ylim=c(-0.1,1), xlim=c(0,1000000), xlab="Position(bp)", ylab="FST", main="100 kbp windows")
    # plot for each population pairwise comparison
    for(a in 1:ncol(all_combinations)) {
      a_rep <- fst_window_100kbp[fst_window_100kbp$population1 == all_combinations[1,a] & fst_window_100kbp$population2 ==   all_combinations[2,a],]
      lines((a_rep$end - 50000), a_rep$fst, col=plot_colors[a])
    }

### Plot D<sub>XY</sub>

    plot(c(1,1), col="white", ylim=c(0,0.7), xlim=c(0,1000000), xlab="Position(bp)", ylab="DXY", main="25 kbp windows")
    # plot for each population pairwise comparison
    for(a in 1:ncol(all_combinations)) {
      a_rep <- dxy_window_25kbp[dxy_window_25kbp$population1 == all_combinations[1,a] & dxy_window_25kbp$population2 == all_combinations[2,a],]
      lines((a_rep$end - 12500), a_rep$dxy, col=plot_colors[a])
    }

    plot(c(1,1), col="white", ylim=c(0,0.7), xlim=c(0,1000000), xlab="Position(bp)", ylab="DXY", main="100 kbp windows")
    # plot for each population pairwise comparison
    for(a in 1:ncol(all_combinations)) {
      a_rep <- dxy_window_100kbp[dxy_window_100kbp$population1 == all_combinations[1,a] & dxy_window_100kbp$population2 == all_combinations[2,a],]
      lines((a_rep$end - 50000), a_rep$dxy, col=plot_colors[a])
    }

&nbsp;

If you ever want to reset plotting dimensions in R or RStudio, just run the following command:
    
    dev.off()

&nbsp;

## 6. Estimating Genetic Structure from SNPs

Next we will be working with a subset of genomic SNPs (n = 5000) from the Brown Creeper populations (the certhia_contact.stru file). 
Here is a map of the subspecies again as well as a sampling map of all the populations used for the SNP dataset we are using today. 
In Panel B, green areas indicate places with dense vegetation, likely different types of forest.

![distribution](https://github.com/jdmanthey/MolEcol2019/blob/master/06_genetic_structure/sampling.png)

Open RStudio and set the working directory to the same one we have been using and download the .stru file from this 
directory on GitHub. 

We'll need to install a package for this exercises, do that here. Packages only need to be installed once. Once they are installed, you 
won't need to run the 'install.packages' command for that package again.

    install.packages("adegenet")
    
Load the libraries:
    
    library("adegenet")
    
### Using PCA

We will use a method called discriminant analysis of principal components (DAPC). This is an extension of the principle
components analysis (PCA) that we discussed in class. If you want more info about this method, a link to the paper is here:
https://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-11-94

To get started we load a structure-formatted file into R. If you are interested in what that looks like, you can open the 
file in a text editor and check it out.

    x <- read.structure("certhia_contact.stru",onerowperind=F,n.ind=24,n.loc=5000,ask=F,sep="\t")

If you already looked at the text input file, you know we have three individuals per population from 8 localities:
    1. Utah (individuals 1-3)
    2. Pinal Mountains (4-6)
    3. Pinaleno Mountains (7-9)
    4. Santa Catalina Mountains (10-12)
    5. Chiricahua Mountains (13-15)
    6. Santa Rita Mountains (16-18)
    7. Huachuca Mountains (19-21)
    8. Central Mexico (22-24)

&nbsp;

Next, we'll use the DAPC program to identify the number of potential clusters in the genetic data. Here, we are setting the
maximum number of clusters to 6 (number of populations) and running for 1e5 iterations. Here, when the program asks you how 
many PCs to retain, choose 20. This is an amount that maintains most of the variation in the data and is less than the 
number of individuals we have sampled. Then, the program will show you a Bayesian Information Criterion (BIC) plot. The BIC
value will be _lowest_ where the number of genetic clusters is most supported. Choose that number of K (hopefully = 2). 

    grp <- find.clusters(x,max.n.clust=6,n.iter=1e5)
    
Next, we will choose the number of PCs again (choose 20 again) as well as the number of discriminant factors to include, which
has to be less than K. Choose the highest number you can based on the value of K you chose.
    
    dapc1 <- dapc(x,grp$grp)

Now, we'll plot the DAPC results in two ways. The first will be the principal components themselves, colored to the group that
each individual (here each point) belongs to:
    
    plot(dapc1$tab[grp$grp==1,1:2], xlim=c(min(dapc1$tab[,1]), max(dapc1$tab[,1])), ylim=c(min(dapc1$tab[,2]), max(dapc1$tab[,2])), col="blue", pch=19, xlab="PC1", ylab="PC2")
    points(dapc1$tab[grp$grp==2,1:2], col="orange2", pch=19)

We can also look at a STRUCTURE-like type of plot showing the assignment of each individual to the two genetic clusters:
    
    compoplot(dapc1, col=c("blue", "orange2"))



