## Microbiome Community Abundance and Diversity - November 2025

Through the next two classes, we will be analyzing some microbiome data (16S rRNA gene fragment) from the guts of _Sitta carolinensis_ birds
from the sky islands of southern Arizona (similar sampling localities to the _Certhia_ datasets we previously analyzed). 
These samples are from four populations: (1) Chiricahua, (2) Huachuca, (3) Pinaleño, and (4) Santa Catalina Mountain Ranges. 
If you need a refresher of where these locations are, check out a previous lab page (toward the bottom): 
https://github.com/jdmanthey/MolEcol2025/tree/main/05_genetic_structre

This list of activities is long, and that is why we are budgeting 2 classes for them all. Go at your own pace.

### 1. Download data

Download the zip file from this link: [link](https://drive.google.com/file/d/1FXaJmU4JY4f9ibGp4hykT31W7NcN5OS_/view?usp=drive_link), unzip it,
and set that directory (class_subsample, not raw_data) as your working directory in RStudio.

### 2. Install necessary packages

    # We are using the DADA2 package (https://benjjneb.github.io/dada2/)
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("dada2")

    # We need the DECIPHER package for some functions
    BiocManager::install("DECIPHER")
    
    # And phyloseq for diversity and plotting of the community data
    BiocManager::install("phyloseq")
    
    # We need a couple other packages that are not included in the above installs
    install.packages("reshape")
    install.packages("picante")
    install.packages("viridis")
    install.packages("phangorn")


### 3. Load packages one at a time to make sure they work

    library(dada2)
    library(phangorn)
    library(DECIPHER)
    library(phyloseq)
    library(ggplot2)
    library(reshape)
    library(RColorBrewer)
    library(plyr)
    library(picante)
    library(viridis)

### 4. Set up of analysis

The data we have today is a subset of the full data collected from the birds. We have 5000 paired-end Illumina sequencing
reads for each individual. These reads target one of the variable regions in the 16S rRNA gene in microorganisms. 

    # set the path to the location of the sequencing reads
    path <- "raw_data"
    list.files(path)
    
The files should have been listed with the second command. If not, the directories are incorrect.

    # sort the order of the forward and reverse reads
    fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
    fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

    # extract sample names from files
    sample.names <- paste(sapply(strsplit(fnFs, "_"), `[`, 1), sapply(strsplit(fnFs, "_"), `[`, 4), sapply(strsplit(fnFs, "_"), `[`, 5), sep="_")
    sample.names

<details>
  <summary>Click to reveal the answer</summary>
  <p></p>
  <p>The answer is 4.</p>
  <p></p>
</details>

You can see we have file names for 11 individuals and that these names include the locality information.

    # Specify the full path to the fnFs and fnRs
    fnFs <- file.path(path, fnFs)
    fnRs <- file.path(path, fnRs)
    
### 5. Error profiling

    # Look at a few sequences and check out their error profiles
    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])
    
    # file paths for putting the filtered reads
    filt_path <- file.path(path, "filtered")
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

    # filter and trim the samples
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)
    head(out)
    
    # learn the error rates for the sequencing
    errF <- learnErrors(filtFs, multithread=FALSE)
    errR <- learnErrors(filtRs, multithread=FALSE)

    # look at plots of errors
    plotErrors(errF, nominalQ=TRUE)
    plotErrors(errR, nominalQ=TRUE)

<details>
  <summary>Click to show quality summaries</summary>
  <p></p>
  <p>The answer is 4.</p>
  <p></p>
</details>

<details>
  <summary>Click to show error plots</summary>
  <p></p>
  <p>The answer is 4.</p>
  <p></p>
</details>

### 6. Dereplicate and call sequence variants

    # dereplicate all reads
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    # rename the dereplicated reads files
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

    # run the main file to call all of the unique sequences
    dadaFs <- dada(derepFs, err=errF, pool=T,multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, pool=T,multithread=TRUE)

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 7. Merge forward/reverse reads, and remove chimeras

    # merge the forward and reverse sequences
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    
    # make a table of all sequences
    seqtab <- makeSequenceTable(mergers)
    # Inspect distribution of sequence lengths
    table(nchar(getSequences(seqtab)))
    
    # keep all mergers with length near the mode (253)
    seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(252,254)]
    
    # remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    # proportion of sequences not chimeric
    sum(seqtab.nochim)/sum(seqtab)

### 8. Summarize filtering

    # summarize the filtering
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
    rownames(track) <- paste0(sapply(strsplit(sample.names, "_"), "[[", 1), "_", sapply(strsplit(sample.names, "_"), "[[", 3))
    track

### 9. Assign taxonomy

Here, we will use the GreenGenes database formatted for DADA2 to classify the 16S sequences we have here. Note that the
classification is only accurate to the family level, and any inferences about genus or species may or may not be accurate.

    taxa <- assignTaxonomy(seqtab.nochim, "gg_13_8_train_set_97.fa.gz", multithread=TRUE)
    unname(head(taxa))

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 10. Create phylogenetic tree of microbial data

    seqs <- getSequences(seqtab.nochim)
    names(seqs) <- seqs
    alignment <- AlignSeqs(DNAStringSet(seqs))
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm)
    plot(treeNJ, show.tip.label=F)

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 11. Make a phyloseq object with the sample data, phylogeny, and sequence variant table

    # Make a data.frame holding the sample data
    samples.out <- rownames(seqtab.nochim)
    samp.number <- sapply(strsplit(samples.out, "_"), `[`, 1)
    species.id <- sapply(strsplit(samples.out, "_"), `[`, 2)
    location <- sapply(strsplit(samples.out, "_"), `[`, 3)
    species.location <- paste(species.id, location, sep=".")
    micro.df <- data.frame(Sample.number=samp.number, Species.ID=species.id, Location=location, Species.Location=species.location)
    rownames(micro.df) <- samples.out


    # construct a phyloseq object
    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),
    		sample_data(micro.df),
    		tax_table(taxa),
    		phy_tree(treeNJ))
    ps
    # remove cyanobacteria
    ps <- subset_taxa(ps, Phylum != "p__Cyanobacteria")

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 12. Composition

Look at a summary of the numbers of unique sequence variants:
    
    summary(sample_sums(ps))

Find the unique phyla, classes, and families in the dataset:

    sort(get_taxa_unique(ps, "Phylum"))
    sort(get_taxa_unique(ps, "Class"))
    sort(get_taxa_unique(ps, "Family"))

Sample a random family:

	sample(sort(get_taxa_unique(ps, "Family")), 1)
	
The following lines of code will summarize the phyla and classes per sample, and summaries across all samples, followed by
writing those results to two tables each (copy and paste the whole thing):

    i_phylum <- c()
    reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Phylum"]))
    reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Phylum"]
    for(a in 1:length(reps)) {
    	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
    	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
    	i_phylum <- cbind(i_phylum, a.rep)
    }
    colnames(i_phylum) <- reps
    inds_names <- rownames(i_phylum)
    i_class <- c()
    reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Class"]))
    reps <- na.omit(reps)
    reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Class"]
    reps_matrix <- na.omit(reps_matrix)
    for(a in 1:length(reps)) {
    	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
    	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
    	i_class <- cbind(i_class, a.rep)
    }
    colnames(i_class) <- reps

    # phylum output					 
    i_tmp <- c()
    for(a in 1:nrow(i_phylum)) {
    	a.rep <- i_phylum[a,]
    	a.rep <- a.rep / sum(a.rep) * 100
    	i_tmp <- rbind(i_tmp, a.rep)
    }
    i_phylum <- i_tmp
    i_phylum_output <- cbind(colnames(i_phylum), apply(i_phylum, 2, mean), apply(i_phylum, 2, min), apply(i_phylum, 2, max))
    i_phylum_output <- data.frame(Phylum=as.vector(sapply(strsplit(i_phylum_output[,1], "__"), "[[", 2)), Mean=(as.numeric(i_phylum_output[,2])),Min=(as.numeric(i_phylum_output[,3])), Max=(as.numeric(i_phylum_output[,4])))

    # class output
    i_tmp <- c()
    for(a in 1:nrow(i_class)) {
    	a.rep <- i_class[a,]
    	a.rep <- a.rep / sum(a.rep) * 100
    	i_tmp <- rbind(i_tmp, a.rep)
    }
    i_class <- i_tmp
    i_class_output <- cbind(colnames(i_class), apply(i_class, 2, mean), apply(i_class, 2, min), apply(i_class, 2, max))
    i_class_output <- data.frame(Class=substr(i_class_output[,1], 4, nchar(i_class_output[,1])), Mean=(as.numeric(i_class_output[,2])),Min=(as.numeric(i_class_output[,3])), Max=(as.numeric(i_class_output[,4])))

    write.table(i_phylum_output, file="i_phylum_output.txt", sep="\t", quote=F, row.names=F)
    write.table(i_class_output, file="i_class_output.txt", sep="\t", quote=F, row.names=F)
    write.table(cbind(inds_names, i_phylum), file="i_phylum_output2.txt", sep="\t", quote=F, row.names=F)
    write.table(cbind(inds_names, i_class), file="i_class_output2.txt", sep="\t", quote=F, row.names=F)

You can look at what the summary tables look like here:
  
    i_phylum_output
    i_class_output

Now we can plot the phyla that are found in the samples:

    sample_number <- as.numeric(sample_data(ps)$Sample.number)
    phyla <- melt(as.data.frame(cbind(sample_number, i_phylum)), id=c("sample_number"), value.name="phyla")
    phyla$sample_number <- as.character(phyla$sample_number)
    ggplot(data=phyla, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(phyla$variable)))))

Do the same plot, but for classes:

    class <- melt(as.data.frame(cbind(sample_number, i_class)), id=c("sample_number"), value.name="class")
    class$sample_number <- as.character(class$sample_number)
    ggplot(data=class, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(class$variable)))))

### 13. Alpha diversity

Here, we'll take a look at a couple ways you can investigate alpha diversity. There is a base function in the phyloseq 
package to look at many types of alpha diversity:

	estimate_richness(ps)

However, we'll just be looking at a few to keep things simple. We'll save them to the object named 'alpha.'
	
	alpha <- cbind(estimate_richness(ps)[,c(1,6,7)], pd((otu_table(ps)@.Data), ps@phy_tree, include.root=F)$PD)
	colnames(alpha) <- c("Observed_SVs", "Shannon", "Simpson", "PD")
	alpha

Remember that sampling depth can influence estimates of diversity, especially if they are incompletely sampled communities.
Let's plot some rarefactions curves to take one look at this factor:

	par(mar=c(4.5,4.5,3,3)) # plotting margins
	rarecurve_table <- otu_table(ps)
  	class(rarecurve_table) <- "matrix"
  	rarecurve(rarecurve_table, step=50, cex=0.5, label=F, xlim=c(0,10000))


If we had the full datasets, the picture may look more completely sampled. Anyways, there are ways of using the full dataset
and also ways of rarefying prior to estimating alpha diversity and beta diversity. We will not get into those methods in these
activities.

Let's plot a couple types of alpha diversity:

	plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"))

We can also plot the relationships among each of the estimates of alpha diversity. Here, Observed_SVs = number of observed
sequence variants and PD = phylogenetic diversity.

	plot(alpha, pch=19, cex=1)

### 14. Beta diversity

Now we'll measure beta diversity among communities in a few different ways. 

First, we'll make a plot of the Bray-Curtis distance:
	
	ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
	wu.df <- data.frame(
		x = as.vector(ord.nmds.bray$points[,1]),
		y = as.vector(ord.nmds.bray$points[,2]),
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.nmds.bray, color="Species.Location", shape="Species.Location", title="Bray NMDS") + geom_point(size=3)

And unweighted Unifrac distance:

	ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)
	wu.df <- data.frame(
		x = ord.unweighted.unifrac$vectors[,1],
		y = ord.unweighted.unifrac$vectors[,2],
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.unweighted.unifrac, color="Species.Location", shape="Species.Location",title="Unweighted Unifrac PCoA") + geom_point(size=3)
	
And weighted Unifrac distance:

	ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)
	wu.df <- data.frame(
		x = ord.weighted.unifrac$vectors[,1],
		y = ord.weighted.unifrac$vectors[,2],
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.weighted.unifrac, color="Species.Location", shape="Species.Location",title="Weighted Unifrac PCoA") + geom_point(size=3)
