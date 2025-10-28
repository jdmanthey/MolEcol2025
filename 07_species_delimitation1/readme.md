## Species delimitation in _Piranga_ tanagers

Today we will be loading phylogenies, plotting them in different ways, and running quick preliminary species delimitation
analyses. We will be working with some phylogenetic data results from a previous publication: 
[https://academic.oup.com/sysbio/article/65/4/640/1753369]


#### Download the data and set up the working directory.

Download the trees.zip file to your computer, unzip the file, and set that as your working directory for these exercises.


#### Check that you are in the right directory

    list.files()

You should see all the files from the folder you unzipped. If not, you are in the wrong directory.


#### Install packages

    install.packages("ape")
    install.packages("phangorn")
    install.packages("phytools")
    
    
#### Load packages

    library(ape)
    library(phangorn)
    library(phytools)


#### Read in the trees file and plot

    # read in file
    ptrees <- read.nexus("piranga.trees")

Since there are too many trees for some of these computers to work with quickly, we will randomly sample 100:

    ptrees <- ptrees[sample(seq(from=2, to=10001, by=1), 100)]

This file contains 10,001 samples of phylogenetic tree estimation from the same dataset. We can plot one of the trees like so:

    plot(ptrees[[1]], cex=0.5)

## Answer Question #1

In the trees, most of the individuals are from the genus _Piranga_ and the outgroup is the Northern Cardinal, (_Cardinalis cardinalis_).

There are other ways to visualize trees as well. We can plot an unrooted tree with the following code. This is the exact same phylogeny
as the previous plot.

    plot(unroot(ptrees[[1]]),type="unrooted",cex=0.6, use.edge.length=FALSE,lab4ut="axial", no.margin=TRUE)

Or we can plot a circular style tree:

    plotTree(ptrees[[1]],type="fan",fsize=0.7,lwd=1, ftype="i", cex=0.5)

Since there are 100 phylogenies that we sampled, we can also plot all the trees at once. Here, each tree will be partially translucent so
that we can see all the trees simultaneously:

    densiTree(ptrees, cex=0.5, consensus=ptrees[[1]])

## Answer Question #2

Now we'll work with a subset of the phylogenies, excluding the outgroup and a clade with all long branches.

    tips_to_drop <- c("C_cardinalis", "P_leucoptera", "P_erythrocephala", "P_rubriceps")
    ptrees_pruned <- lapply(ptrees, drop.tip, tips_to_drop)
    class(ptrees_pruned) <- "multiPhylo"

Let's plot the pruned trees:

    densiTree(ptrees_pruned, consensus=ptrees_pruned[[1]])

## Answer Question #3

## Species Delimitation

We will be using some scripts from Noah Reid to run a species delimitation method called bGMYC: [http://nreid.github.io/software/]
We'll load these scripts here:

    source("bgmyc.r")
    
#### Set base plotting scheme

Some of the plotting functions in the bGMYC scripts will change the base plotting settings. Here, we will save the original 
settings:

    old_par <- par(no.readonly = TRUE)
    
If you need to reset the plotting functions later (if plots look like they have strange sizes or margins), use:

    par(old_par)

#### Running bGMYC

We will run the species delimitation MCMC for the pruned trees. It is a short run for time purposes. If you wanted to do a real analysis
for publication, you'd likely run this much longer (1-2 orders of magnitude longer).

    result_multi_piranga <- bgmyc.multiphylo(ptrees_pruned, mcmc=10000, burnin=4000, thinning=200, t1=2, t2=12, start=c(1,1,11))
    
Now, we will check out the results of delimiting the tips of the phylogeny into different assigned groups:

    result_probmat_piranga <- spec.probmat(result_multi_piranga)
    
And plot:

    plot(result_probmat_piranga, ptrees_pruned[[1]])

## Answer Question #4


