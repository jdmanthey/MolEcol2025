# BLAST TUTORIAL - 28 August 2025

Continuing with the theme of DNA barcoding, today will include activities related to searching databases for matches to 
sequences of unknown affinities. We will be using some DNA sequences that are unknown to you in combination with a search
algorithm (BLAST) and a database of known sequences (available on [NCBI's website](https://www.ncbi.nlm.nih.gov/)). Below, I
will cover a tutorial of using one of the BLAST algorithms with the online database, in the hopes that this information will
allow you to complete a similar exercise with practice sequence(s) and questions for this week's computational activity. 
Once you have gone through and read the tutorial, you will be given an exercise sheet to work through. 

&nbsp;

First some definitions:

BLAST: Basic Local Alignment Search Tool; this tool attempts to find identity between a query sequence and a database of 
known sequences

NCBI: National Center for Biotechnology Information; NCBI has a website that includes many browser-based features, and also
has tools (including BLAST) that can be downloaded for offline and custom applications.

FASTA: pronounced "fast"-"AY" (like the letter A). FASTA is a format for keeping DNA, RNA, or protein sequence data in a text
format. A FASTA file can have one or more sequences and has a minimum of two lines per sequence. The first line includes a 
carrot symbol ">" and then the name of the sequence, and the second (and any subsequent lines prior to another >) contain the
sequence information. See below:

    >Test_sequence
    ATGAACCCCCAAGCAAAACTGATTTTCGCCCTTAGCCTAACTCTAGGCACAACCCTCACAATCTCAAGCA
    ACCATTGAATTATAGCCTGAACCGGACTTGAAATCAATACCCTGGCTATCCTCCCACTAATCTCAAAATC
    TCATCACCCCCGCGCCATTGAAGCTGCAACCAAGTACTTCTTAGTTCAAGCAGCTGCCTCTGCCTTAGTC
    CTGTTCTCAAGCATAACCAATGCTTGACACACTGGCCAATGGGATATTACCCAACTATCGCACCCAACTT
    CATGCCTGATCTTAACTGTAGCTCTTGCAATAAAACTAGGCCTAGTTCCATTCCACTTCTGATTCCCAGA
    AGTCCTCCAAGGTTCTTCCCTAATTGTAGGCCTCCTCCTATCCACTCTCATAAAATTCCCACCAATTACC


## The BLAST Algorithm

BLAST is designed to align a query sequence against a subject database. We will not get into the math here. BLAST does not 
attempt to globally align sequences (i.e., try to match an entire sequence against the entire length of another sequence), but
rather tries to locally match sequences in areas where they are very similar. There are several different implementations of 
the BLAST algorithm, each slightly varying based on the type of input (query) used and the type of data in the database
(subject). Below is a table of three major types of BLAST algorithms you may run into (although there are others!). 

| Search Type | Input      | Database   |
|-------------|------------|------------|
| blastn      | nucleotide | nucleotide |
| blastp      | protein    | protein    |
| blastx      | nucleotide | protein    |

## An empirical example

Here, I will lay out an example of a BLAST search as well as the interpretation of the different methods and components 
involved. First, you need to have a sequence in hand. I will not give the sequence here, but for the exercises, there will 
be sequences available for you in FASTA format. Below is a step by step example of how to do a simple BLAST search, using your 
sequence of interest with the NCBI nucleotide database (a lot of sequencing data from years of scientific study). 

1. Go to NCBI's [website](https://www.ncbi.nlm.nih.gov/). Here you will click on the popular application link for BLAST on 
the right indicated with a red arrow in the image below. 
![NCBI page](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_1.png)
&nbsp;

2. Next, you will be on a page where you can choose some of the most popular BLAST options (the large link buttons). For this
example, we will be choosing a nucleotide BLAST to a nucleotide database, aka blastn. 
![BLAST page](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_2.png)
&nbsp;

3. Once we are on the BLAST webpage for blastn, we are able to input a sequence, choose which default blast applications we
want to use, and even set custom search settings if we so choose. Here, we will only use the default settings, but we will 
change the "Program Selection" section to optimize for blastn rather than megablast. The megablast option is for quicker 
searches when limiting matches to be at least ~95% identical to your query sequence. blastn allows you to obtain matches that
are more distantly related, and potentially gives you more matches, a chance to check all results, and, if you have an unknown 
query sequence, more opportunity to find a match in the BLAST database. Here, we also want to make sure that the database 
we are searching against is the "core nucleotide database (core_nt)."
![blastn option](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_3.png)
&nbsp;

4. As mentioned above, you can also change the BLAST algorithm's settings to your liking at the bottom of the page. Check out
the various settings you could potentially change. We won't change anything here, but I just wanted to make sure you were 
aware of the possibility. 
![settings](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_4.png)
&nbsp;

5. Now is the time to do the search. We can enter one or more sequences in FASTA format in the search box (top arrow). Once we
are satisfied that we have the input and settings we want, click on the BLAST oval button (bottom arrow). 
![BLAST!](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_5.png)
&nbsp;

6. After you click the "BLAST" button, you will have to wait for your results. If you are signed in to NCBI with an account
and have it set up with your contact, you can get emails when the job(s) finish. But usually, if you just have one or two 
searches, you can simply wait for the results. 
![waiting](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_6.png)
&nbsp;

7. Once the search is completed, you will be redirected to a page with the results. Here, you are able to see the details
about matches, as well as ways to filter the results. This following screenshot shows the summary of the results, with four
tabs below that give details about the matches. On the top, the major points of interest are: (1) RID: or the run ID, so that
you can revisit results (before they expire) if you did not save the search; (2) the program used; (3) the database you
BLASTed your query against; (4) the type of molecule you designated the search for (here = dna); (5) the length of the 
sequence you entered (Query Length). 
![top of BLAST results](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_7.png)
&nbsp;

8. The screenshot below shows the "Descriptions" tab of the results. Here, you get a line-by-line summary of results about
the sequences in the database that your query matched, as well as some statistics to tell you how well it matched. The 
columns of this tab are as follows: (1) Description: the sequence's name and some attributes uploaded by the scientist that 
submitted this sequence to the NCBI GenBank database. The amount of information here will vary, but hopefully it will at
least tell you which organism the sequence came from, and the locus / gene that was uploaded; (2) Max Score: the best BLAST
alignment score for the query to subject match. This is based on the search algorithm's math and the details about the match
between the query and the subject. The higher the score the better; (3) Total Score: this is the sum of all the matches for 
your sequence to a subject sequence in the database. Sometimes BLAST results will have one simple match without gaps, and 
sometimes there will be multiple matches with gaps. If the former, the max and total score will be the same. If there are 
multiple matches with gaps, this column will be the sum of the scores of all the matches together; (4) Query Cover: this
represents the % of your query sequence that was covered in the subject sequence match. This ranges from 0-100% and the 
higher the better; (5) E Value: this is the "expectation value" and is a quantitative score of the quality of the match using
multiple metrics from the BLAST algorithm. This value represents the probability of a match by random given the structure of 
your query sequence and the size and structure of the database you searched against. The closer to zero the E value is, the 
better, but remember that this value can change depending on the characteristics of both your search query and the database
used; (6) Per. Ident: the percentage of the query sequence that is identical to the subject sequence in the portion that was
covered by the BLAST search match (i.e., the Query Cover). A higher proportion of matched overlap is generally a better match;
(7) Accession: the database accession of the database subject sequence that the query matched. You can click on this and the
link will take you to a page that gives you all the information about this match. 
![description tab](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_8.png)
&nbsp;

9. A graphical portrayal of the descriptions tab. Each colored line represents a BLAST search match. The color of the line
indicates the max score of the BLAST match. Red = a match greater than a 200 value for the BLAST score. The length of the 
line shows the query coverage of the match. In this case, the query coverage for each match is the entire sequence, and 
therefore the bar is 100% of the horizontal pane in the figure. Each of these lines should be a hyperlink and allow you 
to click on the bar to reach the alignment (next tab) of the match. 
![graphic tab](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_9.png)
&nbsp;

10. The alignments tab of the results shows you the specific alignment of the query sequence vs. different database subject
sequences. Here you can see at the topic general information about the description and accession of the subject sequence. If
you click on the accession (arrow), you will be taken to the webpage for that sequence from the database you searched against
(step 12 here). Below the general information, you see the alignment. This wraps across multiple rows, with each row = 60 base
pairs (bp) of alignment. When the sequences are identical, the Query and Sbjct rows should be the same bases with a line
connecting them. If the sequences have mismatches, there will not be a line connecting the two. 
![alignment tab](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_10.png)
&nbsp;

11. The last tab of the BLAST results is about the taxonomy of the subject match in the database. Here, you will find info
about which taxonomic group / species your query matched, as well as links to explore further. 
![taxonomy tab](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_11.png)
&nbsp;

12. If you click on one of the accession links in the BLAST results tabs, you will be directed to a page that looks like
this screenshot or something similar. This page will give you all of the available descriptions about the sequence, 
including which organism it is from, who submitted it, if it is linked to any scientific publications, and the sequence
itself. If the sequence makes a protein, the amino acids it encodes should also be present. The two arrows below are 
indicating the hyperlink to a linked publication (may or may not exist), as well as a hyperlink to information about the 
taxonomic group that the sequence is from (should almost always exist). 
![accession_page](https://github.com/jdmanthey/MolEcol2025/blob/main/01_blast/Fig_12.png)
&nbsp;

### That concludes the basic BLAST tutorial, and you should move on to the exercise sheet and answer the questions with the unknown sequence(s) you have been assigned.


