# Starting this script on July 4, 2020. 
#I am using the callahan et al 2016 paper still but now I'm going to use the PAIR merged files instead of just the forward reads 
# The latest previous version of this was 200622_MP_PAIRED_practice, so I'm going to copy/paste and edit that script 
#the objective is to make some nice figures, and ideally a PCoA plot! 
#Re-running through this on July 20th.. For the millionith time 
#Also using this tutorial for some parts:https://benjjneb.github.io/dada2/tutorial_1_8.html
#this tutorial was used to help make the figures post creating a phyloseq object: http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# I have a lot of samples, so keeping this organized might be a challenge
## = output or potential output 
#Re-running entire code starting 201108

#Load the packages####
#usually have to run this chunk of code twice
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE)) #needed to do this step for updated version of R 4.0.3
  install.packages("BiocManager")

BiocManager::install("BiocStyle")
library("BiocStyle")
library("knitr")
library("BiocStyle")


#opts_chunk$set(cache = FALSE,fig.path="dadafigure/")
#read_chunk(file.path("src", "bioinformatics.R")) #this didnt work 

if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])}
.inst <- bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)}
#run this loop if the packages don't install: 

BiocManager::install("DECIPHER")
BiocManager::install("phangorn") 
BiocManager::install("phyloseq")
install.packages("gridExtra")
BiocManager::install("dada2", version = "3.12") 


cran_packages <- c("ggplot2", "gridExtra")
bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- cran_packages %in% installed.packages()
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

#This is the output I should see, if any say false, something hasn't installed correctly:
##ggplot2 gridExtra     dada2  phyloseq  DECIPHER  phangorn 
##TRUE      TRUE      TRUE      TRUE      TRUE      TRUE 

###### Data processing (Filtering, Trimming, Errors, taxa assignment etc.) from lines 58-325

#QUALITY CHECK, FILTER, TRIM ######

#path to the files being used 

AB_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files AB")
ABK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files ABK")
Fe_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files Fe")
FeK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged FeK")
S_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files S")
SK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files SK")
NL_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files NL")
NK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged NK")
all_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/ALL_PEAR")


#path to a new file where I'll put the filtered 
AB_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/AB_filt2")
ABK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/ABK_filt")
Fe_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/Fe_filt")
FK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/FK_filt")
S_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/S_filt")
SK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/SK_filt")
NL_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/NL_filt")
NK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/NK_filt")

#Start with the live 

fns_AB = sort(list.files(AB_path, full.names = TRUE)) #done 
fns_Fe = sort(list.files(Fe_path, full.names = TRUE)) 
fns_S = sort(list.files(S_path, full.names = TRUE))
fns_NL = sort(list.files(NL_path, full.names = TRUE))

######## AB3-T4 is missing? ... D'Arcy said that it must've been a crappy sequence, because it's not in the raw data 
#check quality before 
ABi <- sample(length(fns_AB), 17) 
for(i in ABi) { print(plotQualityProfile(fns_AB[i]) + ggtitle("AB paired")) }
#AB1_t4 has poor quality in the middle of the sequence, what to do? 
#it also seems like the begging of the reverse reads are bad quality, not sure what to do there either?


#NEED TO REPEAT THIS STEP FOR ALL. 


#filtering & trimming
#AB
filt_AB <- file.path(AB_filt_path, basename(fns_AB))
for(i in seq_along(fns_AB)) {
  filterAndTrim(fns_AB, 
                filt_AB, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}
#Fe
if(!file_test("-d", Fe_filt_path)) dir.create(Fe_filt_path)
filt_Fe <- file.path(Fe_filt_path, basename(fns_Fe))
for(i in seq_along(fns_Fe)) {
  filterAndTrim(fns_Fe, 
                filt_Fe, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}
#S
if(!file_test("-d", S_filt_path)) dir.create(S_filt_path)
filt_S <- file.path(S_filt_path, basename(fns_S))
for(i in seq_along(fns_S)) {
  filterAndTrim(fns_S, 
                filt_S, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}
#NL
if(!file_test("-d", NL_filt_path)) dir.create(NL_filt_path)
filt_NL <- file.path(NL_filt_path, basename(fns_NL))
for(i in seq_along(fns_NL)) {
  filterAndTrim(fns_NL, filt_NL, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}

  #CHECK TRIMMED QUALITY#
#ii <- sample(length(fns_all), 158) # 158 is the number of samples
#for(i in ii) { print(plotQualityProfile(fns[i]) + ggtitle("paired")) }
### I didn't run this previous bit of code because trying to look at 158 samples plotted at once sounds absurd, instead, I'm going to split them up and look at them
##Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale [this message comes up for every plot]


Fe_filt_i <- sample(length(filt_Fe), 18) 
for(i in Fe_filt_i) {print(plotQualityProfile(filt_Fe[i]) + ggtitle("Fe filt paired")) }
##Error: Input/Output
##no input files found
##dirPath: /Users/judymalas/Desktop/Landfill sequencing data processing /PEAR_merged/Fe_filt/Fe2_t5_PEAR.assembled.fastq
##pattern: character(0)

NL_filt = sample(length(filt_NL), 30)
for(i in NL_filt) {print(plotQualityProfile(filt_NL[i])+ ggtitle("NL filt"))}

#Trim off 10 from each end for ALL
#if(!file_test("-d", all_filt_path)) dir.create(all_filt_path)
#filt_all <- file.path(all_filt_path, basename(fns_all))
#for(i in seq_along(fns_all)) { filterAndTrim(fns_all, filt_all, rev = NULL,  filt.rev =NULL,  trimLeft=c(10), truncLen=c(290), maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}

#what does the dark grey bar at the top mean? e.g. in ABK3-t7
##### so I stopped here last time and R broke, we're gonna try again
# Today is July 20th, I tried this a couple days ago and R shut down. But when I look at  the "all_filt" folder 
#in my PEAR_merged directory, there are files in there.. so I'm not sure if it worked or not. I'm gonna try 
  #so trimming all of the files in one run is not an option
  #today (july 20th) I am trimming each of the live experiments seperately, and then maybe I can put them all in one file together.






#######DEREPLICATION####
#The sequence data is imported into R from demultiplexed fastq files (i.e. one fastq for each sample) 
#and simultaneously dereplicated to remove redundancy. (Callahan, 5) 
#We name the resulting derep-class objects by their sample name. (Callahan, 5)
#derep <- derepFastq(filt_all)
##Error in derepFastq(filt_all) : Not all provided files exist.

#going to try with AB filtered 
derep_AB = derepFastq(filt_AB)
#that seemed to work so I'm not sure what the issue is
sam.names <- sapply(strsplit(basename(filt_AB), "_"), `[`, 1)
names(derep_AB) <- sam.names
##Convergence after  9  rounds

#derep_All = derepFastq(filt_all, verbose=TRUE)
##Error in derepFastq(filt_all, verbose = TRUE) : 
  ##Not all provided files exist.
#still not sure what the issue is 

#https://github.com/benjjneb/dada2/issues/711 
#check if the file exists, and only use the ones that do 
#exists <- file.exists(filt_all)
#print(exists) #so all are false 

#exists_ab= file.exists(filt_AB)
#print(exists_ab)
#all are true.. so it looks like "filt_all" is not valid
#ok, might have to re-do some of the filtering and trimming, can't do it all in one shot

#derep_Fe = derepFastq(filt_Fe)
#that didn't seem to work so I'm not sure what the issue is
 
#exists_Fe = file.exists(filt_Fe)
#print(exists_Fe)#two don't exist.. I'm not sure what that means? Maybe the quality was bad. (#9 and #11)
#only derep the ones that exist 
#derep_Fe = derepFastq(filt_Fe[exists_Fe], verbose= TRUE)
#names(derep_Fe) = sample.names[exists_Fe]
##Error in sample.names[exists_Fe] : 
##object of type 'closure' is not subsettable
# i need to be able to change the names 
# i re-ran the filter and trim, and now they all exist. 

#exists_S = file.exists(filt_S)
#print(exists_S) #all exist

#exists_NL = file.exists(filt_NL)
#print(exists_NL) #they all exist 

#derep Fe, S, NL
derep_Fe = derepFastq(filt_Fe)
sam.names_Fe <- sapply(strsplit(basename(filt_Fe), "_"), `[`, 1)
names(derep_Fe) <- sam.names_Fe

derep_S = derepFastq(filt_S)
sam.names_S <- sapply(strsplit(basename(filt_S), "_"), `[`, 1)
names(derep_S) <- sam.names_S

derep_NL = derepFastq(filt_NL)
sam.names_NL<- sapply(strsplit(basename(filt_NL), "_"), `[`, 1)
names(derep_NL) <- sam.names_NL


#####LEARNING ERRORS #####
#The DADA2 method relies on a parameterized model of substitution errors to distinguish 
#sequencing errors from real biological variation.
#Because error rates can (and often do) vary substantially between sequencing runs and PCR protocols, 
#the model parameters can be discovered from the data itself using a form of unsupervised learning in which sample inference is alternated with parameter estimation until both are jointly consistent.
#Parameter learning is computationally intensive, as it requires multiple iterations of the sequence inference algorithm, and therefore it is often useful to estimate the error rates from a (sufficiently large) subset of the data. (Callahan, 5)

#dd_AB = dada(derep_AB[1:17], err=NULL, selfConsist=TRUE) #1:17 is the number of samples, change as needed 
#dd_Fe = dada(derep_Fe[1:18], err=NULL, selfConsist=TRUE) ## Self-consistency loop terminated before convergence.
#dd_S = dada(derep_S[1:18], err=NULL, selfConsist=TRUE) ## Self-consistency loop terminated before convergence.
#dd_NL = dada(derep_NL[1:30], err=NULL, selfConsist=TRUE)
#data objects were created for dd_Fe and dd_S, but I'm not sure if it means if the loop was terminated before convergence

#then inspect the fit b/w the observed error rates (black points)and the fitted error rates (black lines)
#plotErrors(dd_AB)
#plotErrors(dd_Fe)
#plotErrors(dd_S)
#plotErrors(dd_NL)

#I'm going to try this another way: https://benjjneb.github.io/dada2/tutorial_1_8.html
# in this tutorial, they learn the error rates before dereplication for some reason.. but I'm going to try after and see wat happens

err_AB = learnErrors(derep_AB, multithread = TRUE) 
err_Fe = learnErrors(derep_Fe, multithread = TRUE)
err_S = learnErrors(derep_Fe, multithread = TRUE)
err_NL = learnErrors(derep_Fe, multithread = TRUE)

plotErrors(err_AB, nominalQ = TRUE) #plot of the error rates 

#####DADA2 SEQUENCE INFERENCE METHOD#####
dada_AB <- dada(derep_AB, err=err_AB, pool=TRUE, multithread = TRUE)
dada_Fe <- dada(derep_Fe, err=err_Fe, pool=TRUE, multithread = TRUE)
dada_S <- dada(derep_S, err=err_S, pool=TRUE, multithread = TRUE)
dada_NL <- dada(derep_NL, err=err_NL, pool=TRUE, multithread = TRUE)

##2 samples were pooled: 38189 reads in 23861 unique sequences.

#CONSTRUCT SEQUENCE TABLE####
#seqtab.all <- makeSequenceTable(dadas[!grepl("Mock", names(dadas))])
#seqtab <- makeSequenceTable(dada_AB, dada_Fe)

##Error in FUN(X[[i]], ...) : 
##Unrecognized format: Requires named integer vector, fastq filename, dada-class, derep-class, sequence matrix, or a data.frame with $sequence and $abundance columns.
#The issue is I'm not sure if I can run these seperately, I might have to go back and do all the samples together in the same run 
 
#actually I may as well try to do them seperately first, and then try doing it all together

seqtab.AB = makeSequenceTable(dada_AB)
seqtab.Fe = makeSequenceTable(dada_Fe)
seqtab.S = makeSequenceTable(dada_S)
seqtab.NL = makeSequenceTable(dada_NL)


#https://benjjneb.github.io/dada2/bigdata.html
#this tutorial is for big data so I may try to follow it to go back and look @ my samples 

#REMOVE CHIMERAS####
#remove chimeras by comparing each inferred sequence to the others in the table, and removing those that can be reproduced by stitching together two more abundant sequences.

nochim_seqtab.AB <- removeBimeraDenovo(seqtab.AB, method="consensus", multithread=TRUE, verbose=TRUE)
nochim_seqtab.Fe <- removeBimeraDenovo(seqtab.Fe, method="consensus", multithread=TRUE, verbose=TRUE)
nochim_seqtab.S <- removeBimeraDenovo(seqtab.S, method="consensus", multithread=TRUE, verbose=TRUE)
nochim_seqtab.NL <- removeBimeraDenovo(seqtab.NL, method="consensus", multithread=TRUE, verbose=TRUE)

dim(nochim_seqtab.AB)
## 17 381
sum(nochim_seqtab.AB)/sum(seqtab.AB)
##0.8655657

#CLASSIFY TAXONOMY####
#The dada2 package implements the naive Bayesian classifier method for this purpose. This classifier compares sequence variants to a training set of classified sequences (Callahan, 8).
#I can use the silva dataset instead of whatever they used here 

ref_fasta = "~/Desktop/Landfill sequencing data processing /silva_nr_v138_train_set.fa.gz"
print(ref_fasta)

taxtab.AB = assignTaxonomy(nochim_seqtab.AB, refFasta = ref_fasta)
colnames(taxtab.AB) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

taxtab.Fe = assignTaxonomy(nochim_seqtab.Fe, refFasta = ref_fasta)
colnames(taxtab.Fe) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

taxtab.S = assignTaxonomy(nochim_seqtab.S, refFasta = ref_fasta)
colnames(taxtab.S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

taxtab.NL = assignTaxonomy(nochim_seqtab.NL, refFasta = ref_fasta)
colnames(taxtab.NL) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#wohooo that worked. And it looks much like it was able to classify up to the genus level so we're good to move on 

 
#LOAD METADATA #####

#The phyloseq package organizes and synthesizes the different data types from a typical amplicon sequencing experiment
#into a single data object that can be easily manipulated. 
#The last bit of information needed is the sample data contained in a .csv file. <--- think that's the metadata

AB_metadata = read.csv("~/Desktop/Landfill sequencing data processing /PEAR_merged/AB_metadata.csv")
Fe_metadata = read.csv("~/Desktop/Landfill sequencing data processing /PEAR_merged/Fe_metadata.csv")
S_metadata = read.csv("~/Desktop/Landfill sequencing data processing /PEAR_merged/S_metadata.csv")
NL_metadata = read.csv("~/Desktop/Landfill sequencing data processing /PEAR_merged/NL_metadata.csv")

#small interlude: 
#shit, so I may have fudged something real bad. It looks like The AB time points merged together into one bottle.. Looking at the seqtab, AB1t1-t7 are exactly the same 
#i'm going to try re-running everything on just the AB bottles to see if I get a different result
#if i do I'm going to have to re-run everything again
#actually the other seq.tables don't look the same as AB, they look like they have unique values, so maybe I just need to redo AB 
#Alright we're back. I only had to redo AB 

#This part is about metadata files match the seqtabs 

#row names in nochim seq tab and in the metadata file need to be identitical
all(rownames(seqtab.AB) %in% AB_metadata$sample)
#this needs to be true in order to continue, the sample names in the metadata table need to be the same as in the seq table 
rownames(AB_metadata) = AB_metadata$sample #this added a row name that was the same as the seqtab row name
rownames(nochim_seqtab.AB) = c("AB1","AB1.1","AB1.2","AB1.3","AB1.4","AB1.5", "AB2",   "AB2.1", "AB2.2", "AB2.3", "AB2.4" ,"AB2.5", "AB3"  , "AB3.1", "AB3.2", "AB3.3" ,"AB3.4")
rownames(nochim_seqtab.AB)
identical(rownames(nochim_seqtab.AB), rownames(AB_metadata)) 
##TRUE

# NL 
all(rownames(seqtab.NL) %in% NL_metadata$sample)
#this needs to be true in order to continue, the sample names in the metadata table need to be the same as in the seq table 
rownames(NL_metadata) = NL_metadata$sample #this added a row name that was the same as the seqtab row name
print(NL_metadata$sample)
rownames(nochim_seqtab.NL) =  c("NL1", "NL1.1", "NL1.2", "NL1.3", "NL2", "NL2.1", "NL2.2", "NL2.3", "NL3", "NL3.1", "NL3.2", "NL3.3", "NL4", "NL4.1", "NL4.2", "NL4.3", "NL4.4", "NL4.5", "NL5", "NL5.1", "NL5.2", "NL5.3", "NL5.4", "NL5.5", "NL6", "NL6.1", "NL6.2", "NL6.3", "NL6.4", "NL6.5")
#this step is because otherwise rownames appear like this: 
##> rownames(nochim_seqtab.NL)
##[1] "NL1" "NL1" "NL1" "NL1" "NL2" "NL2" "NL2" "NL2" "NL3" "NL3" "NL3" "NL3" "NL4" "NL4" "NL4" "NL4"
##[17] "NL4" "NL4" "NL5" "NL5" "NL5" "NL5" "NL5" "NL5" "NL6" "NL6" "NL6" "NL6" "NL6" "NL6"
rownames(nochim_seqtab.NL)
rownames(NL_metadata)
identical(rownames(nochim_seqtab.NL), rownames(NL_metadata)) 

#S
all(rownames(seqtab.S) %in% S_metadata$sample)
#this needs to be true in order to continue, the sample names in the metadata table need to be the same as in the seq table 
rownames(S_metadata) = S_metadata$sample #this added a row name that was the same as the seqtab row name
print(S_metadata$sample)
rownames(nochim_seqtab.S)
rownames(nochim_seqtab.S) =  c ("S1", "S1.1", "S1.2", "S1.3", "S1.4", "S1.5", "S2", "S2.1", "S2.2", "S2.3", "S2.4", "S2.5", "S3", "S3.1", "S3.2", "S3.3", "S3.4", "S3.5")
rownames(nochim_seqtab.S)
rownames(S_metadata)
identical(rownames(nochim_seqtab.S), rownames(S_metadata)) 


#Fe
all(rownames(seqtab.Fe) %in% Fe_metadata$sample)
#this needs to be true in order to continue, the sample names in the metadata table need to be the same as in the seq table 
##FALSE 
#Fe is fudged up with the row names, so I have to change them in the metadata file to match seqtab.Fe

rownames(Fe_metadata) = Fe_metadata$sample #this added a row name that was the same as the seqtab row name
print(Fe_metadata$sample)
rownames(nochim_seqtab.Fe)
rownames(nochim_seqtab.Fe) =  c ("Fe1", "Fe1.1", "Fe1.2", "Fe1.3", "Fe1.4", "Fe1.5", "Fe2", "Fe2.1", "Fe2.2", "fe2", "Fe2.3", "Fe2.4", "Fe3", "Fe3.1", "Fe3.2", "Fe3.3", "Fe3.4", "Fe3.5")
rownames(nochim_seqtab.Fe)
rownames(Fe_metadata)
identical(rownames(nochim_seqtab.Fe), rownames(Fe_metadata)) 



#CONSTRUCT PHYLOGENTIC TREES, the tree will be added to the phyloseq object#######
#We begin by performing a multiple-alignment using the DECIPHER R package11.

seqs_AB <- getSequences(nochim_seqtab.AB)
names(seqs_AB) <- seqs_AB # This propagates to the tip labels of the tree
alignment_AB <- AlignSeqs(DNAStringSet(seqs_AB), anchor=NA)

#The phangorn R package is then used to construct a phylogenetic tree. 
#Here we first construct a neighbor-joining tree, 
#and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

phang.align <- phyDat(as(alignment_AB, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


#NL tree construction
seqs_NL <- getSequences(nochim_seqtab.NL)
names(seqs_NL) <- seqs_NL # This propagates to the tip labels of the tree
alignment_NL <- AlignSeqs(DNAStringSet(seqs_NL), anchor=NA)
phang.align <- phyDat(as(alignment_NL, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR_NL <- update(fit, k=4, inv=0.2)
fitGTR_NL <- optim.pml(fitGTR_NL, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


#MAKE PHYLOSEQ OBJECTS######
#AB
##This is what I used before I knew I had to add a tree to the ojbect: 
#going tree-less this time ps_AB = phyloseq(otu_table(nochim_seqtab.AB, taxa_are_rows = FALSE), sample_data(AB_metadata), tax_table(taxtab.AB), phy_tree(fitGTR$tree))
ps3_AB = subset_samples(ps_AB, sample!="AB1.3") #excluding AB1.3 from future analses because the read count is too low

#NL
#use this if you want to include the tree: ps_NL = phyloseq(otu_table(nochim_seqtab.NL, taxa_are_rows = FALSE), sample_data(NL_metadata), tax_table(taxtab.NL), phy_tree(fitGTR_NL$tree))
ps_NL = phyloseq(otu_table(nochim_seqtab.NL, taxa_are_rows = FALSE), sample_data(NL_metadata), tax_table(taxtab.NL))
#Can I filter the ps obeject to exlude NL 1-3?
#this chunk removes NL 1-3 from the phyloseq object "new_ps_NL"
new_ps_NL = subset_samples(ps_NL, sample!="NL1" & sample!="NL1.1" & sample!="NL1.2" & sample!="NL1.3"& sample!="NL2"& sample!="NL2.1"& sample!="NL2.2"& sample!="NL2.3"& sample!="NL3"& sample!="NL3.1"& sample!="NL3.2"& sample!="NL3.3")
new_ps_NL

#S
ps_S = phyloseq(otu_table(nochim_seqtab.S, taxa_are_rows = FALSE), sample_data(S_metadata), tax_table(taxtab.S))

#Fe
ps_Fe = phyloseq(otu_table(nochim_seqtab.Fe, taxa_are_rows = FALSE), sample_data(Fe_metadata), tax_table(taxtab.Fe))





##### MAKING FIGURES #######
#BAR PLOTS TRIAL ONE: IGNORE THIS ONE#####
#Make Bar plot#
top20 <- names(sort(taxa_sums(ps_AB), decreasing=TRUE))[1:20] #top 20 sequences
ps.top20 <- transform_sample_counts(ps_AB, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Genus") + facet_wrap(~Bottle, scales="free_x")
#i think I'm going to need to do some filtering because these bar plots look like shit
#Plus I don't like these bar plots in general, there must be a different way


###FILTERING TRAIL ONE: IGNORE THIS ONE#### 

rank_names(ps_AB) #show available ranks in the dataset 
table(tax_table(ps_AB)[, "Phylum"], exclude = NULL) #create table, number of features for each phyla 

#This shows a few phyla for which only one feature was observed. 
#Those may be worth filtering, and we’ll check that next.
#First, notice that in this case, six features were annotated with a Phylum of NA. (in the AB table, only 3 were uncharacterized)
#These features are probably artifacts in a dataset like this, and should be removed. 
#The following ensures that features with ambiguous phylum annotation are also removed. 
#Note the flexibility in defining strings that should be considered ambiguous annotation
#Callahan et al.(2016) pg 10 

ps_AB_filt <- subset_taxa(ps_AB, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
table(tax_table(ps_AB_filt)[, "Phylum"], exclude = NULL) # in ps_AB_filt, none are classified as "NA"

#A useful next step is to explore feature prevalence in the dataset, 
#which we will define here as the number of samples in which a taxa appears at least once.

# Compute prevalence of each feature, store as data.frame:
prevdf_AB = apply(X = otu_table(ps_AB_filt),
               MARGIN = ifelse(taxa_are_rows(ps_AB_filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame:
prevdf_AB = data.frame(Prevalence = prevdf_AB,
                    TotalAbundance = taxa_sums(ps_AB_filt),
                    tax_table(ps_AB_filt))

#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf_AB, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#I'm not sure I understand how this is working
#column "1" is the mean of prevalence and column "2" is the sum of prevelance

############Phylum         1     2
#1  Acidobacteriota  3.000000    3
#2 Actinobacteriota 11.400000   57
#3     Bacteroidota  3.000000   12
#4    Cyanobacteria  7.000000    7
#5     Deinococcota  4.000000    4
#6       Firmicutes 12.634578 6431
#7   Proteobacteria  8.666667  286

#so then if I take the sum/mean, I would get how many  times it appeared in 1 sample, maybe.. 
#No, actually, it means that of the 554 taxa, Firmicutes appears 6431 times amoung samples
#In the "prevdf_AB" dataset, there's a column for prevalance, that's how many samples it apperaed in. 
#ugh I don't know what column one is.. 
#I know column two is how many times total that taxa appeared in the samples
#I got it!! 
#so column one is the mean of prevalance, how many samples it appears in/taxa 
#column two is the total number of samples it apperas in across all taxa
#so Bacteriodota is the phylum of 4 taxa, and it appears 12 times (column 2) and then 12/4 = 3 
# based on this, I need to fiter out: 
#Cyanobacteria, it appears in 7 samples, has one classifcation, and it's even specific enough to reach order/genus, plus it's low abundance
#Deinococcota, it's low abundance and low number of samples 
#I think I can leave all the other ones alone for now

filterphyla_AB = c("Cyanobacteria", "Deinococcota") #Define the phyla to filter 
#filter enteries with defined phyla 

ps1_AB = subset_taxa(ps_AB_filt, !Phylum %in% filterphyla_AB)

#Prevelance Filter (Callahan p. 11)
#the next tep is completely "unsupervised" because it relys only on the data in this experiment (rather than the taxanomic ref database)
#this means that this filtering step is divorced from the taxomic classification, so it can be apllied when the annotation may be unrealiable

#First, explore relationship of prevalnce and total read count for each feature 
#this may reveal outliers that need to be removed, and provides insight into the ranges of either feature that might be useful 

# Subset to the remaining phyla (that we didn't filter out in the previous step)
prevdf_AB_1 = subset(prevdf_AB, Phylum %in% get_taxa_unique(ps1_AB, "Phylum"))
ggplot(prevdf_AB_1, aes(TotalAbundance, Prevalence / nsamples(ps_AB_filt),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#I guess the idea of this is to define a prevelance threshold, i.e. 5% of samples
#I'm not sure if I wanna do this right now, but I'll come back to it if I do (Callahan p. 12)

###########PLOTTING TREES#####

AB_tree = plot_tree(ps_AB_1, method = "treeonly")
#it says there is no tree in the object I have provided
# I need to add the tree into the phyloseq object 
AB_tree  #yikes, so this is why we need to agglomerate the data to make a decent looking tree 

# How many genera would be present after filtering?
length(get_taxa_unique(ps_AB_1, taxonomic.rank = "Genus"))
##58... that doesn't seem like too many but appearntly it's a lot 


#Combine features that descend from the same genus 
ps_AB_agg = tax_glom(ps_AB_1, "Genus", NArm = TRUE)
AB_tree1 = plot_tree(ps1_AB, method="sampledodge")
AB_tree1
#The "sampledodge" option results in points drawn next to leaves if individuals from that taxa were observed, and a separate point is drawn for each sample.

#use "taxonomy free" based alternative to combine data features based on specific tree height corresponding to the phylogentic deistance b/w features 
h1 = 0.4
ps2_AB = tip_glom(ps_AB_1, h = h1)
AB_tree2 = plot_tree(ps2_AB, method="treeonly")
AB_tree2 #tree 1 looks much better than this

#okay these trees kinda suck for now. 

table(tax_table(ps_AB_agg)[, "Phylum"], exclude = NULL) #list number of features in each phylum
#mkae a plot of just Firmicutes 
AB_firm_agg= subset_taxa(ps_AB_agg, Phylum =="Firmicutes")
plot_tree(AB_firm_agg, color = "time", shape = "Family", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)
#still way too much information on there 


#FILTERING#####
#is there a way to classify the others as "Other"
#maybe I need to do it in the taxatab 
#I lost track of which AB phyloseq object I should be using... 
#so PS3 AB excluses AB1.3 but it is not filtered yet
# I think I need to run the filtering step on ps3_AB

ps3_AB = subset_samples(ps_AB, sample!="AB1.3") #excluding AB1.3 from future analYses because the read count is too low
rank_names(ps3_AB) #show available ranks in the dataset
ps3_AB =  subset_taxa(ps3_AB, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #removes any that are characterized as NA
table(tax_table(ps3_AB)[, "Phylum"], exclude = NULL) #check to make sure that none are "NA"

# Compute prevalence of each feature, store as data.frame:
prevdf_AB = apply(X = otu_table(ps3_AB),
                  MARGIN = ifelse(taxa_are_rows(ps3_AB), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame:
prevdf_AB = data.frame(Prevalence = prevdf_AB,
                       TotalAbundance = taxa_sums(ps3_AB),
                       tax_table(ps_AB_filt))

#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf_AB, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#I'm not sure I understand how this is working
#column "1" is the mean of prevalence and column "2" is the sum of prevelance
#so column one is the mean of prevalance, how many samples it appears in/taxa 
#column two is the total number of samples it apperas in across all taxa
#so Bacteriodota is the phylum of 4 taxa, and it appears 12 times (column 2) and then 12/4 = 3 
# based on this, I need to fiter out: 
#Cyanobacteria, it appears in 7 samples, has one classifcation, and it's even specific enough to reach order/genus, plus it's low abundance
#Deinococcota, it's low abundance and low number of samples 
#I think I can leave all the other ones alone for now

filterphyla_AB = c("Cyanobacteria", "Deinococcota") #Define the phyla to filter 
#filter enteries with defined phyla 

AB = subset_taxa(ps3_AB, !Phylum %in% filterphyla_AB)  #to make things simplier, I'm going to call this processed object just "AB"

##Filter ps_NL 
table(tax_table(ps_NL)[, "Phylum"], exclude = NULL) #3 are NA
ps_NL =  subset_taxa(ps_NL, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #removes any that are characterized as NA
table(tax_table(ps_NL)[, "Phylum"], exclude = NULL) #check to make sure that none are "NA"
# Compute prevalence of each feature, store as data.frame:
prevdf_NL = apply(X = otu_table(ps_NL),
                  MARGIN = ifelse(taxa_are_rows(ps_NL), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame:
prevdf_NL = data.frame(Prevalence = prevdf_NL,
                       TotalAbundance = taxa_sums(ps_NL),
                       tax_table(ps_NL))

#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf_NL, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#I'm not sure I understand how this is working
#column "1" is the mean of prevalence and column "2" is the sum of prevelance
#so column one is the mean of prevalance, how many samples it appears in/taxa 
#column two is the total number of samples it apperas in across all taxa
#going to filter out Cyanobacteria, Campilobacteria, Crenacrchaeota, Desulfobacterota 


filterphyla_NL = c("Cyanobacteria", "Campilobacterota", "Crenarchaeota", "Desulfobacterota" ) #Define the phyla to filter 
#filter enteries with defined phyla 

NL = subset_taxa(ps_NL, !Phylum %in% filterphyla_NL)
pre_NL  = subset_samples(NL, sample!="NL1" & sample!="NL1.1" & sample!="NL1.2" & sample!="NL1.3"& sample!="NL2"& sample!="NL2.1"& sample!="NL2.2"& sample!="NL2.3"& sample!="NL3"& sample!="NL3.1"& sample!="NL3.2"& sample!="NL3.3")


#I'm going to remove the "NA" phyla from S and Fe, I'm not going to do the other filtering steps right now, but I will come back and do them later
rank_names(ps_Fe) #show available ranks in the dataset
Fe =  subset_taxa(ps_Fe, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #removes any that are characterized as NA
table(tax_table(Fe)[, "Phylum"], exclude = NULL)

rank_names(ps_S) #show available ranks in the dataset
S =  subset_taxa(ps_S, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #removes any that are characterized as NA
table(tax_table(S)[, "Phylum"], exclude = NULL)



#####MERGE PHYLOSEQ OBJECTS####
#Make objects with both of these objects
AB_NL = merge_phyloseq(AB, NL) #Merging AB and NL(all of NL)

AB_NL1 = merge_phyloseq(AB, pre_NL)

AB_NL1_S = merge_phyloseq(AB_NL1, S)#combine S with NL and AB 

AB_NL1_S_Fe = merge_phyloseq(AB_NL1_S, Fe) #combine all experiments into one phyloseq object

all= AB_NL1_S_Fe
experiments= merge_phyloseq(AB, Fe)
experiments=merge_phyloseq(experiments, S)

experiments
######ABUNDANCE PLOTS#### 

plot_bar(ps_AB_agg, fill="Genus")#this is the one with the agglomerated taxa, looks decent
plot_bar(ps1_AB, fill="Genus")#okay yeah this looks much worse 

ps2_AB = transform_sample_counts(ps1_AB, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
ps4_AB = transform_sample_counts(ps3_AB, function(x) {x/sum(x)}) #this one excludes sample AB1.3

#I need to somehow filter taxa to make a better abundance plot
top30_AB = names(sort(taxa_sums(ps4_AB), decreasing=TRUE))[1:30] #sort by the top30, the "decreasing=TRUE" part is important, otherwise the top20 are just arbitrary in the list, this sorts it in decending order


#make a plot of top 30, seperate in a matrix by time point
ps4_AB.top30 <- prune_taxa(top30_AB, ps4_AB) 
plot_bar(ps4_AB.top30, x="Bottle.Number", fill="Genus") + facet_wrap(~time, scales="free_x")

library(ggplot2)
#add ggplot layer to remove the OTU seperation lines
p = plot_bar(ps4_AB.top30, x="Bottle.Number", fill="Genus") #top 30 
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x")


#top100
top100_AB = names(sort(taxa_sums(ps4_AB), decreasing=TRUE))[1:100]
ps4_AB.top100 <- prune_taxa(top100_AB, ps4_AB)
p = plot_bar(ps4_AB.top100, x="Bottle.Number", fill="Genus")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x")

#all
p = plot_bar(ps4_AB, x="Bottle.Number", fill="Genus")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x")
#too crowded 

#top150
top150_AB = names(sort(taxa_sums(ps4_AB), decreasing=TRUE))[1:150]
ps4_AB.top150 <- prune_taxa(top150_AB, ps4_AB)
p = plot_bar(ps4_AB.top150, x="Bottle.Number", fill="Genus")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#AB + NL, total abundance 
top150_ABNL = names(sort(taxa_sums(AB_NL), decreasing=TRUE))[1:150]
top150_ABNL <- prune_taxa(top150_ABNL, AB_NL)
p = plot_bar(top150_ABNL, x="Bottle", fill="Genus",title="Total Abundance AB and Control, top 150")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#total abundance top200
top200_ABNL = names(sort(taxa_sums(AB_NL), decreasing=TRUE))[1:200]
top200_ABNL <- prune_taxa(top200_ABNL, AB_NL)
p = plot_bar(top200_ABNL, x="Bottle", fill="Genus",title="Total Abundance AB and Control, top 200")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#total abundance, all 
p = plot_bar(AB_NL, x="Bottle", fill="Genus",title="Total Abundance AB and Control")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#AB + NL, relative abundance
AB_NL_rel = transform_sample_counts(AB_NL, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL = names(sort(taxa_sums(AB_NL_rel), decreasing=TRUE))[1:200]
top200_ABNL <- prune_taxa(top200_ABNL, AB_NL_rel)
p = plot_bar(top200_ABNL, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, top 200")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 



###I'm going to try subsetting top 200 and then merging the phyloseq objects. 

top200_AB = names(sort(taxa_sums(AB), decreasing=TRUE))[1:200]
top200_AB <- prune_taxa(top200_AB, AB)
top200_NL = names(sort(taxa_sums(NL), decreasing=TRUE))[1:200]
top200_NL <- prune_taxa(top200_NL, NL)

top200_ABNL_2 = merge_phyloseq(top200_AB, top200_NL)
p = plot_bar(top200_ABNL_2, x="Bottle", fill="Genus",title="Total Abundance AB and Control, top 200")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#It doesn't look like it makes a huge difference if I cut off the low abudance ones before or after merging. 


#Trying to see if I can somehow rename the ones that dont appear in the top 200 as "other" 

taxa_sums(AB_NL_rel)
sample_sums(AB_NL_rel)

AB_NL_new = data.frame(RelativeAbundance = taxa_sums(AB_NL_rel),
                       tax_table(AB_NL))

AB_new = data.frame(TotalAbundance= taxa_sums(AB), 
  RelativeAbundance = taxa_sums(AB_rel), tax_table(AB))
#not quite 


#I want to make a relative abundance plot with the agglomerated taxa
AB_agg_rel = transform_sample_counts(AB_agg, function(x) {x/sum(x)}) #transform AB_agg to relative abundance
AB_agg_rel_plot = plot_bar(AB_agg_rel, x="Bottle", fill="Genus", title = "Relative Abundance Antibiotics Experiment") #plot AB_agg_rel
AB_agg_rel_plot + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x")  #add ggplot layer to remove the OTU seperation lines




###### NMDS PLOTS#####
library(colorspace) #gives more options for colors

#AB
ps.prop.AB <- transform_sample_counts(AB, function(otu) otu/sum(otu))
ord.nmds.bray.AB <- ordinate(ps.prop.AB, method="NMDS", distance="bray")
AB_NMDS = plot_ordination(ps.prop.AB, ord.nmds.bray.AB, color="time", title="Antibiotics NMDS")
AB_NMDS

#NL
ps.prop_NL <- transform_sample_counts(NL, function(otu) otu/sum(otu))
ord.nmds.bray_NL <- ordinate(ps.prop_NL, method="NMDS", distance="bray")
NL_NMDS = plot_ordination(ps.prop_NL, ord.nmds.bray_NL, color="time", title="Control Experiment Bray NMDS")
NL_NMDS

AB_NL = merge_phyloseq(AB, NL) #Merging AB and NL(all of NL)

#Control V. Antibiotics experiment NMDS 
ps.prop_ABNL <- transform_sample_counts(AB_NL, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL <- ordinate(ps.prop_ABNL, method="NMDS", distance="bray")
ABNL_NMDS = plot_ordination(ps.prop_ABNL, ord.nmds.bray_ABNL, color="time", shape="Spike", title="Control v. Antibiotics Experiment Bray NMDS")
ABNL_NMDS 

#1
#NL(without NL 1-3)
new_ps.prop_NL <- transform_sample_counts(pre_NL, function(otu) otu/sum(otu))
new_ord.nmds.bray_NL <- ordinate(new_ps.prop_NL, method="NMDS", distance="bray")
new_NL_NMDS = plot_ordination(new_ps.prop_NL, new_ord.nmds.bray_NL, color="Spike", shape="time" )
NLcol= "gray50"
new_NL_NMDS + scale_color_manual(values= NLcol) + geom_point(size=3.5)

AB_NL1 = merge_phyloseq(AB, pre_NL) #merge NL (only 1-4) + AB

#2
#Control V. Antibiotics experiment NMDS, but without the extra 3 bottles of NL1-3 
#THIS IS ONE THAT I WILL USE
ps.prop_ABNL1 <- transform_sample_counts(AB_NL1, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1 <- ordinate(ps.prop_ABNL1, method="NMDS", distance="bray")
ABNLcols= c("blue", "gray50")
ABNL1_NMDS = plot_ordination(ps.prop_ABNL1, ord.nmds.bray_ABNL1, color="Spike", shape="time", title="Control v. Antibiotics Experiment Bray NMDS") 
ABNL1_NMDS + scale_color_manual(values= ABNLcols) + geom_point(size = 3.5) 

#combine S with NL and AB 
AB_NL1_S = merge_phyloseq(AB_NL1, ps_S)

#3 
#make an NMDS plot with S, NL, and AB 
ps.prop_ABNL1S <- transform_sample_counts(AB_NL1_S, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1S <- ordinate(ps.prop_ABNL1S, method="NMDS", distance="bray")
ABNLScols= c("blue", "gray50", "gold")
ABNL1S_NMDS = plot_ordination(ps.prop_ABNL1S, ord.nmds.bray_ABNL1S, color="Spike", shape="time", title="Bray NMDS") 
ABNL1S_NMDS + scale_color_manual(values= ABNLScols) + geom_point(size = 3.5) 

#S (this is still unfilerted) "ps_S" is the unfiltered version, once I filter it, I will name the phyloseq obect "S"
ps.prop_S <- transform_sample_counts(ps_S, function(otu) otu/sum(otu))
ord.nmds.bray_S <- ordinate(ps.prop_S, method="NMDS", distance="bray")
S_NMDS = plot_ordination(ps.prop_S, ord.nmds.bray_S, color="time", title="Sulfur Experiment Bray NMDS")
cols_6 = divergingx_hcl(6, palette = "Zissou 1")
S_NMDS + scale_colour_manual(values=cols_6) 


#All 3 experiments plus control 
AB_NL1_S_Fe = merge_phyloseq(AB_NL1_S, ps_Fe)

#4
#make an NMDS plot
ps.prop_ABNL1SFe <- transform_sample_counts(AB_NL1_S_Fe, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1SFe <- ordinate(ps.prop_ABNL1SFe, method="NMDS", distance="bray")
ABNLSFecols= c("blue", "gray50", "gold", "firebrick1", "forestgreen")
ABNL1SFe_NMDS = plot_ordination(ps.prop_ABNL1SFe, ord.nmds.bray_ABNL1SFe, color="Spike", shape="time") 
ABNL1SFe_NMDS + scale_color_manual(values= ABNLSFecols) + geom_point(size = 3.5) 


#Fe,the unfiltered version, once I filter it, I wil name the object "Fe"
ps.prop_Fe<- transform_sample_counts(ps_Fe, function(otu) otu/sum(otu))
ord.nmds.bray_Fe<- ordinate(ps.prop_Fe, method="NMDS", distance="bray")
Fe_NMDS = plot_ordination(ps.prop_Fe, ord.nmds.bray_Fe, color="time", title="Iron Experiment Bray NMDS")
cols_6 = divergingx_hcl(6, palette = "Zissou 1")
Fe_NMDS + scale_colour_manual(values=cols_6) 


#####this creates a split plot of taxa and samples.. althought I'm not sure exactly what the taxa plot represents
NL_NMDS = plot_ordination(new_ps.prop_NL, type="split", new_ord.nmds.bray_NL, color="time") 
NL_NMDS
#try it again but with the ABNL1 plot
ABNL1_NMDS_split = plot_ordination(ps.prop_ABNL1, ord.nmds.bray_ABNL1, type="split", color="Spike", title="Control v. Antibiotics Experiment Bray NMDS") 
ABNL1_NMDS_split



###ALPHA DIVERSITY#### 

ABNLSFecols= c("blue", "gray50", "gold", "firebrick1", "forestgreen")

cols_6 = divergingx_hcl(6, palette = "Zissou 1")
colors= c("red", "orange", "yellow", "green", "blue", "purple")

NL_alpha = prune_species(speciesSums(pre_NL)>0, pre_NL)
NL_bottles = merge_samples(NL_alpha, "time")
plot_richness(NL_bottles, measures="Chao1")+ ylim(450,700)


AB_alpha = prune_species(speciesSums(AB)>0, AB)
AB_bottles = merge_samples(AB_alpha, "time") #this merges AB1, AB2, and AB3 into one at each time point 
AB_alpha = plot_richness(AB_bottles, measures="Chao1") + ylim(450,700)
AB_alpha

S_alpha = prune_species(speciesSums(ps_S)>0, ps_S)
S_bottles = merge_samples(S_alpha, "time")
plot_richness(S_bottles, measures="Chao1") + ylim(450,700)

Fe_alpha = prune_species(speciesSums(ps_Fe)>0, ps_Fe)
Fe_bottles = merge_samples(Fe_alpha, "time")
plot_richness(Fe_bottles, measures="Chao1") + ylim(450,700)

#### trying to make a diversity plot with all of them on the same plot ### 

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

 
all_alpha = prune_species(speciesSums(AB_NL1_S_Fe)>0, AB_NL1_S_Fe) #it's normal to get an error message
plot_richness(all_alpha, measures="Chao1")


####### AGGlOMERATING TAXA #####
AB_agg = tax_glom(AB, "Genus", NArm = TRUE)


AB_NL1 = merge_phyloseq(AB, pre_NL) #this is merging AB and NL but without NL 1-3

###### NMDS PLOTS TO USE IN GSA PRESENTATION #########

#1
#NL(without NL 1-3)
new_ps.prop_NL <- transform_sample_counts(pre_NL, function(otu) otu/sum(otu))
new_ord.nmds.bray_NL <- ordinate(new_ps.prop_NL, method="NMDS", distance="bray")
new_NL_NMDS = plot_ordination(new_ps.prop_NL, new_ord.nmds.bray_NL, color="Spike", shape="time", title = "NMDS plot" )
NLcol= "gray50"
new_NL_NMDS + scale_color_manual(values= NLcol) + geom_point(size=3.5)  + ylim(-1,1.5) + xlim(-2,2)

#2
#Control V. Antibiotics experiment NMDS, but without the extra 3 bottles of NL1-3 
#THIS IS ONE THAT I WILL USE
ps.prop_ABNL1 <- transform_sample_counts(AB_NL1, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1 <- ordinate(ps.prop_ABNL1, method="NMDS", distance="bray")
ABNLcols= c("blue", "gray50")
ABNL1_NMDS = plot_ordination(ps.prop_ABNL1, ord.nmds.bray_ABNL1, color="Spike", shape="time", title="NMDS plot") 
ABNL1_NMDS + scale_color_manual(values= ABNLcols) + geom_point(size = 3.5) + ylim(-1,1.5) + xlim(-2,2)


#3 
#make an NMDS plot with S, NL, and AB 
ps.prop_ABNL1S <- transform_sample_counts(AB_NL1_S, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1S <- ordinate(ps.prop_ABNL1S, method="NMDS", distance="bray")
ABNLScols= c("blue", "gray50", "gold")
ABNL1S_NMDS = plot_ordination(ps.prop_ABNL1S, ord.nmds.bray_ABNL1S, color="Spike", shape="time", title="NMDS plot") 
ABNL1S_NMDS + scale_color_manual(values= ABNLScols) + geom_point(size = 3.5) + ylim(-1,1.5) + xlim(-2,2)

#4
#make an NMDS plot
ps.prop_ABNL1SFe <- transform_sample_counts(AB_NL1_S_Fe, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1SFe <- ordinate(ps.prop_ABNL1SFe, method="NMDS", distance="bray")
ABNLSFecols= c("blue", "gray50", "gold", "firebrick1", "forestgreen")
ABNL1SFe_NMDS = plot_ordination(ps.prop_ABNL1SFe, ord.nmds.bray_ABNL1SFe, color="Spike", shape="time", title="NMDS plot") 
ABNL1SFe_NMDS + scale_color_manual(values= ABNLSFecols) + geom_point(size = 3.5) + ylim(-1,1.5) + xlim(-2,2)

######ABUNDANCE PLOTS FOR GSA#####


#AB + NL, relative abundance, top 200 
AB_NL1_rel = transform_sample_counts(AB_NL1, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL1 = names(sort(taxa_sums(AB_NL1_rel), decreasing=TRUE))[1:200]
top200_ABNL1 <- prune_taxa(top200_ABNL1, AB_NL1_rel)
AB_NL1_rel_plot = plot_bar(top200_ABNL1, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot+ geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip()



#another one 
AB_NL1_rel_plot1 = plot_bar(top200_ABNL1, x="time", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot1 + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~Bottle, scales="free_x") + coord_flip()+ scale_x_discrete(limits = rev(levels("time")))

#another one 
AB_NL1_rel_plot2 = plot_bar(top200_ABNL1, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot2 + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

#another one, this time at the family level
AB_NL1_rel_plot3 = plot_bar(top200_ABNL1, x="Bottle", fill="Family",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot3 + geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+
  facet_wrap(~time, scales="free_x") + labs(y= "Relative Abundance", x = " ") +
  theme(axis.text.x = element_text(angle = 360), axis.ticks =element_blank())


#Last one, back to genus
AB_NL1_rel_plot4 = plot_bar(top200_ABNL1, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot4 + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+
  facet_wrap(~time, scales="free_x") + labs(y= "Relative Abundance", x = " ") +
  theme(axis.text.x = element_text(angle = 360), axis.ticks =element_blank())






#####Figures to copy#####
#AB + NL, relative abundance
AB_NL_rel = transform_sample_counts(AB_NL, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL = names(sort(taxa_sums(AB_NL_rel), decreasing=TRUE))[1:200]
top200_ABNL <- prune_taxa(top200_ABNL, AB_NL_rel)
p = plot_bar(top200_ABNL, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, top 200")
p + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") 

AB_NL1_rel = transform_sample_counts(AB_NL1, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL1 = names(sort(taxa_sums(AB_NL1_rel), decreasing=TRUE))[1:200]
top200_ABNL1 <- prune_taxa(top200_ABNL1, AB_NL1_rel)
AB_NL1_rel_plot = plot_bar(top200_ABNL1, x="Bottle", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot+ geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip()
#####Plots for AGU2020#####

#NMDS plot with all experiments
ps.prop_ABNL1SFe <- transform_sample_counts(AB_NL1_S_Fe, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1SFe <- ordinate(ps.prop_ABNL1SFe, method="NMDS", distance="bray")
ABNLSFecols= c("blue", "gray50", "gold", "firebrick1", "forestgreen")
ABNL1SFe_NMDS = plot_ordination(ps.prop_ABNL1SFe, ord.nmds.bray_ABNL1SFe, color="Spike", shape="time", title="NMDS plot") 
ABNL1SFe_NMDS + scale_color_manual(values= ABNLSFecols) + geom_point(size = 3.5) + ylim(-1,1.5) + xlim(-2,2)


#Transform into relative counts 
NL_rel = transform_sample_counts(NL, function(x) {x/sum(x)}) 
AB_rel = transform_sample_counts(AB, function(x) {x/sum(x)})
S_rel = transform_sample_counts(S, function(x) {x/sum(x)}) 
Fe_rel = transform_sample_counts(Fe, function(x) {x/sum(x)}) 
AB_NL1_S_Fe_rel = transform_sample_counts(AB_NL1_S_Fe, function(x) {x/sum(x)}) 
all_rel = transform_sample_counts(AB_NL1_S_Fe, function(x) {x/sum(x)}) # this is the same as the previous one, just a shorter name
exp_rel = transform_sample_counts(experiments, function(x) {x/sum(x)})

#top 150 taxa of NL, family 
#+ coord_flip() #this flips the x and y axes 
top150_NL = names(sort(taxa_sums(NL_rel), decreasing=TRUE))[1:150]
top150_NL <- prune_taxa(top150_NL, NL_rel)
p1 = plot_bar(top150_NL, x="Bottle", fill="Family",title="Relative Abundance Control, top 150")
p1 + geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip() 
#I like this one. 

#top 150 taxa of AB, family 
top150_AB = names(sort(taxa_sums(AB_rel), decreasing=TRUE))[1:150]
top150_AB <- prune_taxa(top150_AB, AB_rel)
p2 = plot_bar(top150_AB, x="Bottle", fill="Family",title="Relative Abundance Antibiotics Spike, top 150")
p2 + geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip()
# I don't think it makes sense to split up the figures anymore, the colors will change in this figure unless I manually assignn them 
# i don't want to manually assign colors RN 


####Top 100, all experiments######
#all experiments, top 100, order#####
AB_NL1_S_Fe
top100_all = names(sort(taxa_sums(all_rel), decreasing = TRUE)) [1:100]
top100_all = prune_taxa(top100_all, all_rel)
p3 = plot_bar(top100_all, x="Bottle", fill="Order",title="Relative Abundance")
p3 + geom_bar(aes(color=Order, fill=Order), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip()
#I want to reoder samples so that NL is on top not S 


#experiments only, top 100, order
top100_exp = names(sort(taxa_sums(exp_rel), decreasing = TRUE)) [1:100]
top100_exp = prune_taxa(top100_exp, exp_rel)
p4 = plot_bar(top100_exp, x="Bottle", fill="Order",title="Relative Abundance")
p4 + geom_bar(aes(color=Order, fill=Order), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip()

#control only, top 100, order
top100_NL = names(sort(taxa_sums(NL_rel), decreasing=TRUE))[1:100]
top100_NL <- prune_taxa(top100_NL, NL_rel)
p1 = plot_bar(top100_NL, x="Bottle", fill="Order",title="Relative Abundance Control, top 100")
p1 + geom_bar(aes(color=Order, fill=Order), position="stack", stat="identity")+ facet_wrap(~time, scales="free_x") + coord_flip() 
#I like this one. 







#####Top 150, all experiments; by Spike ######

#Genus
top150_all = names(sort(taxa_sums(all_rel), decreasing=TRUE))[1:150]
top150_all <- prune_taxa(top150_all, all_rel)
p5 = plot_bar(top150_all, x="time", fill="Genus",title="Relative Abundance, top 150")
p5 + geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~Spike, scales="free_x") + coord_flip() 

#Family
top150_all = names(sort(taxa_sums(all_rel), decreasing=TRUE))[1:150]
top150_all <- prune_taxa(top150_all, all_rel)
p5 = plot_bar(top150_all, x="time", fill="Family",title="Relative Abundance, top 150")
p5 + geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+ facet_wrap(~Spike) + coord_flip() 


#Figure for EREF report######

top200_all = names(sort(taxa_sums(all_rel), decreasing=TRUE))[1:200]
top200_all <- prune_taxa(top200_all, all_rel)
p = plot_bar(top200_all, x="time", fill="Family",title="Relative Abundance, Top 200 Taxa")
p + geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+ facet_wrap(~Bottle, scales="free_x") 

AB_NL1_rel = transform_sample_counts(AB_NL1, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL1 = names(sort(taxa_sums(AB_NL1_rel), decreasing=TRUE))[1:200]
top200_ABNL1 <- prune_taxa(top200_ABNL1, AB_NL1_rel)
AB_NL1_rel_plot = plot_bar(top200_ABNL1, x="time", fill="Genus",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot+ geom_bar(aes(color=Genus, fill=Genus), position="stack", stat="identity")+ facet_wrap(~Bottle, scales="free_x") 

AB_NL1_rel = transform_sample_counts(AB_NL1, function(x) {x/sum(x)}) #transforms sample counts into relative sample counts
top200_ABNL1 = names(sort(taxa_sums(AB_NL1_rel), decreasing=TRUE))[1:200]
top200_ABNL1 <- prune_taxa(top200_ABNL1, AB_NL1_rel)
AB_NL1_rel_plot = plot_bar(top200_ABNL1, x="time", fill="Family",title="Relative Abundance AB and Control, Top 200 Taxa")
AB_NL1_rel_plot+ geom_bar(aes(color=Family, fill=Family), position="stack", stat="identity")+ facet_wrap(~Bottle, scales="free_x") 

ps.prop_ABNL1SFe <- transform_sample_counts(AB_NL1_S_Fe, function(otu) otu/sum(otu))
ord.nmds.bray_ABNL1SFe <- ordinate(ps.prop_ABNL1SFe, method="NMDS", distance="bray")
ABNL1SFe_NMDS = plot_ordination(ps.prop_ABNL1SFe, ord.nmds.bray_ABNL1SFe, color="Spike", shape="time", title="NMDS plot") 
ABNL1SFe_NMDS + ylim(-1,1.5) + xlim(-2,2)














#unused commands from "landfill paper script. Rmd"######


I didn't end up using these commands:
ABK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files ABK")
FeK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged FeK")
NK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged NK")
SK_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/PEAR merged files SK")


 fns_NK = sort(list.files(NK_path, full.names = TRUE))
fns_FeK = sort(list.files(FeK_path, full.names = TRUE)) 
fns_SK = sort(list.files(SK_path, full.names = TRUE))
fns_ABK = sort(list.files(ABK_path, full.names = TRUE))


ABK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/ABK_filt")
FK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/FK_filt")
SK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/SK_filt")
NK_filt_path = file.path("~/Desktop/Landfill sequencing data processing /PEAR_merged/NK_filt")


if(!file_test("-d", NK_path)) dir.create(NK_filt_path)
filt_NK <- file.path(NK_filt_path, basename(fns_NK))
for(i in seq_along(fns_NK)) {filterAndTrim(fns_NK, NK_filt_path, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}



if(!file_test("-d", ABK_path)) dir.create(ABK_filt_path)
filt_AB <- file.path(killed_filt_path, basename(fns_AB))
for(i in seq_along(fns_AB)) {filterAndTrim(fns_AB, AB_filt_path, 
                rev = NULL, 
                filt.rev =NULL, 
                trimLeft=c(10), truncLen=c(290),
                maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)}


NKi <- sample(length(fns_NK), 33) 
for(i in NKi) { print(plotQualityProfile(fns_NK[i]) + ggtitle("NK paired")) }

SKi <- sample(length(fns_SK), 21) 
for(i in SKi) { print(plotQualityProfile(fns_SK[i]) + ggtitle("SK paired")) }

FeKi <- sample(length(fns_FeK), 21) 
for(i in FeKi) { print(plotQualityProfile(fns_FeK[i]) + ggtitle("FeK paired")) }

ABKi <- sample(length(fns_ABK), 21) 
for(i in ABKi) { print(plotQualityProfile(fns_ABK[i]) + ggtitle("ABK paired")) }





















#space to try things

