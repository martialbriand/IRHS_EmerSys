######################################################################
#Chesneau et al., 2020, "Temporal dynamics of bacterial assemblages during seed development and maturation". 
#Data collected in 2016-2017
#Authors: Guillaume Chesneau, IRHS, guillaume.chesneau <at> inrae.fr
######################################################################


# Before you start
# This script has been made on R 3.6.1 (and Rstudio 1.2.5019)
R.version
RStudio.Version()$version
#Source package richness.R and graphical_methods.R are needed 
# The following packages (and their dependencies) are needed to run the whole analysis
# ape
# dada2
# decipher
# dendextend
# dplyr
# ggplot2
# gridExtra
# phangorn
# pheatMap
# phyloseq
# picante
# RColorBrewer
# UpsetR
# vegan
# viridis

####

#Set working direction
setwd("~/These2/Script R/F2S-diversite/F2S-year-1/Script final_200405/")

# A] SEQUENCES to ASV----

#....1)cutadapt v1.15(Remove primers)----
#for i in `cat group`; do cutadapt --discard-untrimmed -o $i.gyrB.R1.fq -p $i.gyrB.R2.fq -g MGNCCNGSNATGTAYATHGG -G ACNCCRTGNARDCCDCCNGA -e 0.1  -O 20 $i*L001_R1_001.fastq.gz $i*_L001_R2_001.fastq.gz; done
#for i in `cat group`; do cutadapt --discard-untrimmed -o $i.16S.R1.fq -p $i.16S.R2.fq -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -e 0 -O 19 $i*_L001_R1_001.fastq.gz $i*_L001_R2_001.fastq.gz; done

#....2)Dada2 (fq -> ASV)----
library(dada2); packageVersion("dada2") #'1.14.0'

#Adapt from https://benjjneb.github.io/dada2/tutorial.html
#Start with 16S rRNA gene (v4)
path <- "~/These2/Script R/F2S-diversite/F2S-year-1/manual4/fastq/16S/" 
list.files(path)
fnFs <- sort(list.files(path, pattern="16S.R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="16S.R2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:8]) #select 200
plotQualityProfile(fnRs[1:8]) #select 150
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150),
                     maxN=0, maxEE=c(1,1), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # 
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #without pooling or pseudo-pooling (no need to detect rare ASV)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
dim(seqtab2)

#Detect/Remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

#Summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#IDtaxa taxonomic affiliation
library(DECIPHER); packageVersion("DECIPHER") # '2.14.0'
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/Documents/DB/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa2 <- taxid
taxa.print2 <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print2) <- NULL
head(taxa.print2)
saveRDS(seqtab.nochim, "~/Documents/Papers/Chesneau_2019/analyses_191118/16S_ASV.rds")
saveRDS(taxa2, "~/Documents/Papers/Chesneau_2019/analyses_191118/16S_taxo.rds")

#Continue with gyrB
path <- "~/These2/Script R/F2S-diversite/F2S-year-1/manual4/fastq/gyrB/" 
list.files(path)
fnFs <- sort(list.files(path, pattern="gyrB.R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="gyrB.R2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:8]) #select 200
plotQualityProfile(fnRs[1:8]) #select 150
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150),
                     maxN=0, maxEE=c(1,1), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # 
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #without pooling or pseudo-pooling (no need to detect rare ASV)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Since gyrB is a protein-coding genes only triplets should be conserved (244-247-250-253-256-259-262-265-268)
seqtab244 <- seqtab[,nchar(colnames(seqtab)) %in% 244]
seqtab247 <- seqtab[,nchar(colnames(seqtab)) %in% 247]
seqtab250 <- seqtab[,nchar(colnames(seqtab)) %in% 250]
seqtab253 <- seqtab[,nchar(colnames(seqtab)) %in% 253]
seqtab256 <- seqtab[,nchar(colnames(seqtab)) %in% 256]
seqtab259 <- seqtab[,nchar(colnames(seqtab)) %in% 259]
seqtab262 <- seqtab[,nchar(colnames(seqtab)) %in% 262]
seqtab265 <- seqtab[,nchar(colnames(seqtab)) %in% 265]
seqtab268 <- seqtab[,nchar(colnames(seqtab)) %in% 268]

#Merge all files
seq.final <- cbind(seqtab244, seqtab247, seqtab250, seqtab253, seqtab256, seqtab259, seqtab262, seqtab265, seqtab268) 
dim(seq.final)
sum(seq.final)/sum(seqtab)

#Detect/Remove chimera
seqtab.nochim <- removeBimeraDenovo(seq.final, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seq.final)

#Summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
saveRDS(seqtab.nochim, "~/Documents/Papers/Chesneau_2019/analyses_191118/gyrB_ASV.rds")

#RDP taxonomic affiliation
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/DB/train_set_gyrB_v4.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
saveRDS(taxa, "~/Documents/Papers/Chesneau_2019/analyses_191118/gyrB_taxo.rds")

# B] PRE-PROCESSING reads----

#Set working direction
setwd("~/These2/Script R/F2S-diversite/F2S-year-1/Script final_200405/")

## Load required libraries
library("data.table"); packageVersion("data.table") #'1.12.8'
library("DECIPHER"); packageVersion("DECIPHER") # '2.14.0'
library("dplyr"); packageVersion("dplyr") #'0.8.3'
library("ggplot2"); packageVersion("ggplot2") # '3.2.1'
library("gridExtra"); packageVersion("gridExtra") # '2.3'
library("phangorn"); packageVersion("phangorn") # "2.5.5'
library("phyloseq"); packageVersion("phyloseq") # '1.30.0'
library("UpSetR"); packageVersion("UpsetR") #'1.4.0'
library("vegan"); packageVersion("vegan") # '2.5.6'

#....1) Create phyloseq objects----
#Read files
asvgyrB <- readRDS("gyrB_ASV.rds")
taxgyrB <- readRDS("gyrB_taxo.rds")
asv16S <- readRDS("16S_ASV.rds")
tax16S <- readRDS("16S_taxo.rds")
tax16S[5, 7] = "pseudomonas"
design <- read.csv("summary_design_file.csv", sep = ";", row.names = 1, check.names=FALSE)

#Create Phyloseq Objects
psgyrB <- phyloseq(tax_table(taxgyrB), sample_data(design),
               otu_table(asvgyrB, taxa_are_rows = FALSE)) #1444 taxa, 56 samples
ps16S <- phyloseq(tax_table(tax16S), sample_data(design),
                  otu_table(asv16S, taxa_are_rows = FALSE)) # 260 taxa, 56 samples

#Remove all ASVs not affiliated at the kingdom level and parE
psgyrB0 <- subset_taxa(psgyrB, !is.na(Kingdom) & !Kingdom %in% c("parE") & !is.na(Phylum))
ps16S0 <- subset_taxa(ps16S, !is.na(domain) & !is.na(phylum) & !order %in% c("Chloroplast") & !family %in% c("Mitochondria"))

#Remove mock and nc
psgyrB1 <- subset_samples(psgyrB0, sample_names(psgyrB0)  != "mock" & sample_names(psgyrB0)  != "nc")
ps16S1 <- subset_samples(ps16S0, sample_names(ps16S0)  != "mock" & sample_names(ps16S0)  != "nc")

#Rename ASV
dna.gyrB <- Biostrings::DNAStringSet(taxa_names(psgyrB1))
names(dna.gyrB) <- taxa_names(psgyrB1)
psgyrB2 <- merge_phyloseq(psgyrB1, dna.gyrB)
taxa_names(psgyrB2) <- paste0("ASV", seq(ntaxa(psgyrB2)))
dna.16S <- Biostrings::DNAStringSet(taxa_names(ps16S1))
names(dna.16S) <- taxa_names(ps16S1)
ps16S2 <- merge_phyloseq(ps16S1, dna.16S)
taxa_names(ps16S2) <- paste0("ASV", seq(ntaxa(ps16S2)))

#Phylogenetic tree
seqs.gyrB <- refseq(psgyrB2)
alignment.gyrB <- AlignSeqs(DNAStringSet(seqs.gyrB), anchor=NA) 
phang.align.gyrB <- phyDat(as(alignment.gyrB, "matrix"), type="DNA")
dm.gyrB <- dist.ml(phang.align.gyrB)
treeNJ.gyrB <- NJ(dm.gyrB)# Note, tip order != sequence order 
saveRDS(treeNJ.gyrB, "gyrB_NJ.rds")
seqs.16S <- refseq(ps16S2)
alignment.16S <- AlignSeqs(DNAStringSet(seqs.16S), anchor=NA) 
phang.align.16S <- phyDat(as(alignment.16S, "matrix"), type="DNA")
dm.16S <- dist.ml(phang.align.16S)
treeNJ.16S <- NJ(dm.16S)# Note, tip order != sequence order 
saveRDS(treeNJ.16S, "16S_NJ.rds")

treeNJ.gyrB <- readRDS("gyrB_NJ.rds")
dgyrB <- merge_phyloseq(psgyrB2, treeNJ.gyrB)
d16S <- merge_phyloseq(ps16S2, treeNJ.16S)
saveRDS(dgyrB, "gyrB_PS.rds")
saveRDS(d16S, "16S_PS.rds")

# ....2) Representation of phylum abundance/prevalence----
#Read phyloseq object
dgyrB <- readRDS("gyrB_PS.rds") #1155 taxa, 54 samples
d16S <- readRDS("16S_PS.rds") #232 taxa, 54 samples

#Compute prevalence of each feature, store as data.frame
prevdf16S <- apply(X = otu_table(d16S),
                   MARGIN = ifelse(taxa_are_rows(d16S), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
prevdfgyrB <- apply(X = otu_table(dgyrB),
                    MARGIN = ifelse(taxa_are_rows(dgyrB), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf16S <- data.frame(Prevalence = prevdf16S,
                        TotalAbundance = taxa_sums(d16S),
                        tax_table(d16S))
prevdfgyrB <- data.frame(Prevalence = prevdfgyrB,
                         TotalAbundance = taxa_sums(dgyrB),
                         tax_table(dgyrB))

# ....3) Graphical representation of phylum abundance/prevalence----
theme_set(theme_bw())
phyla_16S <-ggplot(prevdf16S, aes(TotalAbundance, Prevalence / nsamples(d16S),color=phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + labs(title="A") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum, ncol=5) + 
  scale_color_manual(values=c("#FA8258", "#FAAC58", "#8F7434", "#A93B35", "#1D8D47", "#A1A4A2", "#CF49B8")) +
  scale_fill_manual(values=c("#FA8258", "#FAAC58", "#8F7434", "#A93B35", "#1D8D47", "#A1A4A2", "#CF49B8")) +
  theme(legend.position="none",
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14))

phyla_gyrB <- ggplot(prevdfgyrB, aes(TotalAbundance, Prevalence / nsamples(dgyrB),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  labs(title="B") + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum, ncol=5) + 
  scale_color_manual(values=c("#FA8258", "#FAAC58", "#8F7434", "#3E57C6", "#1D8D47", "#4FCAC3", "#CFDB3A", "#9DE0A0", "#1D8D47", "#CF49B8", "#414246")) +
  scale_fill_manual(values=c("#FA8258", "#FAAC58", "#8F7434", "#3E57C6","#1D8D47", "#4FCAC3", "#CFDB3A", "#9DE0A0", "#1D8D47", "#CF49B8", "#414246")) +
  theme(legend.position="none",
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14))


# ....4) Assess sequencing depth----
sdt16S <- data.table(as(sample_data(d16S), "data.frame"),
                     TotalReads = sample_sums(d16S), keep.rownames = TRUE)
sdtgyrB <- data.table(as(sample_data(dgyrB), "data.frame"),
                      TotalReads = sample_sums(dgyrB), keep.rownames = TRUE)
sdt16S$logreads <-log(sdt16S$TotalReads, 10)
sdtgyrB$logreads <-log(sdtgyrB$TotalReads, 10)

# ....5) Plot bacterial CFU and number of bacterial reads (16S rRNA gene and gyrB)----
#Reorder legend
sdt16S$PlantStructure <- factor(sdt16S$PlantStructure, 
                                 levels = c("Flower_Bud","Open_Flower","Fruits","Seeds"))
sdtgyrB$PlantStructure <- factor(sdtgyrB$PlantStructure, 
                                 levels = c("Flower_Bud","Open_Flower","Fruits","Seeds"))

#Plot CFU
plotCFU <- ggplot(data=sdtgyrB, aes_string(x='TimePoint_days',y='LogCFU_per_g', color = 'PlantStructure')) + 
  stat_summary(geom = "errorbar", fun.data = mean_se) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  theme(legend.position = "right") + 
  theme(legend.title = element_text(face = "bold", size = 9)) + 
  ylab("Log CFU per gram") + 
  theme(axis.title.y = element_text(face="bold", size = 10)) + 
  theme(axis.title.x = element_blank()) + 
  facet_grid (~Plant, scales = "free") +
  ggtitle("A")

#Plot read 16S
plotread16S <- ggplot(data=sdt16S, aes_string(x='TimePoint_days',y='logreads', color = 'PlantStructure')) + 
  stat_summary(geom = "errorbar", fun.data = mean_se) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  theme(legend.position = "right") + 
  theme(legend.title = element_text(face = "bold", size = 9)) + 
  ylab("Log reads") + 
  theme(axis.title.y = element_text(face="bold", size = 10)) + 
  theme(axis.title.x = element_blank()) + 
  facet_grid (~Plant, scales = "free") +
  ylim(1,6) +
  ggtitle("B")

#Plot read gyrB
plotreadgyrB <- ggplot(data=sdtgyrB, aes_string(x='TimePoint_days',y='logreads', color = 'PlantStructure')) + 
  stat_summary(geom = "errorbar", fun.data = mean_se) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  theme(legend.position = "right") + 
  theme(legend.title = element_text(face = "bold", size = 9)) + 
  ylab("Log reads") + 
  theme(axis.title.y = element_text(face="bold", size = 10)) + 
  theme(axis.title.x = element_blank()) + 
  facet_grid (~Plant, scales = "free") +
  ylim(1,6) +
  ggtitle("C")


FigCFUreads <- grid.arrange(
  plotCFU,
  plotread16S,
  plotreadgyrB,
  ncol=1,
  nrow=3)

#Remove B14 and R18 because of the low number of CFU/sequencing coverage
dgyrBf <- subset_samples(dgyrB, !TimePoint_days %in% c("B14", "R18"))
d16Sf <- subset_samples(d16S, !TimePoint_days %in% c("B14", "R18"))

#Remove empty ASV 
dgyrBf <- filter_taxa(dgyrBf, function(x) sum(x) > 0, TRUE)
d16Sf <- filter_taxa(d16Sf, function(x) sum(x) > 0, TRUE)

# ....6) Perform rarefaction curves----
#Load source file
source("Rscripts/richness.R")

#Perfrom rarefaction curves
rcurve16 <- ggrare(d16Sf, 
                   step = 100, 
                   color = "TimePoint_days", 
                   se = FALSE) + 
  ylim (0, 150) + 
  xlim (0, 125000) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14)) +
  labs (title = "C")

#Most samples reach asymptope around 5,000 reads but lot of samples are removed at this threshold
rcurvegyrB <- ggrare(dgyrBf, 
                     step = 100, 
                     color = "TimePoint_days", 
                     se = FALSE)+ 
  ylim (0, 150) + 
  xlim (0, 125000) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14)) +
  labs (title = "D")

#Figure S2----

tiff("Figures/FigS2.tiff", width=15, height=11, units="in", res=600)
FigS2 <- grid.arrange(
  phyla_16S,
  phyla_gyrB,
  rcurve16,
  rcurvegyrB,
  ncol=2,
  nrow=2)  
dev.off()

saveRDS(dgyrBf, "gyrBf_PS.rds")
saveRDS(d16Sf, "16Sf_PS.rds")

#Final dataset, use for all subsequent analyses----
dgyrB <- readRDS("gyrBf_PS.rds") # 1082 taxa, 48 samples

#Create a new variable (floral)
Floral = c("Flower_Bud", "Open_Flower", "Fruits")
sample_data(dgyrB)$floral <- get_variable(dgyrB, "PlantStructure") %in% Floral
dgyrB.b <- subset_samples(dgyrB, Plant=="Bean")
dgyrB.b <- filter_taxa(dgyrB.b, function(x) sum(x) > 0, TRUE)
dgyrB.r <- subset_samples(dgyrB, Plant=="Radish")
dgyrB.r <- filter_taxa(dgyrB.r, function(x) sum(x) > 0, TRUE)

#Merge by floral or seeds
merged_b <- merge_samples(dgyrB.b, "floral")
merged_r <- merge_samples(dgyrB.r, "floral")

#Extract count datasets from phyloseq
count_b <- t(otu_table(merged_b))
count_r <- t(otu_table(merged_r))

#Transform in binary matrix
count_b[count_b> 0] <- 1
count_r[count_r > 0] <- 1

#Convert to dataframe
df_gyrB_bean <- as.data.frame(count_b)
df_gyrB_radish <- as.data.frame(count_r)

#Upset

Upset_gyrBb <- upset(df_gyrB_bean,
                     nintersects = 50, 
                     order.by = "freq", 
                     decreasing ="TRUE")
Upset_gyrBr <- upset(df_gyrB_radish,
                     nintersects = 50, 
                     order.by = "freq", 
                     decreasing ="TRUE")

# C] ALPHA-DIVERSITY----

## Load required libraries
library(FSA) ; packageVersion("FSA") # '0.8.26'
library(ggplot2); packageVersion("ggplot2") # '3.2.1'
library(picante); packageVersion("picante") # '1.8'
library(rcompanion) ; packageVersion("rcompanion") # '2.3.21'

#Subset seeds samples
dgyrB.s <- subset_samples(dgyrB, PlantStructure=="Seeds")
dgyrB.s <- filter_taxa(dgyrB.s, function(x) sum(x) > 0, TRUE)

# ....1) Rarefy datasets----
rgyrB <- rarefy_even_depth(dgyrB.s, sample.size = 44000, rngseed = 400)

# ....2) Data processing----
#Table of alpha-diversity estimators
table_rgyrB <- estimate_richness(rgyrB, split = TRUE, measures = NULL)

#Add Pielou evenness H/log(obs)
table_rgyrB$Evenness <- (table_rgyrB$Shannon)/(log(table_rgyrB$Observed))

#Add Faith PD
rgyrBPD <- as.data.frame.matrix(t(otu_table(rgyrB)))

FaithPDgyrB <- pd(samp = t(rgyrBPD),tree = phy_tree(rgyrB), include.root = F)

#Bind design + table of alpha_div
sdrgyrB <- sample_data(dgyrB.s)
datargyrB <- cbind(sample_data(sdrgyrB), table_rgyrB) 
datargyrBPD <- merge(datargyrB, FaithPDgyrB, by="row.names", all=TRUE)

#Subset by Plant species
datargyrBPDb <- subset(datargyrBPD, !Plant %in% c("Radish"))
datargyrBPDr <- subset(datargyrBPD, !Plant %in% c("Bean"))

#Summary
summary(datargyrBPDb)
summary(datargyrBPDr)

# ....3) Significant tests---- 
#Because of the low number of replicates per samples (3 replicates), we can't use any parametric tests.
#In fact, the power of the normality tests are to low for samples with less than 30 values (Razali et Wah, 2011)
#Kruskall-Wallis, non-parametric test
#Follow by Dunn post-hoc

##Observed
#Bean
kruskal.test(Observed ~ TimePoint_days, datargyrBPDb) #not significant 

#Radish
kruskal.test(Observed ~ TimePoint_days, datargyrBPDr) #not significant 

##Evenness
#Bean
kruskal.test(Evenness ~ TimePoint_days, datargyrBPDb) #not significative

#Radish
kruskal.test(Evenness ~ TimePoint_days, datargyrBPDr) #not significative 

##Faith PD 
#Bean
kruskal.test(PD ~ TimePoint_days, datargyrBPDb) #significant P = 0.04344

#Dunn test (post-hoc)
PT <- dunnTest(PD ~ TimePoint_days, datargyrBPDb)
PT2 <- PT$res

groupgyrBPDb <- cldList(comparison = PT2$Comparison, 
                        p.value    = PT2$P.adj,
                        threshold  = 0.05) # Significative groups

#Radish
kruskal.test(PD ~ TimePoint_days, datargyrBPDr) # not significative 

#Dunn test
PT <- dunnTest(PD ~ TimePoint_days, datargyrBPDr)
PT2 <- PT$res

groupgyrBPDr <- cldList(comparison = PT2$Comparison, 
                        p.value    = PT2$P.adj,
                        threshold  = 0.5) # ns
groupgyrBPDr$Letter <- c("a","a","a","a","a","a")
#Bind PD stats
groupgyrBPD <- rbind(groupgyrBPDb, groupgyrBPDr)

#Label groups on plot
names(groupgyrBPD)[1]<-paste("TimePoint_days")
groupgyrBPD$TimePoint_days      <- c("B22","B28","B35","B42",
                                     "R25","R32","R39","R46","R53","R67")
groupgyrBPD$Plant      <- c("Bean","Bean","Bean","Bean",
                            "Radish","Radish","Radish","Radish","Radish","Radish")

# ....4) Plot alpha diversity index----

#Reorder legend
datargyrBPD$PlantStructure <- factor(datargyrBPD$PlantStructure, 
                                    levels = c("Flower_Bud","Open_Flower","Fruits","Seeds"))

#Plot
gyrBObs <- ggplot(datargyrBPD, aes(x = TimePoint_days , y = Observed, color = PlantStructure)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, size = 1) +
  stat_summary(fun.y = mean, geom = "point", size = 4) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14)) +
  facet_wrap(~Plant, scales="free_x") +
  scale_color_manual(values=c("#bf80ff"))+
  ggtitle("A")

gyrBEvenness <- ggplot(datargyrBPD, aes(x = TimePoint_days, y = Evenness, color = PlantStructure)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, size = 1) +
  stat_summary(fun.y = mean, geom = "point", size = 4) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14)) +
  guides(fill = guide_legend( c("Flower_Bud","Open_Flower","Fruits","Seeds"))) +
  facet_wrap(~Plant, scales="free_x") +
  scale_color_manual(values=c("#bf80ff"))+
  ggtitle("B")

gyrBPD <- ggplot(datargyrBPD, aes(x = TimePoint_days, y = PD, color = PlantStructure)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, size = 1) +
  stat_summary(fun.y = mean, geom = "point", size = 4) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14)) +
  geom_text(data = groupgyrBPD, aes(label = Letter, y = 8), colour="black", size=4, vjust =0.5) +
  facet_wrap(~Plant, scales="free_x") +
  scale_color_manual(values=c("#bf80ff"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  ggtitle("A")

#D] TAXONOMIC COMPOSITION----

## Load required libraries and source
library(reshape2); packageVersion("reshape2") # '1.4.3'
source("Rscripts/graphical_methods.R")

#Subset plants samples
rgyrB.sb <- subset_samples(rgyrB, Plant=="Bean")
rgyrB.sb <- filter_taxa(rgyrB.sb, function(x) sum(x) > 0, TRUE)
rgyrB.sr <- subset_samples(rgyrB, Plant=="Radish")
rgyrB.sr <- filter_taxa(rgyrB.sr, function(x) sum(x) > 0, TRUE)

#Plot composition (Order-level)
taxgyrBb <- plot_composition(rgyrB.sb, "Kingdom", "gyrB", "Order", numberOfTaxa=7, fill="Order") + 
  facet_wrap(~TimePoint_days, scales="free_x", nrow=1) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("C") +
  scale_color_manual(values=c("#5882FA", "darkorchid1", "pink3", "limegreen", "darkgreen", "darkorange", "darkred", "grey", "black")) +
  scale_fill_manual(values=c("#5882FA", "darkorchid1", "pink3", "limegreen", "darkgreen", "darkorange", "darkred", "grey", "black")) 

taxgyrBr <- plot_composition(rgyrB.sr, "Kingdom", "gyrB", "Order", numberOfTaxa=7, fill="Order") + 
  facet_wrap(~TimePoint_days, scales="free_x", nrow=1) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  ggtitle("D") + 
  scale_color_manual(values=c("#5882FA", "darkorchid1", "pink3", "darkblue", "limegreen", "darkgreen", "darkorange", "grey", "black")) +
  scale_fill_manual(values=c("#5882FA", "darkorchid1", "pink3", "darkblue", "limegreen", "darkgreen", "darkorange", "grey", "black"))

#E] BETA-DIVERSITY----

## Load required libraries
library(betapart); packageVersion("betapart") # '1.5.1'
library(dplyr); packageVersion("dplyr") # '0.8.3'
library(ggplot2); packageVersion("ggplot2") # '3.2.1'
library(gridExtra); packageVersion("gridExtra") # '2.3'
library(vegan); packageVersion("vegan") # '2.5.6'

# ....1) Pcoa / capscale----
rgyrB.sbPcoa <-transform_sample_counts(rgyrB.sb,function(x)1E6*x/sum(x))
uunifa.pcoa.gyrB_sb <-ordinate(rgyrB.sbPcoa,"PCoA","uunifrac")
wunifrac.pcoa.gyrB_sb<-ordinate(rgyrB.sbPcoa,"PCoA","wunifrac")

rgyrB.srPcoa <-transform_sample_counts(rgyrB.sr,function(x)1E6*x/sum(x))
uunifa.pcoa.gyrB_sr <-ordinate(rgyrB.srPcoa,"PCoA","uunifrac")
wunifrac.pcoa.gyrB_sr<-ordinate(rgyrB.srPcoa,"PCoA","wunifrac")

i.uunifrac.pcoa.gyrBb<-plot_ordination(rgyrB.sb,
                                  uunifa.pcoa.gyrB_sb,
                                  type="samples",
                                  color="TimePoint_days")+geom_point(size=3)+ggtitle("G")

i.wunifrac.pcoa.gyrBb<-plot_ordination(rgyrB.sb,
                                  wunifrac.pcoa.gyrB_sb,
                                  type="samples",
                                  color="TimePoint_days")+geom_point(size=3)+ggtitle("G") +
  scale_color_manual(values=c("#E3CEF6", "#BE81F7", "#AC58FA", "#9A2EFE")) + 
  theme(axis.text.x=element_text(angle = 70, hjust = 1)) +
  ggtitle("E")

i.uunifrac.pcoa.gyrBr<-plot_ordination(rgyrB.sr,
                                  uunifa.pcoa.gyrB_sr,
                                  type="samples",
                                  color="TimePoint_days")+geom_point(size=3)+ggtitle("G")

i.wunifrac.pcoa.gyrBr<-plot_ordination(rgyrB.sr,
                                  wunifrac.pcoa.gyrB_sr,
                                  type="samples",
                                  color="TimePoint_days")+geom_point(size=3)+ggtitle("G")+
  scale_color_manual(values=c("#E3CEF6", "#D0A9F5", "#BE81F7", "#AC58FA", "#9A2EFE", "#8904B1")) + 
  theme(axis.text.x=element_text(angle = 70, hjust = 1)) +
  ggtitle("F")

#Capscale

## convert sample_data to data.frame
metadatagyrBb <- as(sample_data(rgyrB.sb), "data.frame") 
metadatagyrBr <- as(sample_data(rgyrB.sr), "data.frame")

#Bean
#Capscale on Jaccard
dist.jac.gyrBb <- phyloseq::distance(rgyrB.sb, method = "jaccard")
cap.jac.gyrBb <- capscale(dist.jac.gyrBb ~ TimePoint_days,
                          data = metadatagyrBb)
anova.jac.gyrBb  <- anova(cap.jac.gyrBb, permutations = 999)
print(anova.jac.gyrBb )# not significant

#Capscale on unweighted Unifrac
dist.uUF.gyrBb <- phyloseq::distance(rgyrB.sb, method = "uunifrac")
cap.uUF.gyrBb <- capscale(dist.uUF.gyrBb ~ TimePoint_days,
                              data = metadatagyrBb)
anova.uUF.gyrBb  <- anova(cap.uUF.gyrBb, permutations = 999)
print(anova.uUF.gyrBb )# not significant

#Capscale on Bray-Curtis
dist.bray.gyrBb <- phyloseq::distance(rgyrB.sb, method = "bray")
cap.bray.gyrBb <- capscale(dist.bray.gyrBb ~ TimePoint_days,
                          data = metadatagyrBb)
anova.bray.gyrBb  <- anova(cap.bray.gyrBb, permutations = 999)
print(anova.bray.gyrBb )# not significant

#Capscale on weighted Unifrac
dist.wUF.gyrBb <- phyloseq::distance(rgyrB.sb, method = "wunifrac")
cap.wUF.gyrBb <- capscale(dist.wUF.gyrBb ~ TimePoint_days,
                              data = metadatagyrBb)
anova.wUF.gyrBb  <- anova(cap.wUF.gyrBb, permutations = 999)
print(anova.wUF.gyrBb )# not significant

#Radish
#Capscale on Jaccard
dist.jac.gyrBr <- phyloseq::distance(rgyrB.sr, method = "jaccard")
cap.jac.gyrBr <- capscale(dist.jac.gyrBr ~ TimePoint_days,
                          data = metadatagyrBr)
anova.jac.gyrBr  <- anova(cap.jac.gyrBr, permutations = 999)
print(anova.jac.gyrBr )# not significant

#Capscale on unweighted Unifrac
dist.uUF.gyrBr <- phyloseq::distance(rgyrB.sr, method = "uunifrac")
cap.uUF.gyrBr <- capscale(dist.uUF.gyrBr ~ TimePoint_days,
                          data = metadatagyrBr)
anova.uUF.gyrBr  <- anova(cap.uUF.gyrBr, permutations = 999)
print(anova.uUF.gyrBr )# significant (P=0.003)
adonis.uUF.gyrBr  <- adonis(dist.uUF.gyrBr ~ TimePoint_days, data = metadatagyrBr, perm = 999)
print(adonis.uUF.gyrBr)# significant (P=0.006, 39%)

#Capscale on Bray-Curtis
dist.bray.gyrBr <- phyloseq::distance(rgyrB.sr, method = "bray")
cap.bray.gyrBr <- capscale(dist.bray.gyrBr ~ TimePoint_days,
                           data = metadatagyrBr)
anova.bray.gyrBr  <- anova(cap.bray.gyrBr, permutations = 999)
print(anova.bray.gyrBr)# not significant

#Capscale on weighted Unifrac
dist.wUF.gyrBr <- phyloseq::distance(rgyrB.sr, method = "wunifrac")
cap.wUF.gyrBr <- capscale(dist.wUF.gyrBr ~ TimePoint_days,
                          data = metadatagyrBr)
anova.wUF.gyrBr  <- anova(cap.wUF.gyrBr, permutations = 999)
print(anova.wUF.gyrBr)# not significant

#Figure S3----
FigS3a <- grid.arrange(
  taxgyrBb,
  taxgyrBr,
  ncol=2,
  nrow = 1)

FigS3b <- grid.arrange(
  i.wunifrac.pcoa.gyrBb,
  i.wunifrac.pcoa.gyrBr,
  ncol=2,
  nrow = 1)

tiff("Figures/FigS3.tiff", width=11, height=15, units="in", res=600)
FigS3 <- grid.arrange(
  gyrBObs,
  gyrBEvenness,
  FigS3a,
  FigS3b,
  ncol=1,
  nrow = 4)
dev.off()

tiff(file="saving_plot3.tiff",
     width=6, height=4, units="in", res=100)
hist(Temperature, col="steelblue")
dev.off()

# ....2) Mutlivariate dispersion between stages (Bean + Radish)----
dist_rgyrB_seed <- phyloseq::distance(physeq = rgyrB, method="wunifrac", type="samples")
var_rgyrB_seed <- sample_data(rgyrB)
stage_rgyrB_seed <- as.character(var_rgyrB_seed$TimePoint_days)
bd_rgyrB_seed <- betadisper(dist_rgyrB_seed, group=stage_rgyrB_seed, type="centroid")
permutest(bd_rgyrB_seed)
var_rgyrB_seed$dispersion <- bd_rgyrB_seed$distances

#Calculte observed means

var_rgyrB_seed_mean <- var_rgyrB_seed %>% 
  group_by(TimePoint_days) %>% 
  summarise(
    dispersion  = mean(dispersion))

var_rgyrB_seed_mean$Plant <- c("Bean","Bean","Bean","Bean","Radish","Radish","Radish","Radish", "Radish", "Radish")
var_rgyrB_seed_mean$TimePoint_days <- c("B22","B28","B35","B42","R25","R32","R39", "R46", "R53", "R67")
var_rgyrB_seed_mean$PlantStructure <- c("Seeds","Seeds","Seeds","Seeds","Seeds","Seeds","Seeds", "Seeds", "Seeds", "Seeds")

#Plot
disp.rgyrB_seed <- ggplot(data=var_rgyrB_seed, 
                    aes_string(x='TimePoint_days',y='dispersion', color = 'PlantStructure')) + 
  geom_point(data = var_rgyrB_seed_mean,size=4) +
  stat_summary(geom = "errorbar", fun.data = mean_se, size=1) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x=element_text(angle = 70, hjust = 1, size = 14),
        axis.title.x = element_blank()) + 
  scale_color_manual(values=c("#bf80ff")) +
  facet_wrap(~Plant, scales = "free_x") +
  ggtitle("B")

# ....3) Bacteria Partitioning Beta Diversity ----
#Merge by habitat
merged_b <- merge_samples(rgyrB.sb, "TimePoint_days")
merged_r <- merge_samples(rgyrB.sr, "TimePoint_days")

#Extract count datasets from phyloseq
count_b <- t(otu_table(merged_b))
count_r <- t(otu_table(merged_r))

#Transform in binary matrix
count_b[count_b> 0] <- 1
count_r[count_r > 0] <- 1

#Convert to dataframe
df_gyrB_bean <- as.data.frame(t(count_b))
df_gyrB_radish <- as.data.frame(t(count_r))

#Betapart

phy_tree(merged_b) <-root(phy_tree(merged_b), sample(taxa_names(merged_b),1), resolve.root = TRUE)
phylo.bean <- phy_tree(merged_b)
phylo.beta.pair(df_gyrB_bean, phylo.bean, index.family="jaccard")
phy_tree(merged_r) <-root(phy_tree(merged_r), sample(taxa_names(merged_r),1), resolve.root = TRUE)
phylo.radish <- phy_tree(merged_r)
phylo.beta.pair(df_gyrB_radish, phylo.radish, index.family="jaccard")

#Plot beta components
dbetacomp <- read.csv("Betapart_components_bean_radish.csv", sep = ";", dec = ",",check.names=FALSE)

beta_comp_plot_gyrB <-ggplot(dbetacomp, aes(X, values,fill=Index)) +
  geom_col(position= "fill") +ggtitle("B") + ylab("Index") + 
  theme(axis.title.x=element_blank()) + 
  scale_color_manual(values=c("#e6ccff", "#ce99ff")) +
  scale_fill_manual(values=c("#e6ccff", "#ce99ff")) +
  facet_grid(~Plant, scales = "free")+
  ggtitle ("C") +
  theme( axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 70, hjust = 1, size=14),
         axis.title.y = element_text(size=14),
         axis.text.y = element_text( size=14))

#Figure 1----
tiff("Figures/Fig1.tiff", width=11, height=11, units="in", res=600)
Fig1 <- grid.arrange(
  gyrBPD,
  disp.rgyrB_seed,
  beta_comp_plot_gyrB,
  ncol=1,
  nrow = 3)
dev.off()

#F) HEATMAP----

## Load required libraries
library(pheatmap); packageVersion("pheatmap") #'1.0.12'
library(RColorBrewer); packageVersion("RColorBrewer") #'1.1.2'
library(viridis); packageVersion("viridis") #'0.5.1'
library(dendextend); packageVersion("dendextend") #'1.13.2'
library(vegan); packageVersion("vegan") #'2.5.6'

##....1) Data processing----
#Transform data into relative abundance
dgyrB2 <- transform_sample_counts(dgyrB, function(x) x / sum(x) ) 
dgyrB3 <- transform_sample_counts(dgyrB2,function(x) ifelse(x>0.001,x,0))
dgyrB4 <- filter_taxa(dgyrB3, function(x) sum(x > 0.001) > (0.025*length(x)), TRUE) #Keep taxa with relative abundance higher than 1??? and found across at least 2 samples

#Execute prevalence filter, using `prune_taxa()` function
keepTaxa.gyrB <- colnames(otu_table(dgyrB4))
dgyrBf <- prune_taxa(keepTaxa.gyrB, dgyrB)

#Subset and by plant species
dgyrB.b <- subset_samples(dgyrBf, Plant %in% c("Bean"))
dgyrB.b <- filter_taxa(dgyrB.b,function(x) mean (x)>0, TRUE)
dgyrB.r <- subset_samples(dgyrBf, Plant %in% c("Radish"))
dgyrB.r <- filter_taxa(dgyrB.r,function(x) mean (x)>0, TRUE)

#ASV Table
dgyrB.b.count <- t(otu_table(dgyrB.b))
dgyrB.r.count <- t(otu_table(dgyrB.r))

#....2) Heatmap plots----
#Data normalization (Log10 + 1)
dgyrB.b.count.norm <-log10(dgyrB.b.count + 1)
dgyrB.r.count.norm <-log10(dgyrB.r.count + 1)

#Cosine similarity
norm.b <- apply(as.matrix(dgyrB.b.count.norm),1,function(x) norm(as.matrix(x),"f")) 
norm.r <- apply(as.matrix(dgyrB.r.count.norm),1,function(x) norm(as.matrix(x),"f")) 
S.b <- dgyrB.b.count.norm %*% t(dgyrB.b.count.norm)
S.r <- dgyrB.r.count.norm %*% t(dgyrB.r.count.norm)
divide_one_norm.b <- S.b/norm.b
divide_one_norm.r <- S.r/norm.r
cosine.b <- t(divide_one_norm.b)/norm.b
cosine.r <- t(divide_one_norm.r)/norm.r

hclust_gyrBb <- hclust(dist(cosine.b), method = "complete")
hclust_gyrBr <- hclust(dist(cosine.r), method = "complete")

#Heatmap Bean
#Columns labels
my_sample_col_bean <- data.frame(sample = rep(c("Flower_Bud","Open_Flower","Fruits","Seeds"), c(3,3,3,12)))
row.names(my_sample_col_bean) <- colnames(dgyrB.b.count.norm)

#Columns colours
my_colour_bean = list(
  sample = c(Flower_Bud = "brown2", 
             Open_Flower = "chartreuse4",
             Fruits = "deepskyblue2",
             Seeds = "#bf80ff"),
  cluster = c("cluster 1" = "#bf80ff", "cluster 2" = "#FE9A2E", "cluster 3" = "#a50000", "cluster 4" = "#F781F3"))

#Cluster annotation
as.dendrogram(hclust_gyrBb) %>%
  plot(horiz = TRUE)
my_gene_col_gyrBb <- cutree(tree = as.dendrogram(hclust_gyrBb), k = 4) #Can be changed according to the number of clusters wanted
my_gene_col_gyrBb <- data.frame(my_gene_col_gyrBb)
colnames(my_gene_col_gyrBb) <- c("cluster")

#Generate cluster names
dict = c("1" = 'cluster 1', "2" = "cluster 2", "3" = "cluster 3", "4" = "cluster 4")

#Add cluster name for legend
for (i in 1:4){my_gene_col_gyrBb <- replace(my_gene_col_gyrBb, my_gene_col_gyrBb == names(dict[i]), dict[i])}

#Figure 2----
tiff("Figures/Fig2.tiff", width=11, height=11, units="in", res=600)
Heatmap_gyrBb <- pheatmap(mat = dgyrB.b.count.norm,
                          gaps_col=c(6,9,21),
                          annotation_col = my_sample_col_bean,
                          annotation_row = my_gene_col_gyrBb,
                          annotation_colors = my_colour_bean,
                          cluster_cols =   FALSE,
                          clustering_distance_rows = dist(cosine.b), method = "euclidean",
                          border_color      = NA,
                          show_colnames     = TRUE,
                          show_rownames     = TRUE,
                          drop_levels       = TRUE,
                          fontsize          = 14,
                          main              = "",
                          width = 18, height = 15)

tiff("Figures/Fig2.tiff", width=11, height=11, units="in", res=600)
Heatmap_gyrBb
dev.off()

#Heatmap Radish
#Columns labels
my_sample_col_radish <- data.frame(sample = rep(c("Flower_Bud","Open_Flower","Fruits","Seeds"), c(3,3,3,18)))
row.names(my_sample_col_radish) <- colnames(dgyrB.r.count.norm)

#Columns colours
my_colour_radish = list(
  sample = c(Flower_Bud = "brown2", 
             Open_Flower = "chartreuse4",
             Fruits = "deepskyblue2",
             Seeds = "#bf80ff"),
  cluster = c("cluster 1" = "#bf80ff", "cluster 2" = "#FE9A2E", "cluster 3" = "#a50000", "cluster 4" = "#F781F3"))

#Cluster annotation
as.dendrogram(hclust_gyrBr) %>%
  plot(horiz = TRUE)
my_gene_col_gyrBr <- cutree(tree = as.dendrogram(hclust_gyrBr), k = 4) #Can be changed according to the number of clusters wanted
my_gene_col_gyrBr <- data.frame(my_gene_col_gyrBr)
colnames(my_gene_col_gyrBr) <- c("cluster")

#Generate cluster names
dict = c("1" = 'cluster 1', "2" = "cluster 2", "3" = "cluster 3", "4" = "cluster 4")

#Add cluster name for legend
for (i in 1:4){my_gene_col_gyrBr <- replace(my_gene_col_gyrBr, my_gene_col_gyrBr == names(dict[i]), dict[i])}

#Figure S4----
Heatmap_gyrBr <- pheatmap(mat = dgyrB.r.count.norm,
                          gaps_col=c(6,9),
                          annotation_col = my_sample_col_radish,
                          annotation_row = my_gene_col_gyrBr,
                          annotation_colors = my_colour_radish,
                          cluster_cols =   FALSE,
                          clustering_distance_rows = dist(cosine.r), method = "euclidean",
                          border_color      = NA,
                          show_colnames     = TRUE,
                          show_rownames     = TRUE,
                          drop_levels       = TRUE,
                          fontsize          = 14,
                          main              = "",
                          width = 18, height = 15)

tiff("Figures/FigS4.tiff", width=11, height=15, units="in", res=600)
Heatmap_gyrBr
dev.off()
#G] DATA FOR ITOL (Interactive Tree of Life)----

## Load required libraries
library(ape); packageVersion("ape") #'5.3'

#....1) Extract data for Itol----
#Extract phylogenetic tree with relative abundance higher than 1??? and found at least in 2 samples
tree_data <- rbind( as.matrix(my_gene_col_gyrBb),as.matrix (my_gene_col_gyrBr))
tree_data <- row.names(tree_data)
tree_data <- unique(tree_data)

pruned.tree<-drop.tip(phy_tree(dgyrB),phy_tree(dgyrB)$tip.label[-match(tree_data, phy_tree(dgyrB)$tip.label)])

#....2) Export tree for Figure S5--------
write.tree(phy = phy_tree(pruned.tree), file = "tree_extract_FigS5.tre")


######################################################################
#END OF SCRIPT
#Chesneau et al., 2020, "Temporal dynamics of bacterial assemblages during seed development and maturation". 
#Data collected in 2016-2017
#Authors: Guillaume Chesneau, IRHS, guillaume.chesneau <at> inra.fr
######################################################################