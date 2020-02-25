### February 2019
### Gloria Torres-Cortes & Alberto Prado 

### Bee2Seed 2017, 2018 16S data on Brassica napus seed microbiomes 

### cleaning memeory
rm(list=ls())

## setting working directory
setwd("~/Desktop/Data_Bee2Seed")

## loading packages
library(vegan)
library(RColorBrewer) 
library(phytools)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(grid)
library(gridExtra)
library(scales)
library(phyloseq)
library(gridBase)
library(devtools)
devtools::install_github("gauravsk/ranacapa")
library(ranacapa)

#Load a ggplot theme
theme_set(theme_bw())

#Load custom R functions for graphics
source("~/Desktop/Data_Bee2Seed/R_scripts/graphical_methods.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/tree_methods.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/plot_merged_trees.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/specificity_methods.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/ternary_plot.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/richness.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/edgePCA.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/copy_number_correction.R")
source("~/Desktop/Data_Bee2Seed/R_scripts/prevalence.R")


### Index ####
## Section 1 Loading data  
## Section 2 Sample selection
## Section 3 Rarefaction curves and alpha diversity
## Section 4 Ordination
## Section 5 Redundancy Analysis
## Section 6 Beta dispersion
## Section 7 Community composition

#########################################################################################
######### Section 1 loading data 
#########################################################################################

##########
#Load manually curated/generated files 
##########

asv16S <- as.matrix(read.csv2("asv_count.csv", row.names = 1, check.names = F))
asv16Staxo <- read.csv2("asv_tax.csv")
asv16Staxo <- apply(asv16Staxo, 2, as.character)
str(asv16Staxo)
tax16S <- as.data.frame(replace(asv16Staxo, is.na(asv16Staxo), "Unclassified"))
rownames(tax16S) <- tax16S$X 
tax16S <- tax16S[,-1]
tax16S <- as.matrix(tax16S)
tree16S <- read.newick("asv_seq2.tre")

sampledata <- read.csv2("design_mb1.csv")
sampledata$Year <- as.factor(sampledata$Year)
row.names(sampledata) <- sampledata$X
sampledata <- sampledata[,-c(1,9)]
sampledata$pline <- ifelse(sampledata$Experiment=="Sterile","Sterile","Fertile")
str(sampledata)

##########
#Load PhyloSeq and sequencing depth
##########
asv16S <- phyloseq(otu_table(asv16S, taxa_are_rows = TRUE),
                  tax_table(tax16S),
                  sample_data(sampledata),
                  phy_tree(tree16S))

### Cleaning unwanted taxa 
#removing chloroplasts 
asv16Sc <- subset_taxa(asv16S, Class!="Chloroplast")
#removing NA
asv16Su <- subset_taxa(asv16Sc, Phylum!="Unclassified")
#removing Archaea beacuse primers should not target this kingdom.
asv16Sa <- subset_taxa(asv16Su, Kingdom!="Archaea")
#Remove ASV with zero count
asv16Sf <- filter_taxa(asv16Sa, function(x) sum(x) >= 1, TRUE)

## Number of taxa (2764) and number of samples (198) 
nsamples(asv16S) # 198
nsamples(asv16Sf) # 198
ntaxa(asv16S) # 4439
ntaxa(asv16Sc) #4316
ntaxa(asv16Su) # 3976
ntaxa(asv16Sa) # 3974
ntaxa(asv16Sf) # here we can see that due to all the samples we have previously excluded (i.e. contaminated samples, seedlings, nectar samples after the exp, etc.)
################ their are lots of taxa at 0 counts 

# Create table, number of features for each phyla
table(tax_table(asv16Sf)[, "Phylum"], exclude = NULL)
## most abundant phyla Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes 

# Compute prevalence of each feature, store as data.frame
prevd16S <- apply(X = otu_table(asv16Sf),
                   MARGIN = ifelse(taxa_are_rows(asv16Sf), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevd16S <- data.frame(Prevalence = prevd16S ,
                        TotalAbundance = taxa_sums(asv16Sf),
                        tax_table(asv16Sf))

Phylum_prev <- plyr::ddply(prevd16S, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

ggplot(prevd16S, aes(TotalAbundance, Prevalence / nsamples(asv16Sf),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#Sequencing deep
sdt_all <- data.table(as(sample_data(asv16Sf), "data.frame"),
                        TotalReads = sample_sums(asv16Sf), keep.rownames = TRUE)

#Total number of reads per year
totalreads <- colSums (sdt_all[ , 12]) # 6,337,890
totalreads/198 ### 32,009



#########################################################################################
######### Section 2 Sample selection
#########################################################################################
##########
##Selecting new phyloseq objects: seeds, bee, nectar and pollen 
##########

#Subsetting seeds samples
seeds <- subset_samples(asv16Sf, Material== "Seed")
seeds <- filter_taxa(seeds, function(x) sum(x) >= 1, TRUE)
nsamples(seeds) # 151
ntaxa(seeds) # 2099
sdt_seeds <- data.table(as(sample_data(seeds), "data.frame"),
                      TotalReads = sample_sums(seeds), keep.rownames = TRUE)
totalreadsseeds <- colSums (sdt_seeds[ , 12]) # 5,687,939
totalreadsseeds/151 # 37,668

#Subsetting bee, nectar, seedling and pollen samples
bee <- subset_samples(asv16Sf, Material== "Bee")
bee <- filter_taxa(bee, function(x) sum(x) >= 1, TRUE)
ntaxa(bee)
nsamples(bee)
nectar <- subset_samples(asv16Sf, Material== "Nectar")
nectar <- filter_taxa(nectar, function(x) sum(x) >= 1, TRUE)
ntaxa(nectar)
nsamples(nectar)
pollen <- subset_samples(asv16Sf, Material== "Pollen")
pollen <- filter_taxa(pollen, function(x) sum(x) >= 1, TRUE)
ntaxa(pollen)
nsamples(pollen)


##########
##Selecting new phyloseq objects for PCoA, Redundancy analysis and Beta dispersion:
#Apis2018, Osmia, Sterile2018, Sterile2017 and Apis2017
##########
#Fertile
Apis <- subset_samples(asv16Sf, Experiment== "Honey")
Apis <- filter_taxa(Apis, function(x) sum(x) >= 1, TRUE)
nsamples(Apis) #95
ntaxa(Apis)
Apis2018 <- subset_samples(Apis, Year== "2018")
Apis2018 <- filter_taxa(Apis2018, function(x) sum(x) >= 1, TRUE)
nsamples(Apis2018)#75
ntaxa(Apis2018)
Apis2017 <- subset_samples(Apis, Year== "2017")
Apis2017 <- filter_taxa(Apis2017, function(x) sum(x) >= 1, TRUE)
nsamples(Apis2017)#20
ntaxa(Apis2017)
OM <- subset_samples(asv16Sf, Experiment== "Osmia")
OM <- filter_taxa(OM, function(x) sum(x) >= 1, TRUE)
nsamples(OM)
ntaxa(OM)
Tunnels <- subset_samples(asv16Sf, Year== "2018")
Tunnels <- subset_samples(Tunnels, Material!= "Bee")
Tunnels <- subset_samples(Tunnels, pline!= "Sterile")
Tunnels <- filter_taxa(Tunnels, function(x) sum(x) >= 1, TRUE)
nsamples(Tunnels) #134
ntaxa(Tunnels)

### Apis: merging samples by plant
ApisSeeds <- subset_samples(Apis2018, Material== "Seed")
ApisSeeds<- filter_taxa(ApisSeeds, function(x) sum(x) >= 1, TRUE)
nsamples(ApisSeeds)#54
ntaxa(ApisSeeds)
HB_merged <- merge_samples(ApisSeeds,"plantID")
head(sample_data(HB_merged))
nsamples(HB_merged) #18
ntaxa(HB_merged)
## Fixing metadata
sample_data(HB_merged)$Project <- levels(sample_data(ApisSeeds)$Project)
sample_data(HB_merged)$Experiment <- levels(sample_data(ApisSeeds)$Experiment)
sample_data(HB_merged)$Material <- levels(sample_data(ApisSeeds)$Material)
sample_data(HB_merged)$plantID <- levels(sample_data(ApisSeeds)$plantID) 
sample_data(HB_merged)$Plant <- levels(sample_data(ApisSeeds)$Plant)
sample_data(HB_merged)$Treatment <- levels(sample_data(ApisSeeds)$Treatment)

### Osmia: merging samples by plant
OMSeeds <- subset_samples(seeds, Experiment== "Osmia")
OMSeeds<- filter_taxa(OMSeeds, function(x) sum(x) >= 1, TRUE)
ntaxa(OMSeeds)
nsamples(OMSeeds)
om_merged <- merge_samples(OMSeeds,"plantID")
head(sample_data(om_merged))
nsamples(om_merged) #20
ntaxa(om_merged)
## Fixing metadata
sample_data(om_merged)$Project <- levels(sample_data(OMSeeds)$Project)
sample_data(om_merged)$Experiment <- levels(sample_data(OMSeeds)$Experiment)
sample_data(om_merged)$Material <- levels(sample_data(OMSeeds)$Material)
sample_data(om_merged)$plantID <- levels(sample_data(OMSeeds)$plantID) 
sample_data(om_merged)$Plant <- levels(sample_data(OMSeeds)$Plant)
sample_data(om_merged)$Treatment <- levels(sample_data(OMSeeds)$Treatment)

###Sterile line
Sterile <- subset_samples(seeds, Experiment== "Sterile")
Sterile <- filter_taxa(Sterile, function(x) sum(x) >= 1, TRUE)
nsamples(Sterile) #32
ntaxa(Sterile)
Sterile2018 <- subset_samples(Sterile, Year== "2018")
Sterile2018 <- filter_taxa(Sterile2018, function(x) sum(x) >= 1, TRUE)
nsamples(Sterile2018) #20
ntaxa(Sterile2018)
Sterile2017 <- subset_samples(Sterile, Year== "2017")
Sterile2017 <- filter_taxa(Sterile2017, function(x) sum(x) >= 1, TRUE)
nsamples(Sterile2017) #12
ntaxa(Sterile2017)

#2017
ApisSeeds2017 <- subset_samples(Apis2017, Material== "Seed")
ApisSeeds2017 <- filter_taxa(ApisSeeds2017, function(x) sum(x) >= 1, TRUE)
nsamples(ApisSeeds2017)#10
ntaxa(ApisSeeds2017)

##########
##Selecting new phyloseq objects for Redundancy analysis (section 5): 
#HB2017plant, HB2018plant
##########

HB2017plant <- subset_samples(Apis2017, Material== "Seed"| Material== "Pollen"| Material== "Nectar")
HB2017plant <- filter_taxa(HB2017plant, function(x) sum(x) >= 1, TRUE)
nsamples(HB2017plant)
ntaxa(HB2017plant)

HB2018plant <- subset_samples(Apis2018, Material== "Seed"| Material== "Pollen"| Material== "Nectar")
HB2018plant <- subset_samples(HB2018plant, Treatment!= "Pollenball")
HB2018plant <- filter_taxa(HB2018plant, function(x) sum(x) >= 1, TRUE)
nsamples(HB2018plant)
ntaxa(HB2018plant)

#########################################################################################
######### Section 3 Rarefaction curves and alpha diversity 
#########################################################################################

#Rarefaction curves with ggplot
#svg("rarefraction.svg")
seedsplot = ggrare(seeds, step = 10, label = NULL, color = "Treatment",
       plot = TRUE, parallel = FALSE, se = FALSE) + xlim (0, 40000) + geom_vline(xintercept = 12000) + labs(subtitle="Rarefaction curves") + theme(axis.text.x = element_text(size=7) + scale_colour_manual(values = c("yelow", "blue", "green")))
#dev.off()

#Rarefy datasets at 12000
rseeds <- rarefy_even_depth(seeds, sample.size = 12000, rngseed = 712)

#Remove samples that have been lost due to rarefaction
seedsf <- subset_samples(seeds, 
                         sample_names(seeds)  != "MS10-R1-D12_S188" & 
                           sample_names(seeds)  != "OB10-P2-G07_S151" &
                           sample_names(seeds)  != "OB7-P3-G06_S143" &
                           sample_names(seeds)  != "HB42-B2-F05_S38")
nsamples(seedsf) # 147

#Table of alpha-diversity estimators
table_rseeds <- estimate_richness(rseeds, split = TRUE, measures = NULL)

#Add Pielou evenness H/log(obs)----
table_rseeds$Evenness <- (table_rseeds$Shannon)/(log(table_rseeds$Observed))

#Add Faith PD----
#Convert phyloseq otu_table into data.frame
rseedPD <- as.data.frame.matrix(otu_table(rseeds))
library(picante)
FaithPDseeds <- pd(samp = t(rseedPD),tree = phy_tree(rseeds), include.root = F)

#Bind sampledata + table of alpha_div
datarseeds <- cbind(sample_data(seedsf), table_rseeds) 
datarseeds <- merge(datarseeds, FaithPDseeds, by="row.names", all=TRUE)

#Subset by year
datarseeds17 <- subset(datarseeds, !Year %in% c("2018"))
datarseeds18 <- subset(datarseeds, !Year %in% c("2017"))


############################################
#Graphical representation of alpha diversity
#########################################

#2018 dataset
p1 <- ggplot(data=datarseeds18, aes_string(x='Experiment',y='Observed', fill='Treatment')) +
  geom_boxplot(lwd=c(0.2)) + 
  ggtitle("A") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10)) + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1", "coral2"), aesthetics = c("colour", "fill")) + scale_x_discrete(labels=c("Honey" = "Apis", "Osmia" = "Osmia",
                                                                                                                                                                                                                      "Sterile" = "Apis")) +facet_grid(~pline, scales="free_x")

p2 <- ggplot(data=datarseeds18, aes_string(x='Experiment',y='Evenness', fill='Treatment')) +
  geom_boxplot(lwd=c(0.2)) + 
  ggtitle("B") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10)) + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1", "coral2"), aesthetics = c("colour", "fill"))+ scale_x_discrete(labels=c("Honey" = "Apis", "Osmia" = "Osmia",
                                                                                                                                                                                                                     "Sterile" = "Apis")) +facet_grid(~pline, scales="free_x")

p3 <- ggplot(data=datarseeds18, aes_string(x='Experiment',y='Shannon', fill='Treatment')) +
  geom_boxplot(lwd=c(0.2)) + 
  ggtitle("C") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10)) + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1", "coral2"), aesthetics = c("colour", "fill"))+ scale_x_discrete(labels=c("Honey" = "Apis", "Osmia" = "Osmia",
                                                                                                                                                                                                                     "Sterile" = "Apis")) +facet_grid(~pline, scales="free_x")
p4 <- ggplot(data=datarseeds18, aes_string(x='Experiment',y='PD', fill='Treatment')) +
  geom_boxplot(lwd=c(0.2)) + 
  ggtitle("D") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10)) + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1", "coral2"), aesthetics = c("colour", "fill"))+ scale_x_discrete(labels=c("Honey" = "Apis", "Osmia" = "Osmia",
                                                                                                                                                                                                                     "Sterile" = "Apis")) +facet_grid(~pline, scales="free_x")
#svg("alphadiv.svg")
alphadiv2018 <- grid.arrange(p1,p2,p3,p4, ncol=2,nrow=2, top="alpha diversity 2018") 
#dev.off()

#2017 dataset
p5 <- ggplot(data=datarseeds17, aes_string(x='Experiment',y='Observed', fill='Treatment')) +
  geom_boxplot() + 
  ggtitle("A") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10))

p6 <- ggplot(data=datarseeds17, aes_string(x='Experiment',y='Evenness', fill='Treatment')) +
  geom_boxplot() + 
  ggtitle("B") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10))

p7 <- ggplot(data=datarseeds17, aes_string(x='Experiment',y='Shannon', fill='Treatment')) +
  geom_boxplot() + 
  ggtitle("C") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10))

p8 <- ggplot(data=datarseeds17, aes_string(x='Experiment',y='PD', fill='Treatment')) +
  geom_boxplot() + 
  ggtitle("D") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10))


alphadiv2017 <- grid.arrange(p5,p6,p7,p8, ncol=2,nrow=2, top="alpha diversity 2017")

#####################
## alpha diversity statistics 
##########

#####2018
#Shapiro-Wilk normality test
library(nortest)
hist(datarseeds18$Observed) # not normally distributed -> some outliers
with(datarseeds18, shapiro.test(Observed)) # indeed confirm by Shapiro test
hist(datarseeds18$Evenness) # not normally distributed -> some outliers
with(datarseeds18, shapiro.test(Evenness)) # indeed confirm by Shapiro test
hist(datarseeds18$Shannon) # not normally distributed -> some outliers
with(datarseeds18, shapiro.test(Shannon)) # indeed confirm by Shapiro test
hist(datarseeds18$PD) # normally distributed 
with(datarseeds18, shapiro.test(PD)) # confirm by Shapiro test; however for sake of clarity we will perform non-parametric tests for all indexes

#Kruskall-Wallis
kruskal.test(Observed ~ ExpTreatment, datarseeds18) # p<0.05
kruskal.test(Evenness ~ ExpTreatment, datarseeds18) # p<0.05
kruskal.test(Shannon ~ ExpTreatment, datarseeds18) # p<0.05
kruskal.test(PD ~ ExpTreatment, datarseeds18) # p<0.05

#Dunn's test
library(FSA)
Dunnobs <- dunnTest(Observed ~ ExpTreatment, datarseeds18)
Dunnobs2 <- Dunnobs$res
library(rcompanion)
cldList(comparison = Dunnobs2$Comparison, 
        p.value    = Dunnobs2$P.adj,
        threshold  = 0.05) # Autogamous vs BP significant for Honeybee

DunnEv <- dunnTest(Evenness ~ ExpTreatment, datarseeds18)
DunnEv2 <- DunnEv$res
cldList(comparison = DunnEv2$Comparison, 
        p.value    = DunnEv2$P.adj,
        threshold  = 0.05) # ns

DunnSha <- dunnTest(Shannon ~ ExpTreatment, datarseeds18)
DunnSha2 <- DunnSha$res
cldList(comparison = DunnSha2$Comparison, 
        p.value    = DunnSha2$P.adj,
        threshold  = 0.05) # Autogamous vs BP significant for Honeybee and Osmia

DunnPD <- dunnTest(PD ~ ExpTreatment, datarseeds18)
DunnPD2 <- DunnPD$res
cldList(comparison = DunnPD2$Comparison, 
        p.value    = DunnPD2$P.adj,
        threshold  = 0.05) # Autogamous vs BP significant for Honeybee

#Wilcox test in the Fertile line for Apis
wilcox.test(Observed ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Honey",]) # Statistically significant
wilcox.test(Evenness ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Honey",]) # Not statistically significant
wilcox.test(Shannon ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Honey",]) # Not statistically significant
wilcox.test(PD ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Honey",]) # Not statistically significant

#Wilcox test in the Sterile line for Apis
wilcox.test(Observed ~ Treatment, datarseeds18[datarseeds18$pline=="Sterile",]) # Statistically significant
wilcox.test(Evenness ~ Treatment, datarseeds18[datarseeds18$pline=="Sterile",]) # Not statistically significant
wilcox.test(Shannon ~ Treatment, datarseeds18[datarseeds18$pline=="Sterile",]) # Not statistically significant
wilcox.test(PD ~ Treatment, datarseeds18[datarseeds18$pline=="Sterile",]) # Not statistically significant

#Wilcox test in the Fertile line for Osmia
wilcox.test(Observed ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Osmia",]) # Not statistically significant
wilcox.test(Evenness ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Osmia",]) # Not statistically significant
wilcox.test(Shannon ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Osmia",]) # Not statistically significant
wilcox.test(PD ~ Treatment, datarseeds18[datarseeds18$pline=="Fertile"&datarseeds18$Experiment=="Osmia",]) # Not statistically significant



######################################################################
#####2017
#Shapiro-Wilk normality test
library(nortest)
hist(datarseeds17$Observed) # not normally distributed 
hist(datarseeds17$Evenness) # not normally distributed 
hist(datarseeds17$Shannon) # not normally distributed 
hist(datarseeds17$PD) # not normally distributed 
#Not enough samples for performing Shapiro-Wilk test

#Subset by fertile/sterile plant

datarseeds17f <- subset(datarseeds17, Experiment %in% c("Honey"))
datarseeds17s <- subset(datarseeds17, Experiment %in% c("Sterile"))

# t test
ric_2017f <- t.test(Observed ~ Treatment, datarseeds17f)
ric_2017f # Not statistically significant p = 0.3103
ric_2017s <- t.test(Observed ~ Treatment, datarseeds17s)
ric_2017s # Not statistically significant p = 0.186

eve_2017f <- t.test(Evenness ~ Treatment, datarseeds17f)
eve_2017f # Not statistically significant p = 0.2586
eve_2017s <- t.test(Evenness ~ Treatment, datarseeds17s)
eve_2017s # Not statistically significant p = 0.1893 

Sha_2017f <- t.test(Shannon ~ Treatment, datarseeds17f)
Sha_2017f # Not statistically significant p = 0.4291 
Sha_2017s <- t.test(Shannon ~ Treatment, datarseeds17s)
Sha_2017s # Not statistically significant p = 0.375

pd_2017f <- t.test(PD ~ Treatment, datarseeds17f)
pd_2017f # Not statistically significant p = 0.1753
pd_2017s <- t.test(PD ~ Treatment, datarseeds17s)
pd_2017s # Not statistically significant p = 0.1954  

#Non parametric tests

#Kruskall-Wallis
kruskal.test(Observed ~ ExpTreatment, datarseeds17) # Not statistically significant
kruskal.test(Evenness ~ ExpTreatment, datarseeds17) # Not statistically significant
kruskal.test(Shannon ~ ExpTreatment, datarseeds17) # Not statistically significant
kruskal.test(PD ~ ExpTreatment, datarseeds17) # Not statistically significant

#Wilcox test
wilcox.test(Observed ~ Treatment, datarseeds17f) # Not statistically significant
wilcox.test(Evenness ~ Treatment, datarseeds17f) # Not statistically significant
wilcox.test(Shannon ~ Treatment, datarseeds17f) # Not statistically significant
wilcox.test(PD ~ Treatment, datarseeds17f) # Not statistically significant

wilcox.test(Observed ~ Treatment, datarseeds17s) # In wilcox.test.default(x = c(74, 117, 83, 91, 69, 90), y = c(98,  :
#impossible de calculer la p-value exacte avec des ex-aequos
wilcox.test(Evenness ~ Treatment, datarseeds17s) # Not statistically significant
wilcox.test(Shannon ~ Treatment, datarseeds17s) # Not statistically significant
wilcox.test(PD ~ Treatment, datarseeds17s) # Not statistically significant

#############################################################################################
#########  Section 4 Ordination plots 
######### 2018: Apis, Apis merged samples, Osmia, Osmia merged samples, sterile line
########  2017: Apis, sterile line
############################################################################################



#PCoA using unweighted uniFrac distance
######################
### 2018 fertile line

#Figure 1A Apis
Apis.pcoa.16S <- ordinate(Apis2018, "PCoA", "unifrac")

Apis.PCoA <- plot_ordination(Apis2018,
                               Apis.pcoa.16S,
                               type="samples",
                               color="Treatment",
                               shape="Material"
) + geom_point(size=4) + ggtitle("PCoA Apis all samples") + scale_shape_manual(values=c(15,17,4,19,1))

#svg("Apis.PCoA.svg")
Apis.PCoA <- Apis.PCoA + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1", "darkorange","coral2","deepskyblue2","darkgoldenrod","coral2"), aesthetics = c("colour", "fill")) 
#dev.off()
Apis.PCoA

#Figure 1B Osmia
FertileOM.pcoa.16S <- ordinate(OM, "PCoA", "unifrac")

Osmia.PCoA <- plot_ordination(OM,
                                  FertileOM.pcoa.16S,
                                  type="samples",
                                  color="Treatment",
                                  shape="Material"
) + geom_point(size=4) + ggtitle("PCoA Osmia all samples") + scale_shape_manual(values=c(15,17,4,19,1))


Osmia.PCoA <- Osmia.PCoA + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1","deepskyblue2", "darkorange","coral2","deepskyblue2"), aesthetics = c("colour", "fill"))  + scale_shape_manual(values=c(15, 17, 16))
Osmia.PCoA

#Figure 1C Apis grouped samples

Apismerged.pcoa.16S <- ordinate(HB_merged, "PCoA", "unifrac")

Apismerged.PCoA <- plot_ordination(HB_merged,
                               Apismerged.pcoa.16S,
                               type="samples",
                               color="Treatment",
                               shape="Material",
                               label="plantID"
) + geom_point(size=4) + ggtitle("PCoA Apis merged samples") + scale_shape_manual(values=c(16))

Apismerged.PCoA <- Apismerged.PCoA + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1"))
Apismerged.PCoA

#Figure 1D Osmia grouped samples

OMmerged.pcoa.16S <- ordinate(om_merged, "PCoA", "unifrac")

OMmerged.PCoA <- plot_ordination(om_merged,
                                 OMmerged.pcoa.16S,
                                 type="samples",
                                 color="Treatment",
                                 shape="Material",
                                 label="plantID"
) + geom_point(size=4) + ggtitle("PCoA Osmia mergeed samples") + scale_shape_manual(values=c(16))


OMmerged.PCoA <- OMmerged.PCoA + scale_colour_manual(values = c("chartreuse4", "darkgoldenrod1"))
OMmerged.PCoA
#svg("orination_2018.svg")
OrdinationPlots <- grid.arrange(Apis.PCoA, Osmia.PCoA, Apismerged.PCoA, OMmerged.PCoA, ncol=2,nrow=2, top="Fig 3. Ordination Plots 2018")
#dev.off()

### Additional Ordination plots
#### 2018 Sterile line


ST2018.pcoa.16S <- ordinate(Sterile2018, "PCoA", "unifrac")

ST2018.PCoA <- plot_ordination(Sterile2018,
                               ST2018.pcoa.16S,
                               type="samples",
                               color="Treatment",
                               shape="Material"
) + geom_point(size=4) + ggtitle("PCoA Apis all samples") + scale_shape_manual(values=c(15,17,4,19,1))

ST2018.PCoA


################
## 2017 Apis

HB2017.pcoa.16S <- ordinate(Apis2017, "PCoA", "unifrac")

HB2017.PCoA <- plot_ordination(Apis2017,
                      HB2017.pcoa.16S,
                      type="samples",
                      color="Treatment",
                      shape="Material"
) + geom_point(size=2) + ggtitle("PCoA Apis all samples") + scale_shape_manual(values=c(15,17,4,19,1))

HB2017.PCoA


#2017 Sterile line

ST2017.pcoa.16S <- ordinate(Sterile2017, "PCoA", "unifrac")

ST2017.PCoA <- plot_ordination(Sterile2017,
                               ST2017.pcoa.16S,
                               type="samples",
                               color="Treatment",
                               shape="Material"
) + geom_point(size=2) + ggtitle("PCoA Apis all samples") + scale_shape_manual(values=c(15,17,4,19,1))

ST2017.PCoA


############################################################################
#############################################################################################
######### Section 5 Redundancy analysis using capscale
##### Year effect
##### 2018: Tunnel effect, Plant tissue/material, Apis plantID, Apis Treatment, Osmia plantID, Osmia Treatment, ST
##### 2017: Plant tissue/material, plantID, Treatment, ST 
#############################################################################################

HB2018plant <- filter_taxa(HB2018plant, function(x) sum(x) >= 1, TRUE)

###############################################
## Year effect, Capscale on year using all data
metadata_abu1 <- as(sample_data(asv16Sf), "data.frame") ## convert sample_data to data.frame
dist.uf.abu1 <- phyloseq::distance(asv16Sf,method = "unifrac")
cap.uf.abu1 <- capscale(dist.uf.abu1  ~ Year,
                       data = metadata_abu1)
cap.uf.abu1
#The year accounts for 5.7% of the variance
anova.uf.abu1 <- anova(cap.uf.abu1, permutations = 999)
print(anova.uf.abu1)
#The effect is significant p=0.001
adonis.uf.abu1 <- adonis(dist.uf.abu1 ~ Year, data = metadata_abu1, perm = 9999)
print(adonis.uf.abu1)
#F= 11.837, p=0.0001

##################################
# 2018
#Capscale on tunnel effect

metadata_abu0 <- as(sample_data(Tunnels), "data.frame") ## convert sample_data to data.frame
dist.uf.abu0 <- phyloseq::distance(Tunnels,method = "unifrac")
cap.uf.abu0 <- capscale(dist.uf.abu0  ~ Experiment, data = metadata_abu0)
cap.uf.abu0
#Tunnels accounts for 4.5% of the variance
anova.uf.abu0 <- anova(cap.uf.abu0, permutations = 999)
print(anova.uf.abu0)
#The effect is significant p=0.001
adonis.uf.abu0 <- adonis(dist.uf.abu0 ~ Experiment, data = metadata_abu0, perm = 9999)
print(adonis.uf.abu0)
#F=6.22, p=0.0001

#Capscale on plant material
metadata_abu2 <- as(sample_data(HB2018plant), "data.frame") ## convert sample_data to data.frame
dist.uf.abu2 <- phyloseq::distance(HB2018plant,method = "unifrac")
cap.uf.abu2 <- capscale(dist.uf.abu2  ~ Material, data = metadata_abu2)
cap.uf.abu2
#plant material accounts for 9.29% of the variance
anova.uf.abu2 <- anova(cap.uf.abu2, permutations = 999)
print(anova.uf.abu2)
#The effect is significant p=0.001
adonis.uf.abu2 <- adonis(dist.uf.abu2 ~ Material, data = metadata_abu2, perm = 9999)
print(adonis.uf.abu2)
#F=3.28, p=0.0001

#Capscale on plantID 2018 HB seeds 
metadata_abu3 <- as(sample_data(HB_merged), "data.frame") ## convert sample_data to data.frame
dist.uf.abu3 <- phyloseq::distance(HB_merged,method = "unifrac")
cap.uf.abu3 <- capscale(dist.uf.abu3  ~ Plant,
                        data = metadata_abu3)
cap.uf.abu3
#Plant ID accounts for 44.84% of the variance in the seed microbial community
anova.uf.abu3 <- anova(cap.uf.abu3, permutations = 999)
print(anova.uf.abu3)
#The effect is NOT significant p=0.858
adonis.uf.abu3 <- adonis(dist.uf.abu3 ~ Plant, data = metadata_abu3, perm = 9999)
print(adonis.uf.abu3)
#F=0,92, p=0.86 

#Capscale on treatment 2018 HB seeds 
metadata_abu4 <- as(sample_data(HB_merged), "data.frame") ## convert sample_data to data.frame
dist.uf.abu4 <- phyloseq::distance(HB_merged,method = "unifrac")
cap.uf.abu4 <- capscale(dist.uf.abu4  ~  Treatment,
                       data = metadata_abu4)
cap.uf.abu4
#Pollination treatment accounts for 12.25% of the variance in the seed microbial community
anova.uf.abu4 <- anova(cap.uf.abu4, permutations = 999)
print(anova.uf.abu4)
#The effect is significant p=0.001
adonis.uf.abu4 <- adonis(dist.uf.abu4 ~  Treatment, data = metadata_abu4, perm = 9999)
print(adonis.uf.abu4)
#F=2.233, p=0.001
plot(cap.uf.abu4, choices=c(1,2), display = c("wa","bp"), main="Apis fertile line 2018")
points(cap.uf.abu4, choices=c(1,2), display = c("wa"), pch=19, col=as.factor(metadata_abu3$Treatment))

## plot with ggplot
smry <- summary(cap.uf.abu4)
df1  <- data.frame(smry$sites[,1:2])       # CAP1 and MDS1
df2  <- data.frame(smry$species[,1:2])     # loadings for CAP1 and MDS1
rda.plot <- ggplot(df1, aes(x=CAP1, y=MDS1)) + 
  geom_text(aes(label=rownames(df1)),size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()
rda.plot <- rda.plot + geom_point(aes(x=CAP1, y=MDS1, col=metadata_abu4$Treatment))
rda.plot

#ggsave(file="RDA.svg", plot=rda.plot, width=6, height=6)

#Capscale on treatmetn 2018 ST seeds 
metadata_abu5 <- as(sample_data(Sterile2018),"data.frame") ## convert sample_data to data.frame
dist.uf.abu5 <- phyloseq::distance(Sterile2018, method = "unifrac")
cap.uf.abu5 <- capscale(dist.uf.abu5  ~ Treatment,
                        data = metadata_abu5)
cap.uf.abu5
#6.4% 
anova.uf.abu5 <- anova(cap.uf.abu5, permutations = 999)
print(anova.uf.abu5)
#The efect is NOT saignificant
adonis.uf.abu5 <- adonis(dist.uf.abu5 ~ Treatment, data = metadata_abu5, perm = 9999)
print(adonis.uf.abu5)
#F=1.16, p=0.1836

#Capscale on plant ID 2018 OM
metadata_abu6 <- as(sample_data(FeromSeeds),"data.frame")
dist.uf.abu6 <- phyloseq::distance(FeromSeeds, method = "unifrac")
cap.uf.abu6 <- capscale(dist.uf.abu6  ~ Plant,
                        data = metadata_abu6)
cap.uf.abu6
#16.71% of the variance is explained by the treatment
anova.uf.abu6 <- anova(cap.uf.abu6, permutations = 999)
print(anova.uf.abu6)
#The efect is NOT saignificant p=0.436
adonis.uf.abu6 <- adonis(dist.uf.abu6 ~ Treatment, data = metadata_abu6, perm = 9999)
print(adonis.uf.abu6)
#F=1.33, p=0.0235


#Capscale on treatment 2018 OM
metadata_abu7 <- as(sample_data(om_merged),"data.frame")
dist.uf.abu7 <- phyloseq::distance(om_merged, method = "unifrac")
cap.uf.abu7 <- capscale(dist.uf.abu7  ~ Treatment,
                       data = metadata_abu7)
cap.uf.abu7
#5.97 % of the variance is explained by the treatment
anova.uf.abu7 <- anova(cap.uf.abu7, permutations = 999)
print(anova.uf.abu7)
#The efect is NOT saignificant
adonis.uf.abu7 <- adonis(dist.uf.abu7 ~ Treatment, data = metadata_abu7, perm = 9999)
print(adonis.uf.abu7)
#F=1.14, p=0.08


########################################
# 2017
#Capscale on plant material
metadata_abu8 <- as(sample_data(HB2017plant), "data.frame") ## convert sample_data to data.frame
dist.uf.abu8 <- phyloseq::distance(HB2017plant,method = "unifrac")
cap.uf.abu8 <- capscale(dist.uf.abu8  ~ Material, data = metadata_abu8)
cap.uf.abu8
#plant material accounts for 35.54% of the variance
anova.uf.abu8 <- anova(cap.uf.abu8, permutations = 999)
print(anova.uf.abu8)
#The effect is significant p=0.001
adonis.uf.abu8 <- adonis(dist.uf.abu8 ~ Material, data = metadata_abu8, perm = 9999)
print(adonis.uf.abu8)
#F=3.584, p=0.0002

#Capscale on treatment 2017 Apis
metadata_abu9 <- as(sample_data(ApisSeeds2017), "data.frame") ## convert sample_data to data.frame
dist.uf.abu9 <- phyloseq::distance(ApisSeeds2017,method = "unifrac")
cap.uf.abu9 <- capscale(dist.uf.abu9  ~ Treatment, data = metadata_abu9)
cap.uf.abu9
# Treatment accounts for 12.97% of the variance
anova.uf.abu9 <- anova(cap.uf.abu9, permutations = 999)
print(anova.uf.abu9)
#The effect is NOT significant
adonis.uf.abu9 <- adonis(dist.uf.abu9 ~ Treatment, data = metadata_abu9, perm = 9999)
print(adonis.uf.abu9)
#F=1.1919, p=0.128

#Capscale on treatment 2017 Apis Steriel line
metadata_abu10 <- as(sample_data(Sterile2017), "data.frame") ## convert sample_data to data.frame
dist.uf.abu10 <- phyloseq::distance(Sterile2017,method = "unifrac")
cap.uf.abu10 <- capscale(dist.uf.abu10  ~ Treatment, data = metadata_abu10)
cap.uf.abu10
# Treatment accounts for 9.688% of the variance
anova.uf.abu10 <- anova(cap.uf.abu10, permutations = 999)
print(anova.uf.abu10)
#The effect is NOT significant p=0.273
adonis.uf.abu10 <- adonis(dist.uf.abu10 ~ Treatment, data = metadata_abu10, perm = 9999)
print(adonis.uf.abu9)
#F=1.0727, p=0.2857

########################################################################################################################
#####################################################################################################################
######### Section 6 Beta dispersion
#########################################################################
## 2018 HBseeds, HB_merged, OMseeds, OM_merged
##############

############### Dispersion

############################################################################################
### variation amongst plants

dist_m1 <- distance(physeq=HB_merged, method="uniFrac", type="samples")
hist(dist_m1)
shapiro.test(dist_m1)

HBseeds <- as.data.frame(sample_data(ApisSeeds))
Treat <- rep(c(1,2),9)
Treat2 <- rep(c(1,1,1,2,2,2),9)
HBseeds_m <- sample_data(HB_merged)

bd1 <- betadisper(dist_m1, group=HBseeds_m$Treatment, type="centroid")
bd1
plot(bd1, axes=c(1,2), cex=1, label=TRUE, label.cex=0.5, col=Treat)
anova(bd1)
boxplot(bd1, ylab = "Distance to centroid", col=c("black","red"))

dist_c1 <- unlist(bd1[3])
HBseeds_m$dist_c1 <- dist_c1
boxplot(HBseeds_m$dist_c1~HBseeds_m$Treatment, notch=TRUE, ylab="Distance to centroid")
a1 <- aov(HBseeds_m$dist_c1~HBseeds_m$Treatment)
summary(a1)
wilcox.test(HBseeds_m$dist_c1~HBseeds_m$Treatment)


distance_cent1 <-HBseeds_m[,c(6,12)]
distance_cent1$type <- "Variation between plants"
names(distance_cent1) <- c("Treatment","distance", "type")


dp1 <- ggplot(data=distance_cent1, aes_string(x='Treatment',y='distance', fill=NULL)) +
  geom_boxplot(outlier.size=-1) + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none") + theme(axis.title.x=element_blank()) + ylab("Distance to centroid") + ylim(0.3,0.6)

dp1
### Bee pollination reduces variance within plants!!!

###Osmia

dist_m2 <- distance(physeq=om_merged, method="uniFrac", type="samples")
hist(dist_m2)
shapiro.test(dist_m2)
OMseeds <- as.data.frame(sample_data(om_merged))
Treat3 <- rep(c(1,2),10)
dist_c2 <- unlist(bd2[3])
OMseeds$dist_c2 <- dist_c2


bd2 <- betadisper(dist_m2, group=OMseeds$Treatment, type="centroid")
bd2
plot(bd2, axes=c(1,2), cex=1, label=TRUE, label.cex=0.5, col=Treat3)
anova(bd2)
boxplot(bd2, ylab = "Distance to centroid", col=c("black","red"))
## same trand as Apis but not significant
#NO differences


dp2 <- ggplot(data=OMseeds, aes_string(x='Treatment',y='dist_c2', fill=NULL)) +
  geom_boxplot(outlier.size=-1) + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none") + theme(axis.title.x=element_blank()) + ylab("Distance to centroid") + ylim(0.3,0.6)

dp2



#####################################
###STERILE LINE

dist_m3 <- distance(physeq=Sterile2018, method="uniFrac", type="samples")
hist(dist_m3)
shapiro.test(dist_m3)
ST_2018 <- sample_data(Sterile2018)
bd3 <- betadisper(dist_m3, group=ST_2018$Treatment, type="centroid")
bd3
Treat4 <- factor(ST_2018$Treatment, labels=c(1,2))
plot(bd3, axes=c(1,2), cex=1, label=TRUE, label.cex=0.5, col=Treat4)
anova(bd3)
boxplot(bd5, ylab = "Distance to centroid", col=c("black","red"))
### No difference

STseeds_m <- sample_data(Sterile2018)
dist_c3 <- unlist(bd3[3])
STseeds_m$dist_c3 <- dist_c3
distance_cent3 <-STseeds_m[,c(6,12)]
distance_cent3$type <- "Variation between plants"
names(distance_cent3) <- c("Treatment","distance", "type")


dp3 <- ggplot(data=distance_cent3, aes_string(x='Treatment',y='distance', fill=NULL)) +
  geom_boxplot(outlier.size=-1) + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none") + theme(axis.title.x=element_blank()) + ylab("Distance to centroid") + ylim(0.3,0.6)

dp3

#Fig 4
#svg("betadispersion.svg")
grid.arrange(dp1,dp3, dp2,ncol=3,nrow=1, top="Fig 4. Beta-dispersion 2018")
#dev.off()


####################
## 2017

dist_m5 <- distance(physeq=ApisSeeds2017, method="uniFrac", type="samples")
hist(dist_m5)
shapiro.test(dist_m5)

HB_2017seeds <- as.data.frame(sample_data(ApisSeeds2017))
Treat6 <- factor(HB_2017seeds$Treatment, label=c(1,2))

bd6 <- betadisper(dist_m5, group=HB_2017seeds$Treatment, type="centroid")
bd6
plot(bd6, axes=c(1,2), cex=1, label=TRUE, label.cex=0.5, col=Treat6)
scores(bd6, display=c("sites","centroids"))
anova(bd6)
## Not significant


#########################################################################################
######### Section 7 Community composition
#########################################################################################

#Abundance

p2017 = plot_bar(HB2017seeds, fill = "Class")
p2018HB = plot_bar(HB_merged, fill = "Class")
p2019OM = plot_bar(OM_merged, fill= "Class")


tax_HB2017_order <- plot_composition(HB2017seeds, "Kingdom", "Bacteria", "Order", numberOfTaxa=5, fill="Order") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_HB2017_order
tax_HB2017_genus <- plot_composition(HB2017seeds, "Kingdom", "Bacteria", "Genus", numberOfTaxa=5, fill="Genus") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_HB2017_genus

tax_HB2018_order <- plot_composition(HB_merged, "Kingdom", "Bacteria", "Order", numberOfTaxa=10, fill="Order") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_HB2018_order
tax_HB2018_genusHB <- plot_composition(HB_merged, "Kingdom", "Bacteria", "Genus", numberOfTaxa=10, fill="Genus") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_HB2018_genusHB

tax_OM2018_order <- plot_composition(OM_merged, "Kingdom", "Bacteria", "Order", numberOfTaxa=10, fill="Order") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_OM2018_order
tax_OM2018_genusOM <- plot_composition(OM_merged, "Kingdom", "Bacteria", "Genus", numberOfTaxa=10, fill="Genus") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_OM2018_genusOM


Taxa_genus <- grid.arrange(tax_HB2017_genus, tax_HB2018_genusHB, tax_OM2018_genusOM, ncol=3,nrow=1, top="Taxa composition")


#Acinetobacter dominating year 2018



tax_seed_asv <- plot_composition(mergedHB, "Kingdom", "Bacteria", "Species", numberOfTaxa=10, fill="Species") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_seed_asv

plot_heatmap(mergedHB, taxa.label = "Genus")



tax_aci_genus <- plot_composition(MergedHB_aci, "Kingdom", "Bacteria", "Genus", numberOfTaxa=10, fill="Genus") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_aci_genus
tax_aci_genusHB <- plot_composition(HB_aci, "Kingdom", "Bacteria", "Genus", numberOfTaxa=10, fill="Genus") + facet_wrap(~Treatment, scales="free_x", nrow=1)
tax_aci_genusHB

###############################################################################
###############################################################################
### Comparing bee associated taxa

### grouping OTUs by Genus
seed_genera <- tax_glom(seeds, taxrank="Genus")
seed_genera <- transform_sample_counts(seed_genera, function(x) x / sum(x) * 100)
OT1 <- as.data.frame(otu_table(seed_genera))
TT1 <- as.data.frame(tax_table(seed_genera))
OT1 <- merge(OT1,TT1, by=0)
tail(OT1)
OT2 <- OT1

bee_asso_gen <- c("Arsenophonus","Bombella","Gilliamella","Frischella","Moraxella","Lactobacillus","Snodgrassella","Spiroplasma")
OT1 <- OT1[OT1$Genus%in%bee_asso_gen,-c(1,154:158)]
OT1 <- melt(OT1,id.vars = c("Genus"))
names(OT1) <- c("Genus","SmpName","ReadsPercentage")
OT1 <- merge(OT1, sampledata, by.x="SmpName",by.y=0)

## 2017
OT2017hb <- OT1[OT1$Year=="2017"&OT1$Experiment=="Honey",]
OT2017ST <- OT1[OT1$Year=="2017"&OT1$Experiment=="Sterile",]

bp1 <- ggplot(data=OT2017hb, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) 
bp1 + facet_grid(~Genus, scales="free_y")

bp2 <- ggplot(data=OT2017ST, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) + facet_grid(~Genus)
bp2

OT2017hb <-droplevels(OT2017hb)

## Figure 4A insert bee associated genera 2017
#svg("genera2017.svg")
par(mfrow = c(1,4))
boxplot(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Arsenophonus",])
boxplot(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Frischella",])
boxplot(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Moraxella",])
boxplot(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Bombella",])
#dev.off()
par(mfrow = c(1,1))

## 2018
OT2018hb <- OT1[OT1$Year=="2018"&OT1$Experiment=="Honey",]
OT2018ST <- OT1[OT1$Year=="2018"&OT1$Experiment=="Sterile",]
bp3 <- ggplot(data=OT2018hb, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) + facet_grid(~Genus)
bp3

bp4 <- ggplot(data=OT2018ST, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) + facet_grid(~Genus)
bp4

## Osmia
OTOsmia <- OT1[OT1$Year=="2018"&OT1$Experiment=="Osmia",]

bp5 <- ggplot(data=OTOsmia, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) + facet_grid(~Genus)
bp5

## most abundant genera 2018
most_gen <- c("Acinetobacter","Pantoea","Pseudomonas","Sphingoboium","Paracoccus","Streptococcus")
OT2 <- OT2[OT2$Genus%in%most_gen,-c(1,154:158)]
OT2 <- melt(OT2,id.vars = c("Genus"))
names(OT2) <- c("Genus","SmpName","ReadsPercentage")
OT2 <- merge(OT2, sampledata, by.x="SmpName",by.y=0)
OT2 <- OT2[OT2$Year=="2018"&OT2$Experiment=="Honey",]
OT2 <- droplevels(OT2)

bp6 <- ggplot(data=OT2, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) + facet_grid(~Genus)
bp6

## Figure 4B most abundabt genera 2018
#svg("genera2018.svg")
par(mfrow = c(1,4))
boxplot(ReadsPercentage~Treatment, data=OT2[OT2$Genus=="Acinetobacter",])
boxplot(ReadsPercentage~Treatment, data=OT2[OT2$Genus=="Pantoea",])
boxplot(ReadsPercentage~Treatment, data=OT2[OT2$Genus=="Paracoccus",])
boxplot(ReadsPercentage~Treatment, data=OT2[OT2$Genus=="Pseudomonas",])
#dev.off()
par(mfrow = c(1,1))

####################################################
#### STATS
# 2017
#HB
#Arsenophonus
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Arsenophonus",])
# difference in read abundance

#Bombella
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Bombella",])
# difference in read abundance

#Frischella
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Frischella",])
# NO difference in read abundance

#Lactobacillus
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Lactobacillus",])
# Higher abundance of Lactobacillus in Autogamous

#Gilliamella
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Gilliamella",])
# NO difference in read abundance

#Snodgrassella
wilcox.test(ReadsPercentage~Treatment, data=OT2017hb[OT2017hb$Genus=="Snodgrassella",])
# NO difference in read abundance

#########################################
#Sterile line 2107
#Arsenophonus
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Arsenophonus",])
# NO difference in read abundance

#Bombella
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Bombella",])
# NO difference in read abundance

#Frischella
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Frischella",])
# NO difference in read abundance

#Lactobacillus
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Lactobacillus",])
# NO difference in read abundance

#Gilliamella
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Gilliamella",])
# NO difference in read abundance

#Moraxella
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Moraxella",])
# NO difference in read abundance

#Snodgrassella
wilcox.test(ReadsPercentage~Treatment, data=OT2017ST[OT2017ST$Genus=="Snodgrassella",])
# NO difference in read abundance
################################
# 2018 #####
# HB
# No Arsenophonus
# No Bombella
# No Frischella

#Lactobacillus
wilcox.test(ReadsPercentage~Treatment, data=OT2018hb[OT2018hb$Genus=="Lactobacillus",])
# Difference in read abundance

#Gilliamella
wilcox.test(ReadsPercentage~Treatment, data=OT2018hb[OT2018hb$Genus=="Gilliamella",])
# NO difference in read abundance

#Moraxella
wilcox.test(ReadsPercentage~Treatment, data=OT2018hb[OT2018hb$Genus=="Moraxella",])
# NO difference in read abundance

#Snodgrassella
wilcox.test(ReadsPercentage~Treatment, data=OT2018hb[OT2018hb$Genus=="Snodgrassella",])
# NO difference in read abundance
#########################################
#Sterile line 2107

#Lactobacillus
wilcox.test(ReadsPercentage~Treatment, data=OT2018ST[OT2018ST$Genus=="Lactobacillus",])
# NO difference in read abundance

#Moraxella
wilcox.test(ReadsPercentage~Treatment, data=OT2018ST[OT2018ST$Genus=="Moraxella",])
# Difference in read abundance

#Snodgrassella
wilcox.test(ReadsPercentage~Treatment, data=OT2018ST[OT2018ST$Genus=="Snodgrassella",])
# NO difference in read abundance

# 2018 OSMIA #####
# Osmia
# No Arsenophonus
# No Bombella
# No Frischella
# No Gilliamella

#Lactobacillus
wilcox.test(ReadsPercentage~Treatment, data=OTOsmia[OTOsmia$Genus=="Lactobacillus",])
# Difference in read abundance

#Moraxella
wilcox.test(ReadsPercentage~Treatment, data=OTOsmia[OTOsmia$Genus=="Moraxella",])
# NO difference in read abundance

#Snodgrassella
wilcox.test(ReadsPercentage~Treatment, data=OTOsmia[OTOsmia$Genus=="Snodgrassella",])
# NO difference in read abundance


##############################################################################
##################################
#####
#### Plotting graphs based on lefese results by ASV


OT3 <- as.data.frame(otu_table(ApisSeeds))
OT3 <- t(OT3)
OT3 <- decostand(OT3, method = "total")
apply(OT3, 1, sum)
OT3 <- as.data.frame(t(OT3))

TT3 <- as.data.frame(tax_table(ApisSeeds))
OT3 <- merge(OT3,TT3, by=0)
str(OT3)

lefse_ASV <- c("ASV0049","ASV0044","ASV0040","ASV0038","ASV0026","ASV0022","ASV0021","ASV0015")
OT3 <- OT3[OT3$Row.names%in%lefse_ASV,c(1:55,61)]
OT3 <- melt(OT3,id.vars = c("Row.names","Genus"))
names(OT3) <- c("ASV","Genus","SmpName","ReadsPercentage")
OT3 <- merge(OT3, sampledata, by.x="SmpName",by.y=0)
OT3$ASV <- paste(OT3$ASV,OT3$Genus) 

bp11 <- ggplot(data=OT3, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) 
bp11 + facet_wrap(Genus ~ ., scales="free_y")

OT3mean <- ddply(OT3, c("ASV","Treatment"), summarise,
                 mean= mean(ReadsPercentage),
                 sd=  sd(ReadsPercentage),
                 N= length(ReadsPercentage),
                 se=  sd /sqrt(N))
bp1 <- ggplot(OT3mean, aes(x = Treatment, y = mean, fill=Treatment))
bp1 <- bp1 + geom_bar(stat='identity', color="black") + ylab("Read Percentage")
bp1 <- bp1 + geom_errorbar(aes(ymin=mean, ymax=mean+se), size=0.5,   
                      width=.25,position=position_dodge(.9)) 
#svg("ASV_lefse2018.svg")
bp1 + theme_bw() + facet_wrap(ASV ~ ., scales="free")
#dev.off()

##############################################################################

OT4 <- as.data.frame(otu_table(ApisSeeds2017))
OT4 <- t(OT4)
OT4 <- decostand(OT4, method = "total")
apply(OT4, 1, sum)
OT4 <- as.data.frame(t(OT4))


TT4 <- as.data.frame(tax_table(ApisSeeds2017))
OT4 <- merge(OT4,TT4, by=0)

lefse_ASV2017 <- c("ASV0004","ASV0009","ASV0023","ASV0058","ASV0075","ASV0100","ASV0226","ASV0194")
OT4 <- OT4[OT4$Row.names%in%lefse_ASV2017,c(1:11,17)]
OT4 <- melt(OT4,id.vars = c("Row.names","Genus"))
names(OT4) <- c("ASV","Genus","SmpName","ReadsPercentage")
OT4 <- merge(OT4, sampledata, by.x="SmpName",by.y=0)
OT4$ASV <-paste(OT4$ASV,OT4$Genus)

bp12 <- ggplot(data=OT4, aes_string(x='Treatment',y='ReadsPercentage', fill=NULL)) +
  geom_boxplot() + geom_point(size=2, alpha=1/5) +
  theme(legend.position = "none")  + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7)) 
bp12 + facet_wrap(ASV ~ ., scales="free_y")

OT4mean <- ddply(OT4, c("ASV","Treatment"), summarise,
                 mean= mean(ReadsPercentage),
                 sd=  sd(ReadsPercentage),
                 N= length(ReadsPercentage),
                 se=  sd /sqrt(N))
bp12 = ggplot(OT4mean, aes(x = Treatment, y = mean, fill=Treatment))
bp12 = bp12 + geom_bar(stat='identity', color="black") + ylab("Read Percentage")
bp12 = bp12 + geom_errorbar(aes(ymin=mean, ymax=mean+se), size=0.5,   
                          width=.25,position=position_dodge(.9)) 

#svg("ASV_lefse2017.svg")
bp12 + theme_bw() + facet_wrap(ASV ~ ., scales="free")
#dev.off()




