########Torres-Cortes et al 2018
#Functional Microbial Features Driving Community Assembly During Seed Germination and Emergence
###############################################################################################

### cleaning memory
rm(list=ls())

## setting working directory
setwd("~/Documents/Torres-Cortes_etal_2018")

#Load the different packages.
library("phytools")
library("ggplot2")
library(plyr)
library(reshape2)
library(grid)
library(gridExtra)
library(scales)
library(ampvis)
library("DESeq2")
library("phyloseq")
library("beanplot")

#Load a ggplot theme
theme_set(theme_bw())

#Load custom R functions for graphics
source("Functions/graphical_methods.R")
source("Functions/tree_methods.R")
source("Functions/plot_merged_trees.R")
source("Functions/specificity_methods.R")
source("Functions/ternary_plot.R")
source("Functions/richness.R")
source("Functions/edgePCA.R")
source("Functions/copy_number_correction.R")
source("Functions/import_frogs.R")
source("Functions/prevalence.R")

##########################################################################################################
                                      #Taxonomic profiling#
#####################

#Load abundant bacterial genus dataset (all abundant genus (>0.1%))
genus_abund <- as.matrix(read.csv("Input_tables/genus_abund_ba.txt", sep = "\t", row.names = 1, check.names=FALSE))
genus_taxabund <- as.matrix(read.csv("Input_tables/genus_abund_tax_ba.txt", sep = "\t", row.names = 1))
sampledata <- read.csv("Input_tables/design.txt", sep = "\t", row.names = 1, check.names=FALSE)
abu <- phyloseq(otu_table(genus_abund, taxa_are_rows = TRUE),
                     tax_table(genus_taxabund),
                     sample_data(sampledata))

#Change the order of factors
sample_data(abu)$Stage= factor(sample_data(abu)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(abu)$Stage)

#Split phyloseq object in bean + radish
bean <- subset_samples(abu, Plant=="Bean")
radish <- subset_samples(abu, Plant=="Radish")

###################################################################################
#############################
#Analysis of alpha-diversity_Figure1c#
#############################
# 1- Rarefy datasets at 100,000 reads per sample
rbean <- rarefy_even_depth(bean, sample.size = 100000, rngseed = 722)
rrad <- rarefy_even_depth(radish, sample.size = 100000, rngseed = 722)

# 2- Graphical overview of alpha diversity
pbean <- plot_richness(rbean, color="Stage", measures=c("Observed", "InvSimpson"), x ="Stage") + geom_point(size=3) 
prad <- plot_richness(rrad, color="Stage", measures=c("Observed", "InvSimpson"), x ="Stage") + geom_point(size=3)
alpha <- grid.arrange(pbean, prad, nrow=2, ncol=1)
#decline of genus richness/diversity in bean and radish during emergence.

# 3-Table of alpha-diversity estimators
tbean <- estimate_richness(rbean, split = TRUE, measures = NULL)
tradish <- estimate_richness(rrad, split = TRUE, measures = NULL)

# 4- Bind sampledata + table of alpha_div
dbbean <- cbind(sampledata, tbean)
dbrad <- cbind(sampledata, tradish)

# 5- ANOVA + HSD on Richness
ric_bean <- aov(Observed ~ Stage, dbbean)
summary(ric_bean)
TukeyHSD(ric_bean)

ric_rad <- aov(Observed ~ Stage, dbrad)
summary(ric_rad)
TukeyHSD(ric_rad)


# 6- ANOVA + HSD on Diversity
div_bean <- aov(InvSimpson ~ Stage, dbbean)
summary(div_bean)
TukeyHSD(div_bean)

div_rad <- aov(InvSimpson ~ Stage, dbrad)
summary(div_rad)
TukeyHSD(div_rad)
#statistically significant decrease of richness/diversity between seeds, germinating seeds and seelings.

#############################
##Taxonomic composition_Figure1a#
#############################

# 1- Overview of taxonomic composition (Class-level)
ptaxbean <- plot_composition(bean, "Kingdom", "Bacteria", "Class", numberOfTaxa=7, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)
ptaxrad <- plot_composition(radish, "Kingdom", "Bacteria", "Class", numberOfTaxa=7, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)
graph_tax_man <- grid.arrange(ptaxbean, ptaxrad, nrow=2, ncol=1)
#Bean seed bacterial assemblages are mainly composed of Gammaproteobacteria
#Radish seed bacterial assemblages are dominated by Actinobacteria, Gammaproteobacteria and Alphaproteobacteria
#Gammaproteobacteria is the most abundant taxa in seedlings

#genera
ptaxbeang <- plot_composition(bean, "Kingdom", "Bacteria", "Genus", numberOfTaxa=7, fill="Genus") + facet_wrap(~Stage, scales="free_x", nrow=1)
ptaxradg <- plot_composition(radish, "Kingdom", "Bacteria", "Genus", numberOfTaxa=7, fill="Genus") + facet_wrap(~Stage, scales="free_x", nrow=1)
graph_tax_man <- grid.arrange(ptaxbeang, ptaxradg, nrow=2, ncol=1)

#####Boxplot_Figure1d####

# 2- Convert read to percentage
per.bean <- transform_sample_counts(bean, function(x) x / sum(x) * 100)
per.rad <- transform_sample_counts(radish, function(x) x / sum(x) * 100)

#Change the order of factors
sample_data(per.bean)$Stage= factor(sample_data(per.bean)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(per.bean)$Stage)

sample_data(per.rad)$Stage= factor(sample_data(per.rad)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(per.rad)$Stage)

#3- Box plot on top 5 genus
gen.bean <-amp_rabund(data = per.bean, 
                      order.group = c("Seed", "Germinating", "Seedling"),
                      tax.aggregate = "Genus",
                      tax.show = 5,
                      scale.seq = 100,
                      group = "Stage")

gen.rad <-amp_rabund(data = per.rad,
                     order.group = c("Seed", "Germinating", "Seedling"),
                     tax.aggregate = "Genus",
                     tax.show = 5,
                     scale.seq = 100,
                     group = "Stage")

graph_gen <- grid.arrange(gen.bean, gen.rad, nrow=1, ncol=2)
# Pantoea, Pseudomonas, Enterobacter, Erwinia and Unknown are the top 5 taxa in bean
# Pseudomonas, Pantoea, Streptomyces, Arthrobacter and Bacillus are the top 5 taxa in radish

#4- Box plot on top 7 Class
gen.bean.c <-amp_rabund(data = per.bean, 
                      order.group = c("Seed", "Germinating", "Seedling"),
                      tax.aggregate = "Class",
                      tax.show = 7,
                      scale.seq = 100,
                      group = "Stage")

gen.rad.c <-amp_rabund(data = per.rad,
                     order.group = c("Seed", "Germinating", "Seedling"),
                     tax.aggregate = "Class",
                     tax.show = 7,
                     scale.seq = 100,
                     group = "Stage")

graph_gen <- grid.arrange(gen.bean.c, gen.rad.c, nrow=1, ncol=2)


#############################
##Performe DESeq2 for assessing differences in RA of taxa#
############################

#1. Bean dataset
beandds <- phyloseq_to_deseq2(bean, ~ Stage)
beandds <- DESeq(beandds, test="LRT", reduced= ~ 1)
beanres = results(beandds, contrast=c("Stage","Seedling","Seed"))

alpha = 0.01

sigtab = beanres[which(beanres$padj < alpha), ]
sigtabup = sigtab[which(sigtab$log2FoldChange>=3), ]
sigtabdown = sigtab[which(sigtab$log2FoldChange<=-3), ]
sigtabup2 = cbind(as(sigtabup, "data.frame"), as(tax_table(bean)[rownames(sigtabup), ], "matrix"))
sigtabdown2 = cbind(as(sigtabdown, "data.frame"), as(tax_table(bean)[rownames(sigtabdown), ], "matrix"))

head(sigtabup2)
head(sigtabdown2)
dim(sigtabup2)
dim(sigtabdown2)
sigtabf = rbind(sigtabup2, sigtabdown2)
head(sigtabf)
dim(sigtabf)

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Family order
x = tapply(sigtabf$log2FoldChange, sigtabf$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Family = factor(as.character(sigtabf$Family), levels=names(x))

beandf = ggplot(sigtabf, aes(x=Family, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandf)



# Class
x = tapply(sigtabf$log2FoldChange, sigtabf$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Class = factor(as.character(sigtabf$Class), levels=names(x))

beandf = ggplot(sigtabf, aes(x=Class, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandf)


#  order
x = tapply(sigtabf$log2FoldChange, sigtabf$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Order = factor(as.character(sigtabf$Order), levels=names(x))

beandf = ggplot(sigtabf, aes(x=Order, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandf)

#  genus
x = tapply(sigtabf$log2FoldChange, sigtabf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Genus = factor(as.character(sigtabf$Genus), levels=names(x))

beandfgenus = ggplot(sigtabf, aes(x=Genus, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandfgenus)

#  genus enriched
x = tapply(sigtabup2$log2FoldChange, sigtabup2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabup2$Genus = factor(as.character(sigtabup2$Genus), levels=names(x))

beandfgenus = ggplot(sigtabup2, aes(x=Genus, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandfgenus)


#2. Radish dataset
radds <- phyloseq_to_deseq2(radish, ~ Stage)
radds <- DESeq(radds, test="LRT", reduced= ~ 1)
radres = results(radds, contrast=c("Stage","Seedling","Seed"))

alpha = 0.01

sigtab = radres[which(radres$padj < alpha), ]
sigtabup = sigtab[which(sigtab$log2FoldChange>=3), ]
sigtabdown = sigtab[which(sigtab$log2FoldChange<=-3), ]
sigtabup2 = cbind(as(sigtabup, "data.frame"), as(tax_table(radish)[rownames(sigtabup), ], "matrix"))
sigtabdown2 = cbind(as(sigtabdown, "data.frame"), as(tax_table(radish)[rownames(sigtabdown), ], "matrix"))

head(sigtabup2)
head(sigtabdown2)
dim(sigtabup2)
dim(sigtabdown2)
sigtabf = rbind(sigtabup2, sigtabdown2)
head(sigtabf)
dim(sigtabf)

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Family order
x = tapply(sigtabf$log2FoldChange, sigtabf$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Family = factor(as.character(sigtabf$Family), levels=names(x))

raddf = ggplot(sigtabf, aes(x=Family, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddf)

# Class
x = tapply(sigtabf$log2FoldChange, sigtabf$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Class = factor(as.character(sigtabf$Class), levels=names(x))

raddf = ggplot(sigtabf, aes(x=Class, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddf)


# Order order

x = tapply(sigtabf$log2FoldChange, sigtabf$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Family = factor(as.character(sigtabf$Order), levels=names(x))

raddf = ggplot(sigtabf, aes(x=Order, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddf)

#  genus
x = tapply(sigtabf$log2FoldChange, sigtabf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabf$Genus = factor(as.character(sigtabf$Genus), levels=names(x))

raddfgenus = ggplot(sigtabf, aes(x=Genus, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddfgenus)

#  genus enriched
x = tapply(sigtabup2$log2FoldChange, sigtabup2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabup2$Genus = factor(as.character(sigtabup2$Genus), levels=names(x))

raddfgenus = ggplot(sigtabup2, aes(x=Genus, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddfgenus)


############################
##Beta-Diversity analyses_Figure1b#
############################

# 1-Transform to even sampling depth
propbean <- transform_sample_counts(bean, function(x) 1E6 * x/sum(x))
proprad <- transform_sample_counts(radish, function(x) 1E6 * x/sum(x))

# 2- Calculate Jaccard distance
jac.pcoa.bean <- ordinate(propbean, "PCoA", "jaccard")
jac.pcoa.rad <- ordinate(proprad, "PCoA", "jaccard")

p.jac.pcoa.bean <- plot_ordination(propbean, jac.pcoa.bean, type="samples", color="Stage", title="bean_Jaccard") + geom_point(size=3) + ggtitle("A")
print(p.jac.pcoa.bean)
p.jac.pcoa.rad <- plot_ordination(proprad, jac.pcoa.rad, type="samples", color="Stage", title="radish_Jaccard") + geom_point(size=3) + ggtitle("B")
print(p.jac.pcoa.rad)

# 3- Calculate Bray-Curtis distance
bray.pcoa.bean <- ordinate(propbean, "PCoA", "bray")
bray.pcoa.rad <- ordinate(proprad, "PCoA", "bray")


p.bray.pcoa.bean <- plot_ordination(propbean, bray.pcoa.bean, type="samples", color="Stage", title="bean_BrayCurtis") + geom_point(size=3) + ggtitle("BrayCurtis Bean")
print(p.bray.pcoa.bean)
p.bray.pcoa.rad <- plot_ordination(proprad, bray.pcoa.rad, type="samples", color="Stage", title="rad_BrayCurtis") + geom_point(size=3) + ggtitle("BrayCurtis Radish")
print(p.bray.pcoa.rad)

#Graphical representation bray
graph_bdiv_allman <- grid.arrange(p.bray.pcoa.bean, p.bray.pcoa.rad)

#Graphical representation all
#graph_bdiv_allman <- grid.arrange(p.jac.pcoa.bean, p.jac.pcoa.rad, p.bray.pcoa.bean, p.bray.pcoa.rad , nrow=2, ncol=2)

#CAPscale on Bray-Curtis in bean
metadata_abu <- as(sample_data(bean), "data.frame") ## convert sample_data to data.frame
dist.bc.abu <- phyloseq::distance(bean, method = "bray")
cap.bc.abu <- capscale(dist.bc.abu  ~ Stage,
                       data = metadata_abu)
anova.bc.abu <- anova(cap.bc.abu, permutations = 999)
print(anova.bc.abu)
adonis.bc.abu <- adonis(dist.bc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.bc.abu)

# 8- CAPscale on Jaccard in bean
dist.jac.abu <- phyloseq::distance(bean, method = "jaccard")
cap.jac.abu <- capscale(dist.jac.abu ~ Stage,
                        data = metadata_abu)
anova.jac.abu <- anova(cap.jac.abu, permutations = 999)
print(anova.jac.abu)
adonis.jac.abu <- adonis(dist.jac.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.jac.abu)

#CAPscale on Bray-Curtis in radish
metadata_abu <- as(sample_data(radish), "data.frame") ## convert sample_data to data.frame
dist.bc.abu <- phyloseq::distance(radish, method = "bray")
cap.bc.abu <- capscale(dist.bc.abu  ~ Stage,
                       data = metadata_abu)
anova.bc.abu <- anova(cap.bc.abu, permutations = 999)
print(anova.bc.abu)
adonis.bc.abu <- adonis(dist.bc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.bc.abu)

# 8- CAPscale on Jaccard in radish
dist.jac.abu <- phyloseq::distance(radish, method = "jaccard")
cap.jac.abu <- capscale(dist.jac.abu ~ Stage,
                        data = metadata_abu)
anova.jac.abu <- anova(cap.jac.abu, permutations = 999)
print(anova.jac.abu)
adonis.jac.abu <- adonis(dist.jac.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.jac.abu)


#########################################################################################################
                                   #Functional profiling#
######################
################################################################################
#########################
#Load COGs####
COGs <- as.matrix(read.csv("Input_tables/COG_all.csv", sep = "\t", row.names = 1, check.names=FALSE))
Class <- as.matrix(read.csv("Input_tables/cog_tax.csv", sep = "\t", row.names = 1))
sampledata <- read.csv("Input_tables/design.txt", sep = "\t", row.names = 1, check.names=FALSE)
test<- phyloseq(otu_table(COGs, taxa_are_rows = TRUE),
                tax_table(Class),
                sample_data(sampledata))


#Change the order of factors
sample_data(test)$Stage= factor(sample_data(test)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(test)$Stage)

#Split phyloseq object in bean + radish
bean2 <- subset_samples(test, Plant=="Bean")
radish2 <- subset_samples(test, Plant=="Radish")


#Rarefy datasets at 10,000 reads per sample
rbean2 <- rarefy_even_depth(bean2, sample.size = 10000, rngseed = 722)
rrad2 <- rarefy_even_depth(radish2, sample.size = 10000, rngseed = 722)

#########################
#Graphical overview of alpha diversity####
pbean2 <- plot_richness(rbean2, color="Stage", measures=c("Observed", "InvSimpson"), x ="Stage") + geom_point(size=3) 
prad2 <- plot_richness(rrad2, color="Stage", measures=c("Observed", "InvSimpson"), x ="Stage") + geom_point(size=3)
alpha <- grid.arrange(pbean2, prad2, nrow=2, ncol=1)

#Table of alpha-diversity estimators
tbean2 <- estimate_richness(rbean2, split = TRUE, measures = NULL)
tradish2 <- estimate_richness(rrad2, split = TRUE, measures = NULL)

# 4- Bind sampledata + table of alpha_div
dbbean2 <- cbind(sampledata, tbean2)
dbrad2 <- cbind(sampledata, tradish2)

# 5- ANOVA + HSD on Richness
ric_bean2 <- aov(Observed ~ Stage, dbbean2)
summary(ric_bean2)
TukeyHSD(ric_bean2)

ric_rad2 <- aov(Observed ~ Stage, dbrad2)
summary(ric_rad2)
TukeyHSD(ric_rad2)


# 6- ANOVA + HSD on Diversity
div_bean2 <- aov(InvSimpson ~ Stage, dbbean2)
summary(div_bean2)
TukeyHSD(div_bean2)

div_rad2 <- aov(InvSimpson ~ Stage, dbrad2)
summary(div_rad2)
TukeyHSD(div_rad2)

#########################
#General classification. Figure2A####
#Overview of broad functional categories composition (COG)
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

ptaxbean2 <- plot_composition(bean2, "Category", "Broad", "Class", numberOfTaxa=15, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)+ scale_fill_manual(values = phylum_colors)
(ptaxbean2)
ptaxrad2 <- plot_composition(radish2, "Category", "Broad", "Class", numberOfTaxa=15, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)+ scale_fill_manual(values = phylum_colors)
(ptaxrad2)
graph_func_man <- grid.arrange(ptaxbean2, ptaxrad2, nrow=2, ncol=1)


ptaxbean2 <- plot_composition(bean2, "Category", "Broad", "Class", numberOfTaxa=15, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)
(ptaxbean2)
ptaxrad2 <- plot_composition(radish2, "Category", "Broad", "Class", numberOfTaxa=15, fill="Class") + facet_wrap(~Stage, scales="free_x", nrow=1)
(ptaxrad2)
graph_func_man <- grid.arrange(ptaxbean2, ptaxrad2, nrow=2, ncol=1)

###Overview mobilome

mobilome <- read.table("~/Documents/Phyloseq/Matthieu/mobilome.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


mob = ggplot(mobilome, aes(x=ID, y=Percentage, fill=Stage, group= interaction(Plant, Stage))) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none")

#########################
#Box plot on top ten COG####

#Convert read to percentage
per.bean2 <- transform_sample_counts(bean2, function(x) x / sum(x) * 100)
#write.csv(otu_table(per.bean2), file="~/Documents/toto.csv")
per.rad2 <- transform_sample_counts(radish2, function(x) x / sum(x) * 100)
#write.csv(otu_table(per.rad2), file="~/Documents/toto.csv")
#write.csv(otu_table(radish2), file="test3")

#Change the order of factors
sample_data(per.bean2)$Stage= factor(sample_data(per.bean2)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(per.bean2)$Stage)

sample_data(per.rad2)$Stage= factor(sample_data(per.rad2)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(per.rad2)$Stage)


#3- Box plot on top ten COG
gen.bean2 <-amp_rabund(data = per.bean2,
                      order.group = c("Seed", "Germinating", "Seedling"),
                      tax.aggregate = "Class",
                      tax.show = 10,
                      scale.seq = 100,
                      group = "Stage")
gen.rad2 <-amp_rabund(data = per.rad2,
                     order.group = c("Seed", "Germinating", "Seedling"),
                     tax.aggregate = "Class",
                     tax.show = 10,
                     scale.seq = 100,
                     group = "Stage")
graph_gen <- grid.arrange(gen.bean2, gen.rad2, nrow=1, ncol=2)

#Save the plot

ggsave("top10-cog.pdf", plot = graph_gen, width = 11, height = 8)

#Box plot top COG without S and R cathegories
gen.bean3 <-amp_rabund(data = per.bean2,
                       order.group = c("Seed", "Germinating", "Seedling"),
                       tax.aggregate = "Class",
                       tax.show = c("K","G","M", "E", "T", "P", "L", "C", "N"),
                       scale.seq = 100,
                       group = "Stage")

gen.rad3 <-amp_rabund(data = per.rad2,
                      order.group = c("Seed", "Germinating", "Seedling"),
                      tax.aggregate = "Class",
                      tax.show = c("U", "K", "G", "O", "T", "L", "J", "E", "N"),
                      scale.seq = 100,
                      group = "Stage")


graph_gen <- grid.arrange(gen.bean3, gen.rad3, nrow=1, ncol=2)


#########################
#Beta-Diversity analyses. Figure2C#
#########################

#Transform to even sampling depth
propbean2 <- transform_sample_counts(bean2, function(x) 1E6 * x/sum(x))
proprad2 <- transform_sample_counts(radish2, function(x) 1E6 * x/sum(x))

#Calculate Jaccard distance
jac.pcoa.bean2 <- ordinate(propbean2, "PCoA", "jaccard")
jac.pcoa.rad2 <- ordinate(proprad2, "PCoA", "jaccard")
p.jac.pcoa.bean2 <- plot_ordination(propbean2, jac.pcoa.bean2, type="samples", color="Stage", title="bean_Jaccard") + geom_point(size=3) + ggtitle("A")
print(p.jac.pcoa.bean2)
p.jac.pcoa.rad2 <- plot_ordination(proprad2, jac.pcoa.rad2, type="samples", color="Stage", title="radish_Jaccard") + geom_point(size=3) + ggtitle("B")
print(p.jac.pcoa.rad2)

#Calculate Bray-Curtis distance
bray.pcoa.bean2 <- ordinate(propbean2, "PCoA", "bray")
bray.pcoa.rad2 <- ordinate(proprad2, "PCoA", "bray")
p.bray.pcoa.bean2 <- plot_ordination(propbean2, bray.pcoa.bean2, type="samples", color="Stage", title="bean_BrayCurtis") + geom_point(size=3) + ggtitle("C")
print(p.bray.pcoa.bean2)
p.bray.pcoa.rad2 <- plot_ordination(proprad2, bray.pcoa.rad2, type="samples", color="Stage", title="rad_BrayCurtis") + geom_point(size=3) + ggtitle("D")
print(p.bray.pcoa.rad2)

#Graphical representation
graph_bdiv_allman <- grid.arrange(p.jac.pcoa.bean2, p.jac.pcoa.rad2, p.bray.pcoa.bean2, p.bray.pcoa.rad2 , nrow=2, ncol=2)
graph_bdiv_allman <- grid.arrange( p.bray.pcoa.bean2, p.bray.pcoa.rad2)

#Save the plot

ggsave("cog_bray.pdf", plot = graph_bdiv_allman, width = 11, height = 8)


#CAPscale on Bray-Curtis in bean
metadata_abu <- as(sample_data(bean2), "data.frame") ## convert sample_data to data.frame
dist.bc.abu <- phyloseq::distance(bean2, method = "bray")
cap.bc.abu <- capscale(dist.bc.abu  ~ Stage,
                       data = metadata_abu)
anova.bc.abu <- anova(cap.bc.abu, permutations = 999)
print(anova.bc.abu)
adonis.bc.abu <- adonis(dist.bc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.bc.abu)


#CAPscale on Bray-Curtis in radish
metadata_abu <- as(sample_data(radish2), "data.frame") ## convert sample_data to data.frame
dist.bc.abu <- phyloseq::distance(radish2, method = "bray")
cap.bc.abu <- capscale(dist.bc.abu  ~ Stage,
                       data = metadata_abu)
anova.bc.abu <- anova(cap.bc.abu, permutations = 999)
print(anova.bc.abu)
adonis.bc.abu <- adonis(dist.bc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.bc.abu)


#CAPscale on Jaccard in bean
metadata_abu <- as(sample_data(bean2), "data.frame") ## convert sample_data to data.frame
dist.jc.abu <- phyloseq::distance(bean2, method = "jaccard")
cap.jc.abu <- capscale(dist.jc.abu  ~ Stage,
                       data = metadata_abu)
anova.jc.abu <- anova(cap.jc.abu, permutations = 999)
print(anova.jc.abu)
adonis.jc.abu <- adonis(dist.jc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.jc.abu)


#CAPscale on Jaccard in radish
metadata_abu <- as(sample_data(radish2), "data.frame") ## convert sample_data to data.frame
dist.jc.abu <- phyloseq::distance(radish2, method = "jaccard")
cap.jc.abu <- capscale(dist.jc.abu  ~ Stage,
                       data = metadata_abu)
anova.jc.abu <- anova(cap.jc.abu, permutations = 999)
print(anova.jc.abu)
adonis.jc.abu <- adonis(dist.jc.abu ~ Stage, data = metadata_abu, perm = 9999)
print(adonis.jc.abu)

#########################
#Perform DESeq2 for assessing differences in RA of COG#
#########################

#Bean dataset
beandds2 <- phyloseq_to_deseq2(bean2, ~ Stage)
beandds2 <- DESeq(beandds2, test="LRT", reduced= ~ 1)
beanres2 = results(beandds2, contrast=c("Stage","Seedling","Seed"))

alpha = 0.01

sigtab2 = beanres2[which(beanres2$padj < alpha), ]
sigtabup3 = sigtab2[which(sigtab2$log2FoldChange>=5), ]
sigtabdown3 = sigtab2[which(sigtab2$log2FoldChange<=-5), ]
sigtabup32 = cbind(as(sigtabup3, "data.frame"), as(tax_table(bean2)[rownames(sigtabup3), ], "matrix"))
sigtabdown32 = cbind(as(sigtabdown3, "data.frame"), as(tax_table(bean2)[rownames(sigtabdown3), ], "matrix"))

#List of COG with foldchange>=5
write.table(sigtabup32, "upBean_deseq.csv")

head(sigtabup32)
head(sigtabdown32)
dim(sigtabup32)
dim(sigtabdown32)
sigtabf2 = rbind(sigtabup32, sigtabdown32)
head(sigtabf2)
dim(sigtabf2)

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
#Class order
x = tapply(sigtabf2$log2FoldChange, sigtabf2$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabf2$Class = factor(as.character(sigtabf2$Class), levels=names(x))

beandf2 = ggplot(sigtabf2, aes(x=Class, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandf2)

#COG order
x = tapply(sigtabf2$log2FoldChange, sigtabf2$COG, function(x) max(x))
x = sort(x, TRUE)
sigtabf2$COG = factor(as.character(sigtabf2$COG), levels=names(x))

beandf2 = ggplot(sigtabf2, aes(x=COG, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(beandf2)


#Radish
radds2 <- phyloseq_to_deseq2(radish2, ~ Stage)
radds2 <- DESeq(radds2, test="LRT", reduced= ~ 1)
radres2 = results(radds2, contrast=c("Stage","Seedling","Seed"))

alpha = 0.01

sigtab2 = radres2[which(radres2$padj < alpha), ]
sigtabup2 = sigtab2[which(sigtab2$log2FoldChange>=5), ]
sigtabdown2 = sigtab2[which(sigtab2$log2FoldChange<=-5), ]
sigtabup22 = cbind(as(sigtabup2, "data.frame"), as(tax_table(radish2)[rownames(sigtabup2), ], "matrix"))
sigtabdown22 = cbind(as(sigtabdown2, "data.frame"), as(tax_table(radish2)[rownames(sigtabdown2), ], "matrix"))

head(sigtabup22)
head(sigtabdown22)
dim(sigtabup22)
dim(sigtabdown22)
sigtabf2 = rbind(sigtabup22, sigtabdown22)
head(sigtabf2)
dim(sigtabf2)

#List of COG with foldchange>=5
write.table(sigtabup22, "upradish_deseq.csv")

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Class order
x = tapply(sigtabf2$log2FoldChange, sigtabf2$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabf2$Class = factor(as.character(sigtabf2$Class), levels=names(x))

raddf = ggplot(sigtabf2, aes(x=Class, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddf)

#COG order
x = tapply(sigtabf2$log2FoldChange, sigtabf2$COG, function(x) max(x))
x = sort(x, TRUE)
sigtabf2$COG = factor(as.character(sigtabf2$COG), levels=names(x))

raddf = ggplot(sigtabf2, aes(x=COG, y=log2FoldChange, color=log2FoldChange)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) +
  scale_color_gradient2(limits= c(-10, 10), low="blue", mid= "grey", high="red")
(raddf)

###############################
#Perform Anova in general COG cathegories
##############################

# Create a factor corresponding to the General COG cathegories
genfac = factor(tax_table(test)[, "Class"])


# Tabulate the counts for each Class in each sample
gentab = apply(otu_table(test), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(gentab)[, 1:10]

#Create new phyloseq

gentabtaxo <- as.matrix(read.csv("~/Documents/Phyloseq/Matthieu/gentab_tax.txt", sep = "\t", row.names = 1))
gentabp <- phyloseq(otu_table(gentab, taxa_are_rows = TRUE),
                    tax_table(gentabtaxo),
                    sample_data(sampledata))

#Change the order of factors
sample_data(gentabp)$Stage= factor(sample_data(gentabp)$Stage, levels= c("Seed", "Germinating", "Seedling"))
levels(sample_data(gentabp)$Stage)

#Split phyloseq object in bean + radish
beangeneral2 <- subset_samples(gentabp, Plant=="Bean")
radishgeneral2 <- subset_samples(gentabp, Plant=="Radish")

per.beangeneral2 <- transform_sample_counts(beangeneral2, function(x) x / sum(x) * 100)
per.radishgeneral2 <- transform_sample_counts(radishgeneral2, function(x) x / sum(x) * 100)

#Extract raw data

per.beangeneral2data =as(otu_table(per.beangeneral2), "matrix")
per.radishgeneral2data = as(otu_table(per.radishgeneral2), "matrix")

write.csv(otu_table(per.beangeneral2), file="~/Documents/bean.csv")
write.csv(otu_table(per.radishgeneral2), file="~/Documents/radish.csv")

#Prepare for anova

tperbean = t(per.bengeneral2data)
tperradish = t(per.radishgeneral2data)
  
tperbean = as.data.frame(tperbean)
tperradish = as.data.frame(tperradish)
  
tperbean$stage <- c("Germinating", "Germinating", "Germinating", "Seedling", "Seedling", "Seedling", "Seeds", "Seeds", "Seeds")
tperradish$stage <- c("Germinating", "Germinating", "Germinating", "Seedling", "Seedling", "Seedling", "Seeds", "Seeds", "Seeds")

#Statistics in bean samples

Bean_K = aov(K ~ stage, data = tperbean)
Bean_G = aov(G ~ stage, data = tperbean )
Bean_M = aov(M ~ stage, data = tperbean )
Bean_E = aov(E ~ stage, data = tperbean )
Bean_T = aov(T ~ stage, data = tperbean )
Bean_P = aov(P ~ stage, data = tperbean )
Bean_L = aov(L ~ stage, data = tperbean )
Bean_C = aov(C ~ stage, data = tperbean )

summary.aov(Bean_K)
summary.aov(Bean_G)
summary.aov(Bean_M)
summary.aov(Bean_E)
summary.aov(Bean_T)
summary.aov(Bean_P)
summary.aov(Bean_L)
summary.aov(Bean_C)

#Only significant G (0.012)

TukeyHSD(Bean_G)

#Statistics in radish samples


Radish_U = aov(U ~ stage, data = tperradish)
Radish_K = aov(K ~ stage, data = tperradish)
Radish_G = aov(G ~ stage, data = tperradish)
Radish_O = aov(O ~ stage, data = tperradish)
Radish_T = aov(T ~ stage, data = tperradish)
Radish_L = aov(L ~ stage, data = tperradish)
Radish_E = aov(E ~ stage, data = tperradish)
Radish_J = aov(J ~ stage, data = tperradish)

summary.aov(Radish_U)
summary.aov(Radish_K)
summary.aov(Radish_G)
summary.aov(Radish_O)
summary.aov(Radish_T)
summary.aov(Radish_L)
summary.aov(Radish_E)
summary.aov(Radish_J)

#Significant U, K, G, O, T, L, E (P<0.5)

TukeyHSD(Radish_U)
TukeyHSD(Radish_K)
TukeyHSD(Radish_G)
TukeyHSD(Radish_O)
TukeyHSD(Radish_T)
TukeyHSD(Radish_L)
TukeyHSD(Radish_E)

#############################################################
##Multifunctional redundancy
#############################################################

library("ggplot2")
library("plotly")

#Load the data

Mbean_Fradish <- read.table("Input_tables/ML_shannon_bean_radish.txt", row.names =1, header = TRUE, sep="\t", na.strings=c("NA", "-", "?"))

#Plot and regresion line

plotbeanradish <- ggplot(Mbean_Fradish, aes(x=Genus, y=MF, color=Plant)) +
  scale_color_manual(values=c("Black", "Red")) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line

plotbeanradish

##################################################################################################
                                 ##Copiotrophic functional traits
###################################################################################################
#####################
#rrn copy number
####################

#Load the different tables

B = read.table("Input_tables/rrn_factor_bean.txt", header=TRUE)
R = read.table("Input_tables/rrn_factor_radish_2.txt", header=TRUE)

#Create Boxplots

bean_rrn = ggplot(B, aes(x=as.factor(stage), y=diversity)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("rrn copy number")

radish_rrn = ggplot(R, aes(x=as.factor(stage), y=diversity)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("rrn copy number")

grid.arrange(bean_rrn, radish_rrn, nrow=1, ncol=2)


#Save the plots

ggsave("bean_rrn.pdf", plot = bean, width = 11, height = 8)
ggsave("radish_rrn.pdf", plot = radish, width = 11, height = 8)

#ANOVA
#load data
B = read.table("~/Documents/Metagenomics_paper/GT_06_2017/rrn_factor_bean.txt", header=TRUE)
R = read.table("~/Documents/Metagenomics_paper/GT_06_2017/rrn_factor_radish.txt", header=TRUE)

#Anova
Ba= aov(diversity ~ stage, data = B )
Ra= aov(diversity ~ stage, data = R )
Ra2= aov(Diversity ~ Stage, data = R2 )

#Summary
SumBa = summary.aov(Ba)
SumRa = summary.aov(Ra)
SumRa2 = summary.aov(Ra2)

#####################
#Classified reads
####################

#Load the different tables

Bc = read.table("Input_tables/class_bean.txt", header=TRUE)
Rc = read.table("Input_tables/class_radish.txt", header=TRUE)

#Create Boxplots

beanc = ggplot(Bc, aes(x=as.factor(stage), y=classified)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("Read Classified (%)")

radishc = ggplot(Rc, aes(x=as.factor(stage), y=classified)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("Read Classified (%)")

grid.arrange(beanc, radishc, nrow=1, ncol=2)

#Save the plots

ggsave("bean_class.pdf", plot = beanc, width = 11, height = 8)
ggsave("radish_class.pdf", plot = radishc, width = 11, height = 8)


#Anova
Bca = aov(classified ~ stage, data = Bc )
Rca= aov(classified ~ stage, data = Rc )

#Summary
SumBa = summary.aov(Bca)
SumRa = summary.aov(Rca)

# Tukey multiple comparisons of means 95% family-wise confidence level
Tuka = TukeyHSD(Bca)
TukaR = TukeyHSD(Rca)

####################
#Average genome size
####################

#Load the different tables

Bags = read.table("Input_tables/bean_ags.txt", header=TRUE)
Rags = read.table("Input_tables/radish_ags.txt", header=TRUE)


#Create Boxplots


bean_plot_ags = ggplot(Bags, aes(x=as.factor(stage), y=AGS)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("AGS")

radish_plot_ags = ggplot(Rags, aes(x=as.factor(stage), y=AGS)) + 
  geom_boxplot(fill="black", alpha=0.2) + 
  ylab("AGS")

grid.arrange(bean_plot_ags, radish_plot_ags, nrow=1, ncol=2)


plot(bean_plot_ags)
plot(radish_plot_ags)

#Save the plots

ggsave("bean_ags.pdf", plot = bean_plot_ags, width = 11, height = 8)
ggsave("radish_ags.pdf", plot = radish_plot_ags, width = 11, height = 8)


#Anova
Bags_an = aov(AGS ~ stage, data = Bags )
Rags_an = aov(AGS ~ stage, data = Rags )

#Summary
SumBa = summary.aov(Bags_an)
SumRa = summary.aov(Rags_an)

# Tukey multiple comparisons of means 95% family-wise confidence level
Tuka = TukeyHSD(Bags_an)
TukaR = TukeyHSD(Rags_an)

#####################################################################################################
                                   ####Grotwh curves analysis####
#####################################################################################################
#load packages

library(growthcurver) 
library(ggplot2)
library(gridExtra)
library(dplyr)
library(multcomp)

# load samples

dbean <- read.table("Input_tables/growth_curve_beana.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dradish <- read.table("Input_tables/growth_curve_radisha.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Grotwhcurve

gc_out_bean <- SummarizeGrowthByPlate(dbean, bg_correct = "blank")
gc_out_radish <- SummarizeGrowthByPlate(dradish, bg_correct = "blank")

write.table(gc_out_bean, "gc_out_bean", sep="\t") 
write.table(gc_out_radish, "gc_out_radish", sep="\t") 

# Generate plots for all of the growth curves 


gc_out_beangraph <- SummarizeGrowthByPlate(dbean, plot_fit = TRUE, 
                                           plot_file = "gc_plots_bean_blank.pdf")

gc_out_radishgraph <- SummarizeGrowthByPlate(dradish, plot_fit = TRUE, 
                                             plot_file = "gc_plots_radish_blank.pdf")

#Plot auc

# Remove last row

gc_out_bean2 = gc_out_bean[-c(25), ] 
gc_out_radish2 = gc_out_radish[-c(25), ] 

# Add id replicate number and oligotroph/copiotroph colum

gc_out_bean2$id <- c("R04", "R04", "R04", "R10", "R10", "R10", "R19", "R19", "R19", "R27", "R27", "R27", "H3", "H3", "H3", "H19", "H19", "H19","H13", "H13", "H13", "ASV1", "ASV1", "ASV1" )

gc_out_bean2$growth <- c("copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph","copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph","copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph","oligotroph", "oligotroph", "oligotroph")

gc_out_radish2$id <- c("R04", "R04", "R04", "R10", "R10", "R10", "R19", "R19", "R19", "R27", "R27", "R27", "H3", "H3", "H3", "H19", "H19", "H19","H13", "H13", "H13", "ASV1", "ASV1", "ASV1" )

gc_out_radish2$growth <- c("copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph","copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph","copiotroph", "copiotroph", "copiotroph", "copiotroph", "copiotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph", "oligotroph","oligotroph", "oligotroph", "oligotroph")


#Plot it
bean = ggplot(gc_out_bean2, aes(x=id, y=auc_e, fill=growth)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none")

radish = ggplot(gc_out_radish2, aes(x=id, y=auc_e, fill=growth)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none")

graph_grothw <- grid.arrange(bean, radish , nrow=1, ncol=2)


#Anova

gc_out_bean2$id = as.factor(gc_out_bean2$id)
gc_out_radish2$id = as.factor(gc_out_radish2$id)

Ba_id= aov(auc_e ~ id, data = gc_out_bean2 )
Ra_id= aov(auc_e ~ id, data = gc_out_radish2 )

#Summary

SumBaid = summary.aov(Ba_id)
SumRaid = summary.aov(Ra_id)


# Tukey multiple comparisons of means 95% family-wise confidence level

Tukai = TukeyHSD(Ba_id)
TukaRi = TukeyHSD(Ra_id)

#Plot Tukey comparisons
#Specify all pair-wise comparisons amon levels of variable id

tukbean = glht(Ba_id, linfct = mcp(id= "Tukey"))

tukradish = glht(Ra_id, linfct = mcp(id="Tukey"))

# Extract information

tukbean.cld = cld(tukbean)
tukradish.cld = cld(tukradish)

# Use sufficiently large upper margin

old.par <- par(mai=c(1,1,1.25,1), no.readonly=TRUE)
 ### plot
plot(tukbean.cld) 
  par(old.par)

  old.par <- par(mai=c(1,1,2,1), no.readonly=TRUE)
plot(tukradish.cld) 
  par(old.par)



#DeSeq data for Figure 5

Beandeseq <- read.csv("Input_tables/sigtabf.txt", sep = "\t", row.names = 1, check.names=FALSE)

#Plot

#  genus

beandfgenus = ggplot(Beandeseq, aes(x=ID, y=log2FoldChange, color=Plant)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ylim(-10, +10) + geom_hline(yintercept = 0) + scale_color_manual(values=c("darkgreen", "red"))

(beandfgenus)

#AUC and Deseq data together

grid.arrange(beandfgenus, graph_grothw, nrow=2, ncol=1)


