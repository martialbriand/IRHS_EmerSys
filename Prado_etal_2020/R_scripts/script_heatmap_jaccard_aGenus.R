#Draw a heatmap with a clique file genererated with Clark (adapted from http://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/)

#Open R or Rstudio

#load the gplots package for heatmap.2
library(gplots)  
#load the vegan package for hierachical clustering.
library(vegan)
#load the RColorBrewer package for better colour options
library(RColorBrewer) 

#load data
all.data = read.table("Z:/1 Ongoing projects/Bee2Seed/Data_2017/Data_Gloria_2017/R_scripts/aGenera.txt", header=TRUE)

#open the dataframe
View(all.data)

#Convert the group ids to row names
row.names(all.data) <- all.data$Group
all.data <- all.data[, -1]

#Transform the raw counts of reads to proportions within a sample
data.prop = (all.data/rowSums(all.data))*100000

#Transform to log 10
data.prop.1 = (data.prop + 1)
data.prop.1.log = log(data.prop.1, 10)

#Determine the maximum relative abundance for each column
maxab = apply(data.prop.1.log, 2, max)
head(maxab)

#Remove the clique that have less than 1% as their maximum relative abundance
n1 = names(which(maxab < 3))
data.prop.selected = data.prop.1.log[, -which(names(data.prop) %in% n1)]

#Calculate the Bray dissimilarity matrix on the full dataset:
data.dist <- vegdist(data.prop.1.log, method = "bray")

#Perform average linkage hierarchical clustering. 
row.clus <- hclust(data.dist, "aver")

#Cluster the clique by co-occurence pattern. DON'T FORGET TO TRANSPOSE THE DATASET.
data.dist.otu = vegdist(t(data.prop.selected), method = "jaccard")
col.clus = hclust(data.dist.otu, "aver")


#Load heatmap 
heatmap.2(as.matrix(data.prop.selected), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = rev(heat.colors(256)), margins = c(7, 3), trace = "none", density.info = "none", main = "aGenus bray", lhei = c(2, 8))

#Most abundant without clustering
heatmap(as.matrix(data.prop.selected), Rowv = NA, Colv = NA, col = rev(heat.colors(256)), margins = c(10, 2))

# make the heatmap with Rowv = as.dendrogram(row.clus)

heatmap(as.matrix(data.prop.selected), Rowv = as.dendrogram(row.clus), Colv = NA, col = rev(heat.colors(256)), margins = c(10, 3))

#Alternatively, we may want to change the colors of the heatmap.
#Creates a colour palette
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(15, 7), trace = "none", density.info = "none", main = "gyrB", lhei = c(2, 8)) 

#Highlighted the different contrast with different colors
#Load the design file
factors = read.table("~/partage_windows7/EMERSYS/metaSEED_Axe_3/Analysis-2013-2014/heatmap/ITS.inoc.design", header=TRUE)
v1 = as.vector(factors$group)
v2 = as.vector(factors$color)
test = cbind(v1, v2)
heatmap.2(as.matrix(data.prop.selected), 
          Rowv = as.dendrogram(row.clus), 
          Colv = as.dendrogram(col.clus), 
          col = rev(heat.colors(256)), 
          margins = c(6, 1), 
          trace = "none", 
          density.info = "none", 
          lhei = c(1, 4), 
          cexCol = 1.2,
          cexRow =0.5,
          RowSideColors = test[,2]) 
