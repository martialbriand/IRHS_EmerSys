rm(list = ls())

setwd("~/Desktop/Data_Bee2Seed")
library(metacoder)
library(readr) # Loads the readr package so we can use `read_tsv`
library(dplyr)

otu_data <- read.csv2("asv_count.csv", check.names = F) # You might need to change the path to the file
names(otu_data)[1] <- "OTU ID"
otu_data <- as_tibble(otu_data)
print(otu_data) # You can also enter just `otu_data` to print it

#### Sample data

sample_data <- read.csv2("design_mb1.csv") # each "c" means a column of "character"
sample_data$Year <- factor(sample_data$Year, labels=c("Y2017","Y2018"))
str(sample_data)
names(sample_data)[1] <-"Sample"
sample_data <- sample_data[,-2]
sample_data <- as_tibble(sample_data, col_types = "ccccccccccc") # You can also enter `sample_data` to print it
print(sample_data)

sample_data <- filter(sample_data, Year == "Y2018")
sample_data <- filter(sample_data, Experiment == "Honey")
keep <- c("OTU ID", as.character(sample_data$Sample))
otu_data <- otu_data[,names(otu_data)%in%keep]

tax_data <- read.csv2("asv_tax.csv")
names(tax_data)[1] <- "OTU ID" 
print(tax_data) # You can also enter `tax_data` to print it
tax_data$taxonomy <- ifelse(is.na(tax_data$Kingdom),"Unassigned", 
                            paste("Root;k__", tax_data$Kingdom,";p__",tax_data$Phylum,
                                  ";c__",tax_data$Class,";o__",tax_data$Order,";f__",
                                  tax_data$Family,";g__",tax_data$Genus,sep=""))

print(tax_data) 

tax_data$`OTU ID` <- as.character(tax_data$`OTU ID`) # Must be same type for join to work
otu_data <- left_join(otu_data, tax_data,
                      by = "OTU ID") # identifies cols with shared IDs
print(otu_data)

tail(colnames(otu_data), n = 10) # `tail` returns the last n elements

head(otu_data$taxonomy, 10)

##### MAking tax-map object

library(taxa)
obj <- parse_tax_data(otu_data,
                      class_cols = "taxonomy", # The column in the input table
                      class_sep = ";") # What each taxon is seperated by
print(obj)

## Obj is an r object that has many tibbles (dataframes)

print(obj$data$tax_data)

obj <- parse_tax_data(otu_data,
                      class_cols = "taxonomy",
                      class_sep = ";",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
print(obj)

head(taxon_names(obj))

head(taxon_ranks(obj))

obj$data$class_data <- NULL

names(obj$data) <- "otu_counts"
print(obj)

print(sample_data)

obj <- filter_taxa(obj, taxon_names != "")
print(obj)

head(taxon_names(obj))

head(all_names(obj), 20)

obj <- filter_taxa(obj, taxon_names == "Bacteria", subtaxa = TRUE)

obj <- filter_taxa(obj, !taxon_names == "Cyanobacteria/Chloroplast", subtaxa = FALSE, supertaxa = FALSE,
                   drop_obs = TRUE)
obj <- filter_taxa(obj, !taxon_names == "Chloroplast", subtaxa = FALSE, supertaxa = FALSE,
                   drop_obs = TRUE)
obj <- filter_taxa(obj, !taxon_names == "Streptophyta", subtaxa = FALSE, supertaxa = FALSE,
                   drop_obs = TRUE)
obj <- filter_taxa(obj, !taxon_names == "NA", subtaxa = FALSE, supertaxa = TRUE,
                   drop_obs = TRUE)
print(obj)

## removing taxa with 0 observations

obj$data$otu_counts <- obj$data$otu_counts[c("taxon_id", "OTU ID", as.character(sample_data$Sample))]
tail(obj$data$otu_counts)
print(obj)

n_obs(obj)

has_no_reads <- rowSums(obj$data$otu_counts[, as.character(sample_data$Sample)]) == 0 
sum(has_no_reads)
obj <- filter_obs(obj, "otu_counts", ! has_no_reads, drop_taxa = TRUE) # note the ! negation operator

print(obj)

### plotting 
library(metacoder)

heat_tree(obj,
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs)

map_data(obj, taxon_names, n_obs)

obj %>% 
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")

obj %>% 
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")

obj %>%
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs
            #,output_file = "plot_all.svg"
            )

set.seed(1) # Each number will produce a slightly different result for some layouts
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford"
            #,output_file = "heat_tree_seeds_allHB.svg"
            )

###############################
##### By Material
#### 
### grouping data by sample type
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_counts",
                                       cols = as.character(sample_data$Sample),
                                       groups = as.character(sample_data$Material))


set.seed(1) # Each number will produce a slightly different result for some layouts
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = NA,
            node_size = n_obs,
            node_color = Bee, 
            initial_layout = "re", layout = "da",
            output_file = "heat_tree_Bee_2018HB.svg"
            )
set.seed(1) # Each number will produce a slightly different result for some layouts
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = NA,
            node_size = n_obs,
            node_color = Pollen, 
            initial_layout = "re", layout = "da",
            output_file = "heat_tree_Pollen_2018HB.svg"
  )

set.seed(1) # Each number will produce a slightly different result for some layouts
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = NA,
            node_size = n_obs,
            node_color = Nectar, 
            initial_layout = "re", layout = "da",
            output_file = "heat_tree_Nectar_2018HB.svg"
  )
