library(readr)
healthy_taxa_metadata_non_smokers <- read_csv("healthy_taxa_metadata_non-smokers.csv")
View(healthy_taxa_metadata_non_smokers)

#Set reference group for analysis
healthy_taxa_metadata_non_smokers$Med_group
temp_levels <- factor(healthy_taxa_metadata_non_smokers$Med_group,
                      levels=c("0", "1", "2"))
temp_levels
healthy_taxa_metadata_non_smokers$Med_group <- temp_levels
row.names(healthy_taxa_metadata_non_smokers) <- healthy_taxa_metadata_non_smokers$Sample


##METADATA ANALYSIS

#Descriptive analysis of covariates - sex, age, BMI
shapiro.test(healthy_taxa_metadata_non_smokers$Age)
shapiro.test(healthy_taxa_metadata_non_smokers$BMI)

aggregate(Age ~ Med_group, data = healthy_taxa_metadata_non_smokers, summary)
aggregate(BMI ~ Med_group, data = healthy_taxa_metadata_non_smokers, summary)

Sex_Table <- CrossTable(healthy_taxa_metadata_non_smokers$Sex , 
                        healthy_taxa_metadata_non_smokers$Med_group, 
                        prop.r=TRUE, chisq = TRUE)
kruskal.test(healthy_taxa_metadata_non_smokers$Age, 
             healthy_taxa_metadata_non_smokers$Med_group)
library(PMCMR)
posthoc.kruskal.dunn.test(healthy_taxa_metadata_non_smokers$Age, 
          as.factor(healthy_taxa_metadata_non_smokers$Med_group), "bonferroni")


##TAXA ANALYSIS

#preprep data for analysis.

#create taxa table (remove any extra columns & set row names)
taxa_table <- as.data.frame(healthy_taxa_metadata_non_smokers[,-c(1:12)])
row.names(taxa_table) <- row.names(healthy_taxa_metadata_non_smokers)

#filter samples to a specified sequencing depth (5000)
rowSums(taxa_table)
which(rowSums(taxa_table) < 5000) 
less_5K <- which(rowSums(taxa_table) < 5000) 
length(less_5K)
depth_filtered <- taxa_table[-less_5K,]

remove_rare_features <- function( table , cutoff_pro=0.1, parallel=1 ) {
  if(cutoff_pro==0){
    message("No filtering will be done due to cutoff_pro set to 0")
    return(table)
  }
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  if(parallel <= 1){
    for ( i in 1:nrow(table) ) {
      row_nonzero <- length( which( table[ i , ]  > 0 ) )
      if ( row_nonzero > cutoff ) {
        row2keep <- c( row2keep , i)
      }
    }
    return( table [ row2keep , , drop=F ])}}

#flip table so taxa are rows, samples columns
flip_depth_filtered <- t(depth_filtered)
#the cutoff proportion used to identify rare features
rare_filter_table <- remove_rare_features(flip_depth_filtered,
                                          cutoff_pro = 0.1)
dim(rare_filter_table)
dim(flip_depth_filtered)

#create metadata table with only variables of interest
#(remove extra columns, set row names, filter metadata 
#to match samples in the depth filtered taxa table).
metadata <- as.data.frame(healthy_taxa_metadata_non_smokers[ ,c(1:12)])
row.names(metadata) <- metadata$Sample
filtered_metadata <- metadata[colnames(rare_filter_table),]


#Corncob differential abundance (DA) analysis
library(corncob)
library(phyloseq)
#create the object containing the metadata and taxa to be tested.
otu_tab <- otu_table(rare_filter_table,taxa_are_rows = TRUE)
sam_data <- sample_data(filtered_metadata)
phylo <- phyloseq(otu_tab, sam_data)

#runs DA analysis, returns plot and, DA results
results <- differentialTest(formula = ~Med_group,
                            formula_null = ~1,
                            phi.formula = ~Med_group,
                            phi.formula_null = ~Med_group,
                            data = phylo,
                            fdr_cutoff = 0.1,
                            test = "Wald")
plot(results)
results$p_fdr
results$significant_taxa
results$significant_models

#runs DA analysis while controlling for covariates, returns plot and, DA results
results2 <- differentialTest(formula = ~Med_group + Sex + BMI + Age,
                            formula_null = ~1 + Sex + BMI + Age,
                            phi.formula = ~Med_group + Sex + BMI + Age,
                            phi.formula_null = ~Med_group + Sex + BMI + Age,
                            data = phylo,
                            fdr_cutoff = 0.1,
                            test = "Wald")
plot(results2)
results2$p_fdr
results2$significant_taxa
results2$significant_models


#ALDEx2 DA analysis
#install.packages("BiocManager")
#BiocManager::install("ALDEx2")
library(ALDEx2)

#create model containing the variables to be tested.
matrixmodel <- model.matrix(~Age+Sex+BMI+Med_group, filtered_metadata) 
View(matrixmodel)
#for uncorrected, remove covariates.

#generate Monte Carlo samples of the Dirichlet distribution, performs centred log-ratio transformation.
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
#calculates the expected values for each coefficient of a glm model on the data returned by aldex.clr
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)


#Maaslin2 DA analysis
#BiocManager::install("Maaslin2")
library(Maaslin2)

#finds multivariable associations between taxa and metadata based on GLM 
results <- Maaslin2(rare_filter_table, filtered_metadata, "test", 
        transform = "AST", fixed_effects = c("Age", "Sex", "BMI", "Med_group" ))
View(results[[1]])
Med_group_results <- results[[1]][which(results[[1]]$name=="Med_group1"),]
View(Med_group_results)
Med_group_results2 <- results[[1]][which(results[[1]]$name=="Med_group2"),]
View(Med_group_results2)


#ANCOM2 DA analysis
library(tidyr)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(exactRankTests)
deps = c("exactRankTests", "nlme", "dplyr", "ggplot2", "compositions")
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(dep)
  }
  library(dep, character.only = TRUE)
}
source("~/Desktop/DECLERCQ_et_al_2021_supporting_files/ANCOM2.R")

#preprocessing step 
preprocess <- feature_table_pre_process(feature_table = rare_filter_table, 
                      meta_data = filtered_metadata, sample_var = "Sample", 
                      out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000)
#run main ANCOM function with preprocessed data, adjusting p-values for multiple comparisons
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
          preprocess$structure_zeros,"Med_group", "BH", 0.1, "Age + Sex + BMI")
# to remove covariates, replace with NULL
View(rez[[1]])


# Beta Diversity
#use of previously generated Bray Curtis dissimilarity matrix or Weighted UniFrac distance matrix.
# Bray Curtis
bray_curtis_distance <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/bray_curtis_distance_matrix.tsv",
                      sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
#match metadata samples to bray-curtis samples 
intersect(rownames(healthy_taxa_metadata_non_smokers), 
          rownames(bray_curtis_distance))
samples_to_keep <- intersect(rownames(healthy_taxa_metadata_non_smokers),
                             rownames(bray_curtis_distance))
filtered_bray_curtis_distance <- bray_curtis_distance[samples_to_keep, samples_to_keep]
dim(filtered_bray_curtis_distance)
filtered_metadata <- healthy_taxa_metadata_non_smokers[samples_to_keep,]
row.names(filtered_metadata) <- filtered_metadata$Sample

#adonis test - permutational ANOVA of dissimilarities
set.seed(23)
adonis2(filtered_bray_curtis_distance ~ filtered_metadata$Med_group,
        permutations = 1000, by="margin")

set.seed(23)
adonis2(filtered_bray_curtis_distance ~ filtered_metadata$Med_group
        + filtered_metadata$BMI
        + filtered_metadata$Sex
        + filtered_metadata$Age, 
        permutations = 1000, by="margin")

#plotting Bray Curtis dissimilarity
library(vegan)
library(devtools)
library(ggord)
library(ggplot2)
bray_curtis_pcoa <- cmdscale(filtered_bray_curtis_distance, k=2, eig = TRUE)
barplot(bray_curtis_pcoa$eig[1:10])
component1 <- bray_curtis_pcoa$eig[1]/sum(bray_curtis_pcoa$eig)
component2 <- bray_curtis_pcoa$eig[2]/sum(bray_curtis_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        class=filtered_metadata$Med_group)
bray_plot <- ggplot(plot_data, aes(x=pc1, y=pc2, color=class))+ geom_point()+
  stat_ellipse() + scale_color_manual(values=c("Orange", "Dark Grey", "Blue"), 
                                      name="", labels=c("None", "Single", "Multi"))+
  theme_bw() + xlab("PC1(23.9%)") + ylab("PC2(6.5%)") + theme(legend.position = c(0.078,0.925)) +
  theme(text = element_text(size = 18)) +theme(legend.title = element_blank())
bray_plot


# Weighted Unifrac
weighted_unifrac <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/weighted_unifrac_distance_matrix.tsv",
                      sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(weighted_unifrac))
samples_to_keep <- intersect(rownames(healthy_taxa_metadata_non_smokers), 
                             rownames(weighted_unifrac))
filtered_weighted_unifrac <- weighted_unifrac[samples_to_keep, samples_to_keep]
dim(filtered_weighted_unifrac)
unifrac__metadata <- healthy_taxa_metadata_non_smokers[samples_to_keep,]

set.seed(23)
adonis2(filtered_weighted_unifrac ~ unifrac__metadata$Med_group, 
        permutations = 1000, by="margin")

set.seed(23)
adonis2(filtered_weighted_unifrac ~ unifrac__metadata$Med_group
        + unifrac__metadata$BMI
        + unifrac__metadata$Sex
        + unifrac__metadata$Age, 
        permutations = 1000, by="margin")

#plotting weighted unifrac
weighted_unifrac_pcoa <- cmdscale(filtered_weighted_unifrac, k=2, eig = TRUE)
barplot(weighted_unifrac_pcoa$eig[1:10])
component1 <- weighted_unifrac_pcoa$eig[1]/sum(weighted_unifrac_pcoa$eig)
component2 <- weighted_unifrac_pcoa$eig[2]/sum(weighted_unifrac_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        class=filtered_metadata$Med_group)
unifrac_plot <- ggplot(plot_data, aes(x=pc1, y=pc2, color=class))+ geom_point()+
  stat_ellipse() + scale_color_manual(values=c("Orange", "Dark Grey", "Blue"), 
                                      name="", labels=c("None", "Single", "Multi"))+
  theme_bw() + xlab("PC1(46.9%)") + ylab("PC2(18.0%)") + theme(legend.position = c(0.078,0.925)) +
  theme(text = element_text(size = 18)) +theme(legend.title = element_blank())
unifrac_plot


#Alph Diversity
#Shannon Diversity Index - taxa diversity in a community
shannon <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/shannon.tsv", sep="\t", 
                     header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(shannon))
samples_to_keep1 <- intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(shannon))
filtered_shannon <- data.frame(shannon[samples_to_keep1,])
filter_metadata1 <- healthy_taxa_metadata_non_smokers[samples_to_keep1,]
row.names(filter_metadata1) <- filter_metadata1$Sample

boxplot(filtered_shannon$shannon.samples_to_keep1... ~ filter_metadata1$Med_group, las=1,
        ylab= "Alpha Diversity Measure", xlab ="Medication Group", 
        ylim=c(1,10),
        col=c("Orange", "Dark Grey", "Blue"),names = c("None", "Single", "Multi"),
        par(cex.lab=1.5), par(cex.axis=1.5), par(cex.main=1.5),
        main="Shannon")
shapiro.test(filtered_shannon$shannon.samples_to_keep1...)
kruskal.test(filtered_shannon$shannon.samples_to_keep1...,
          as.factor(filter_metadata1$Med_group))
posthoc.kruskal.dunn.test(filtered_shannon$shannon.samples_to_keep1...,
                          as.factor(filter_metadata1$Med_group), "bonferroni")

#evenness - the uniformity between taxa in a community
evenness <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/evenness.tsv", sep="\t", 
                      header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(evenness))
samples_to_keep1 <- intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(evenness))
filtered_evenness <- data.frame(evenness[samples_to_keep1,])
filter_metadata2 <- healthy_taxa_metadata_non_smokers[samples_to_keep1,]

boxplot(filtered_evenness$evenness.samples_to_keep1... ~ filter_metadata2$Med_group, las=1,
        ylab= "Alpha Diversity Measure", xlab ="Medication Group", 
        ylim=c(0.15,1.0),
        col=c("Orange", "Dark Grey", "Blue"), names = c("None", "Single", "Multi"),
        par(cex.lab=1.5), par(cex.axis=1.5), par(cex.main=1.5),
        main="Evenness")
shapiro.test(filtered_evenness$evenness.samples_to_keep1...)
kruskal.test(filtered_evenness$evenness.samples_to_keep1...,
             as.factor(filter_metadata2$Med_group))
posthoc.kruskal.dunn.test(filtered_evenness$evenness.samples_to_keep1...,
                          as.factor(filter_metadata2$Med_group), "bonferroni")

#richness - number of different taxa in a community
richness <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/richness.tsv", sep="\t", 
                      header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(richness))
samples_to_keep1 <- intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(richness))
filtered_richness <- data.frame(richness[samples_to_keep1,])
filter_metadata3 <- healthy_taxa_metadata_non_smokers[samples_to_keep1,]

boxplot(filtered_richness$richness.samples_to_keep1... ~ filter_metadata3$Med_group, las=1,
        ylab= "Alpha Diversity Measure", xlab ="Medication Group", 
        ylim=c(0,650),
        col=c("Orange", "Dark Grey", "Blue"), names = c("None", "Single", "Multi"),
        par(cex.lab=1.5), par(cex.axis=1.5), par(cex.main=1.5),
        main="Richness")
shapiro.test(filtered_richness$richness.samples_to_keep1...)
kruskal.test(filtered_richness$richness.samples_to_keep1...,
             as.factor(filter_metadata3$Med_group))
posthoc.kruskal.dunn.test(filtered_richness$richness.samples_to_keep1...,
                                       as.factor(filter_metadata3$Med_group), "bonferroni")

#Faith's Phylogenetic diversity - uses phylogentic distances to calculate diversity of community
faiths <-read.table("~/Desktop/DECLERCQ_et_al_2021_supporting_files/faiths_pd.tsv", sep="\t", 
                    header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(faiths))
samples_to_keep1 <- intersect(rownames(healthy_taxa_metadata_non_smokers), rownames(faiths))
filtered_faiths <- data.frame(faiths[samples_to_keep1,])
filter_metadata4 <- healthy_taxa_metadata_non_smokers[samples_to_keep1,]

boxplot(filtered_faiths$faiths.samples_to_keep1... ~ filter_metadata4$Med_group, las=1,
        ylab= "Alpha Diversity Measure", xlab ="Medication Group", 
        ylim=c(0,70),
        col=c("Orange", "Dark Grey", "Blue"), names = c("None", "Single", "Multi"),
        par(cex.lab=1.5), par(cex.axis=1.5), par(cex.main=1.5),
        main="Faiths")
shapiro.test(filtered_faiths$faiths.samples_to_keep1...)
kruskal.test(filtered_faiths$faiths.samples_to_keep1...,
             as.factor(filter_metadata4$Med_group))
posthoc.kruskal.dunn.test(filtered_faiths$faiths.samples_to_keep1...,
                                       as.factor(filter_metadata4$Med_group), "bonferroni")
             
             
             