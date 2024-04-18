
#Install the packages, IF YOU NEED TO :)
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# 
##############################################

getwd()
###Set your working directory path
setwd("/Users/Sandra.Gomez/Desktop/Microbiome_analysis/Tar_spot_overtime/Single_reads/R/qiime_results")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("metadata.tsv")
#summary of metadata
str(metadata) 
levels(metadata$`sample-time`)
colnames(metadata)[3] <- "sample.time"
metadata$sample.time.ord = factor(metadata$sample.time, c("Early", "Middle", "Late"))
str(metadata)

#two ways to assign row names
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

jaccard_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")

my_colors <- c(
  '#31a354','#fec44f', '#e34a33','#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
  '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928', 
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD"
)

jaccard_meta <- jaccard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
my_column <- "sample.time.ord"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=my_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

##same graph than before but with centroids

centroids_jc <- aggregate(cbind(PC1,PC2)~get(my_column),jaccard_meta,mean)
colnames(centroids_jc)[1] <- "sample.time.ord"

ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  theme(axis.text.x = element_text(size = 10,),
        axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  labs(y=(paste0("PC2 (", round(100*jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")), 
       x = (paste0("PC1 (", round(100*jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")), 
       title = "PCoA Jaccard")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
scale_color_manual(values=my_colors, name = my_column)
ggsave(paste0("output/Jaccard_", my_column,".pdf"), height=3, width=4.5, device="pdf")

## SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids_Wuni <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)
colnames(centroids_Wuni)[1] <- "sample.time.ord"

ggplot(Uwuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  theme(axis.text.x = element_text(size = 10,),
        axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  labs(y=(paste0("PC2 (", round(100*Uwuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")), 
       x = (paste0("PC1 (", round(100*Uwuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")), 
       title = "PCoA Unweighted Unifract")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  scale_color_manual(values=my_colors, name = my_column)
ggsave(paste0("output/Unweighted_", my_column,".pdf"), height=3, width=4.5, device="pdf")


#Same PERMANOVA as in qiime (beta significance)
weighted_dist_mat<-read_qza("core-metrics-results/weighted_unifrac_distance_matrix.qza")
weighted_dm <- as.matrix(weighted_dist_mat$data) 
rownames(weighted_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub_weighted <- metadata[match(rownames(weighted_dm),metadata$SampleID),]
rownames(weighted_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out_weighted <- adonis2(weighted_dm ~ sample.time.ord, data = metadata_sub_weighted)

write.table(PERMANOVA_out_weighted,"output/Sample_time_Adonis_overall_weighted.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

Sample_time_Pair_weighted <- pairwise.adonis2(weighted_dm ~ sample.time.ord, data = metadata_sub_weighted)
write.table(Sample_time_Pair_weighted,"output/sample_time_Adonis_pairwise_weighted.csv",sep=",", row.names = TRUE) 

