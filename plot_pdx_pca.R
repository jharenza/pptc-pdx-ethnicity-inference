dir <- "~/snp-array-files/pdx_ethnicity_analysis/"

setwd(dir)
library(ggplot2)

mds <- read.table(paste(dir, "PCA.plink.eigenvec", sep=""),
                  stringsAsFactors = F, head=T)
mds$C1 <- mds$PC1
mds$C2 <- mds$PC2
mds$C3 <- mds$PC3

clin <- read.table(paste(dir, "2019-02-09-pdx-clinical-final-for-paper.txt", sep=""),
                    stringsAsFactors = F, head=T, sep="\t")

# Subset results to models in PPTC
rownames(clin) <- clin$snp.array.sample.ID
mds$In_PPTC <- clin[mds$IID, "Model.Part.of.PPTC"]
mds[grep("^NA", mds$IID), "In_PPTC"] <- "HapMap"
mds <- subset(mds, In_PPTC == "yes" | In_PPTC == "HapMap")  

# Add histology
mds$Histology_detailed <- clin[mds$IID, "Histology.Detailed"]
mds[grep("^NA", mds$IID), "Histology_detailed"] <- "HapMap" 



# Add model
mds$Model <- clin[mds$IID, "Model"]
mds[grep("^NA", mds$IID), "Model"] <- "HapMap" 

# Add ethnicity for HapMap individuals
hapmap_matchups <- read.table(paste(dir, "hapmap_ancestry_info.txt", sep=""), sep='\t', stringsAsFactors = F)
rownames(hapmap_matchups) <- hapmap_matchups$V2
mds$Ethnicity <- hapmap_matchups[mds$IID,]$V7
mds[grep("^PPTC|COG", mds$IID), "Ethnicity"] <- "PDX"

# Add reported ethnicity for PDXs
mds$"Reported_Ethnicity" <- clin[mds$IID, "Ethnicity"]
mds[grep("^NA", mds$IID), "Reported_Ethnicity"] <- "HapMap" 

# Subset to 9 ethnicities (include all African populations, but exclude the highly overlapping East Asian populations)
mds_hapmap9 <- subset(mds, Ethnicity %in% c("CEU", "YRI", "ASW", "CHD", "GIH", "MEX", "TSI", "PDX", "LWK", "MKK"))
mds_hapmap9[mds_hapmap9$Ethnicity != "PDX", "Ethnicity"] <- paste(" HapMap: ", mds_hapmap9[mds_hapmap9$Ethnicity != "PDX", "Ethnicity"], sep="")
  # inserted a space in front of HapMap to help with legend plotting
ethcolors <- c("red4", "orangered", "yellow", "green4", "cyan", "darkblue", "purple", "maroon3", "black", "gray60")

# Combine into 4 general ethnicity groups
mds_hapmap4 <- mds
mds_hapmap4[mds_hapmap4$Ethnicity %in% c("ASW", "LWK", "MKK", "YRI"), "Ethnicity"] <- " HapMap: African"
mds_hapmap4[mds_hapmap4$Ethnicity %in% c("CEU", "TSI"), "Ethnicity"] <- " HapMap: European"
mds_hapmap4[mds_hapmap4$Ethnicity %in% c("CHD", "CHB", "JPT"), "Ethnicity"] <- " HapMap: East Asian"
mds_hapmap4[mds_hapmap4$Ethnicity %in% c("GIH", "MEX"), "Ethnicity"] <- " HapMap: South Asian or Hispanic"

# Read in standard histology color codes
histcolors_data <- read.table(paste(dir, "2019-02-09-all-hist-colors.txt", sep=""), sep='\t', stringsAsFactors = F, comment.char = "")
histcolors <- histcolors_data$V2
names(histcolors) <- histcolors_data$V1

### Plot all HapMap populations and all PDX histologies

mds_hm <- mds_hapmap9[grep("^NA", mds_hapmap9$IID),]
mds_pdx <- mds_hapmap9[grep("^PPTC|COG", mds_hapmap9$IID),]

ethcolors_9 <- c("red4", "orangered", "goldenrod1", "chartreuse2", "green4", "dodgerblue", "cyan", "darkblue", "purple")
names(ethcolors_9) <- c(" HapMap: ASW", " HapMap: CEU", " HapMap: CHD", " HapMap: GIH", " HapMap: LWK", " HapMap: MEX",
                      " HapMap: MKK", " HapMap: TSI", " HapMap: YRI")
colors_combined <- c(histcolors, ethcolors_9)

shape_hm <- 2
shape_pdx <- 1
size_hm <- 1
size_pdx <- 2
alpha_hm <- 0.3
alpha_pdx <- 1

p <- ggplot(data=NULL, aes(x=C1, y=C2)) +
  geom_point(data=mds_hm, aes(color=Ethnicity), shape=shape_hm, size=size_hm, alpha=alpha_hm) +
  geom_point(data=mds_pdx, aes(color=Histology_detailed), shape=shape_pdx, size=size_pdx, alpha=alpha_pdx) +
  scale_color_manual(values=colors_combined)  +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(shape = c(rep(shape_hm, 9), rep(shape_pdx, 29)),
                                                  size = c(rep(size_hm, 9), rep(size_pdx, 29)),
                                                  alpha = c(rep(alpha_hm, 9), rep(alpha_pdx, 29))))) +
  theme(legend.title=element_blank()) + xlab("Principal Component 1") + ylab("Principal Component 2")

ggsave(p, file=paste(dir, "plots/PDX_pca_40kSNPs_252_detailed-populations.pdf", sep=""), width=10, height=6)

### Plot all HapMap populations and reported PDX ethnicities

reportedcolors <- c("red", "orange", "yellow", "green", "blue", "purple", "gray60")
names(reportedcolors) <- c("African American", "European", "Hispanic or Latino", "Mixed", "Non-Hispanic", "Other", "Unknown")
# ethcolors_9 <- c("red4", "orangered", "goldenrod1", "chartreuse2", "green4", "dodgerblue", "cyan", "darkblue", "purple")
# names(ethcolors_9) <- c(" HapMap: ASW", " HapMap: CEU", " HapMap: CHD", " HapMap: GIH", " HapMap: LWK", " HapMap: MEX",
#                         " HapMap: MKK", " HapMap: TSI", " HapMap: YRI")
colors_combined <- c(reportedcolors, ethcolors_9)

shape_hm <- 2
shape_pdx <- 20
size_hm <- 1
size_pdx <- 2
alpha_hm <- 0.3
alpha_pdx <- 1

p <- ggplot(data=NULL, aes(x=C1, y=C2)) +
  geom_point(data=mds_hm, aes(color=Ethnicity), shape=shape_hm, size=size_hm, alpha=alpha_hm) +
  geom_point(data=mds_pdx, aes(color=Reported_Ethnicity), shape=shape_pdx, size=size_pdx, alpha=alpha_pdx) +
  scale_color_manual(values=colors_combined)  +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(shape = c(rep(shape_hm, 9), rep(shape_pdx, 7)),
                                                  size = c(rep(size_hm, 9), rep(size_pdx, 7)),
                                                  alpha = c(rep(alpha_hm, 9), rep(alpha_pdx, 7))))) +
  theme(legend.title=element_blank()) + xlab("Principal Component 1") + ylab("Principal Component 2")

ggsave(p, file=paste(dir, "plots/PDX_pca_40kSNPs_252_reported-ethnicities.pdf", sep=""), width=10, height=6)


### Plot HapMap grouped into 4 general populations (with boxes) and all PDX histologies

ethnicity_coords = read.table("ethnicity_coordinates_40kSNPs.txt", sep='\t', head=T, stringsAsFactors = F)

mds_hm <- mds_hapmap4[grep("^NA", mds_hapmap4$IID),]
mds_pdx <- mds_hapmap4[grep("^PPTC|COG", mds_hapmap4$IID),]

ethcolors_4 <- c("pink", "lightgoldenrod", "lightgreen", "lightskyblue")
names(ethcolors_4) <- c(" HapMap: African", " HapMap: East Asian", " HapMap: European", " HapMap: South Asian or Hispanic")
colors_combined <- c(histcolors, ethcolors_4)

shape_hm <- 2
shape_pdx <- 1
size_hm <- 1
size_pdx <- 2
alpha_hm <- 0.3
alpha_pdx <- 1

p <- ggplot(data=NULL, aes(x=C1, y=C2)) +
  geom_point(data=mds_hm, aes(color=Ethnicity), shape=shape_hm, size=size_hm, alpha=alpha_hm) +
  geom_point(data=mds_pdx, aes(color=Histology_detailed), shape=shape_pdx, size=size_pdx, alpha=alpha_pdx) +  # could change stroke
  scale_color_manual(values=colors_combined)  +
  geom_polygon(ethnicity_coords, mapping=aes(group=Ethnicity), fill=NA, linetype=5, size=0.3,
               color=ethnicity_coords$BoxColor) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(shape = c(rep(shape_hm, 4), rep(shape_pdx, 29)),
                                                  size = c(rep(size_hm, 4), rep(size_pdx, 29)),
                                                  alpha = 1))) +
                                                  #alpha = c(rep(alpha_hm, 4), rep(alpha_pdx, 25))))) +
  theme(legend.title=element_blank()) + xlab("Principal Component 1") + ylab("Principal Component 2")

ggsave(p, file=paste(dir, "plots/PDX_pca_40kSNPs_252_grouped-populations_boxes.pdf", sep=""), width=10, height=6)



### Annotate samples with inferred ethnicity

library(sp)

mds_pdx$Inferred_Ethnicity <- "Mixed_or_Unknown"

annotate_ethnicity <- function(cutoff_coords, pc_results, ethnicity){
  cutoff <- subset(cutoff_coords, Ethnicity==ethnicity)
  res <- point.in.polygon(pc_results$C1, pc_results$C2, cutoff$C1, cutoff$C2)
  pc_results[res==1, "Inferred_Ethnicity"] <- ethnicity
  return(pc_results)
}

mds_pdx <- annotate_ethnicity(ethnicity_coords, mds_pdx, "European")
mds_pdx <- annotate_ethnicity(ethnicity_coords, mds_pdx, "African")
mds_pdx <- annotate_ethnicity(ethnicity_coords, mds_pdx, "EastAsian")
mds_pdx <- annotate_ethnicity(ethnicity_coords, mds_pdx, "SouthAsianOrHispanic")

table(mds_pdx$Inferred_Ethnicity)
table(mds_pdx$Inferred_Ethnicity, mds_pdx$Reported_Ethnicity)

mds_pdx_simple <- mds_pdx[,c("Model", "Histology_detailed", "Reported_Ethnicity", "Inferred_Ethnicity", "PC1", "PC2", "PC3")]
write.table(mds_pdx_simple, "inferred_ethnicities_40kSNPs.txt", quote=F, row.names=F, col.names=T, sep="\t")



### Plot in 3D
library(plotly)
plot_ly(mds_hapmap9, x=~C1, y=~C2, z=~C3, color=~Ethnicity)
