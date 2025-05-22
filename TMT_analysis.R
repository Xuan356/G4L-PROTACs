library(qPLEXanalyzer)
library(gridExtra)
library(pander)
library(readxl)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(ggrepel)

# human annotation
data(human_anno)

# load data
metadata <- read_excel("peptide.xlsx", sheet = "metadata")
intensities <- read_excel("peptide.xlsx", sheet = "abundance")

# Convert intensities and metadata to MSnset
MSnset_data <- convertToMSnset(intensities, metadata = metadata, indExpData = c(13:30), Sequences = 3, Accessions = 8, rmMissing = FALSE)

non_zero_rows <- apply(exprs(MSnset_data), 1, function(row) all(row != 0))
MSnset_data <- MSnset_data[non_zero_rows,]
MSnset_data <- MSnset_data[which(fData(MSnset_data)[,5]==1 & fData(MSnset_data)[,6]==1),]

pdf('results/peptide_intensity_distribution_displot_before_norm.pdf')
intensityPlot(MSnset_data, title = "Peptide intensity distribution")
dev.off()

pdf('results/peptide_intensity_distribution_boxplot_before_norm.pdf')
intensityBoxplot(MSnset_data, title = "Peptide intensity distribution")
dev.off()

pdf('results/relative_peptide_log_intensity_boxplot_before_norm.pdf')
rliPlot(MSnset_data, title = "Relative Peptide intensity")
dev.off()

MSnset_norm <- normalizeScaling(MSnset_data, scalingFunction = median)

pdf('results/peptide_intensity_distribution_displot_normed.pdf')
intensityPlot(MSnset_norm, title = "Peptide intensity distribution")
dev.off()

pdf('results/peptide_intensity_distribution_boxplot_normed.pdf')
intensityBoxplot(MSnset_norm, title = "Peptide intensity distribution")
dev.off()

pdf('results/relative_peptide_log_intensity_boxplot_normed.pdf')
rliPlot(MSnset_norm, title = "Relative Peptide intensity")
dev.off()

MSnset_Pnorm <- summarizeIntensities(MSnset_norm, summarizationFunction = sum, annotation = human_anno)

write.exprs(MSnset_Pnorm, file = "protein_levels_normed.txt", sep = "\t", row.names=TRUE)

pdf('results/protein_correlation_heatmap.pdf')
corrPlot(MSnset_Pnorm)
dev.off()

pdf('results/protein_PCA.pdf')
pcaPlot(MSnset_Pnorm, labelColumn = "TechRep", pointsize = 3, x.nudge=0.6, labelsize = 3)
dev.off()

pdf('results/protein_hclust.pdf')
hierarchicalPlot(MSnset_Pnorm)
dev.off()

contrasts <- c(protac1_vs_DMSO = "PROTAC1 - DMSO",
               protac2_vs_DMSO = "PROTAC2 - DMSO",
               protac3_vs_DMSO = "PROTAC3 - DMSO",
               protac1_vs_NegProtac = "PROTAC1 - NegProtac",
               protac2_vs_NegProtac = "PROTAC2 - NegProtac",
               protac3_vs_NegProtac = "PROTAC3 - NegProtac")

labels <- c("p1d", "p2d", "p3d", "p1n", "p2n", "p3n")

diffstats <- computeDiffStats(MSnset_Pnorm, contrasts = contrasts)

for (i in 1:length(contrasts))
{
    diffexp <- getContrastResults(diffstats, contrast = contrasts[i], writeFile = TRUE)

    opath <- paste(labels[i], "diff_MA.pdf", sep='_')
    pdf(opath, width = 6, height = 4)
    print(maVolPlot(diffstats, contrast = contrasts[i], plotType = "MA", title = contrasts[i]))
    dev.off()
  
    opath <- paste(labels[i], "diff_volcano.pdf", sep='_')
    pdf(opath, width = 6, height = 4)
    print(maVolPlot(diffstats, contrast = contrasts[i], plotType = "Volcano", title = contrasts[i]))
    dev.off()
}
