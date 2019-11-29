source("~/code/single_cell/R/seurat.R", chdir = TRUE)

library(Matrix)
library(Seurat)
library(cellrangerRkit)
library(plyr)

set.seed(0) # Set seed for entire stochastic process.

args = commandArgs(TRUE)
output = args[1]
labelFile = args[2]
feature = args[3]
inputs = args[4:length(args)]

cells = read.table(file = labelFile, header = TRUE, sep = ",")

mat = do.call("cbind", lapply(inputs, loadMatrix))

sMat = seuratPipeline(mat, cells)

df = GetClusters(sMat)
names(df) = c("cell,cluster")

write.csv(df, stdout(), row.names = FALSE, quote = FALSE)

write.csv( sMat@dr$tsne@cell.embeddings
         , file = paste0(output, "_tsne_locations.csv")
         , quote = FALSE
         )

pdf(paste0("./seurat_img/", output, ".pdf"), useDingbats=FALSE)
p = TSNEPlot(object = sMat, do.return = TRUE)
p + theme(aspect.ratio = 1)
dev.off()

pdf(paste0("./seurat_img/label_", output, ".pdf"), useDingbats=FALSE)
p = TSNEPlot(object = sMat, group.by = "label", do.return = TRUE)
p + theme(aspect.ratio = 1)
dev.off()

if (feature != "") {
  pdf(paste0("./seurat_img/feature_", feature, "_", output, ".pdf"), useDingbats=FALSE)
  FeaturePlot( object = sMat
            , features.plot = c(feature)
            , cols.use = c("grey", "red")
            , reduction.use = "tsne"
              )
  dev.off()
}
