source("~/code/single_cell/R/seurat.R", chdir = TRUE)

library(cellrangerRkit)

args = commandArgs(TRUE)
seed = args[1]
n = args[2]
inputLabel = args[3]
labelFile = args[4]
output = args[5]
inputs = args[6:length(args)]

set.seed(seed)

print("Loading files:")
print(inputs)

# Load matrix.
gbm = do.call("cbind", lapply(inputs, loadMatrix))

# Load labels.
cells = read.table(file = labelFile, header = TRUE, sep = ",")

# Subset matrix by celltype.
cellsSubset = cells[cells$label == inputLabel,]
gbmSubset = gbm[,colnames(gbm) %in% cellsSubset$item]

# Subset matrix randomly.
gbmSubset = gbmSubset[,sample(1:ncol(gbmSubset), n, replace = FALSE)]

# Save matrix.
dir.create(output, recursive = TRUE)

# expression matrix
writeMM(gbmSubset, paste0(output, "/matrix.mtx"))

# data frame of genes
outRowDf = data.frame(x = rownames(gbmSubset), y = rownames(gbmSubset))
write.table(outRowDf, file = paste0(output, "/genes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data frame of cell barcodes
write.table(colnames(gbmSubset), file = paste0(output, "/barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
